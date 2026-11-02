#!/usr/bin/env python
"""UniProtKB annotations for all tier-1/2 accessions of an atlas snapshot.

Download (cursor-paginated REST search, JSON with evidences, gzip) per query:
  * tier 1: one query per reference proteome  `proteome:UPxxxxxxxxx`
  * tier 2: `reviewed:true`, partitioned by taxon (REVIEWED_PARTS + a NOT-remainder) for parallel download
  * leftovers (accessions of our set not returned above): `accession:(A OR B ...)` in batches,
    then `sec_acc:` lookups (merged/demerged entries); still-missing accessions are reported.
Raw pages are stored as JSON-lines (one UniProt entry per line) under
  /tmp/pe_big/raw/uniprot/<release>/<query_tag>.jsonl.gz      (untrusted data; parsed only as JSON)
The release is taken from the `X-UniProt-Release` response header and stored in every output.

Outputs (data/annotations/):
  uniprot_proteins.parquet       one row per atlas accession (accession = atlas accession, uniprot_acc = primary)
  uniprot_site_features.parquet  Binding site / Active site / Site / Mutagenesis features, one row each,
                                 with ligand, ChEBI, ECO codes, experimental flag (ECO:0000269/ECO:0007744 only)
  uniprot_meta.json              release, queries, raw-file sha256, coverage

Usage:
  .venv/bin/python -I scripts/annotations/uniprot.py --accessions data/annotations/accessions.parquet
  (add --parse-only to rebuild parquet from existing raw files)
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import gzip
import hashlib
import json
import os
import re
import sys
import time

import pandas as pd
import requests

BASE = "https://rest.uniprot.org/uniprotkb/search"
FIELDS = ("accession,id,gene_primary,protein_name,organism_id,length,sequence,sequence_version,fragment,"
          "protein_existence,annotation_score,reviewed,ec,cc_cofactor,"
          "xref_pfam,xref_panther,xref_interpro,xref_supfam,xref_gene3d,xref_pdb,"
          "ft_binding,ft_act_site,ft_site,ft_mutagen,ft_signal,ft_transmem,ft_region")
EXPERIMENTAL_ECO = {"ECO:0000269", "ECO:0007744"}
SITE_TYPES = {"Binding site", "Active site", "Site", "Mutagenesis"}
RAW_ROOT = "/tmp/pe_big/raw/uniprot"
# disjoint taxon partitions of Swiss-Prot (remainder query = reviewed AND NOT any of these)
REVIEWED_PARTS = {"pseudomonadota": 1224, "bacillota": 1239, "actinomycetota": 201174, "bacteria_other": None,
                  "archaea": 2157, "viruses": 10239, "metazoa": 33208, "fungi": 4751, "viridiplantae": 33090}


def reviewed_partitions():
    named = {k: v for k, v in REVIEWED_PARTS.items() if v is not None}
    phyla = [1224, 1239, 201174]
    out = [(k, f"reviewed:true AND taxonomy_id:{v}") for k, v in named.items()]
    out.append(("bacteria_other", "reviewed:true AND taxonomy_id:2 " + " ".join(f"NOT taxonomy_id:{t}" for t in phyla)))
    out.append(("rest", "reviewed:true " + " ".join(f"NOT taxonomy_id:{t}" for t in [2] + list(named.values()))))
    return out


def session():
    s = requests.Session()
    s.headers["User-Agent"] = "pocket-evolution-annotations/1.0 (ETH Zurich research)"
    s.headers["Accept-Encoding"] = "gzip"
    return s


def get(s, url, params=None, tries=8):
    for i in range(tries):
        try:
            r = s.get(url, params=params, timeout=300)
            if r.status_code == 200:
                return r
            if r.status_code in (400, 404):
                raise RuntimeError(f"{r.status_code} {r.text[:300]}")
        except (requests.ConnectionError, requests.Timeout, requests.exceptions.ChunkedEncodingError) as e:
            print("retry", i, e, file=sys.stderr)
        time.sleep(min(60, 2 ** i))
    raise RuntimeError(f"failed {url} {params}")


def fetch_query(query: str, out_path: str, size=500):
    """Cursor-paginate a search; write entries as JSON lines. Returns (n, release)."""
    if os.path.exists(out_path) and os.path.exists(out_path + ".done"):
        return json.load(open(out_path + ".done"))
    s = session()
    url, params = BASE, {"query": query, "format": "json", "fields": FIELDS, "size": size}
    n, release, total = 0, None, None
    tmp = out_path + ".part"
    with gzip.open(tmp, "wt") as fh:
        while url:
            r = get(s, url, params)
            release = release or r.headers.get("X-UniProt-Release")
            total = total or r.headers.get("X-Total-Results")
            for e in rjson(r)["results"]:
                fh.write(json.dumps(e, separators=(",", ":")) + "\n")
                n += 1
            m = re.search(r'<([^>]+)>;\s*rel="next"', r.headers.get("Link", ""))
            url, params = (m.group(1), None) if m else (None, None)
    os.replace(tmp, out_path)
    info = {"query": query, "n": n, "total": int(total or 0), "release": release}
    json.dump(info, open(out_path + ".done", "w"))
    return info


def rjson(r):
    """Parse a JSON response; UniProt cursor pages are sometimes gzip-compressed twice."""
    b = r.content
    while b[:2] == b"\x1f\x8b":
        b = gzip.decompress(b)
    return json.loads(b)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


# ---------------------------------------------------------------- parsing
def _xref_ids(e, db):
    return [x["id"] for x in e.get("uniProtKBCrossReferences", []) if x["database"] == db]


def _gene3d_supfam(e, db):
    return [x["id"] for x in e.get("uniProtKBCrossReferences", []) if x["database"] == db]


def _protein_name(pd_):
    try:
        return pd_["recommendedName"]["fullName"]["value"]
    except KeyError:
        try:
            return pd_["submissionNames"][0]["fullName"]["value"]
        except (KeyError, IndexError):
            return None


def _ecs(pd_):
    out = []
    for key in ("recommendedName",):
        out += [x["value"] for x in pd_.get(key, {}).get("ecNumbers", [])]
    for alt in pd_.get("alternativeNames", []) + pd_.get("submissionNames", []):
        out += [x["value"] for x in alt.get("ecNumbers", [])]
    return sorted(set(out))


def _segments(features, ftype, desc=None):
    segs = []
    for f in features:
        if f["type"] == ftype and (desc is None or f.get("description") == desc):
            st, en = f["location"]["start"].get("value"), f["location"]["end"].get("value")
            if st is not None and en is not None:
                segs.append(f"{st}-{en}")
    return ";".join(segs)


def parse_entry(e):
    pdesc = e.get("proteinDescription", {})
    feats = e.get("features", [])
    seq = e.get("sequence", {})
    cof = []
    for c in e.get("comments", []):
        if c.get("commentType") == "COFACTOR":
            for x in c.get("cofactors", []):
                chebi = (x.get("cofactorCrossReference") or {}).get("id")
                cof.append(f"{x.get('name')}|{chebi}")
    pdb = _xref_ids(e, "PDB")
    prot = {
        "uniprot_acc": e["primaryAccession"],
        "secondary_accs": ";".join(e.get("secondaryAccessions", [])),
        "uniprot_id": e.get("uniProtkbId"),
        "reviewed": "Swiss-Prot" in e.get("entryType", ""),
        "gene_primary": (e.get("genes") or [{}])[0].get("geneName", {}).get("value"),
        "protein_name": _protein_name(pdesc),
        "organism_id": e.get("organism", {}).get("taxonId"),
        "length": seq.get("length"),
        "seq_crc64": seq.get("crc64"),
        "seq_md5": seq.get("md5"),
        "sequence_version": e.get("entryAudit", {}).get("sequenceVersion"),
        "fragment": pdesc.get("flag") in ("Fragment", "Fragments", "Precursor;Fragment", "Fragment,Precursor")
                    or "Fragment" in str(pdesc.get("flag", "")),
        "protein_existence": e.get("proteinExistence"),
        "annotation_score": e.get("annotationScore"),
        "ec": ";".join(_ecs(pdesc)),
        "cofactors": ";".join(cof),
        "pfam": ";".join(_xref_ids(e, "Pfam")),
        "panther": ";".join(_xref_ids(e, "PANTHER")),
        "interpro": ";".join(_xref_ids(e, "InterPro")),
        "supfam": ";".join(_xref_ids(e, "SUPFAM")),
        "gene3d": ";".join(_xref_ids(e, "Gene3D")),
        "n_pdb": len(pdb),
        "pdb_ids": ";".join(pdb),
        "signal": _segments(feats, "Signal"),
        "transmem": _segments(feats, "Transmembrane"),
        "disordered": _segments(feats, "Region", "Disordered"),
    }
    sites = []
    for f in feats:
        if f["type"] not in SITE_TYPES:
            continue
        loc = f["location"]
        ecos = sorted({ev.get("evidenceCode") for ev in f.get("evidences", []) if ev.get("evidenceCode")})
        lig = f.get("ligand") or {}
        chebi = lig.get("id")
        if chebi and chebi.startswith("ChEBI:"):
            chebi = chebi[len("ChEBI:"):]
        sites.append({
            "uniprot_acc": e["primaryAccession"],
            "type": f["type"],
            "start": loc["start"].get("value"),
            "end": loc["end"].get("value"),
            "start_modifier": loc["start"].get("modifier"),
            "end_modifier": loc["end"].get("modifier"),
            "description": f.get("description") or None,
            "ligand": lig.get("name"),
            "chebi": chebi,
            "ligand_label": lig.get("label"),
            "ligand_note": lig.get("note"),
            "ligand_part": (f.get("ligandPart") or {}).get("name"),
            "alt_seq": ";".join(f.get("alternativeSequence", {}).get("alternativeSequences", [])) or None,
            "eco_codes": ";".join(ecos),
            "evidence_refs": ";".join(sorted({f"{ev.get('source')}:{ev.get('id')}" for ev in f.get("evidences", [])
                                              if ev.get("source")})),
            "experimental": bool(EXPERIMENTAL_ECO & set(ecos)),
        })
    return prot, sites


def iter_raw(path):
    with gzip.open(path, "rt") as fh:
        for line in fh:
            yield json.loads(line)


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--accessions", default="data/annotations/accessions.parquet")
    ap.add_argument("--outdir", default="data/annotations")
    ap.add_argument("--raw-root", default=RAW_ROOT)
    ap.add_argument("--seq-out", default="/tmp/pe_big/work/annotations/uniprot_sequences.parquet",
                    help="full sequences (kept outside the repo; used by biolip.py for residue mapping)")
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--parse-only", action="store_true")
    a = ap.parse_args()

    acc = pd.read_parquet(a.accessions)
    wanted = set(acc.accession)
    s = session()
    r = get(s, BASE, {"query": "accession:P00918", "format": "json", "fields": "accession", "size": 1})
    release = r.headers.get("X-UniProt-Release")
    raw = os.path.join(a.raw_root, release)
    os.makedirs(raw, exist_ok=True)
    print("UniProt release", release, "raw dir", raw)

    jobs = [(f"proteome:{p}", os.path.join(raw, f"proteome_{p}.jsonl.gz"))
            for p in sorted(acc.proteome.dropna().unique())]
    if (acc.tier == 2).any():
        for tag, q in reviewed_partitions():
            jobs.append((q, os.path.join(raw, f"reviewed_{tag}.jsonl.gz")))
    infos = []
    if not a.parse_only:
        with cf.ThreadPoolExecutor(a.workers) as ex:
            futs = {ex.submit(fetch_query, q, p): q for q, p in jobs}
            for f in cf.as_completed(futs):
                info = f.result()
                infos.append(info)
                print(info, flush=True)
    else:
        infos = [json.load(open(p + ".done")) for _, p in jobs]

    # which accessions are covered (primary or secondary)?
    def covered_from(paths):
        got = set()
        for p in paths:
            for e in iter_raw(p):
                for x in [e["primaryAccession"]] + e.get("secondaryAccessions", []):
                    if x in wanted:
                        got.add(x)
        return got

    got = covered_from([p for _, p in jobs])
    missing = sorted(wanted - got)
    print("missing after main queries:", len(missing), flush=True)
    left_path = os.path.join(raw, "leftover_accessions.jsonl.gz")
    if not a.parse_only and missing and not os.path.exists(left_path + ".done"):
        tmp = left_path + ".part"
        n = 0
        with gzip.open(tmp, "wt") as fh:
            for i in range(0, len(missing), 100):
                chunk = missing[i:i + 100]
                for field in ("accession", "sec_acc"):
                    q = f"{field}:(" + " OR ".join(chunk) + ")"
                    url, params = BASE, {"query": q, "format": "json", "fields": FIELDS, "size": 500}
                    while url:
                        rr = get(s, url, params)
                        for e in rjson(rr)["results"]:
                            fh.write(json.dumps(e, separators=(",", ":")) + "\n")
                            n += 1
                        m = re.search(r'<([^>]+)>;\s*rel="next"', rr.headers.get("Link", ""))
                        url, params = (m.group(1), None) if m else (None, None)
        os.replace(tmp, left_path)
        json.dump({"query": "leftover accession/sec_acc batches", "n": n, "release": release},
                  open(left_path + ".done", "w"))
    paths = [p for _, p in jobs] + ([left_path] if os.path.exists(left_path) else [])

    # parse: keep only entries that map to our accessions; map atlas accession -> entry
    prots, sites, seen, seqs = [], [], set(), []
    acc2primary = {}
    for p in paths:
        for e in iter_raw(p):
            ids = [e["primaryAccession"]] + e.get("secondaryAccessions", [])
            hits = [x for x in ids if x in wanted]
            if not hits or e["primaryAccession"] in seen:
                for h in hits:
                    acc2primary.setdefault(h, e["primaryAccession"])
                continue
            seen.add(e["primaryAccession"])
            pr, st = parse_entry(e)
            prots.append(pr)
            seqs.append((e["primaryAccession"], e.get("sequence", {}).get("value")))
            sites.extend(st)
            for h in hits:
                # primary-accession match wins over secondary
                if h == e["primaryAccession"] or h not in acc2primary:
                    acc2primary[h] = e["primaryAccession"]
    P = pd.DataFrame(prots)
    os.makedirs(os.path.dirname(a.seq_out), exist_ok=True)
    pd.DataFrame(seqs, columns=["uniprot_acc", "sequence"]).to_parquet(a.seq_out, index=False, compression="zstd")
    amap = pd.DataFrame({"accession": list(acc2primary), "uniprot_acc": list(acc2primary.values())})
    tiers = acc.groupby("accession").tier.apply(lambda t: ";".join(map(str, sorted(t)))).rename("tiers").reset_index()
    out = tiers.merge(amap, on="accession", how="left").merge(P, on="uniprot_acc", how="left")
    out["uniprot_release"] = release
    out["accession_status"] = "primary"
    out.loc[out.uniprot_acc.isna(), "accession_status"] = "not_found"
    out.loc[out.uniprot_acc.notna() & (out.uniprot_acc != out.accession), "accession_status"] = "secondary"
    for c in ("organism_id", "length", "sequence_version", "n_pdb"):
        out[c] = out[c].astype("Int64")
    out.to_parquet(os.path.join(a.outdir, "uniprot_proteins.parquet"), index=False, compression="zstd")
    S = pd.DataFrame(sites)
    S = S[S.uniprot_acc.isin(set(out.uniprot_acc.dropna()))]
    S["start"] = S.start.astype("Int64")
    S["end"] = S.end.astype("Int64")
    S["uniprot_release"] = release
    S.to_parquet(os.path.join(a.outdir, "uniprot_site_features.parquet"), index=False, compression="zstd")

    meta = {
        "uniprot_release": release,
        "fields": FIELDS,
        "queries": infos,
        "raw_files": {os.path.basename(p): sha256(p) for p in paths},
        "n_accessions": int(len(out)),
        "n_found": int(out.uniprot_acc.notna().sum()),
        "n_secondary": int((out.accession_status == "secondary").sum()),
        "not_found": sorted(out.loc[out.uniprot_acc.isna(), "accession"].tolist())[:5000],
        "site_rows": int(len(S)),
        "site_rows_by_type": S.type.value_counts().to_dict(),
        "site_rows_experimental_by_type": S[S.experimental].type.value_counts().to_dict(),
    }
    json.dump(meta, open(os.path.join(a.outdir, "uniprot_meta.json"), "w"), indent=1)
    print({k: v for k, v in meta.items() if k not in ("queries", "raw_files", "not_found")})


if __name__ == "__main__":
    main()
