#!/usr/bin/env python
"""Map tier-1 atlas accessions to HOG members (tier-1 one-protein-per-gene reference proteome entries).

~7 % of tier-1 atlas accessions (AFDB v6 models) are not in the gene-centric reference proteomes used as
Ancestra input (redundant TrEMBL entries, alternative gene models; mostly zebrafish, maize, rat, soybean).
Each is mapped to the representative proteome member of the same taxon, by the first method that applies:
  direct            accession is itself a proteome member
  identical_seq     a member of the same proteome has the identical sequence (UniProt REST sequence; for
                    deleted/merged entries the UniParc sequence)
  merged_into_member  inactive entry merged (UniProt) into a proteome member
  gene_name+seq     same primary gene name (unique in the proteome) AND DIAMOND identity >= min_identity,
                    query coverage >= min_cov
  seq               DIAMOND best hit in the same proteome, identity >= min_identity, query cov >= min_cov
  gene_name_only    same unique primary gene name, sequence evidence below threshold (flagged; low confidence)
  unmapped          none of the above (excluded from orthology-based analyses, counted)
Sequences / gene names of non-members come from UniProt REST (`uniprotkb/accessions`, current release).

Output: <out_dir>/idgap_mapping.parquet  (atlas_accession, taxid, hog_member_accession, method, pident, qcov,
        gene_primary, n_pockets_kept) for EVERY tier-1 atlas accession, + idgap_summary.tsv.
Usage: .venv/bin/python scripts/orthology/map_idgap.py [--config config/orthology.yaml]
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import os
import re
import subprocess
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor

import duckdb
import pandas as pd

from _common import ANCESTRA_BIN, acc_of_header, atlas_accessions, load_config, log, read_fasta

FIELDS = "accession,gene_primary,organism_id,uniparc_id,length,sequence"


def fetch_batch(accs: list[str]) -> pd.DataFrame:
    url = ("https://rest.uniprot.org/uniprotkb/accessions?" +
           urllib.parse.urlencode({"accessions": ",".join(accs), "fields": FIELDS, "format": "tsv", "size": 500}))
    for i in range(5):
        try:
            with urllib.request.urlopen(url, timeout=300) as r:
                return pd.read_csv(io.StringIO(r.read().decode()), sep="\t", dtype=str)
        except Exception as e:  # noqa: BLE001
            log(f"retry {i}: {e}")
    raise RuntimeError("UniProt fetch failed")


def fetch_inactive(acc: str):
    for i in range(5):
        try:
            with urllib.request.urlopen(f"https://rest.uniprot.org/uniprotkb/{acc}.json", timeout=120) as r:
                d = json.loads(r.read())
            ir = d.get("inactiveReason") or {}
            mt = ir.get("mergeDemergeTo") or []
            return (acc, ir.get("inactiveReasonType"), mt[0] if len(mt) == 1 else None,
                    (d.get("extraAttributes") or {}).get("uniParcId"))
        except Exception as e:  # noqa: BLE001
            log(f"retry {acc} {i}: {e}")
    return (acc, None, None, None)


def fetch_uniparc(upis: list[str]) -> pd.DataFrame:
    q = " OR ".join(f"upi:{u}" for u in upis)
    url = "https://rest.uniprot.org/uniparc/stream?" + urllib.parse.urlencode(
        {"query": q, "fields": "upi,sequence", "format": "tsv"})
    for i in range(5):
        try:
            with urllib.request.urlopen(url, timeout=300) as r:
                return pd.read_csv(io.StringIO(r.read().decode()), sep="\t", dtype=str).rename(
                    columns={"Entry": "upi", "Sequence": "seq"})
        except Exception as e:  # noqa: BLE001
            log(f"retry uniparc {i}: {e}")
    raise RuntimeError("UniParc fetch failed")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config")
    args = ap.parse_args()
    cfg = load_config(args.config)
    W = os.path.join(cfg["work_dir"], "idgap")
    os.makedirs(W, exist_ok=True)
    out_dir = cfg["out_dir"]

    t1 = pd.read_csv(cfg["tier1_proteomes"], sep="\t", header=None,
                     names=["taxid", "proteome", "ptype", "species", "superregnum", "genes", "proteins"])
    # proteome members (accession, taxid, gene, seq md5)
    mem, seqs = [], {}
    for r in t1.itertuples(index=False):
        f = os.path.join(cfg["tier1_fasta_dir"], f"{r.proteome}_{r.taxid}.fasta.gz")
        for h, s in read_fasta(f):
            a = acc_of_header(h)
            g = re.search(r" GN=(\S+)", h)
            mem.append((a, int(r.taxid), g.group(1) if g else None, hashlib.md5(s.encode()).hexdigest()))
            seqs[a] = (int(r.taxid), s)
    mem = pd.DataFrame(mem, columns=["accession", "taxid", "gene", "md5"])
    members = set(mem.accession)

    acc_pq = atlas_accessions(cfg)
    con = duckdb.connect()
    atl = con.execute(f"""SELECT accession atlas_accession, any_value(taxid) taxid, max(n_pockets_kept) n_pockets_kept
                          FROM '{acc_pq}' WHERE tier=1 GROUP BY 1""").df()
    miss = atl[~atl.atlas_accession.isin(members)]
    log(f"tier-1 atlas accessions {len(atl)}; not proteome members {len(miss)}")

    # UniProt info for missing accessions (cached)
    cache = os.path.join(W, "uniprot_missing.tsv")
    if not os.path.exists(cache):
        al = sorted(miss.atlas_accession)
        batches = [al[i:i + 500] for i in range(0, len(al), 500)]
        with ThreadPoolExecutor(4) as ex:
            dfs = list(ex.map(fetch_batch, batches))
        pd.concat(dfs).to_csv(cache, sep="\t", index=False)
    up = pd.read_csv(cache, sep="\t", dtype=str).rename(columns={
        "Entry": "atlas_accession", "Gene Names (primary)": "gene_primary", "Sequence": "seq"})
    up = up.drop_duplicates("atlas_accession")
    log(f"UniProt returned {len(up)} active entries / {len(miss)} missing accessions")
    # inactive (deleted / merged) entries: per-accession JSON -> UniParc id (+ merge target), then UniParc sequence
    cache2 = os.path.join(W, "uniprot_inactive.tsv")
    if not os.path.exists(cache2):
        ina = sorted(set(miss.atlas_accession) - set(up.atlas_accession))
        with ThreadPoolExecutor(8) as ex:
            info = list(ex.map(fetch_inactive, ina))
        inf = pd.DataFrame(info, columns=["atlas_accession", "inactive_reason", "merged_to", "upi"])
        upis = sorted(inf.upi.dropna().unique())
        with ThreadPoolExecutor(4) as ex:
            sq = pd.concat(list(ex.map(fetch_uniparc, [upis[i:i + 100] for i in range(0, len(upis), 100)])))
        inf = inf.merge(sq, on="upi", how="left")
        inf.to_csv(cache2, sep="\t", index=False)
    inf = pd.read_csv(cache2, sep="\t", dtype=str)
    log(f"inactive: {len(inf)}; reasons {inf.inactive_reason.value_counts().to_dict()}; with sequence {inf.seq.notna().sum()}")
    up = pd.concat([up, inf[["atlas_accession", "seq"]].assign(gene_primary=None)], ignore_index=True)
    merged_to = {r.atlas_accession: r.merged_to for r in inf.itertuples() if isinstance(r.merged_to, str)}
    q = miss.merge(up[["atlas_accession", "gene_primary", "seq"]], on="atlas_accession", how="left")
    q["md5"] = q.seq.map(lambda s: hashlib.md5(s.encode()).hexdigest() if isinstance(s, str) else None)

    # 1. identical sequence within the same taxon
    ident = q.merge(mem[["accession", "taxid", "md5"]], on=["taxid", "md5"], how="inner") \
             .sort_values("accession").drop_duplicates("atlas_accession")
    res = {r.atlas_accession: (r.accession, "identical_seq", 100.0, 1.0) for r in ident.itertuples()}
    for a, m in merged_to.items():   # merged into a proteome member (UniProt merge), if not already identical
        if a not in res and m in members:
            res[a] = (m, "merged_into_member", None, None)

    # 2. DIAMOND (per-taxon restriction applied after search: hits outside the taxon are ignored)
    rest = q[~q.atlas_accession.isin(res) & q.seq.notna()]
    qfa, dfa, hits = os.path.join(W, "q.fa"), os.path.join(W, "db.fa"), os.path.join(W, "hits.tsv")
    with open(qfa, "w") as fh:
        for r in rest.itertuples():
            fh.write(f">{r.atlas_accession}\n{r.seq}\n")
    taxa = set(rest.taxid)
    with open(dfa, "w") as fh:
        for a, (t, s) in seqs.items():
            if t in taxa:
                fh.write(f">{a}\n{s}\n")
    dmnd = os.path.join(ANCESTRA_BIN, "diamond")
    subprocess.run([dmnd, "makedb", "--in", dfa, "-d", os.path.join(W, "db"), "--quiet"], check=True)
    subprocess.run([dmnd, "blastp", "-q", qfa, "-d", os.path.join(W, "db"), "-o", hits, "--quiet",
                    "--threads", str(cfg["placement"]["threads"]), "--sensitive", "-k", "25", "-e", "1e-5",
                    "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "qlen", "slen", "evalue", "bitscore",
                    "qstart", "qend"], check=True)
    h = pd.read_csv(hits, sep="\t", header=None,
                    names=["q", "s", "pident", "length", "qlen", "slen", "evalue", "bitscore", "qstart", "qend"])
    h["qcov"] = (h.qend - h.qstart + 1) / h.qlen
    tax_of = {a: t for a, (t, _) in seqs.items()}
    qtax = dict(zip(rest.atlas_accession, rest.taxid))
    h = h[h.s.map(tax_of) == h.q.map(qtax)]
    h = h.sort_values(["q", "bitscore", "s"], ascending=[True, False, True])
    best = h.drop_duplicates("q").set_index("q")
    # gene-name candidates: unique primary gene name within the taxon
    gcount = mem.dropna(subset=["gene"]).groupby(["taxid", "gene"]).accession.agg(list)
    ic = cfg["idgap"]
    for r in rest.itertuples():
        a = r.atlas_accession
        gm = None
        if isinstance(r.gene_primary, str):
            lst = gcount.get((r.taxid, r.gene_primary))
            if lst is not None and len(lst) == 1:
                gm = lst[0]
        if gm is not None:
            hh = h[(h.q == a) & (h.s == gm)]
            if len(hh) and hh.pident.iloc[0] >= ic["min_identity"] and hh.qcov.iloc[0] >= ic["min_cov"]:
                res[a] = (gm, "gene_name+seq", float(hh.pident.iloc[0]), float(hh.qcov.iloc[0]))
                continue
        if a in best.index:
            b = best.loc[a]
            if b.pident >= ic["min_identity"] and b.qcov >= ic["min_cov"]:
                res[a] = (b.s, "seq", float(b.pident), float(b.qcov))
                continue
        if gm is not None:
            hh = h[(h.q == a) & (h.s == gm)]
            res[a] = (gm, "gene_name_only", float(hh.pident.iloc[0]) if len(hh) else None,
                      float(hh.qcov.iloc[0]) if len(hh) else None)

    rows = []
    gene_of = dict(zip(mem.accession, mem.gene))
    qgene = dict(zip(q.atlas_accession, q.gene_primary))
    for r in atl.itertuples():
        a = r.atlas_accession
        if a in members:
            rows.append((a, r.taxid, a, "direct", 100.0, 1.0, gene_of.get(a), r.n_pockets_kept))
        elif a in res:
            m, meth, pi, qc = res[a]
            rows.append((a, r.taxid, m, meth, pi, qc, qgene.get(a), r.n_pockets_kept))
        else:
            rows.append((a, r.taxid, None, "unmapped", None, None, qgene.get(a), r.n_pockets_kept))
    out = pd.DataFrame(rows, columns=["atlas_accession", "taxid", "hog_member_accession", "method", "pident", "qcov",
                                      "gene_primary", "n_pockets_kept"])
    out.to_parquet(os.path.join(out_dir, "idgap_mapping.parquet"), index=False)
    s = out.groupby("method").agg(n=("atlas_accession", "size"),
                                  n_with_pockets=("n_pockets_kept", lambda x: int((x > 0).sum()))).reset_index()
    s.to_csv(os.path.join(out_dir, "idgap_summary.tsv"), sep="\t", index=False)
    pt = out[out.method != "direct"].pivot_table(index="taxid", columns="method", values="atlas_accession",
                                                  aggfunc="count", fill_value=0)
    pt.to_csv(os.path.join(out_dir, "idgap_by_taxon.tsv"), sep="\t")
    many = out[out.method != "direct"].hog_member_accession.value_counts()
    log(s.to_string())
    log(json.dumps({"members_receiving_>1_atlas_acc": int((many > 1).sum()), "max": int(many.max()) if len(many) else 0}))


if __name__ == "__main__":
    main()
