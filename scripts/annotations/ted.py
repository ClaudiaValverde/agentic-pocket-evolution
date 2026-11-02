#!/usr/bin/env python
"""TED (The Encyclopedia of Domains) domain assignments with CATH labels for atlas accessions.

Raw (Zenodo record 13908086, v5, CC BY 4.0; to /tmp/pe_big/raw/ted/):
  ted_365m.domain_summary.cath.globularity.taxid.tsv.gz  (19.9 GB, no header, 21 columns; md5 from Zenodo API)

TED was computed on AFDB v4 models; the atlas uses AFDB v6 (tiers 1-2). Domain sequences are checked against
current UniProt sequences: md5(concatenated segment residues) == md5_domain  ->  seq_matches_uniprot.

Output: data/annotations/ted_domains.parquet
  accession, ted_id, domain_no, consensus_level, chopping, nres_domain, num_segments, plddt, cath_label,
  cath_level (H/T), cath_method, packing_density, norm_rg, globular, seq_matches_uniprot
and ted_meta.json.

Usage: .venv/bin/python -I scripts/annotations/ted.py --accessions data/annotations/accessions.parquet
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess

import duckdb
import pandas as pd

RECORD = "https://zenodo.org/api/records/13908086"
FNAME = "ted_365m.domain_summary.cath.globularity.taxid.tsv.gz"
RAW = "/tmp/pe_big/raw/ted"
COLS = ["ted_id", "md5_domain", "consensus_level", "chopping", "nres_domain", "num_segments", "plddt",
        "num_helix_strand_turn", "num_helix", "num_strand", "num_helix_strand", "num_turn", "proteome_id",
        "cath_label", "cath_level", "cath_method", "packing_density", "norm_rg", "tax_common_name",
        "tax_scientific_name", "tax_lineage"]


def filehash(p, algo):
    h = hashlib.new(algo)
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def domain_md5(seq, chopping):
    if not isinstance(seq, str) or not chopping:
        return None
    parts = []
    for seg in chopping.split("_"):
        s, e = seg.split("-")
        parts.append(seq[int(s) - 1:int(e)])
    return hashlib.md5("".join(parts).encode()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--accessions", default="data/annotations/accessions.parquet")
    ap.add_argument("--sequences", default="/tmp/pe_big/work/annotations/uniprot_sequences.parquet")
    ap.add_argument("--proteins", default="data/annotations/uniprot_proteins.parquet")
    ap.add_argument("--outdir", default="data/annotations")
    ap.add_argument("--threads", type=int, default=16)
    a = ap.parse_args()
    os.makedirs(RAW, exist_ok=True)
    rec_path = os.path.join(RAW, "record_13908086.json")
    if not os.path.exists(rec_path):
        subprocess.run(["curl", "-sSf", "-o", rec_path, RECORD], check=True)
    rec = json.load(open(rec_path))
    finfo = next(f for f in rec["files"] if f["key"] == FNAME)
    raw = os.path.join(RAW, FNAME)
    if not os.path.exists(raw):
        subprocess.run(["curl", "-sSfL", "-o", raw, f"{RECORD}/files/{FNAME}/content"], check=True)
    got_md5 = filehash(raw, "md5")
    assert "md5:" + got_md5 == finfo["checksum"], (got_md5, finfo["checksum"])

    acc = pd.read_parquet(a.accessions, columns=["accession"]).drop_duplicates()
    con = duckdb.connect()
    con.execute(f"set threads={a.threads}")
    con.register("acc", acc)
    colspec = "{" + ",".join(f"'{c}':'VARCHAR'" for c in COLS) + "}"
    cache = os.path.join("/tmp/pe_big/work/annotations", f"ted_subset_{got_md5[:8]}_{len(acc)}.parquet")
    os.makedirs(os.path.dirname(cache), exist_ok=True)
    if os.path.exists(cache):
        T = pd.read_parquet(cache)
    else:
        T = _scan(con, raw, colspec)
        T.to_parquet(cache, index=False)
    T = _post(T, a)
    _finish(T, a, acc, got_md5, raw)


def _scan(con, raw, colspec):
    return con.execute(f"""
        with t as (
          select *, split_part(ted_id, '-', 2) as accession
          from read_csv('{raw}', delim='\t', header=false, columns={colspec}, quote='', escape='')
        )
        select t.accession, ted_id,
               try_cast(regexp_extract(ted_id, 'TED(\\d+)$', 1) as int) as domain_no,
               consensus_level, chopping, nres_domain::int as nres_domain, num_segments::int as num_segments,
               plddt::float as plddt, md5_domain,
               nullif(cath_label, '-') as cath_label, nullif(cath_level, '-') as cath_level,
               nullif(cath_method, '-') as cath_method,
               try_cast(packing_density as float) as packing_density, try_cast(norm_rg as float) as norm_rg,
               split_part(ted_id, '-', 3) as af_fragment
        from t join acc using (accession)""").df()


def _post(T, a):
    T["globular"] = (T.packing_density >= 10.333) & (T.norm_rg < 0.356)
    if os.path.exists(a.sequences):
        S = pd.read_parquet(a.sequences)
        if os.path.exists(a.proteins):
            P = pd.read_parquet(a.proteins, columns=["accession", "uniprot_acc"]).dropna()
            S = P.merge(S, on="uniprot_acc")[["accession", "sequence"]]
        else:
            S = S.rename(columns={"uniprot_acc": "accession"})
        T = T.merge(S, on="accession", how="left")
        T["seq_matches_uniprot"] = [domain_md5(s, c) == m if isinstance(s, str) else None
                                    for s, c, m in zip(T.sequence, T.chopping, T.md5_domain)]
        T = T.drop(columns=["sequence"])
    return T.drop(columns=["md5_domain"])


def _finish(T, a, acc, got_md5, raw):
    T.to_parquet(os.path.join(a.outdir, "ted_domains.parquet"), index=False, compression="zstd")
    meta = {"source": RECORD, "file": FNAME, "zenodo_version": "v5 (2024-10-31)", "license": "CC BY 4.0",
            "md5": got_md5, "sha256": filehash(raw, "sha256"),
            "n_domains": int(len(T)), "n_accessions_with_ted": int(T.accession.nunique()),
            "n_accessions_total": int(len(acc)),
            "n_domains_cath_H": int((T.cath_level == "H").sum()),
            "frac_domains_seq_match": float(pd.Series(T.get("seq_matches_uniprot")).astype(float).mean())
            if "seq_matches_uniprot" in T else None}
    json.dump(meta, open(os.path.join(a.outdir, "ted_meta.json"), "w"), indent=1)
    print(meta)


if __name__ == "__main__":
    main()
