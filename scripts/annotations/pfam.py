#!/usr/bin/env python
"""Pfam clans (Pfam-A.clans.tsv) and Pfam-A domain regions for atlas accessions.

Raw (InterPro FTP, current_release; to /tmp/pe_big/raw/pfam/):
  Pfam-A.clans.tsv.gz    pfamA_acc, clan_acc, clan_id, pfamA_id, description
  Pfam-A.regions.tsv.gz  pfamseq_acc, seq_version, crc64, md5, pfamA_acc, seq_start, seq_end, ali_start, ali_end (~5 GB)
  Pfam.version.gz, md5_checksums (verified)

Outputs (data/annotations/):
  pfam_clans.parquet     pfam -> clan (all Pfam-A families; clan null if unclustered)
  pfam_regions.parquet   accession, pfam, clan, seq_start, seq_end, ali_start, ali_end, crc64,
                         seq_matches_uniprot (Pfam's crc64 == current UniProt crc64, if uniprot_proteins exists)
  pfam_meta.json

Usage: .venv/bin/python -I scripts/annotations/pfam.py --accessions data/annotations/accessions.parquet
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import subprocess

import duckdb
import pandas as pd

FTP = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/"
RAW = "/tmp/pe_big/raw/pfam"
FILES = ["Pfam-A.clans.tsv.gz", "Pfam-A.regions.tsv.gz", "Pfam.version.gz", "md5_checksums"]


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--accessions", default="data/annotations/accessions.parquet")
    ap.add_argument("--proteins", default="data/annotations/uniprot_proteins.parquet")
    ap.add_argument("--outdir", default="data/annotations")
    ap.add_argument("--threads", type=int, default=16)
    a = ap.parse_args()
    os.makedirs(RAW, exist_ok=True)
    for f in FILES:
        p = os.path.join(RAW, f)
        if not os.path.exists(p):
            subprocess.run(["curl", "-sSf", "-o", p, FTP + f], check=True)
    sums = dict(reversed(l.split()) for l in open(os.path.join(RAW, "md5_checksums")) if l.strip())
    for f in ("Pfam-A.clans.tsv.gz", "Pfam-A.regions.tsv.gz"):
        got = md5(os.path.join(RAW, f))
        assert got == sums[f], f"md5 mismatch {f}: {got} != {sums[f]}"
    version = gzip.open(os.path.join(RAW, "Pfam.version.gz"), "rt").read()
    release = next(l.split(":")[1].strip() for l in version.splitlines() if l.startswith("Pfam release"))

    clans = pd.read_csv(os.path.join(RAW, "Pfam-A.clans.tsv.gz"), sep="\t", header=None,
                        names=["pfam", "clan", "clan_id", "pfam_id", "description"], dtype=str)
    clans["pfam_release"] = release
    clans.to_parquet(os.path.join(a.outdir, "pfam_clans.parquet"), index=False, compression="zstd")

    con = duckdb.connect()
    con.execute(f"set threads={a.threads}")
    acc = pd.read_parquet(a.accessions, columns=["accession"]).drop_duplicates()
    con.register("acc", acc)
    con.register("clans", clans[["pfam", "clan"]])
    R = con.execute(f"""
        select r.pfamseq_acc as accession, r.pfamA_acc as pfam, c.clan, r.seq_version,
               r.seq_start::int as seq_start, r.seq_end::int as seq_end,
               r.ali_start::int as ali_start, r.ali_end::int as ali_end, r.crc64
        from read_csv('{os.path.join(RAW, "Pfam-A.regions.tsv.gz")}', delim='\t', header=true,
                      all_varchar=true) r
        join acc on acc.accession = r.pfamseq_acc
        left join clans c on c.pfam = r.pfamA_acc
        order by accession, seq_start""").df()
    if os.path.exists(a.proteins):
        P = pd.read_parquet(a.proteins, columns=["accession", "seq_crc64"])
        R = R.merge(P, on="accession", how="left")
        R["seq_matches_uniprot"] = R.crc64.str.upper() == R.seq_crc64.str.upper()
        R = R.drop(columns=["seq_crc64"])
    R["pfam_release"] = release
    R.to_parquet(os.path.join(a.outdir, "pfam_regions.parquet"), index=False, compression="zstd")
    meta = {"pfam_release": release, "version_file": version,
            "sha256": {f: sha256(os.path.join(RAW, f)) for f in FILES[:2]},
            "n_families": int(len(clans)), "n_families_in_clan": int(clans.clan.notna().sum()),
            "n_regions": int(len(R)), "n_accessions_with_pfam": int(R.accession.nunique()),
            "n_accessions_total": int(len(acc)),
            "frac_regions_seq_match": float(R.seq_matches_uniprot.mean()) if "seq_matches_uniprot" in R else None}
    json.dump(meta, open(os.path.join(a.outdir, "pfam_meta.json"), "w"), indent=1)
    print({k: v for k, v in meta.items() if k != "version_file"})


if __name__ == "__main__":
    main()
