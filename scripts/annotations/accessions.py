#!/usr/bin/env python
"""Build the accession list (tiers 1 and 2) for an atlas snapshot manifest.

Output: data/annotations/accessions.parquet with columns
  accession, tier, taxid, proteome (UPxxxx from the AFDB tar name, tier 1 only), organism, n_res, fragment
One row per (accession, tier); an accession can be in both tiers (Swiss-Prot entries of reference proteomes).

Usage:
  .venv/bin/python scripts/annotations/accessions.py \
      --manifest data/snapshot/20261030T085136Z/manifest.parquet --out data/annotations/accessions.parquet
"""
from __future__ import annotations

import argparse
import json
import os
import re

import duckdb
import pyarrow.parquet as pq

ATLAS_DEFAULT = "/cluster/scratch/gmagnusson/pocketcompare_afdb"


def build(manifest: str, atlas_root: str | None = None, tiers=(1, 2)):
    md = pq.read_metadata(manifest).metadata or {}
    if atlas_root is None:
        try:
            atlas_root = json.loads(md[b"atlas_snapshot"])["atlas_root"]
        except Exception:
            atlas_root = ATLAS_DEFAULT
    m = pq.read_table(manifest, columns=["tier", "structures_path"]).to_pandas()
    m = m[m.tier.isin(tiers)]
    paths = [os.path.join(atlas_root, p) for p in m.structures_path]
    con = duckdb.connect()
    df = con.execute(
        "select accession, tier::int as tier, taxid, organism, source, n_res, fragment "
        "from read_parquet(?, union_by_name=true)", [paths]).df()
    df["proteome"] = df.source.map(lambda s: (re.match(r"(UP\d+)_", s or "") or [None, None])[1])
    df = (df.sort_values(["accession", "tier", "fragment"])
            .drop_duplicates(["accession", "tier"])
            .drop(columns=["source"]).reset_index(drop=True))
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--atlas-root", default=None)
    ap.add_argument("--out", default="data/annotations/accessions.parquet")
    a = ap.parse_args()
    df = build(a.manifest, a.atlas_root)
    df.to_parquet(a.out, index=False)
    print(df.groupby("tier").accession.nunique().to_dict(), "unique:", df.accession.nunique(),
          "proteomes:", df.proteome.nunique())


if __name__ == "__main__":
    main()
