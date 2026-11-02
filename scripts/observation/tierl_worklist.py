#!/usr/bin/env python
"""Build the tier-L (local pockets) work list for the atlas worker (external/pocketcompare/afdb/worker.py).

Input: the context-species list written by the orthology agent (rows with pockets_source == 'local') and the
proteome FASTA files listed in o2_proteomes.tsv. For every proteome accession present in the atlas Foldcomp
entry list (fc_entries_afdb_uniprot_v4.parquet: all AFDB v4 entries minus tier-1/2 accessions) we write

  <root>/fc_entries_tierL.parquet   same schema as the atlas entry list (fc_unit, off, len, key, accession,
                                    fragment, version, taxid), sorted by fc_unit then offset, row groups 100k
  <root>/worklist.parquet           worker-compatible units (source_type foldcomp, source afdb_uniprot_v4,
                                    sources [fc_entries_tierL.parquet], tier = TIER_CODE), <= 4,000 entries each
  <root>/skip_accessions.txt        empty (entries already exclude tier-1/2 accessions)
  <root>/fc_taxa_afdb_uniprot_v4.parquet -> symlink to the atlas file (organism names)
  <root>/tierL_accounting.parquet   per proteome: n_fasta, n_in_fc_entries, n_in_tier12 (skip list), n_missing

Tier code: the atlas and src/pocketevo/sites.py store `tier` as int8, so tier L is stored as TIER_CODE = 4 in the
tables; manifests add tier_label = 'L'. Unit ids start at UNIT0 = 900000 (no collision with atlas units).

  python scripts/observation/tierl_worklist.py --root /tmp/pe_big/tierL
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import time
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

P = Path(__file__).resolve().parents[2]
ATLAS = Path("/cluster/scratch/gmagnusson/pocketcompare_afdb")
TIER_CODE = 4
UNIT0 = 900000
UNIT_MAX = 4000


def fasta_accessions(path: str) -> list[str]:
    out = []
    with gzip.open(path, "rt") as fh:
        for ln in fh:
            if ln.startswith(">"):
                out.append(ln[1:].split("|")[1] if "|" in ln else ln[1:].split()[0])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default="/tmp/pe_big/tierL")
    ap.add_argument("--context", default=str(P / "data/orthology/ancestra_O2/context_species.tsv"))
    ap.add_argument("--proteomes", default=str(P / "data/orthology/ancestra_O2/o2_proteomes.tsv"))
    ap.add_argument("--threads", type=int, default=16)
    a = ap.parse_args()
    t0 = time.time()
    root = Path(a.root)
    root.mkdir(parents=True, exist_ok=True)
    ctx = pd.read_csv(a.context, sep="\t")
    loc = ctx[ctx.pockets_source == "local"].copy()
    prot = pd.read_csv(a.proteomes, sep="\t")
    prot = prot[prot.proteome.isin(loc.proteome)]
    rows = []
    for r in prot.itertuples():
        for acc in fasta_accessions(r.file):
            rows.append((acc, r.proteome, int(r.taxid)))
    acc = pd.DataFrame(rows, columns=["accession", "proteome", "proteome_taxid"]).drop_duplicates("accession")
    skip = pd.read_csv(ATLAS / "skip_accessions.txt", sep="\t", header=None, names=["accession", "first_tier"],
                       dtype={"accession": str, "first_tier": np.int8}, na_filter=False)
    con = duckdb.connect()
    con.execute(f"SET threads={a.threads}")
    con.register("acc", acc)
    e = con.execute(f"""
        SELECT f.off, f.len, f.key, f.accession, f.fragment, f.version, f.taxid, acc.proteome, acc.proteome_taxid
        FROM read_parquet('{ATLAS}/fc_entries_afdb_uniprot_v4.parquet') f JOIN acc USING (accession)
    """).df()
    # one row per accession+fragment; keep fragment 1..n (AFDB v4 has fragments only for human, none here)
    e = e.sort_values(["proteome", "accession", "fragment"]).reset_index(drop=True)
    acc = acc.merge(skip, on="accession", how="left")
    found = set(e.accession)
    acc["in_fc_entries"] = acc.accession.isin(found)
    acc["in_tier12"] = acc.first_tier.isin([1, 2])
    accounting = acc.groupby("proteome").agg(n_fasta=("accession", "size"), n_in_fc_entries=("in_fc_entries", "sum"),
                                             n_in_tier12=("in_tier12", "sum")).reset_index()
    accounting["n_missing"] = accounting.n_fasta - accounting.n_in_fc_entries - accounting.n_in_tier12
    accounting = accounting.merge(loc[["proteome", "taxid", "species"]], on="proteome")
    pq.write_table(pa.Table.from_pandas(accounting, preserve_index=False), root / "tierL_accounting.parquet")
    acc.to_parquet(root / "tierL_accessions.parquet", index=False)
    # units: per proteome, chunks of <= UNIT_MAX entries, in context-list order
    order = {p: i for i, p in enumerate(loc.proteome)}
    e["porder"] = e.proteome.map(order)
    e = e.sort_values(["porder", "accession", "fragment"]).reset_index(drop=True)
    e["chunk"] = e.groupby("proteome").cumcount() // UNIT_MAX
    keys = e[["porder", "chunk"]].drop_duplicates().reset_index(drop=True)
    keys["fc_unit"] = np.arange(len(keys), dtype=np.int32)
    e = e.merge(keys, on=["porder", "chunk"])
    e = e.sort_values(["fc_unit", "off"]).reset_index(drop=True)
    ent = pa.table({"fc_unit": pa.array(e.fc_unit.to_numpy(np.int32)), "off": pa.array(e.off.to_numpy(np.uint64)),
                    "len": pa.array(e.len.to_numpy(np.uint32)), "key": pa.array(e.key.to_numpy(np.uint64)),
                    "accession": pa.array(e.accession.astype(str)), "fragment": pa.array(e.fragment.to_numpy(np.int16)),
                    "version": pa.array(e.version.to_numpy(np.int8)), "taxid": pa.array(e.taxid.to_numpy(np.int32))})
    pq.write_table(ent, root / "fc_entries_tierL.parquet", row_group_size=100_000, compression="zstd",
                   write_statistics=True)
    taxa = pd.read_parquet(ATLAS / "fc_taxa_afdb_uniprot_v4.parquet", columns=["taxid", "tax_name"])
    tn = dict(zip(taxa.taxid, taxa.tax_name))
    wl = []
    for k in keys.itertuples():
        sub = e[e.fc_unit == k.fc_unit]
        lr = loc.iloc[int(k.porder)]
        tax = int(lr.taxid)
        wl.append(dict(unit_id=UNIT0 + int(k.fc_unit), tier=TIER_CODE, source_type="foldcomp", source="afdb_uniprot_v4",
                       sources=["fc_entries_tierL.parquet"], member_start=int(k.fc_unit), member_end=int(k.fc_unit) + 1,
                       byte_start=0, byte_end=0, taxid=tax, code=f"L_{lr.oscode}_{int(k.chunk)}",
                       organism=tn.get(tax, str(lr.species)), lineage=str(lr.lineage), model_version=4,
                       est_structures=len(sub), source_taxids=[tax], source_sizes=[len(sub)],
                       proteome=str(lr.proteome)))
    wl = pd.DataFrame(wl)
    wl.to_parquet(root / "worklist.parquet", index=False)
    (root / "skip_accessions.txt").write_text("")
    lk = root / "fc_taxa_afdb_uniprot_v4.parquet"
    if not lk.exists():
        os.symlink(ATLAS / "fc_taxa_afdb_uniprot_v4.parquet", lk)
    # taxid sanity: Foldcomp taxid vs proteome taxid
    tax_agree = float((e.taxid == e.proteome_taxid).mean())
    summ = dict(n_proteomes=int(len(loc)), n_fasta_accessions=int(len(acc)), n_entries=int(len(e)),
                n_units=int(len(wl)), n_in_tier12=int(acc.in_tier12.sum()),
                n_missing=int(len(acc) - acc.in_fc_entries.sum() - acc.in_tier12.sum()),
                fc_taxid_equals_proteome_taxid=tax_agree, tier_code=TIER_CODE, unit0=UNIT0,
                context_sha=__import__("hashlib").sha256(Path(a.context).read_bytes()).hexdigest(),
                runtime_s=round(time.time() - t0, 1))
    (root / "tierL_worklist.json").write_text(json.dumps(summ, indent=1))
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()
