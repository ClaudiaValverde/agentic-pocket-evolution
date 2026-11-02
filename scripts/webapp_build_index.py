#!/usr/bin/env python
"""Build the web app's accession index from an atlas snapshot manifest.

One row per atlas structure (accession x fragment x tier) with everything the web app needs to
(a) find the unit pocket file of a protein and (b) read its model bytes directly from the source
(AFDB v6 tar member or Foldcomp DB entry) without copying AFDB.

Input : data/snapshot/<snapshot>/manifest.parquet (frozen list of atlas units)
        ATLAS/structures/tier<T>/<unit>.parquet       (accession, fragment, n_res, pLDDT, counts)
        ATLAS/tar_index/<tar>.parquet                 (tier 1/2: member name -> offset_data, size)
        ATLAS/fc_entries_afdb_uniprot_v4.parquet      (tier 3: fc_unit, off, len per entry)
Output: <out>/tier<T>/bucket=<BB>.parquet, BB = f"{crc32(accession) & 0xff:02x}" (256 shards per tier,
        each sorted by accession; a lookup reads one shard per tier). Plus <out>/index.json (provenance).
        Hash-prefix sharding keeps shards balanced (UniProt accession prefixes are very skewed:
        most tier-3 accessions start with "A0A"), so the layout scales to all 214 M AFDB v4 entries
        (~0.8 M rows per shard).

Default output: tiers 1-2 -> results/webapp_index/<snapshot>/ (small, in the repo);
                tier 3   -> /tmp/pe_big/webapp_index/<snapshot>/ (large, regenerable; node-local).

Usage (repo root):
  .venv/bin/python -I scripts/webapp_build_index.py --snapshot 20261030T085136Z [--tiers 1,2,3] [--threads 32]
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import datetime as dt
import hashlib
import json
import os
import sys
import time
import zlib
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

REPO = Path(__file__).resolve().parents[1]
N_BUCKETS = 256

SCHEMA = pa.schema([
    ("accession", pa.string()),
    ("fragment", pa.int16()),
    ("tier", pa.int8()),
    ("unit_id", pa.int32()),
    ("model_version", pa.int8()),
    ("taxid", pa.int32()),
    ("organism", pa.string()),
    ("n_res", pa.int32()),
    ("plddt_mean", pa.float32()),
    ("n_pockets_kept", pa.int32()),
    ("status", pa.string()),
    ("blob_kind", pa.string()),     # "tar_cif_gz" | "foldcomp"
    ("blob_file", pa.string()),     # tar file name (under AFDB v6 dir) or Foldcomp DB name (under ATLAS/foldcomp)
    ("blob_offset", pa.int64()),    # byte offset of the member data / Foldcomp entry
    ("blob_size", pa.int64()),      # byte length
])

STRUCT_COLS = ["accession", "fragment", "model_version", "taxid", "organism", "n_res", "plddt_mean",
               "n_pockets_kept", "status"]


def bucket_of(acc: str) -> int:
    return zlib.crc32(acc.encode()) & (N_BUCKETS - 1)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def read_structures(atlas: Path, row) -> pd.DataFrame:
    f = atlas / row.structures_path
    sch = pq.read_schema(f)
    cols = STRUCT_COLS + [c for c in ("source", "source_file") if c in sch.names]
    df = pq.read_table(f, columns=cols).to_pandas()
    df["unit_id"] = int(row.unit_id)
    df["tier"] = int(row.tier)
    return df


def attach_tar_blobs(df: pd.DataFrame, atlas: Path) -> pd.DataFrame:
    # tier 1: `source` holds the tar name; tier 2: `source_file` holds it (`source` = afdb_swissprot_v6_tar)
    if "source_file" in df and df["source_file"].notna().all() and (df["source_file"] != "").all():
        tar = df["source_file"]
    else:
        tar = df["source"]
    df = df.assign(blob_file=tar.astype(str), blob_kind="tar_cif_gz")
    df["member"] = ("AF-" + df.accession + "-F" + df.fragment.astype(str) + "-model_v"
                    + df.model_version.astype(str) + ".cif.gz")
    out = []
    for tarname, g in df.groupby("blob_file", sort=False):
        ti = pq.read_table(atlas / "tar_index" / f"{tarname}.parquet").to_pandas()
        ti = ti.rename(columns={"name": "member", "offset_data": "blob_offset", "size": "blob_size"})
        out.append(g.merge(ti, on="member", how="left"))
    df = pd.concat(out, ignore_index=True)
    miss = df.blob_offset.isna().sum()
    if miss:
        print(f"  WARNING: {miss} tier-1/2 structures without a tar member", file=sys.stderr)
    return df.drop(columns=["member"])


def attach_fc_blobs(df: pd.DataFrame, atlas: Path, unit_ids: np.ndarray, fc_unit_base: int) -> pd.DataFrame:
    """fc_unit = unit_id - fc_unit_base (worklist member_start of tier-3 Foldcomp units)."""
    fcp = pq.ParquetFile(atlas / "fc_entries_afdb_uniprot_v4.parquet")
    want = set((unit_ids - fc_unit_base).tolist())
    rgs = []
    md = fcp.metadata
    ci = fcp.schema_arrow.get_field_index("fc_unit")
    for i in range(md.num_row_groups):
        st = md.row_group(i).column(ci).statistics
        lo, hi = st.min, st.max
        if any(lo <= u <= hi for u in want):
            rgs.append(i)
    t = fcp.read_row_groups(rgs, columns=["fc_unit", "off", "len", "accession", "fragment"])
    t = t.filter(pc.is_in(t["fc_unit"], pa.array(sorted(want), pa.int32())))
    fc = t.to_pandas()
    fc["unit_id"] = fc.fc_unit.astype(np.int64) + fc_unit_base
    fc = fc.rename(columns={"off": "blob_offset", "len": "blob_size"}).drop(columns=["fc_unit"])
    fc = fc.drop_duplicates(["unit_id", "accession", "fragment"])
    df = df.merge(fc, on=["unit_id", "accession", "fragment"], how="left")
    miss = df.blob_offset.isna().sum()
    if miss:
        print(f"  WARNING: {miss} tier-3 structures without a Foldcomp entry", file=sys.stderr)
    df["blob_kind"] = "foldcomp"
    df["blob_file"] = "afdb_uniprot_v4"
    return df


def write_buckets(df: pd.DataFrame, outdir: Path) -> dict:
    outdir.mkdir(parents=True, exist_ok=True)
    df = df.copy()
    df["_b"] = [bucket_of(a) for a in df.accession]
    sizes = {}
    for b, g in df.groupby("_b"):
        g = g.drop(columns=["_b"]).sort_values(["accession", "tier", "fragment"])
        t = pa.Table.from_pandas(g[[f.name for f in SCHEMA]], schema=SCHEMA, preserve_index=False)
        p = outdir / f"bucket={b:02x}.parquet"
        pq.write_table(t, p, compression="zstd", row_group_size=65536,
                       use_dictionary=["organism", "status", "blob_kind", "blob_file"])
        sizes[p.name] = len(g)
    return sizes


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--snapshot", default="20261030T085136Z")
    ap.add_argument("--tiers", default="1,2,3")
    ap.add_argument("--out12", default=None, help="output dir for tiers 1-2 (default results/webapp_index/<snap>)")
    ap.add_argument("--out3", default=None, help="output dir for tier 3 (default /tmp/pe_big/webapp_index/<snap>)")
    ap.add_argument("--threads", type=int, default=min(32, os.cpu_count() or 4))
    ap.add_argument("--limit-units", type=int, default=0, help="debug: only the first N units per tier")
    a = ap.parse_args(argv)

    snapdir = REPO / "data" / "snapshot" / a.snapshot
    man_path = snapdir / "manifest.parquet"
    man = pd.read_parquet(man_path)
    meta = json.loads(pq.read_schema(man_path).metadata[b"atlas_snapshot"])
    atlas = Path(meta["atlas_root"])
    wl = pd.read_parquet(atlas / "worklist.parquet", columns=["unit_id", "tier", "source_type", "member_start"])
    fcw = wl[(wl.tier == 3) & (wl.source_type == "foldcomp")]
    fc_unit_base = int((fcw.unit_id - fcw.member_start).iloc[0]) if len(fcw) else 0
    assert len(fcw) == 0 or ((fcw.unit_id - fcw.member_start) == fc_unit_base).all()

    out12 = Path(a.out12) if a.out12 else REPO / "results" / "webapp_index" / a.snapshot
    out3 = Path(a.out3) if a.out3 else Path("/tmp/pe_big/webapp_index") / a.snapshot
    tiers = [int(x) for x in a.tiers.split(",")]
    record = {"snapshot": a.snapshot, "manifest_sha256": sha256(man_path), "atlas_root": str(atlas),
              "n_buckets": N_BUCKETS, "bucket_fn": "crc32(accession.encode()) & 0xff -> '%02x'",
              "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(), "tiers": {}}
    for tier in tiers:
        t0 = time.time()
        units = man[man.tier == tier]
        if a.limit_units:
            units = units.head(a.limit_units)
        if units.empty:
            continue
        with cf.ThreadPoolExecutor(a.threads) as ex:
            parts = list(ex.map(lambda r: read_structures(atlas, r), units.itertuples()))
        df = pd.concat(parts, ignore_index=True)
        print(f"tier {tier}: {len(units)} units, {len(df):,} structures read in {time.time()-t0:.1f}s", flush=True)
        if tier in (1, 2):
            df = attach_tar_blobs(df, atlas)
        else:
            df = attach_fc_blobs(df, atlas, units.unit_id.to_numpy(), fc_unit_base)
        out = (out12 if tier in (1, 2) else out3) / f"tier{tier}"
        if out.exists():
            for p in out.glob("bucket=*.parquet"):
                p.unlink()
        sizes = write_buckets(df, out)
        nbytes = sum((out / k).stat().st_size for k in sizes)
        record["tiers"][str(tier)] = {"dir": str(out), "units": int(len(units)), "rows": int(len(df)),
                                      "accessions": int(df.accession.nunique()), "bytes": int(nbytes),
                                      "missing_blob": int(df.blob_offset.isna().sum()),
                                      "seconds": round(time.time() - t0, 1)}
        print(f"tier {tier}: wrote {len(sizes)} shards, {nbytes/1e6:.1f} MB to {out} in {time.time()-t0:.1f}s",
              flush=True)
        # per-tier provenance next to the shards
        (out / "index.json").write_text(json.dumps({**{k: v for k, v in record.items() if k != "tiers"},
                                                    "tier": tier, **record["tiers"][str(tier)]}, indent=1))
    out12.mkdir(parents=True, exist_ok=True)
    prev = out12 / "index.json"
    if prev.exists():  # keep provenance of tiers built in earlier runs
        old = json.loads(prev.read_text())
        if old.get("manifest_sha256") == record["manifest_sha256"]:
            record["tiers"] = {**old.get("tiers", {}), **record["tiers"]}
    (out12 / "index.json").write_text(json.dumps(record, indent=1))
    print(json.dumps(record["tiers"], indent=1))


if __name__ == "__main__":
    main()
