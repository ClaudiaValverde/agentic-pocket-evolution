#!/usr/bin/env python
"""Representative-pocket fingerprints for a subset of tier-3 sites (O2 specificity; unblocks D-spec-O2).

A site's fingerprint is, by definition (scripts/build_sites.py, `site_fp`), the atlas row of its rep_pocket_id:
the 384-d pocketcompare vector stored in ATLAS/pockets/tier<T>/<unit>.parquet. Tier-3 sites were built with
--no-fingerprints, so this script reads those atlas rows for the requested sites only. Nothing is recomputed:
pocketcompare version / method record are those of the atlas file (copied into the output schema metadata, and
checked against the snapshot method_config_sha256). Atlas files are checked against the snapshot (size, mtime).

Selection (default): every P/W member site of the O2_r1 SOG states (rescue overlay, tau 0.2) on a tier-3 model.

    python scripts/partB/extract_site_fp.py --tier 3 --workers 16 \
        --out /tmp/pe_big/sites/<snap>/fingerprints_O2
Output: <out>/tier<T>/<unit:06d>.parquet (site_id, rep_pocket_id, fingerprint halffloat[384]), same layout as
site_fp/, plus <out>/tier<T>.manifest.json.
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys
import time
from multiprocessing import get_context
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.dataset as ds
import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from pocketevo import recon  # noqa: E402


def one_unit(args):
    uid, man_row, atlas, sites_dir, need, out_dir = args
    f_site = Path(sites_dir) / f"{uid:06d}.parquet"
    st = pq.read_table(f_site, columns=["site_id", "rep_pocket_id"],
                       filters=pc.field("site_id").isin(need))
    apath = Path(atlas) / man_row["pockets_path"]
    s = os.stat(apath)
    if s.st_size != man_row["pockets_size"] or abs(s.st_mtime - man_row["pockets_mtime"]) > 1.0:
        raise RuntimeError(f"{apath} changed since the snapshot")
    meta = pq.read_schema(apath).metadata or {}
    method = meta.get(b"pocketcompare", b"").decode()
    fp = pq.read_table(apath, columns=["pocket_id", "fingerprint"],
                       filters=pc.field("pocket_id").isin(st["rep_pocket_id"].to_pylist()))
    pos = {k: i for i, k in enumerate(fp["pocket_id"].to_pylist())}
    idx = pa.array([pos.get(k) for k in st["rep_pocket_id"].to_pylist()], pa.int64())
    t = pa.table({"site_id": st["site_id"].cast(pa.string()), "rep_pocket_id": st["rep_pocket_id"].cast(pa.string()),
                  "fingerprint": fp["fingerprint"].take(idx)})
    t = t.sort_by("site_id")
    n_missing = t.num_rows - int(pc.sum(pc.is_valid(t["fingerprint"])).as_py() or 0)
    t = t.replace_schema_metadata({"pocketcompare": method, "pocketevo_site_fp": json.dumps(
        dict(unit_id=uid, source=str(apath), selection="O2_r1 P/W tier-3 members", script="extract_site_fp.py"))})
    f = Path(out_dir) / f"{uid:06d}.parquet"
    tmp = f.with_suffix(f".tmp{os.getpid()}")
    pq.write_table(t, tmp)
    os.replace(tmp, f)
    return uid, t.num_rows, n_missing, json.loads(method).get("config_sha256") if method else None, \
        json.loads(method).get("pocketcompare_version") if method else None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--snapshot", default=None)
    ap.add_argument("--tier", type=int, default=3)
    ap.add_argument("--states-glob", default=None, help="default /tmp/pe_big/partB/<snap>/O2_r1/states_resc_t0.2")
    ap.add_argument("--out", default=None)
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--limit-units", type=int, default=0, help="timing: first N units only")
    a = ap.parse_args()
    snap = a.snapshot or recon.default_snapshot()
    t0 = time.time()
    sg = a.states_glob or f"/tmp/pe_big/partB/{snap}/O2_r1/states_resc_t0.2/shard=*.parquet"
    st = ds.dataset(sorted(glob.glob(sg)), format="parquet").to_table(
        columns=["site_id"], filter=pc.field("state").isin(["P", "W"]) & (pc.field("model_tier") == a.tier)
        & pc.field("site_id").is_valid())
    need = pc.unique(st["site_id"].cast(pa.string()))
    man = pd.read_parquet(ROOT / f"data/snapshot/{snap}/manifest.parquet")
    mjs = json.loads((ROOT / f"data/snapshot/{snap}/manifest.json").read_text())
    atlas = mjs["atlas_root"]
    man = man[man["tier"] == a.tier].set_index("unit_id")
    sites_dir = f"/tmp/pe_big/sites/{snap}/sites/tier{a.tier}"
    # site -> unit via the site tables (one pass over the small metadata column)
    sd = ds.dataset(sites_dir, format="parquet")
    su = sd.to_table(columns=["site_id", "unit_id"], filter=pc.field("site_id").isin(need))
    by_unit = su.group_by("unit_id").aggregate([("site_id", "list")]).to_pydict()
    units = sorted(zip(by_unit["unit_id"], by_unit["site_id_list"]))
    if a.limit_units:
        units = units[:a.limit_units]
    out = Path(a.out or f"/tmp/pe_big/sites/{snap}/fingerprints_O2") / f"tier{a.tier}"
    out.mkdir(parents=True, exist_ok=True)
    print(f"[{time.time() - t0:.0f}s] {len(need)} needed sites, {su.num_rows} found in site tables, "
          f"{len(units)} units", flush=True)
    jobs = [(int(u), man.loc[int(u)].to_dict(), atlas, sites_dir, pa.array(s, pa.string()), str(out))
            for u, s in units]
    res = []
    with get_context("fork").Pool(a.workers) as pool:
        for i, r in enumerate(pool.imap_unordered(one_unit, jobs, chunksize=1)):
            res.append(r)
            if (i + 1) % 250 == 0:
                print(f"[{time.time() - t0:.0f}s] {i + 1}/{len(jobs)} units, {sum(x[1] for x in res)} sites",
                      flush=True)
    shas = sorted({x[3] for x in res}); vers = sorted({x[4] for x in res})
    rec = dict(snapshot=snap, tier=a.tier, states_glob=sg, n_needed=len(need), n_in_site_tables=su.num_rows,
               n_units=len(res), n_written=sum(x[1] for x in res), n_missing_fp=sum(x[2] for x in res),
               config_sha256=shas, pocketcompare_version=vers,
               snapshot_method_config_sha256=mjs.get("method_config_sha256"),
               wall_s=time.time() - t0, limit_units=a.limit_units)
    rec["method_matches_snapshot"] = shas == sorted(mjs.get("method_config_sha256") or [])
    (out.parent / f"tier{a.tier}.manifest.json").write_text(json.dumps(rec, indent=1))
    print(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
