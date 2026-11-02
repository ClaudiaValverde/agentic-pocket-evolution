#!/usr/bin/env python
"""Freeze a manifest of finished atlas units.

A unit is included only if
  * state/done/<unit_id>.json exists and parses, and
  * pockets/tier<T>/<unit_id:06d>.parquet and structures/tier<T>/<unit_id:06d>.parquet
    both exist and their Parquet footers are readable.
Only Parquet *metadata* is read (row counts, schema, key-value metadata): no data pages,
so this scales to the full ~80k-unit / ~0.8 B-pocket atlas.

Outputs (in <out>/<UTC timestamp>/):
  manifest.parquet   one row per included unit (paths, sizes, mtimes, row counts, json counts, method hash)
  excluded.parquet   units seen on disk but excluded, with the reason
  manifest.json      summary + sha256 of manifest.parquet

Downstream steps must read ONLY the units listed in a manifest. Re-running on the complete atlas is the
same command; a new snapshot id is created each time.

Usage:
  uv run --with duckdb --with pyarrow --with pandas python scripts/atlas_snapshot.py \
      --atlas-root /cluster/scratch/gmagnusson/pocketcompare_afdb --out data/snapshot
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
import re
import sys
from concurrent.futures import ThreadPoolExecutor

import pyarrow as pa
import pyarrow.parquet as pq

UNIT_RE = re.compile(r"^(\d+)\.json$")
PQ_RE = re.compile(r"^(\d{6,})\.parquet$")


def sha256_file(path: str, bufsize: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while chunk := fh.read(bufsize):
            h.update(chunk)
    return h.hexdigest()


def pq_info(path: str) -> dict:
    """Footer-only info for one parquet file; raises if unreadable."""
    st = os.stat(path)
    pf = pq.ParquetFile(path)
    md = pf.metadata
    schema = pf.schema_arrow
    kv = schema.metadata or {}
    method_sha = None
    if b"pocketcompare" in kv:
        try:
            method_sha = json.loads(kv[b"pocketcompare"]).get("config_sha256")
        except Exception:
            method_sha = None
    return dict(
        size=st.st_size,
        mtime=st.st_mtime,
        rows=md.num_rows,
        row_groups=md.num_row_groups,
        has_source=("source" in schema.names),
        method_sha=method_sha,
    )


def scan_unit(atlas: str, unit_id: int, done_path: str) -> dict:
    rec: dict = {"unit_id": unit_id}
    try:
        st = os.stat(done_path)
        with open(done_path) as fh:
            dj = json.load(fh)
    except Exception as e:  # noqa: BLE001
        rec["excluded_reason"] = f"done_json_unreadable: {e!r}"[:300]
        return rec
    tier = int(dj.get("tier", -1))
    rec.update(
        tier=tier,
        done_json=os.path.relpath(done_path, atlas),
        done_json_size=st.st_size,
        done_json_mtime=st.st_mtime,
        done_time=str(dj.get("time", "")),
        json_structures=int(dj.get("structures", -1)),
        json_structures_ok=int(dj.get("structures_ok", -1)),
        json_structures_error=int(dj.get("structures_error", -1)),
        json_pockets_all=int(dj.get("pockets_all", -1)),
        json_pockets_kept=int(dj.get("pockets_kept", -1)),
    )
    for kind in ("pockets", "structures"):
        rel = f"{kind}/tier{tier}/{unit_id:06d}.parquet"
        path = os.path.join(atlas, rel)
        if not os.path.exists(path):
            rec["excluded_reason"] = f"missing_{kind}_parquet"
            return rec
        try:
            info = pq_info(path)
        except Exception as e:  # noqa: BLE001
            rec["excluded_reason"] = f"unreadable_{kind}_parquet: {e!r}"[:300]
            return rec
        rec[f"{kind}_path"] = rel
        rec[f"{kind}_size"] = info["size"]
        rec[f"{kind}_mtime"] = info["mtime"]
        rec[f"{kind}_rows"] = info["rows"]
        rec[f"{kind}_row_groups"] = info["row_groups"]
        rec[f"{kind}_has_source"] = info["has_source"]
        if kind == "pockets":
            rec["method_config_sha256"] = info["method_sha"]
    rec["rows_match_json"] = bool(
        rec["pockets_rows"] == rec["json_pockets_kept"] and rec["structures_rows"] == rec["json_structures"]
    )
    rec["excluded_reason"] = None
    return rec


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--atlas-root", required=True)
    ap.add_argument("--out", required=True, help="parent dir; snapshot goes in <out>/<UTC timestamp>/")
    ap.add_argument("--threads", type=int, default=32)
    ap.add_argument("--snapshot-id", default=None, help="override the UTC-timestamp id (for tests)")
    args = ap.parse_args(argv)

    atlas = os.path.abspath(args.atlas_root)
    t0 = dt.datetime.now(dt.timezone.utc)
    snap_id = args.snapshot_id or t0.strftime("%Y%m%dT%H%M%SZ")
    out_dir = os.path.join(os.path.abspath(args.out), snap_id)
    os.makedirs(out_dir, exist_ok=False)

    # 1. list done jsons (directory listing first = the freeze point)
    done_dir = os.path.join(atlas, "state", "done")
    done = {}
    with os.scandir(done_dir) as it:
        for e in it:
            m = UNIT_RE.match(e.name)
            if m:
                done[int(m.group(1))] = e.path
    # parquet files present on disk without a done json (in progress / orphaned): report only
    orphan = []
    for kind in ("pockets", "structures"):
        kdir = os.path.join(atlas, kind)
        if not os.path.isdir(kdir):
            continue
        for tdir in sorted(os.listdir(kdir)):
            tpath = os.path.join(kdir, tdir)
            if not (tdir.startswith("tier") and os.path.isdir(tpath)):
                continue
            with os.scandir(tpath) as it:
                for e in it:
                    m = PQ_RE.match(e.name)
                    if m and int(m.group(1)) not in done:
                        orphan.append(dict(unit_id=int(m.group(1)), tier=int(tdir[4:]),
                                           excluded_reason=f"no_done_json ({kind} parquet present)"))

    # 2. footers in parallel
    with ThreadPoolExecutor(args.threads) as ex:
        recs = list(ex.map(lambda kv: scan_unit(atlas, kv[0], kv[1]), sorted(done.items())))
    inc = [r for r in recs if r.get("excluded_reason") is None]
    exc = [r for r in recs if r.get("excluded_reason") is not None] + orphan

    cols = ["unit_id", "tier", "done_json", "done_json_size", "done_json_mtime", "done_time",
            "json_structures", "json_structures_ok", "json_structures_error", "json_pockets_all",
            "json_pockets_kept",
            "pockets_path", "pockets_size", "pockets_mtime", "pockets_rows", "pockets_row_groups",
            "pockets_has_source",
            "structures_path", "structures_size", "structures_mtime", "structures_rows",
            "structures_row_groups", "structures_has_source", "method_config_sha256", "rows_match_json"]
    tbl = pa.Table.from_pylist([{c: r.get(c) for c in cols} for r in inc]) if inc else pa.table({c: [] for c in cols})
    tbl = tbl.sort_by("unit_id")
    tbl = tbl.replace_schema_metadata({
        "atlas_snapshot": json.dumps({"atlas_root": atlas, "snapshot_id": snap_id,
                                      "created_utc": t0.isoformat()})})
    man_path = os.path.join(out_dir, "manifest.parquet")
    pq.write_table(tbl, man_path)
    exc_tbl = pa.Table.from_pylist([{"unit_id": r["unit_id"], "tier": r.get("tier"),
                                     "excluded_reason": r["excluded_reason"]} for r in exc],
                                   schema=pa.schema([("unit_id", pa.int64()), ("tier", pa.int64()),
                                                     ("excluded_reason", pa.string())]))
    pq.write_table(exc_tbl, os.path.join(out_dir, "excluded.parquet"))

    # 3. summary
    per_tier = {}
    d = tbl.to_pydict()
    for i, t in enumerate(d["tier"]):
        s = per_tier.setdefault(str(t), dict(units=0, structures=0, pockets=0, pockets_bytes=0,
                                              structures_bytes=0, rows_mismatch=0, missing_source_col=0))
        s["units"] += 1
        s["structures"] += d["structures_rows"][i]
        s["pockets"] += d["pockets_rows"][i]
        s["pockets_bytes"] += d["pockets_size"][i]
        s["structures_bytes"] += d["structures_size"][i]
        s["rows_mismatch"] += int(not d["rows_match_json"][i])
        s["missing_source_col"] += int(not (d["pockets_has_source"][i] and d["structures_has_source"][i]))
    worklist_units = None
    wl = os.path.join(atlas, "worklist.parquet")
    if os.path.exists(wl):
        wt = pq.read_table(wl, columns=["tier"]).column("tier").to_pylist()
        worklist_units = {str(t): wt.count(t) for t in sorted(set(wt))}
        for t, n in worklist_units.items():
            per_tier.setdefault(t, dict(units=0))["worklist_units"] = n
    excl_reasons: dict = {}
    for r in exc:
        key = r["excluded_reason"].split(":")[0]
        excl_reasons[key] = excl_reasons.get(key, 0) + 1
    summary = dict(
        snapshot_id=snap_id,
        created_utc=t0.isoformat(),
        finished_utc=dt.datetime.now(dt.timezone.utc).isoformat(),
        atlas_root=atlas,
        manifest="manifest.parquet",
        manifest_sha256=sha256_file(man_path),
        worklist_sha256=sha256_file(wl) if os.path.exists(wl) else None,
        n_units=tbl.num_rows,
        n_excluded=len(exc),
        excluded_reasons=excl_reasons,
        method_config_sha256=sorted({x for x in d["method_config_sha256"] if x}),
        per_tier=per_tier,
        script=os.path.abspath(__file__),
        argv=sys.argv,
    )
    with open(os.path.join(out_dir, "manifest.json"), "w") as fh:
        json.dump(summary, fh, indent=1)
    print(json.dumps(summary, indent=1))
    print(out_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
