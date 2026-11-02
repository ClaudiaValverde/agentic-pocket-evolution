#!/usr/bin/env python
"""Build the per-protein site table (DESIGN.md §3.1) from a frozen atlas snapshot, one output per unit file.

For every snapshot unit of the requested tiers: read its pocket file, convert residues to UniProt positions,
merge the pockets of each accession with residue Jaccard >= --jaccard (union-find) and write
    <out>/sites/tier<T>/<unit:06d>.parquet       site metadata (no fingerprints)
    <out>/site_fp/tier<T>/<unit:06d>.parquet     site_id, rep_pocket_id, rep fingerprint (float16[384])
                                                 (omitted with --no-fingerprints)
    <out>/state/<unit:06d>.json                  input key, counts, timing
Only one unit's fingerprints are in memory per worker, so this streams over the full atlas.

Accessions whose fragments are spread over several unit files (rare; found by a footer/statistics pre-pass,
cached in <out>/split_accessions.parquet) are processed by the lowest unit_id that contains them; that worker
also reads those accessions' rows from the other units, and the other units drop them.

Resumable: a unit is skipped when its state json has the same input key (sha256 over the input files'
path/size/mtime as recorded in the manifest and checked on disk, the split-accession assignment, --jaccard,
the schema version and the sha256 of pocketevo/sites.py) and both outputs exist.

    python scripts/build_sites.py --manifest data/snapshot/20261030T085136Z/manifest.parquet --tiers 1,2 \
        --out /tmp/pe_big/sites/20261030T085136Z --workers 32
    # timing sample: --sample 50 --seed 0  (deterministic random subset of the selected units)
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pyarrow as pa  # noqa: E402
import pyarrow.compute as pc  # noqa: E402
import pyarrow.parquet as pq  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from pocketevo import sites as S  # noqa: E402


def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def load_manifest(path: Path, tiers: list[int]) -> tuple[pd.DataFrame, dict]:
    mj = json.loads((path.parent / "manifest.json").read_text())
    if sha256_bytes(path.read_bytes()) != mj["manifest_sha256"]:
        raise SystemExit(f"manifest checksum mismatch: {path}")
    man = pd.read_parquet(path)
    man = man[man.tier.isin(tiers)].sort_values("unit_id").reset_index(drop=True)
    if man.unit_id.duplicated().any():
        raise SystemExit("duplicate unit_id in manifest")
    return man, mj


# ------------------------------------------------------------------------------------ split accessions
def _frag_max(path: str) -> int:
    md = pq.ParquetFile(path).metadata
    j = md.schema.to_arrow_schema().get_field_index("fragment")
    if md.num_rows == 0 or j < 0:
        return 1
    mx = 1
    for g in range(md.num_row_groups):
        st = md.row_group(g).column(j).statistics
        if st is None or not st.has_min_max:
            return 1 << 30   # unknown -> must read
        mx = max(mx, st.max)
    return mx


def _acc_range(path: str):
    md = pq.ParquetFile(path).metadata
    j = md.schema.to_arrow_schema().get_field_index("accession")
    if md.num_rows == 0 or j < 0:
        return ("~", "~")   # empty unit file: an empty range that excludes every accession
    lo, hi = None, None
    for g in range(md.num_row_groups):
        st = md.row_group(g).column(j).statistics
        if st is None or not st.has_min_max:
            return None
        lo = st.min if lo is None else min(lo, st.min)
        hi = st.max if hi is None else max(hi, st.max)
    return lo, hi


def find_split_accessions(man: pd.DataFrame, atlas: Path, workers: int) -> pd.DataFrame:
    """accession, unit_id rows for multi-fragment accessions present in > 1 unit file (pocket files)."""
    paths = [str(atlas / p) for p in man.pockets_path]
    with ProcessPoolExecutor(workers) as ex:
        fmax = list(ex.map(_frag_max, paths, chunksize=64))
    multi = set()
    for p, m in zip(paths, fmax):
        if m > 1:
            t = pq.read_table(p, columns=["accession", "fragment"])
            multi |= set(pc.filter(t["accession"], pc.greater(t["fragment"], 1)).to_pylist())
    if not multi:
        return pd.DataFrame({"accession": pd.Series([], dtype=str), "unit_id": pd.Series([], dtype=np.int64)})
    lo_all, hi_all = min(multi), max(multi)
    with ProcessPoolExecutor(workers) as ex:
        ranges = list(ex.map(_acc_range, paths, chunksize=64))
    rows = []
    marr = pa.array(sorted(multi))
    for (uid, p), r in zip(zip(man.unit_id, paths), ranges):
        if r is not None and (r[1] < lo_all or r[0] > hi_all):
            continue
        if r == ("~", "~") or (r is not None and not any(r[0] <= a <= r[1] for a in multi)):
            continue
        acc = pq.read_table(p, columns=["accession"])["accession"].combine_chunks().cast(pa.string())
        for a in set(pc.filter(acc, pc.is_in(acc, marr)).to_pylist()):
            rows.append((a, int(uid)))
    df = pd.DataFrame(rows, columns=["accession", "unit_id"])
    n_units = df.groupby("accession").unit_id.transform("nunique")
    return df[n_units > 1].sort_values(["accession", "unit_id"]).reset_index(drop=True)


# ------------------------------------------------------------------------------------ per unit
def unit_plan(u: dict, split: pd.DataFrame) -> dict:
    """Which split accessions this unit drops (owned elsewhere) and which it owns (+ their other units)."""
    mine = split[split.unit_id == u["unit_id"]]
    owner = split.groupby("accession").unit_id.min()
    drop = sorted(a for a in mine.accession if owner[a] != u["unit_id"])
    own = sorted(a for a in mine.accession if owner[a] == u["unit_id"])
    extra = sorted(set(split[split.accession.isin(own) & (split.unit_id != u["unit_id"])].unit_id.astype(int)))
    return dict(drop=drop, own=own, extra_units=extra)


def input_key(u: dict, plan: dict, extra_rows: list[dict], jaccard: float, code_sha: str) -> str:
    files = [(r["pockets_path"], int(r["pockets_size"]), float(r["pockets_mtime"])) for r in [u] + extra_rows]
    rec = dict(schema=S.SCHEMA_VERSION, jaccard=jaccard, step=S.FRAGMENT_STEP, code=code_sha, files=files,
               drop=plan["drop"], own=plan["own"])
    return sha256_bytes(json.dumps(rec, sort_keys=True).encode())


def check_on_disk(atlas: Path, r: dict):
    st = os.stat(atlas / r["pockets_path"])
    if st.st_size != r["pockets_size"] or abs(st.st_mtime - r["pockets_mtime"]) > 1.0:
        raise RuntimeError(f"{r['pockets_path']} changed since the snapshot (size/mtime)")


def read_pockets(path: Path, accessions=None) -> pa.Table:
    schema = pq.read_schema(path)
    if pq.ParquetFile(path).metadata.num_rows == 0:   # empty unit files carry a truncated schema
        return S.POCKET_SCHEMA.empty_table()
    cols = S.POCKET_COLUMNS + (["source"] if "source" in schema.names else []) + ["fingerprint"]
    filt = pc.field("accession").isin(accessions) if accessions is not None else None
    t = pq.read_table(path, columns=cols, filters=filt)
    if "source" not in t.column_names:
        t = t.append_column("source", pa.nulls(t.num_rows, pa.string()))
    t = t.select(S.POCKET_COLUMNS + ["source", "fingerprint"])
    return t.cast(pa.schema([f.with_type(pa.string()) if pa.types.is_large_string(f.type) else f
                             for f in t.schema]))


def _atomic_write(table: pa.Table, path: Path, meta: dict):
    path.parent.mkdir(parents=True, exist_ok=True)
    md = dict(table.schema.metadata or {})
    md.update({k.encode(): v.encode() for k, v in meta.items()})
    tmp = path.with_suffix(f".tmp{os.getpid()}")
    pq.write_table(table.replace_schema_metadata(md), tmp, compression="zstd")
    os.replace(tmp, path)


def process_unit(args) -> dict:
    u, plan, extra_rows, atlas, out, jaccard, code_sha, force, with_fp = args
    atlas, out = Path(atlas), Path(out)
    uid, tier = int(u["unit_id"]), int(u["tier"])
    f_sites = out / "sites" / f"tier{tier}" / f"{uid:06d}.parquet"
    f_fp = out / "site_fp" / f"tier{tier}" / f"{uid:06d}.parquet"
    f_state = out / "state" / f"{uid:06d}.json"
    key = input_key(u, plan, extra_rows, jaccard, code_sha)
    if not force and f_state.exists() and f_sites.exists() and (f_fp.exists() or not with_fp):
        try:
            if json.loads(f_state.read_text()).get("input_key") == key:
                return dict(unit_id=uid, tier=tier, skipped=True)
        except Exception:  # noqa: BLE001
            pass
    t0 = time.time()
    for r in [u] + extra_rows:
        check_on_disk(atlas, r)
    t = read_pockets(atlas / u["pockets_path"])
    n_in = t.num_rows
    if plan["drop"]:
        t = t.filter(pc.invert(pc.is_in(t["accession"], pa.array(plan["drop"]))))
    n_dropped = n_in - t.num_rows
    extra = [read_pockets(atlas / r["pockets_path"], plan["own"]) for r in extra_rows]
    n_extra = sum(x.num_rows for x in extra)
    if extra:
        t = pa.concat_tables([t] + extra)
    method = (pq.read_schema(atlas / u["pockets_path"]).metadata or {}).get(b"pocketcompare", b"").decode()
    t1 = time.time()
    sites, fps, ps = S.build_sites(t, jaccard=jaccard, with_fingerprints=with_fp)
    # unit_id of a site = the unit file holding it (differs from the rep pocket's unit for split accessions)
    sites = sites.set_column(sites.schema.get_field_index("unit_id"), "unit_id",
                             pa.array(np.full(sites.num_rows, uid, np.int32)))
    t2 = time.time()
    # invariants: every input pocket in exactly one site
    n_members = int(pc.sum(sites["n_pockets"]).as_py() or 0)
    if n_members != t.num_rows or ps.num_rows != t.num_rows:
        raise RuntimeError(f"unit {uid}: member count {n_members} != pockets {t.num_rows}")
    meta = dict(pocketevo_sites=json.dumps(dict(schema=S.SCHEMA_VERSION, jaccard=jaccard, unit_id=uid, tier=tier,
                                                input_key=key, code_sha256=code_sha,
                                                fragment_rule=f"uniprot = author + {S.FRAGMENT_STEP}*(k-1)")),
                pocketcompare=method)
    _atomic_write(sites, f_sites, meta)
    if with_fp:
        _atomic_write(fps, f_fp, meta)
    rec = dict(unit_id=uid, tier=tier, input_key=key, n_pockets_file=n_in, n_pockets_dropped_split=n_dropped,
               n_pockets_extra_split=n_extra, n_pockets=t.num_rows, n_sites=sites.num_rows,
               n_accessions=int(pc.count_distinct(sites["accession"]).as_py()) if sites.num_rows else 0,
               t_read=t1 - t0, t_build=t2 - t1, t_write=time.time() - t2, t_total=time.time() - t0,
               jaccard=jaccard, code_sha256=code_sha, own_split=plan["own"], drop_split=plan["drop"])
    f_state.parent.mkdir(parents=True, exist_ok=True)
    tmp = f_state.with_suffix(f".tmp{os.getpid()}")
    tmp.write_text(json.dumps(rec, indent=1))
    os.replace(tmp, f_state)
    return dict(rec, skipped=False)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True, help="data/snapshot/<id>/manifest.parquet (manifest.json beside)")
    ap.add_argument("--tiers", default="1,2")
    ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--jaccard", type=float, default=S.DEFAULT_JACCARD)
    ap.add_argument("--sample", type=int, default=0, help="process a random subset of N units (timing)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--no-fingerprints", action="store_true",
                    help="skip <out>/site_fp (rep fingerprints = atlas rows of rep_pocket_id; saves ~0.75 KB/site, "
                         "needed for tier 3 at full scale: ~0.8 B sites x 0.77 KB = ~600 GB)")
    a = ap.parse_args()
    tiers = [int(x) for x in a.tiers.split(",")]
    man, mj = load_manifest(Path(a.manifest), tiers)
    atlas = Path(mj["atlas_root"])
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    code_sha = S.code_sha256()
    T0 = time.time()

    # split-accession pre-pass over ALL selected-tier units (cached per manifest + tiers)
    sp_path = out / f"split_accessions.tiers{'-'.join(map(str, tiers))}.parquet"
    sp_key = sha256_bytes(f"{mj['manifest_sha256']}|{tiers}".encode())
    split = None
    if sp_path.exists():
        md = pq.read_schema(sp_path).metadata or {}
        if md.get(b"key", b"").decode() == sp_key:
            split = pd.read_parquet(sp_path)
    if split is None:
        split = find_split_accessions(man, atlas, a.workers)
        tb = pa.Table.from_pandas(split, preserve_index=False)
        pq.write_table(tb.replace_schema_metadata({b"key": sp_key.encode()}), sp_path)
    t_pre = time.time() - T0
    print(f"pre-pass: {split.accession.nunique()} accessions split over units "
          f"{sorted(set(split.unit_id))} ({t_pre:.1f} s)", flush=True)

    units = man
    if a.sample:
        units = man.sample(n=min(a.sample, len(man)), random_state=a.seed).sort_values("unit_id")
    rows = {int(r["unit_id"]): r for r in man.to_dict("records")}
    jobs = []
    for u in units.sort_values("pockets_size", ascending=False).to_dict("records"):
        plan = unit_plan(u, split)
        jobs.append((u, plan, [rows[x] for x in plan["extra_units"]], str(atlas), str(out), a.jaccard,
                     code_sha, a.force, not a.no_fingerprints))
    results = []
    T1 = time.time()
    with ProcessPoolExecutor(a.workers) as ex:
        futs = [ex.submit(process_unit, j) for j in jobs]
        for i, f in enumerate(as_completed(futs), 1):
            results.append(f.result())
            if i % 100 == 0 or i == len(futs):
                print(f"{i}/{len(futs)} units, {time.time() - T1:.0f} s", flush=True)
    wall = time.time() - T1
    df = pd.DataFrame(results)
    done = df[~df.skipped] if len(df) else df
    summ = dict(manifest=str(a.manifest), manifest_sha256=mj["manifest_sha256"], tiers=tiers, out=str(out),
                workers=a.workers, jaccard=a.jaccard, sample=a.sample, seed=a.seed, n_units=len(jobs),
                n_processed=int(len(done)), n_skipped=int(df.skipped.sum()) if len(df) else 0,
                prepass_s=t_pre, wall_s=wall, code_sha256=code_sha)
    if len(done):
        summ.update(pockets=int(done.n_pockets.sum()), sites=int(done.n_sites.sum()),
                    cpu_s_sum=float(done.t_total.sum()), cpu_s_per_unit_mean=float(done.t_total.mean()),
                    s_per_1M_pockets_cpu=float(done.t_total.sum() / done.n_pockets.sum() * 1e6))
    tag = f"run_tiers{'-'.join(map(str, tiers))}" + (f"_sample{a.sample}_seed{a.seed}" if a.sample else "")
    (out / f"{tag}.json").write_text(json.dumps(summ, indent=1))
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()
