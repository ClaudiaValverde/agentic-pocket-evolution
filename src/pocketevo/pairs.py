"""Part A pair universe and exact kNN (DESIGN §4.2, round-2 amendment 1).

Pair unit: an unordered pair of tier-1/2 *sites* on different proteins; each site is represented by the
fingerprint of its highest-burrow-score pocket (`rep_pocket_id`, copied unchanged from the atlas; see
results/sites/<snapshot>/SITES.md). Nothing in this module looks at homology; see homology.py / nulls.py.

Bulk layout (one directory per snapshot, default /tmp/pe_big/partA/<snapshot>/):
  universe.parquet        one row per site, row number = site index `si` used everywhere downstream
  fp.f16.npy              (N, 384) float16 fingerprints in universe order (np.load(mmap_mode='r'))
  knn/I.i32.npy, D.f16.npy  (N, k) neighbour index / inner product, rows in universe order, self included
  knn/blocks/*.done       resumability markers (one per query block)
"""
from __future__ import annotations

import glob
import hashlib
import json
import os
import subprocess
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

REPO = Path(__file__).resolve().parents[2]
PREREG = REPO / "config" / "prereg.yaml"
FROZEN_PREREG_SHA256 = "52f9461890984c210b730af86c40de221307a668e7d4b38098b3316059690d83"
DIM = 384

UNIVERSE_COLS = ["site_id", "accession", "site_index", "taxid", "tier", "unit_id", "source", "n_res", "res_min",
                 "res_max", "max_burrow_score", "plddt_mean_rep", "plddt_frac_lt70_max", "rep_site_volume_A3",
                 "volume_max", "n_fragments", "residues"]


# ----------------------------------------------------------------------------------------------- provenance
def sha256_file(path: str | Path, chunk: int = 1 << 22) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def git_commit() -> str:
    try:
        out = subprocess.run(["git", "-C", str(REPO), "rev-parse", "HEAD"], capture_output=True, text=True,
                             check=True).stdout.strip()
        mine = [f"src/pocketevo/{m}.py" for m in ("pairs", "nulls", "homology", "estimands", "directional")]
        dirty = subprocess.run(["git", "-C", str(REPO), "status", "--porcelain", *mine, "scripts/partA",
                                "tests/test_partA.py", "config/partA_ligand_map.tsv"],
                               capture_output=True, text=True).stdout.strip()
        return out + ("+dirty" if dirty else "")
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def load_prereg() -> dict:
    import yaml  # local import: only scripts need it
    with open(PREREG) as fh:
        return yaml.safe_load(fh)


def prereg_sha256() -> str:
    return sha256_file(PREREG)


def check_prereg(allow_mismatch: bool = False) -> str:
    h = prereg_sha256()
    if h != FROZEN_PREREG_SHA256 and not allow_mismatch:
        raise RuntimeError(f"config/prereg.yaml sha256 {h} != frozen {FROZEN_PREREG_SHA256}")
    return h


def write_sidecar(out_path: str | Path, snapshot_id: str, seed, extra: dict | None = None) -> Path:
    """`<out>.json` with the DESIGN §0 sidecar fields plus code hashes and anything stage-specific."""
    out_path = Path(out_path)
    code = {p.name: sha256_file(p)[:16] for p in sorted((REPO / "src" / "pocketevo").glob("*.py"))
            if p.name in ("pairs.py", "nulls.py", "homology.py", "estimands.py", "directional.py")}
    rec = {"snapshot_id": snapshot_id, "prereg_sha256": prereg_sha256(), "git_commit": git_commit(),
           "seed": seed, "code_sha256_16": code, "written_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    if extra:
        rec.update(extra)
    side = out_path.with_name(out_path.name + ".json") if out_path.suffix != ".json" else out_path
    side.write_text(json.dumps(rec, indent=1, default=_json_default))
    return side


def _json_default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


# ----------------------------------------------------------------------------------------------- universe
def site_files(sites_dir: str | Path, tiers=(1, 2)) -> list[tuple[Path, Path]]:
    out = []
    for t in tiers:
        meta = sorted(glob.glob(str(Path(sites_dir) / "sites" / f"tier{t}" / "*.parquet")))
        for m in meta:
            fp = Path(sites_dir) / "site_fp" / f"tier{t}" / Path(m).name
            if not fp.exists():
                raise FileNotFoundError(fp)
            out.append((Path(m), fp))
    return out


def build_universe(sites_dir: str | Path, out_dir: str | Path, tiers=(1, 2)) -> dict:
    """Concatenate site metadata + rep fingerprints of the given tiers in (tier, unit, row) order."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    files = site_files(sites_dir, tiers)
    n_per = [pq.ParquetFile(m).metadata.num_rows for m, _ in files]
    n = int(sum(n_per))
    fp = np.lib.format.open_memmap(out_dir / "fp.f16.npy.tmp", mode="w+", dtype=np.float16, shape=(n, DIM))
    metas, off, n_bad_norm = [], 0, 0
    for (m, f), k in zip(files, n_per):
        meta = pd.read_parquet(m, columns=UNIVERSE_COLS)
        fpt = pq.read_table(f, columns=["site_id", "fingerprint"])
        if k == 0:
            continue
        if not np.array_equal(fpt.column("site_id").to_numpy(zero_copy_only=False), meta.site_id.to_numpy()):
            raise ValueError(f"site_id order differs between {m} and {f}")
        v = fpt.column("fingerprint").combine_chunks()
        arr = np.asarray(v.values.to_numpy(zero_copy_only=False), dtype=np.float16).reshape(k, DIM)
        nrm = np.linalg.norm(arr.astype(np.float32), axis=1)
        n_bad_norm += int((np.abs(nrm - 1) > 1e-2).sum())
        fp[off:off + k] = arr
        metas.append(meta)
        off += k
    fp.flush()
    del fp
    os.replace(out_dir / "fp.f16.npy.tmp", out_dir / "fp.f16.npy")
    uni = pd.concat(metas, ignore_index=True)
    uni.insert(0, "si", np.arange(len(uni), dtype=np.int64))
    uni.to_parquet(out_dir / "universe.parquet", index=False)
    return {"n_sites": n, "n_files": len(files), "n_bad_norm": n_bad_norm,
            "n_proteins": int(uni.accession.nunique()), "by_tier": uni.tier.value_counts().to_dict()}


# ----------------------------------------------------------------------------------------------- kNN
def knn_exact(fp_path: str | Path, out_dir: str | Path, k: int = 1000, device: int = 0, block: int = 32768,
              fp16_storage: bool = False, log=print) -> dict:
    """Exact (brute-force) inner-product kNN of every site against all sites, FAISS GPU flat index.

    Resumable: each finished query block writes `blocks/<b>.done`; I/D memmaps are opened r+ on restart.
    Self-matches are kept (column 0 normally); H0/self are removed downstream.
    """
    import faiss
    out_dir = Path(out_dir)
    (out_dir / "blocks").mkdir(parents=True, exist_ok=True)
    X = np.load(fp_path, mmap_mode="r")
    n = X.shape[0]
    t0 = time.time()
    res = faiss.StandardGpuResources()
    cfg = faiss.GpuIndexFlatConfig()
    cfg.device = device
    cfg.useFloat16 = fp16_storage
    index = faiss.GpuIndexFlatIP(res, DIM, cfg)
    for s in range(0, n, 500_000):
        index.add(np.ascontiguousarray(X[s:s + 500_000], dtype=np.float32))
    t_add = time.time() - t0
    ip, dp = out_dir / "I.i32.npy", out_dir / "D.f16.npy"
    mode = "r+" if ip.exists() and dp.exists() else "w+"
    I = np.lib.format.open_memmap(ip, mode=mode, dtype=np.int32, shape=(n, k))
    D = np.lib.format.open_memmap(dp, mode=mode, dtype=np.float16, shape=(n, k))
    nb = (n + block - 1) // block
    t_search, done_before = 0.0, 0
    for b in range(nb):
        mark = out_dir / "blocks" / f"{b:05d}.done"
        if mark.exists():
            done_before += 1
            continue
        s, e = b * block, min(n, (b + 1) * block)
        q = np.ascontiguousarray(X[s:e], dtype=np.float32)
        t1 = time.time()
        d, i = index.search(q, k)
        t_search += time.time() - t1
        if (i < 0).any():
            raise RuntimeError(f"block {b}: faiss returned -1 indices")
        I[s:e] = i.astype(np.int32)
        D[s:e] = d.astype(np.float16)
        if b % 10 == 0:
            I.flush(); D.flush()
            log(f"knn block {b + 1}/{nb} rows {s}-{e} search {t_search:.0f}s")
        mark.write_text(json.dumps({"rows": [s, e], "s": time.time() - t1}))
    I.flush(); D.flush()
    return {"n": n, "k": k, "block": block, "fp16_storage": fp16_storage, "t_add_s": t_add,
            "t_search_s_this_run": t_search, "blocks": nb, "blocks_done_before": done_before,
            "wall_s_this_run": time.time() - t0}
