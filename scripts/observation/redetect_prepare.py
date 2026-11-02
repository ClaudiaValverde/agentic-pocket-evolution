#!/usr/bin/env python
"""Re-detection experiment, compute part (bulk in /tmp/pe_big/obs):

  1. extract the AFDB v6 model of every sampled accession from its proteome tar (atlas tar index; seek)
  2. decompress the AFDB v4 Foldcomp entry (afdb_uniprot_v4, offsets from the full .index/.lookup)
  3. sequence check v6 vs v4 (exact sequence equality of the modelled residues; mismatches skipped + counted)
  4. `pocketcompare describe --af-mode --backend fast` on the v4 models (frozen method = atlas method record)
  5. seeds subset: describe_burrow_seed.py with BURROW s0 (control: must reproduce the atlas), s1, s2 on v6

  python scripts/observation/redetect_prepare.py --sample /tmp/pe_big/obs/redetect_sample.parquet --workers 6
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

P = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(P / "src"))
sys.path.insert(0, str(P / "scripts"))
from pocketevo import observation as O  # noqa: E402

AFDB_V6 = Path("/cluster/work/beltrao/shared/alphafold/v6")
FC_DB = "/cluster/scratch/gmagnusson/pocketcompare_afdb/foldcomp/afdb_uniprot_v4"
PY = str(P / ".venv/bin/python")
PCBIN = str(P / ".venv/bin/pocketcompare")


def _decomp(rows):
    out = []
    for name, off, ln, dest in rows:
        try:
            n, pdb = O.foldcomp_read(FC_DB, off, ln)
            if not n.startswith(name):
                raise ValueError(f"entry {n} != {name}")
            Path(dest).write_text(pdb)
            out.append((name, None))
        except Exception as e:  # noqa: BLE001
            out.append((name, f"{type(e).__name__}: {e}"))
    return out


def _seq(path):
    try:
        s, num, pl = O.read_model(path)
        return s, int(num.min()) if len(num) else -1, len(s), float(pl.mean()) if len(pl) else np.nan, None
    except Exception as e:  # noqa: BLE001
        return None, -1, 0, np.nan, f"{type(e).__name__}: {e}"


def chunks(xs, k):
    k = max(1, k)
    return [xs[i::k] for i in range(k) if xs[i::k]]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", default="/tmp/pe_big/obs/redetect_sample.parquet")
    ap.add_argument("--work", default="/tmp/pe_big/obs")
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--skip-seeds", action="store_true")
    a = ap.parse_args()
    env = dict(os.environ, OMP_NUM_THREADS="1", NUMBA_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    W = Path(a.work)
    v6d, v4d = W / "v6", W / "v4"
    v6d.mkdir(parents=True, exist_ok=True)
    v4d.mkdir(parents=True, exist_ok=True)
    S = pd.read_parquet(a.sample)
    T = {}
    t0 = time.time()
    # 1. v6 extraction
    from extract_afdb_member import extract
    for tar, g in S.groupby("source"):
        need = [m for m in g.member if not (v6d / m).exists()]
        if need:
            extract(AFDB_V6 / tar, need, v6d)
    T["extract_v6_s"] = time.time() - t0
    # 2. v4 decompression
    t1 = time.time()
    rows = [(r.fc_name, int(r.off), int(r.len), str(v4d / f"{r.fc_name}.pdb")) for r in S.itertuples()
            if not (v4d / f"{r.fc_name}.pdb").exists()]
    errs = {}
    with ProcessPoolExecutor(min(4, a.workers)) as ex:
        for res in ex.map(_decomp, chunks(rows, 4 * min(4, a.workers))):
            for n, e in res:
                if e:
                    errs[n] = e
    T["decompress_v4_s"] = time.time() - t1
    # 3. sequence check
    t2 = time.time()
    p6 = [str(v6d / m) for m in S.member]
    p4 = [str(v4d / f"{n}.pdb") for n in S.fc_name]
    with ProcessPoolExecutor(a.workers) as ex:
        r6 = list(ex.map(_seq, p6, chunksize=16))
        r4 = list(ex.map(_seq, p4, chunksize=16))
    S["seq_v6_md5"] = [O.seq_md5(r[0]) if r[0] else None for r in r6]
    S["seq_v4_md5"] = [O.seq_md5(r[0]) if r[0] else None for r in r4]
    S["len_v6"] = [r[2] for r in r6]
    S["len_v4"] = [r[2] for r in r4]
    S["start_v6"] = [r[1] for r in r6]
    S["start_v4"] = [r[1] for r in r4]
    S["plddt_v6"] = [r[3] for r in r6]
    S["plddt_v4"] = [r[3] for r in r4]
    S["v4_error"] = [errs.get(n) or r[4] for n, r in zip(S.fc_name, r4)]
    S["seq_identical"] = (S.seq_v6_md5 == S.seq_v4_md5) & S.seq_v6_md5.notna()
    # identity for mismatches (global, simple): fraction of equal residues at equal positions if same length
    ident = []
    for (s6, *_), (s4, *_) in zip(r6, r4):
        if not s6 or not s4:
            ident.append(np.nan)
        elif s6 == s4:
            ident.append(1.0)
        else:
            from difflib import SequenceMatcher
            ident.append(SequenceMatcher(None, s6, s4, autojunk=False).ratio())
    S["seq_identity_v6_v4"] = ident
    T["seqcheck_s"] = time.time() - t2
    S.to_parquet(W / "redetect_models.parquet", index=False)
    # 4. describe v4 (identical sequence only)
    t3 = time.time()
    ok = S[S.seq_identical]
    lst = W / "v4_inputs.txt"
    lst.write_text("\n".join(str(v4d / f"{n}.pdb") for n in ok.fc_name) + "\n")
    out4 = W / "v4_describe.parquet"
    if not out4.exists():
        cmd = [PCBIN, "describe", *[str(v4d / f"{n}.pdb") for n in ok.fc_name], "--af-mode", "--backend", "fast",
               "--workers", str(a.workers), "-o", str(out4)]
        subprocess.run(cmd, check=True, env=env)
    T["describe_v4_s"] = time.time() - t3
    T["n_describe_v4"] = int(len(ok))
    # 5. seeds subset on v6 (s0 control, s1, s2)
    if not a.skip_seeds:
        sub = S[S.in_seed_subset & S.seq_identical]
        paths = [str(v6d / m) for m in sub.member]
        for seed in ("s0", "s1", "s2"):
            t4 = time.time()
            outs = []
            procs = []
            for k, ch in enumerate(chunks(paths, a.workers)):
                o = W / f"seed_{seed}_part{k}.parquet"
                outs.append(o)
                if not o.exists():
                    procs.append(subprocess.Popen([PY, str(P / "scripts/observation/describe_burrow_seed.py"),
                                                   "--seed", seed, "-o", str(o), *ch], env=env,
                                                  stdout=subprocess.DEVNULL, stderr=subprocess.PIPE))
            for p_ in procs:
                _, err = p_.communicate()
                if p_.returncode != 0:
                    raise RuntimeError(f"seed {seed}: {err.decode()[-2000:]}")
            T[f"describe_seed_{seed}_s"] = time.time() - t4
        T["n_seed_subset"] = int(len(sub))
    summ = dict(n_sample=int(len(S)), n_seq_identical=int(S.seq_identical.sum()),
                n_seq_mismatch=int((~S.seq_identical).sum()), n_v4_error=int(S.v4_error.notna().sum()),
                mismatch_len_diff=int((S.len_v6 != S.len_v4).sum()), timings=T,
                wall_s=round(time.time() - t0, 1))
    (W / "redetect_prepare.json").write_text(json.dumps(summ, indent=1, default=float))
    print(json.dumps(summ, indent=1, default=float))


if __name__ == "__main__":
    main()
