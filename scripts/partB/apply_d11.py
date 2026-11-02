#!/usr/bin/env python
"""PI decision D11: an A call additionally requires the member's footprint to span > 10 residues; otherwise U.

Member footprint residues = non-gap positions of that member at the SOG footprint columns (FAMSA MSA, bulk
/tmp/pe_big/sogs/<snap>/O1/msa). Applied to both strata (`state` primary, `state_s70` secondary); D8
small_footprint SOGs are left as they are (excluded downstream). Writes the same shard layout:
/tmp/pe_big/partB/<snap>/O1/states_d11_t<tau>/shard=*.parquet (+ added column n_fp_res) and a counts json.

    python scripts/partB/apply_d11.py --tau 0.2
"""
from __future__ import annotations

import argparse
import glob
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.dataset as ds

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from pocketevo import recon  # noqa: E402

SNAP = recon.cli_value("--snapshot") or recon.default_snapshot()   # --snapshot / $PE_SNAPSHOT / prereg
BULK = None   # set in main from --sog-run
_G = {}


def _shard(path):
    tau, out, min_res = _G["tau"], _G["out"], _G["min_res"]
    st = pd.read_parquet(path)
    ogs = set(st["og_id"])
    fp = _G["fp"]
    msa = _G["msa"]
    need = st["state"].eq("A") | st["state_s70"].eq("A")
    n = np.full(len(st), -1, dtype=np.int32)
    cache = {}
    for i in np.flatnonzero(need.to_numpy()):
        og, sog, acc = st.iat[i, 0], st.iat[i, 1], st.iat[i, 3]
        cols = fp.get(sog)
        aln = msa.get((og, acc))
        if cols is None or aln is None:
            continue
        key = (og, acc)
        if key not in cache:
            cache[key] = np.frombuffer(aln.encode(), dtype=np.uint8) != ord("-")
        a = cache[key]
        c = cols[cols < len(a)]
        n[i] = int(a[c].sum())
    st["n_fp_res"] = n
    cnt = {}
    for col in ("state", "state_s70"):
        m = (st[col] == "A") & (st["n_fp_res"] <= min_res)
        cnt[col] = int(m.sum())
        st.loc[m, col] = "U"
        if col == "state":
            st.loc[m, "reason"] = f"footprint<={min_res}res(D11)"
    st.to_parquet(out / Path(path).name, index=False)
    _ = ogs
    return cnt, int(need.sum()), int((n[need.to_numpy()] < 0).sum())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tau", type=float, default=0.2)
    ap.add_argument("--min-res", type=int, default=10, help="A requires > min_res footprint residues (D11)")
    ap.add_argument("--workers", type=int, default=16)
    recon.add_partb_args(ap, sog_run="O1")
    a = ap.parse_args()
    global BULK
    BULK = Path(f"/tmp/pe_big/sogs/{SNAP}/{a.sog_run}")
    t0 = time.time()
    out = Path(f"/tmp/pe_big/partB/{SNAP}/{a.sog_run}/states_d11_t{a.tau}")
    out.mkdir(parents=True, exist_ok=True)
    fp = ds.dataset(sorted(glob.glob(str(BULK / "sogs/shard=*.parquet"))), format="parquet").to_table(
        columns=["sog_id", "tau", "footprint_columns"]).to_pandas()
    fp = fp[np.isclose(fp["tau"], a.tau)]
    _G["fp"] = {s: np.asarray(c, dtype=int) for s, c in zip(fp["sog_id"], fp["footprint_columns"])}
    msa = ds.dataset(sorted(glob.glob(str(BULK / "msa/shard=*.parquet"))), format="parquet").to_table().to_pandas()
    _G["msa"] = {(o, x): s for o, x, s in zip(msa["og_id"], msa["accession"], msa["aligned"])}
    del msa
    _G.update(tau=a.tau, out=out, min_res=a.min_res)
    files = sorted(glob.glob(str(BULK / f"states_t{a.tau}/shard=*.parquet")))
    tot = {"state": 0, "state_s70": 0}
    n_need = n_miss = 0
    from multiprocessing import get_context
    with ProcessPoolExecutor(a.workers, mp_context=get_context("fork")) as ex:
        for c, nn, nm in ex.map(_shard, files, chunksize=8):
            for k in tot:
                tot[k] += c[k]
            n_need += nn; n_miss += nm
    rec = dict(decision="D11", rule=f"A requires member footprint > {a.min_res} residues, else U", tau=a.tau,
               n_A_rows_checked=n_need, n_without_msa_or_footprint=n_miss, A_to_U=tot, out=str(out),
               runtime_s=time.time() - t0)
    f = ROOT / f"results/partB/{SNAP}/{a.sog_run}/d11_t{a.tau}.json"
    f.parent.mkdir(parents=True, exist_ok=True)
    f.write_text(json.dumps(rec, indent=1))
    recon.write_sidecar(f, SNAP, 0, {"script": "scripts/partB/apply_d11.py", "argv": sys.argv[1:]})
    print(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
