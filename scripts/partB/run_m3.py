#!/usr/bin/env python
"""M3 tree consistency on real SOG states (DESIGN §4.1 M3; prereg part_a.M3) -> results/partB/m3/<snap>/<orth>/
and results/numbers.json["partA_descent"]["M3"].

SOG x species matrix from the member states (strongest state per species: P > A > U; W/P -> 1, A -> 0, U and
no member -> missing). SOGs: footprint > 3 columns (D8) and P in >= 2 species (single-species SOGs are trivially
consistent; the all-SOG value is reported too). M3 = fraction with <= 1 Dollo event on the NCBI guide tree.
Nulls: masked curveball (per-SOG / per-species totals and missing mask fixed; 100 x n(1) swaps; prereg
replicates) and noisy-but-vertical (every informative cell truly present, observed A with the measured per-species
rate P(A | informative, present) = miss q / (miss q + 1 - miss); pre-rescue, matching the input states).

    python scripts/partB/run_m3.py --tau 0.2
"""
from __future__ import annotations

import argparse
import glob
import json
import sys
import time
from multiprocessing import get_context
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.dataset as ds

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from pocketevo import frozen, recon, simulate  # noqa: E402

SNAP = recon.cli_value("--snapshot") or recon.default_snapshot()   # --snapshot / $PE_SNAPSHOT / prereg
_G = {}


def m3_of(X):
    ev, _ = recon.dollo_events(_G["tree"], X)
    d = ev >= 0
    return float(((ev <= 1) & d).sum() / max(1, d.sum()))


def _rep(r):
    Y, done, att = simulate.checkerboard(_G["X"], n_swaps=_G["n_swaps"], seed=_G["seed"] * 1000 + r)
    return r, m3_of(Y), done, att


def update_numbers(key_path: list, value):
    f = ROOT / "results/numbers.json"
    d = json.loads(f.read_text()) if f.exists() else {}
    cur = d
    for k in key_path[:-1]:
        cur = cur.setdefault(k, {})
    cur[key_path[-1]] = value
    f.write_text(json.dumps(d, indent=1, default=float))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tau", type=float, default=0.2)
    ap.add_argument("--state-col", default="state")
    ap.add_argument("--guide", default=str(ROOT / "data/orthology/ancestra_tier1_v0/tier1_guide_tree.nwk"))
    ap.add_argument("--noise-csv", default=str(ROOT / "results/partB/sim/measured/noise_primary_known_gt10.csv"))
    ap.add_argument("--reps", type=int, default=None)
    ap.add_argument("--states-dir", default=None)
    ap.add_argument("--swap-mult", type=int, default=100, help="swaps = mult x min(#1, #0); D22 mixing check 50/200")
    ap.add_argument("--only", default=None, help="run only this SOG set (multi_species | all_with_P)")
    ap.add_argument("--release", "--sog-run", dest="release", default="O1")
    recon.add_partb_args(ap, skip=("sog_run",))
    ap.add_argument("--reps-all", type=int, default=100, help="replicates for the secondary all-SOG variant")
    ap.add_argument("--time-one", action="store_true", help="time one curveball replicate and exit")
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--numbers", action="store_true", help="write results/numbers.json partA_descent.M3")
    a = ap.parse_args()
    t0 = time.time()
    pr = recon.load_prereg()
    seeds = pr["meta"]["seeds"]
    n_rep = a.reps or int(pr["part_a"]["M3"]["nulls"]["curveball"]["n_replicates"]["value"])
    bulk = Path(f"/tmp/pe_big/sogs/{SNAP}/{a.release}")
    st = recon.SpeciesTree.from_file(a.guide)
    sp = st.tree.leaf_names()
    col = {int(st.tree.taxid[v]): i for i, v in enumerate(st.tree.leaves)}
    d = ds.dataset(sorted(glob.glob(str(Path(a.states_dir or bulk / f"states_t{a.tau}") / "shard=*.parquet"))), format="parquet").to_table(
        columns=["sog_id", "taxid", a.state_col]).to_pandas().rename(columns={a.state_col: "state"})
    meta = pd.read_parquet(ROOT / f"results/sogs/{SNAP}/{a.release}/sogs_all_tau.parquet")
    meta = meta[np.isclose(meta["tau"], a.tau)]
    small = set(meta.loc[meta["footprint_size"] <= 3, "sog_id"])
    for flag in ("small_footprint", "empty_footprint"):
        if flag in meta.columns:
            small |= set(meta.loc[meta[flag].astype(bool), "sog_id"])
    n_all = d["sog_id"].nunique()
    d = d[~d["sog_id"].isin(small)]
    rank = d["state"].map({"P": 3, "W": 3, "A": 1, "U": 0, "G": -1}).astype(np.int8)
    d = d.assign(rank=rank, j=d["taxid"].map(col))
    d = d[d["j"].notna()]
    g = d.groupby(["sog_id", "j"])["rank"].max().reset_index()
    ids = pd.Index(sorted(g["sog_id"].unique()))
    X = np.full((len(ids), len(sp)), -1, dtype=np.int8)
    v = g["rank"].to_numpy()
    X[ids.get_indexer(g["sog_id"]), g["j"].astype(int).to_numpy()] = np.where(v == 3, 1, np.where(v == 1, 0, -1))
    nPs = (X == 1).sum(1)
    res = {}
    (ROOT / f"results/partB/m3/{SNAP}/{a.release}/tau{a.tau}").mkdir(parents=True, exist_ok=True)
    rk = json.loads(frozen.path("rescue_k").read_text())   # D9 frozen (snapshot-invariant)
    nz = pd.read_csv(a.noise_csv)
    pa_inf = dict(zip(nz["species"], nz["miss"] * nz["q_present"] / (nz["miss"] * nz["q_present"] + 1 - nz["miss"])))
    fa = np.array([pa_inf.get(s, float(np.mean(list(pa_inf.values())))) for s in sp])[None, :]
    for name, mask in (("multi_species", nPs >= 2), ("all_with_P", nPs >= 1)):
        if a.only and name != a.only:
            continue
        Xm = np.ascontiguousarray(X[mask])
        # INTERP (for PI): prereg n_swaps = 100 x number of 1s; with sparse informative absences (A after D11)
        # a swap needs an informative 0 at both targets, so the minority class bounds mixing: target =
        # 100 x min(#1, #0) successful swaps (each 0-cell is moved ~100 times on average).
        n_sw = int(a.swap_mult * min(int((Xm == 1).sum()), int((Xm == 0).sum())))
        _G.update(tree=st.tree, X=Xm, seed=seeds[0], n_swaps=n_sw)
        print(name, Xm.shape, "ones", int((Xm == 1).sum()), "zeros", int((Xm == 0).sum()), "n_swaps", n_sw, flush=True)
        if a.time_one:
            t1 = time.time(); _rep(0); print(name, Xm.shape, int((Xm == 1).sum()), "s/rep", time.time() - t1)
            continue
        real = m3_of(Xm)
        nr = n_rep if name == "multi_species" else min(n_rep, a.reps_all)
        with get_context("fork").Pool(a.workers) as pool:
            null = sorted(pool.map(_rep, range(nr)))
        nd = np.array([x[1] for x in null])
        rng = np.random.default_rng(seeds[0])
        vert = []
        for _ in range(100):
            Xs = np.where(Xm >= 0, (rng.random(Xm.shape) >= fa).astype(np.int8), -1).astype(np.int8)
            vert.append(m3_of(Xs))
        res[name] = dict(n_sogs=int(mask.sum()), m3=real, curveball_mean=float(nd.mean()),
                         curveball_q025=float(np.quantile(nd, 0.025)), curveball_q975=float(np.quantile(nd, 0.975)),
                         p_null_ge_real=(1 + int((nd >= real).sum())) / (1 + nr), curveball_replicates=nr,
                         swaps_per_replicate_min=int(min(x[2] for x in null)),
                         swaps_target=n_sw, n_zeros=int((Xm == 0).sum()),
                         swaps_rule=f"{a.swap_mult} x min(#1, #0) (INTERP of prereg 100 x #1; see SIM.md)",
                         attempts_per_replicate_max=int(max(x[3] for x in null)),
                         null_method="masked checkerboard swaps (simulate.checkerboard)",
                         n_ones=int((Xm == 1).sum()), vertical_noisy_mean=float(np.mean(vert)),
                         vertical_noisy_sd=float(np.std(vert)))
        pd.DataFrame({"replicate": range(nr), "m3": nd, "swaps": [x[2] for x in null]}).to_csv(
            ROOT / f"results/partB/m3/{SNAP}/{a.release}/tau{a.tau}/curveball_{name}_x{a.swap_mult}.csv", index=False)
        print(f"[{time.time() - t0:.0f}s] {name}: {res[name]}", flush=True)
    if a.time_one:
        return
    out = ROOT / f"results/partB/m3/{SNAP}/{a.release}/tau{a.tau}"
    out.mkdir(parents=True, exist_ok=True)
    rec = dict(snapshot_id=SNAP, orthology=a.orthology, release=a.release, swap_mult=a.swap_mult, tau=a.tau, state_col=a.state_col, n_sogs_total=n_all,
               n_small_footprint_excluded=len(small),
               states="SOG states before rescue (A not rescued; D11 applied)",
               vertical_noisy_rate_source=f"{a.noise_csv} (pre-rescue P(A|informative,present))",
               rescue_k=rk.get("k_min_false_A"), definition=pr["part_a"]["M3"]["definition"], results=res,
               prereg_sha256=recon.prereg_sha256(), runtime_s=time.time() - t0)
    f = out / ("m3.json" if a.swap_mult == 100 and not a.only else f"m3_x{a.swap_mult}{'_' + a.only if a.only else ''}.json")
    f.write_text(json.dumps(rec, indent=1))
    recon.write_sidecar(f, SNAP, seeds, {"script": "scripts/partB/run_m3.py", "argv": sys.argv[1:]})
    if a.numbers and a.swap_mult == 100 and "multi_species" in res:
        update_numbers(["partA_descent", "M3"], {"tau": a.tau, "source": str(f.relative_to(ROOT)),
                                                 "prereg_sha256": rec["prereg_sha256"], **res["multi_species"],
                                                 "all_with_P": res["all_with_P"], "states": rec["states"]})
    print(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
