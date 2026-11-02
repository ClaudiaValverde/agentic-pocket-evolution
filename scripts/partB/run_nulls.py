#!/usr/bin/env python
"""M3 tree consistency + the Part B null-input pipeline on any SOG state table (DESIGN §4.1 M3, §5.4 controls;
prereg part_a.M3, part_b.negative_controls).

M3: SOG x species matrix (strongest member state; P/W -> 1, A -> 0, U -> missing, no member -> missing (G)),
Dollo parsimony on the NCBI guide tree: M3 = fraction of SOGs (with >= 1 present species) explained by <= 1
gain/loss. Nulls: (i) masked curveball (per-SOG and per-species totals and the missing mask preserved;
n_swaps = 100 x number of 1s, 1000 replicates, prereg) -> one-sided p = P(null M3 >= real);
(ii) simulated noisy-but-vertical: every informative cell truly present, observed A with the false-A rate.

Null-input pipeline (same code and thresholds as run_events.py, primary costs; prereg seeds):
  curveball   : the permuted species matrix mapped back to members;
  pseudo_clade: member states shuffled among the leaves of each SOG tree (residues/fingerprints too, for the
                specificity test);
  replicate   : same proteins observed twice (replicate_states/ in the SOG dir, e.g. second fragment window /
                model source / seed) as two copies of the tree; calls on a copy stem are false by construction.
Output: results/partB/nulls/<snapshot>/<orthology>[/<stratum>]/ funnel_real_vs_null.csv, m3.json,
m3_curveball_null.csv, manifest.json (+ sidecars).

    python scripts/partB/run_nulls.py --sogs results/sogs/<snap>/O2 --orthoxml <O2 orthoxml> --snapshot <snap>
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from multiprocessing import get_context
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from pocketevo import events, recon, simulate  # noqa: E402
import run_events as RE  # noqa: E402

_M3 = {}


def m3_value(tree, X) -> tuple:
    ev, _ = recon.dollo_events(tree, X)
    d = ev >= 0
    return float(((ev <= 1) & d).sum() / max(1, d.sum())), int(d.sum())


def _m3_rep(rep):
    X = _M3["X"]
    Y, done, _ = simulate.checkerboard(X, n_swaps=int(_M3["spo"] * min(int((X == 1).sum()), int((X == 0).sum()))),
                                       seed=_M3["seed"] * 1000 + rep)
    v, n = m3_value(_M3["tree"], Y)
    return rep, v, done


def species_order_matrix(assembled, trees, st):
    """SOG x species-tree-leaf matrix in species-tree leaf order."""
    # match species by taxid (the SOG web species_tree.nwk has taxid leaves; gene trees carry UP..._taxid names)
    sp_order = [str(int(st.tree.taxid[v])) for v in st.tree.leaves]
    codes_l, sp_l, ids = [], [], []
    for og, sog_ids, codes, _ in assembled:
        sp = [str(int(x)) for x in trees[og].leaf_taxids()]
        for j, s in enumerate(sog_ids):
            codes_l.append(codes[:, j]); sp_l.append(sp); ids.append(s)
    M = simulate.species_matrix(codes_l, sp_l, sp_order)
    return M, codes_l, sp_l, ids, sp_order


def rebuild(assembled, new_codes_list):
    out, k = [], 0
    for og, sog_ids, codes, conf in assembled:
        c = np.column_stack(new_codes_list[k:k + len(sog_ids)]).astype(np.int8)
        k += len(sog_ids)
        out.append((og, sog_ids, c, conf))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sogs", required=True)
    ap.add_argument("--orthoxml", default=str(RE.O1))
    ap.add_argument("--guide", default=None)
    ap.add_argument("--thresholds", default=str(RE.THRESH))
    ap.add_argument("--snapshot", required=True)
    ap.add_argument("--orthology", default=None)
    ap.add_argument("--stratum", default="primary")
    ap.add_argument("--m3-reps", type=int, default=None, help="default prereg n_replicates (1000)")
    ap.add_argument("--noise-json", default=None, help="measured false-A rates for the vertical-noisy M3 null")
    ap.add_argument("--noise-floor", type=float, default=None)
    ap.add_argument("--workers", type=int, default=recon.n_workers())
    ap.add_argument("--out", default=None)
    ap.add_argument("--states-glob", default=None)
    ap.add_argument("--state-col", default="state")
    ap.add_argument("--sog-meta", default=None)
    ap.add_argument("--tau", type=float, default=None)
    ap.add_argument("--thresholds-le10", default=None)
    ap.add_argument("--thresholds-gt10", default=None)
    ap.add_argument("--prefilter", action="store_true")
    ap.add_argument("--n-seeds", type=int, default=None, help="null seeds 0..n-1 (D24: >= 20); default prereg seeds")
    ap.add_argument("--skip-m3", action="store_true", help="M3 computed by run_m3.py on all SOGs")
    recon.add_partb_args(ap, skip=("snapshot", "orthology"))
    a = ap.parse_args()
    t0 = time.time()
    pr = recon.load_prereg()
    seeds = pr["meta"]["seeds"]
    null_seeds = list(range(a.n_seeds)) if a.n_seeds else seeds
    sog_dir = Path(a.sogs)
    orth = a.orthology or sog_dir.name
    out = Path(a.out) if a.out else ROOT / "results/partB/nulls" / a.snapshot / orth / a.stratum
    out.mkdir(parents=True, exist_ok=True)
    th = RE.load_thresholds(a.thresholds)
    cnt = {}
    states, trees, st, nwk, fpsize = RE.load_inputs(sog_dir, a.orthoxml, a.guide, a.states_glob, a.state_col,
                                                    a.sog_meta, a.tau, prefilter=a.prefilter, counters=cnt)
    assembled = events.assemble(trees, states, cnt)
    del states
    costs = recon.costs_from_prereg(pr)["primary"]
    settings = events.rules_by_cost(th, [costs], pr)
    by = None
    if a.thresholds_le10 and a.thresholds_gt10:
        by = RE.settings_by_class({"le10": a.thresholds_le10, "gt10": a.thresholds_gt10}, [costs], pr)
        th = {"chosen": {c: RE.load_thresholds(p_)["chosen"] for c, p_ in
                         (("le10", a.thresholds_le10), ("gt10", a.thresholds_gt10))}}

    def run(asm, trees_=None, st_=st):
        trees_ = trees_ or trees
        if by is None:
            return events.run_assembled(asm, trees_, settings, st_, workers=a.workers)[0]
        ev_ = []
        for cls, part in RE.split_by_class(asm, fpsize).items():
            ev_ += events.run_assembled(part, trees_, by[cls], st_, workers=a.workers)[0]
        return ev_

    n_chars = sum(len(s) for _, s, _, _ in assembled)

    # ---------------- M3
    M, codes_l, sp_l, ids, sp_order = species_order_matrix(assembled, trees, st)
    X = np.where(M == 1, 1, np.where(M == 0, 0, -1)).astype(np.int8)
    m3_real, n_def = m3_value(st.tree, X)
    if a.skip_m3:
        a.m3_reps = 1
    m3cfg = pr["part_a"]["M3"]["nulls"]["curveball"]
    n_rep = a.m3_reps or int(m3cfg["n_replicates"]["value"])
    spo = 100   # prereg n_swaps "100 x number of 1s"
    _M3.update(X=X, tree=st.tree, seed=seeds[0], spo=spo)
    with get_context("fork").Pool(a.workers) as pool:
        null = sorted(pool.map(_m3_rep, range(n_rep)))
    nd = pd.DataFrame(null, columns=["replicate", "m3", "swaps_done"])
    p_cb = (1 + int((nd["m3"] >= m3_real).sum())) / (1 + n_rep)
    noise = simulate.NoiseModel.from_json(a.noise_json) if a.noise_json else simulate.NoiseModel()
    rng = np.random.default_rng(seeds[0])
    vert = []
    fa = noise.rates(sp_order)["false_A"][None, :]
    for _ in range(100):
        Xs = np.where(X >= 0, (rng.random(X.shape) >= fa).astype(np.int8), -1).astype(np.int8)
        vert.append(m3_value(st.tree, Xs)[0])
    m3 = dict(m3_real=m3_real, n_sogs_defined=n_def, n_sogs=len(ids), curveball_replicates=n_rep,
              curveball_mean=float(nd["m3"].mean()), curveball_q025=float(nd["m3"].quantile(0.025)),
              curveball_q975=float(nd["m3"].quantile(0.975)), p_real_ge_null=p_cb,
              swaps_target=int(spo * int((X == 1).sum())), swaps_done_min=int(nd["swaps_done"].min()),
              vertical_noisy_mean=float(np.mean(vert)), vertical_noisy_sd=float(np.std(vert)),
              vertical_noisy_false_A_source=noise.source,
              interpretation="M3 near the vertical-noisy value and far above curveball = descent with noise")

    print(f"[{time.time() - t0:.0f}s] M3 done", flush=True)
    # ---------------- funnels: real + nulls (primary costs, same thresholds)
    fun = []
    ev_real = run(assembled)
    fr = events.funnel(ev_real, n_chars)
    for typ in ("gain", "loss"):
        fun.append(dict(input="real", seed=None, type=typ, n_sogs=n_chars, **fr[typ]))
    floor = a.noise_floor
    rc = sog_dir / "replicate_cosines.parquet"
    if floor is None and rc.exists():
        floor = events.noise_floor(pd.read_parquet(rc)["cosine"].to_numpy(),
                                   pr["part_b"]["specificity_change"]["noise_floor_quantile"]["value"])
    spec_in = RE.load_spec_inputs(sog_dir) if floor is not None else None
    if spec_in is not None:
        sp_, nt = RE.specificity(assembled, trees, st, spec_in, floor, pr, seeds[0], permutation=False)
        fun.append(RE.spec_funnel_row(sp_, nt, n_chars, input="real", seed=None))
    for sd in null_seeds:
        print(f"[{time.time() - t0:.0f}s] null seed {sd}", flush=True)
        # curveball-permuted P/A
        # masked checkerboard swaps, 100 x min(#1, #0) (same null as run_m3.py; curveball trades stall on the
        # sparse informative matrix)
        Y, done_sw, _ = simulate.checkerboard(X, n_swaps=int(spo * min(int((X == 1).sum()), int((X == 0).sum()))),
                                              seed=sd)
        print(f"  checkerboard seed {sd}: swaps {done_sw}, cells changed {int((Y != X).sum())}", flush=True)
        Mp = np.where(Y == X, M, Y).astype(np.int8)
        new = simulate.remap_members(codes_l, sp_l, sp_order, M, Mp)
        asm = rebuild(assembled, new)
        ev = run(asm)
        fr = events.funnel(ev, n_chars)
        for typ in ("gain", "loss"):
            fun.append(dict(input="null_curveball", seed=sd, type=typ, n_sogs=n_chars, **fr[typ]))
        # pseudo-clades
        rng = np.random.default_rng(sd)
        perm_idx = {}
        asm = []
        for og, sog_ids, codes, conf in assembled:
            c = codes.copy(); cf = conf.copy()
            for j, s in enumerate(sog_ids):
                pi = rng.permutation(codes.shape[0])
                c[:, j] = codes[pi, j]; cf[:, j] = conf[pi, j]
                perm_idx[s] = (og, pi)
            asm.append((og, sog_ids, c, cf))
        ev = run(asm)
        fr = events.funnel(ev, n_chars)
        for typ in ("gain", "loss"):
            fun.append(dict(input="null_pseudo_clade", seed=sd, type=typ, n_sogs=n_chars, **fr[typ]))
        if spec_in is not None:
            resd, fpd, covd = spec_in
            r2, f2, c2 = {}, {}, {}
            for s, (og, pi) in perm_idx.items():
                names = trees[og].leaf_names()
                mp = {names[k]: names[pi[k]] for k in range(len(names))}   # leaf k now carries member pi[k]
                for src, dst in ((resd, r2), (fpd, f2), (covd, c2)):
                    if s in src:
                        dst[s] = {acc: src[s][m] for acc, m in mp.items() if m in src[s]}
            sp_, nt = RE.specificity(asm, trees, st, (r2, f2, c2), floor, pr, sd, permutation=False)
            fun.append(RE.spec_funnel_row(sp_, nt, n_chars, input="null_pseudo_clade", seed=sd))
    # replicate (same-protein) null
    rep_dir = sog_dir / "replicate_states"
    if rep_dir.exists():
        import glob
        rep = pd.concat([pd.read_parquet(f) for f in sorted(glob.glob(f"{rep_dir}/bucket=*.parquet"))],
                        ignore_index=True)
        rep_asm = {og: (s, c) for og, s, c, _ in events.assemble(trees, rep, {})}
        rtrees, asm = {}, []
        for og, sog_ids, codes, conf in assembled:
            if og not in rep_asm or rep_asm[og][0] != sog_ids:
                continue
            tt = simulate.replicate_pair_tree(trees[og])
            rtrees[og] = tt
            asm.append((og, sog_ids, np.vstack([codes, rep_asm[og][1]]), np.vstack([conf, conf])))
        ev = run(asm, rtrees, None)
        n_rep_chars = sum(len(s) for _, s, _, _ in asm)
        fr = events.funnel(ev, n_rep_chars)
        stem = {rtrees[og].tree_id: (1, 1 + trees[og].n_nodes) for og in rtrees}
        for typ in ("gain", "loss"):
            false = sum(1 for e in ev if e["called"] and e["type"] == typ and e["node"] in stem[e["og_id"]])
            fun.append(dict(input="null_replicate", seed=None, type=typ, n_sogs=n_rep_chars, **fr[typ],
                            called_false_by_construction=false))
    else:
        for typ in ("gain", "loss"):
            fun.append(dict(input="null_replicate", seed=None, type=typ, n_sogs=0, note="not run: no replicate_states/"))
    fdf = pd.DataFrame(fun)
    files = []
    f = out / "funnel_real_vs_null.csv"; fdf.to_csv(f, index=False); files.append(f)
    f = out / "m3_curveball_null.csv"; nd.to_csv(f, index=False); files.append(f)
    f = out / "m3.json"; f.write_text(json.dumps(m3, indent=1)); files.append(f)
    man = dict(stage="partB_nulls", snapshot_id=a.snapshot, orthology=orth, stratum=a.stratum,
               prereg_sha256=recon.prereg_sha256(), git_commit=recon.git_commit(), seeds=seeds,
               runtime_s=time.time() - t0, n_in=cnt["sogs_in"], n_out=n_chars,
               n_dropped_by_reason={"og_not_in_orthology": cnt["og_not_in_orthology"]}, member_counters=cnt,
               thresholds_chosen=th["chosen"] if th else None, species_tree=str(nwk), orthoxml=str(a.orthoxml),
               code_sha256=recon.code_sha256(recon, events, simulate, __file__, RE.__file__))
    f = out / "manifest.json"; f.write_text(json.dumps(man, indent=1, default=str)); files.append(f)
    for f in files:
        recon.write_sidecar(f, a.snapshot, seeds, {"script": "scripts/partB/run_nulls.py", "argv": sys.argv[1:]})
    print(fdf.to_string(index=False))
    print(json.dumps(m3, indent=1))


if __name__ == "__main__":
    main()
