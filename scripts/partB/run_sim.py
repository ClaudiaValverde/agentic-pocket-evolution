#!/usr/bin/env python
"""Part B simulation + FDR calibration of event calling (DESIGN §5.1, §8.3; prereg part_b.event_fdr).

On real O1 HOG trees (Ancestra tier-1 accurate OrthoXML; root HOGs with >= 3 species, seeded sample):
  world 'vertical' : SOG present at the root, no change (no-event null) + observation noise
  world 'gain'/'loss': one planted gain/loss on a branch eligible under the prereg rules + the same noise
For a noise grid (false-A / false-P multiplier) and the prereg cost grid (2:1 primary; 1:1, 3:1 at the reference
noise) every character is reconstructed (Sankoff), events extracted and scored with the frozen rules; a grid of
stricter calibrated cuts (min species in the derived clade, min parsimony margin) is evaluated.
FDR of a mixture world with planted fraction pi: FDR = false / all calls; thresholds chosen at the reference noise
= least strict cell with FDR <= 0.10 at pi_ref (max recall). Recall = planted events called on the exact branch.
Also: M3 (fraction of characters with <= 1 Dollo event) in the vertical-noisy worlds.

    PE_WORKERS=16 python scripts/partB/run_sim.py --n-trees 1500 --reps 20 --out results/partB/sim
"""
from __future__ import annotations

import argparse
import itertools
import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from pocketevo import events, recon, simulate  # noqa: E402

ORTHOXML = ROOT / "data/orthology/ancestra_tier1_v0/accurate/hogs.orthoxml.gz"
GUIDE = ROOT / "data/orthology/ancestra_tier1_v0/tier1_guide_tree.nwk"
NOISE_GRID = [0.0, 0.5, 1.0, 2.0, 4.0]
REF_NOISE = 1.0
U_ABSENT_GRID = [0.2, 0.8]
U_ABSENT_REF = 0.5                  # = NoiseModel default          # sensitivity of recall to the U rate of truly absent members
PI_GRID = [0.02, 0.05, 0.10, 0.20]
PI_REF = 0.05
K_GAIN = [3, 4, 5, 6, 8]
K_LOSS = [2, 3, 4, 5, 6, 8]
MARGINS = [1.0, 2.0]               # unambiguous events have margin >= 1 (integer costs)
LINEAGES = [2, 3, 4]               # min outgroup lineages supporting the ancestral state (A for gain, P for loss)
FIELDS = ["type", "node", "unambiguous", "margin", "n_species_P", "n_species_A", "frac_present", "frac_absent",
          "out_A_members", "out_A_lineages", "out_P_lineages", "candidate", "deep", "parent_kind", "n_species",
          "parent_is_root", "n_events_char", "n_clade_events_char", "f1_single_event"]

_G = {}


def _init(trees, rules, noise, seed):
    _G.update(trees=trees, rules=rules, noise=noise, seed=seed)


def _work(job):
    ti, reps = job
    t = _G["trees"][ti]
    rules, noise, seed = _G["rules"], _G["noise"], _G["seed"]
    tix = events.TreeIndex(t)
    sp = t.leaf_species()
    out, chars = [], []
    elig = {typ: simulate.eligible_branches(t, tix, typ, min_species=(3 if typ == "gain" else 2))
            for typ in ("gain", "loss")}
    settings = [(lam, (2.0, 1.0), None) for lam in NOISE_GRID] + \
               [(REF_NOISE, c, None) for c in [(1.0, 1.0), (3.0, 1.0)]] + \
               [(REF_NOISE, (2.0, 1.0), u) for u in U_ABSENT_GRID]
    for li, (lam, (g, l), u_abs) in enumerate(settings):
        nm = noise.scaled(lam)
        if u_abs is not None:
            nm.u_absent = u_abs
            nm.per_species = {sp_: dict(r_, u_absent=u_abs) for sp_, r_ in nm.per_species.items()}
        u_key = round(float(nm.u_absent), 4)
        for wi, world in enumerate(("vertical", "gain", "loss")):
            if world != "vertical" and len(elig[world]) == 0:
                continue
            # same truth + noise draw for every cost setting at a noise level (paired comparison)
            rng = np.random.default_rng([seed, ti, NOISE_GRID.index(lam), wi, int(round(u_key * 100))])
            truth, planted = simulate.simulate_truth(t, reps, world, rng, elig.get(world))
            rates = nm.rates(sp, rng)
            codes = simulate.observe(truth, rates, rng)
            res = recon.sankoff(t, codes, g, l)
            ev = events.extract_events(t, codes, res, rules, tix)
            key = dict(tree=ti, noise=lam, u_absent=u_key, gain_cost=g, loss_cost=l, world=world)
            n_called_type = {}
            for e in ev:
                e = events.apply_rules(e, rules)
                pl = int(planted[e["r"]])
                rec = {k: e[k] for k in FIELDS}
                rec.update(key, r=e["r"], planted=pl,
                           exact=bool(pl >= 0 and e["node"] == pl and e["type"] == world),
                           consistent=bool(pl >= 0 and e["type"] == world and
                                           t.leaf_lo[pl] <= t.leaf_lo[e["node"]] and
                                           t.leaf_hi[e["node"]] <= t.leaf_hi[pl]),
                           lenient=bool(pl >= 0 and e["type"] == world and
                                        set(e["clade_species"].split(";")) ==
                                        set(tix.sp_names[s] for s in np.unique(
                                            tix.leaf_sp[t.leaf_lo[pl]:t.leaf_hi[pl]]))))
                out.append(rec)
            # Dollo M3 on the gene tree (P/W -> 1, A -> 0, U -> missing)
            X = np.where(np.isin(codes, (0, 1)), 1, np.where(codes == 2, 0, -1)).astype(np.int8).T
            dev, _ = recon.dollo_events(t, X)
            chars.append(dict(key, n_chars=reps, n_dollo_defined=int((dev >= 0).sum()),
                              n_dollo_le1=int(((dev >= 0) & (dev <= 1)).sum()),
                              planted_nodes=planted.tolist() if world != "vertical" else []))
            _ = n_called_type
    return out, chars, {"tree": ti, "n_elig_gain": len(elig["gain"]), "n_elig_loss": len(elig["loss"])}


GKEY = ["noise", "u_absent", "gain_cost", "loss_cost"]


def evaluate_grid(ev: pd.DataFrame, chars: pd.DataFrame) -> pd.DataFrame:
    """Per (setting, type, k, margin): calls per character in each world, recall, FDR(pi).

    Truth: recall = planted event called on the exact branch; precision primary = 'consistent' (called clade is
    the planted clade or nested in it, same type: the claimed derived state and the change above it are true);
    strict precision = exact branch only (both reported)."""
    rows = []
    cand = ev[ev["candidate"]]
    nchar = chars.groupby(GKEY + ["world"])["n_chars"].sum()
    for key, sub in cand.groupby(GKEY):
        for typ, ks in (("gain", K_GAIN), ("loss", K_LOSS)):
            st = sub[sub["type"] == typ]
            nsp = st["n_species_P"] if typ == "gain" else st["n_species_A"]
            nlin = st["out_A_lineages"] if typ == "gain" else st["out_P_lineages"]
            for k, m, lin, f1 in itertools.product(ks, MARGINS, LINEAGES, (False, True)):
                if True:
                    sel = st[(nsp >= k) & (st["margin"] >= m - 1e-12) & (nlin >= lin) &
                             (st["f1_single_event"] | (not f1))]
                    rec = dict(zip(GKEY, key), type=typ, k=k, margin=m, lineages=lin, f1=f1)
                    for w in ("vertical", "gain", "loss"):
                        n = nchar.get(tuple(key) + (w,), 0)
                        sw = sel[sel["world"] == w]
                        rec[f"n_chars_{w}"] = int(n)
                        rec[f"calls_{w}"] = float(len(sw) / n) if n else np.nan
                        rec[f"false_{w}"] = float((~sw["consistent"]).sum() / n) if n else np.nan
                        rec[f"falsestrict_{w}"] = float((~sw["exact"]).sum() / n) if n else np.nan
                        rec[f"true_{w}"] = float(sw["exact"].sum() / n) if n else np.nan
                    rec["recall"] = rec[f"true_{typ}"]
                    rec["recall_consistent"] = float(
                        sel[(sel["world"] == typ) & sel["consistent"]][["tree", "r"]].drop_duplicates().shape[0]
                        / rec[f"n_chars_{typ}"]) if rec[f"n_chars_{typ}"] else np.nan
                    for pi in PI_GRID:
                        w0, w1 = 1 - pi, pi / 2
                        tot = w0 * rec["calls_vertical"] + w1 * (rec["calls_gain"] + rec["calls_loss"])
                        fc = w0 * rec["false_vertical"] + w1 * (rec["false_gain"] + rec["false_loss"])
                        fs = w0 * rec["falsestrict_vertical"] + w1 * (rec["falsestrict_gain"] + rec["falsestrict_loss"])
                        rec[f"fdr_pi{pi}"] = float(fc / tot) if tot > 0 else 0.0
                        rec[f"fdrstrict_pi{pi}"] = float(fs / tot) if tot > 0 else 0.0
                    rows.append(rec)
    return pd.DataFrame(rows)


def choose(grid: pd.DataFrame, fdr: float, pi: float, noise: float) -> dict:
    out = {}
    ref = grid[(grid["noise"] == noise) & (np.isclose(grid["u_absent"], U_ABSENT_REF)) & ~grid["f1"]]
    for (g, l), sub in ref.groupby(["gain_cost", "loss_cost"]):
        key = f"{int(g)}:{int(l)}"
        out[key] = {}
        for typ in ("gain", "loss"):
            s = sub[(sub["type"] == typ) & (sub[f"fdr_pi{pi}"] <= fdr) & (sub["recall"] > 0)]
            if s.empty:
                out[key][typ] = None
                continue
            best = s.sort_values(["recall", "k", "lineages", "margin"], ascending=[False, True, True, True]).iloc[0]
            out[key][typ] = {"min_species": int(best["k"]), "min_margin": float(best["margin"]),
                             "min_out_lineages": int(best["lineages"]),
                             "recall": float(best["recall"]), "fdr": float(best[f"fdr_pi{pi}"]),
                             "fdr_all_pi": {str(p): float(best[f"fdr_pi{p}"]) for p in PI_GRID},
                             "fdr_strict": float(best[f"fdrstrict_pi{pi}"]),
                             "recall_consistent": float(best["recall_consistent"]),
                             "false_calls_per_null_sog": float(best["false_vertical"])}
    return out


def plots(grid: pd.DataFrame, chosen: dict, outdir: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    files = []
    prim = grid[(grid["gain_cost"] == 2) & (grid["loss_cost"] == 1) & np.isclose(grid["u_absent"], U_ABSENT_REF) &
                ~grid["f1"]]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    for ax, typ in zip(axes, ("gain", "loss")):
        c = chosen["2:1"][typ]
        k0 = 3 if typ == "gain" else 2
        cells = [("frozen minima", k0, 1.0, 2, "--")]
        if c:
            cells.append(("calibrated", c["min_species"], c["min_margin"], c["min_out_lineages"], "-"))
        for lab, k, m, lin, ls in cells:
            s = prim[(prim["type"] == typ) & (prim["k"] == k) & (prim["margin"] == m) &
                     (prim["lineages"] == lin)].sort_values("noise")
            ax.plot(s["noise"], s["recall"], ls, color="#1f77b4", marker="o", label=f"recall ({lab})")
            ax.plot(s["noise"], 1 - s[f"fdr_pi{PI_REF}"], ls, color="#d62728", marker="s",
                    label=f"precision, pi={PI_REF} ({lab})")
        ax.set_xlabel("noise multiplier (false-A 0.05 x m, false-P 0.02 x m)")
        ax.set_title(f"planted {typ}{'es' if typ == 'loss' else 's'}, cost 2:1, u_absent {U_ABSENT_REF}")
        ax.set_ylim(0, 1.02)
        ax.legend(fontsize=7)
    f = outdir / "precision_recall_vs_noise.png"
    fig.savefig(f, dpi=130)
    plt.close(fig)
    files.append(f)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    ref = prim[prim["noise"] == REF_NOISE]
    cols = {2: "#1b9e77", 3: "#d95f02", 4: "#7570b3"}
    for ax, typ in zip(axes, ("gain", "loss")):
        s = ref[(ref["type"] == typ)]
        for (m, lin), q in s.groupby(["margin", "lineages"]):
            q = q.sort_values("k")
            ax.plot(q["recall"], 1 - q[f"fdr_pi{PI_REF}"], marker="o" if m == 1 else "x", color=cols[int(lin)],
                    ls="-" if m == 1 else ":", label=f"lineages >= {int(lin)}, margin >= {m:g}")
            for _, rr in q.iterrows():
                ax.annotate(str(int(rr["k"])), (rr["recall"], 1 - rr[f"fdr_pi{PI_REF}"]), fontsize=6)
        ax.axhline(0.9, color="grey", lw=0.8, ls=":")
        ax.set_xlabel("recall (exact branch)")
        ax.set_ylabel(f"precision (pi={PI_REF})")
        ax.set_title(f"{typ}: labels = min species k; reference noise")
        ax.legend(fontsize=6)
    f = outdir / "calibration_pr_curves.png"
    fig.savefig(f, dpi=130)
    plt.close(fig)
    files.append(f)
    return files


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--orthoxml", default=str(ORTHOXML))
    ap.add_argument("--n-trees", type=int, default=1500)
    ap.add_argument("--max-leaves", type=int, default=2000)
    ap.add_argument("--reps", type=int, default=20)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--heterogeneity", type=float, default=20.0)
    ap.add_argument("--noise-json", default=None, help="measured rates (results/observation/...); else defaults")
    ap.add_argument("--snapshot", default="20261030T085136Z")
    ap.add_argument("--out", default=str(ROOT / "results/partB/sim"))
    ap.add_argument("--workers", type=int, default=recon.n_workers())
    recon.add_partb_args(ap, skip=("snapshot",))
    a = ap.parse_args()
    t0 = time.time()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    pr = recon.load_prereg()
    rules = events.Rules.from_prereg(pr)
    noise = simulate.NoiseModel.from_json(a.noise_json) if a.noise_json else simulate.NoiseModel()
    noise.heterogeneity = a.heterogeneity
    global U_ABSENT_REF, U_ABSENT_GRID
    U_ABSENT_REF = round(float(noise.u_absent), 4)
    U_ABSENT_GRID = [u for u in (0.2, 0.5, 0.8) if abs(u - U_ABSENT_REF) > 0.05]
    min_sp = pr["sogs"]["min_species_per_root_og"]["value"]
    trees_all = [t for t in recon.iter_orthoxml_trees(a.orthoxml, min_species=min_sp)]
    n_elig_pool = len(trees_all)
    pool_idx = [i for i, t in enumerate(trees_all) if t.n_leaves <= a.max_leaves]
    rng = np.random.default_rng(a.seed)
    pick = sorted(rng.choice(pool_idx, size=min(a.n_trees, len(pool_idx)), replace=False).tolist())
    trees = [trees_all[i] for i in pick]
    del trees_all
    # warm numba caches in the parent
    t_ = trees[0]
    recon.sankoff(t_, np.full((t_.n_leaves, 1), 0, np.int8))
    recon.dollo_events(t_, np.ones((1, t_.n_leaves), np.int8))
    jobs = [(i, a.reps) for i in range(len(trees))]
    ev_rows, ch_rows, tr_rows = [], [], []
    with Pool(a.workers, initializer=_init, initargs=(trees, rules, noise, a.seed)) as pool:
        for evs, chs, tr in pool.imap_unordered(_work, jobs, chunksize=4):
            ev_rows.extend(evs); ch_rows.extend(chs); tr_rows.append(tr)
    ev = pd.DataFrame(ev_rows)
    ch = pd.DataFrame(ch_rows)
    trdf = pd.DataFrame(tr_rows).sort_values("tree")
    trdf["og_id"] = [trees[i].tree_id for i in trdf["tree"]]
    trdf["n_leaves"] = [trees[i].n_leaves for i in trdf["tree"]]
    trdf["n_species"] = [len(set(trees[i].leaf_species())) for i in trdf["tree"]]
    grid = evaluate_grid(ev, ch)
    fdr = pr["part_b"]["event_fdr"]["value"]
    chosen = choose(grid, fdr, PI_REF, REF_NOISE)
    # M3 in vertical-noisy worlds
    m3 = ch[ch["world"] == "vertical"].groupby(GKEY).agg(
        n=("n_dollo_defined", "sum"), le1=("n_dollo_le1", "sum")).reset_index()
    m3["frac_le1"] = m3["le1"] / m3["n"]
    # recall stratified by clade size and node kind at the reference setting with the chosen thresholds
    strat = []
    ref = ev[(ev["noise"] == REF_NOISE) & (ev["gain_cost"] == 2) & (ev["loss_cost"] == 1) &
             np.isclose(ev["u_absent"], U_ABSENT_REF)]
    plant_rows = []
    for _, r in ch[(ch["noise"] == REF_NOISE) & (ch["gain_cost"] == 2) & (ch["loss_cost"] == 1) &
                   np.isclose(ch["u_absent"], U_ABSENT_REF) & (ch["world"] != "vertical")].iterrows():
        t = trees[r["tree"]]
        tix_sp = None
        for rr, pl in enumerate(r["planted_nodes"]):
            plant_rows.append(dict(tree=r["tree"], world=r["world"], r=rr, planted=pl,
                                   parent_kind=recon.KIND_NAME[int(t.kind[t.parent[pl]])],
                                   clade_leaves=int(t.leaf_hi[pl] - t.leaf_lo[pl])))
        _ = tix_sp
    pdf = pd.DataFrame(plant_rows)
    for typ in ("gain", "loss"):
        c = chosen["2:1"][typ]
        if c is None:
            continue
        nsp = "n_species_P" if typ == "gain" else "n_species_A"
        nlin = "out_A_lineages" if typ == "gain" else "out_P_lineages"
        hit = ref[(ref["world"] == typ) & ref["exact"] & ref["candidate"] & (ref[nsp] >= c["min_species"]) &
                  (ref["margin"] >= c["min_margin"] - 1e-12) & (ref[nlin] >= c["min_out_lineages"])][
            ["tree", "r"]].drop_duplicates()
        hit["hit"] = True
        p = pdf[pdf["world"] == typ].merge(hit, on=["tree", "r"], how="left").fillna({"hit": False})
        p["clade_bin"] = pd.cut(p["clade_leaves"], [0, 3, 5, 10, 20, 50, 10**6],
                                labels=["<=3", "4-5", "6-10", "11-20", "21-50", ">50"])
        for (kind, b), s in p.groupby(["parent_kind", "clade_bin"], observed=True):
            strat.append(dict(type=typ, parent_kind=kind, clade_leaves=str(b), n_planted=len(s),
                              recall=float(s["hit"].mean())))
    strat = pd.DataFrame(strat)
    files = []
    for name, df in [("calibration_grid.csv", grid), ("m3_vertical_noisy.csv", m3), ("sim_trees.csv", trdf),
                     ("recall_by_clade.csv", strat)]:
        f = out / name
        df.to_csv(f, index=False)
        files.append(f)
    f = out / "sim_events.parquet"
    ev.to_parquet(f, index=False)
    files.append(f)
    files += plots(grid, chosen, out)
    psha = recon.prereg_sha256()
    thr = {"schema": "partb-thresholds-v1", "prereg_sha256": psha, "snapshot_id": a.snapshot,
           "event_fdr_target": fdr, "pi_ref": PI_REF, "reference_noise_multiplier": REF_NOISE,
           "noise_model": {k: v for k, v in noise.__dict__.items() if k != "per_species"},
           "noise_model_is_placeholder": a.noise_json is None,
           "frozen_rules": rules.__dict__, "chosen": chosen,
           "apply": "Rules.from_prereg(calibrated={'gain_min_species_cal': chosen['2:1']['gain']['min_species'], "
                    "'gain_min_margin': ..., 'loss_min_species_cal': ..., 'loss_min_margin': ...})",
           "real_data_rule": "on real SOGs, re-check FDR with events.empirical_fdr(null calls per SOG from a "
                             "simulation on the same trees with measured rates, n_real_sogs, n_real_calls)"}
    f = out / "thresholds.json"
    f.write_text(json.dumps(thr, indent=1, default=float))
    files.append(f)
    runtime = time.time() - t0
    man = dict(stage="partB_sim", n_root_hogs_ge3_species=n_elig_pool, n_trees=len(trees), reps=a.reps,
               noise_grid=NOISE_GRID, cost_settings=["2:1 all noise", "1:1 and 3:1 at ref noise"],
               n_characters=int(ch["n_chars"].sum()), n_events_stored=len(ev), runtime_s=runtime,
               orthoxml=str(a.orthoxml), orthoxml_sha256=recon.file_sha256(a.orthoxml),
               code_sha256=recon.code_sha256(recon, events, simulate, __file__))
    (out / "manifest.json").write_text(json.dumps(man, indent=1))
    for f in files + [out / "manifest.json"]:
        recon.write_sidecar(f, a.snapshot, a.seed, {"script": "scripts/partB/run_sim.py", "argv": sys.argv[1:]})
    print(json.dumps({"chosen": chosen, "runtime_s": runtime, "n_characters": man["n_characters"]}, indent=1))


if __name__ == "__main__":
    main()
