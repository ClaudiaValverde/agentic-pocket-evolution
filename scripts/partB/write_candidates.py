#!/usr/bin/env python
"""Assemble the post-rescue Part B candidate table and CANDIDATES.md (real O1, snapshot 20261030T085136Z).

Inputs (all produced by scripts/partB/*): pre-rescue events (run_events.py on D11 states), post-rescue events
(run_events.py on <events>/rescue/sog_states), rescue_summary.json, null funnel (run_nulls.py), M3 (run_m3.py),
thresholds (run_sim.py on measured noise), D11 counts (apply_d11.py).
Output: results/partB/<snap>/O1/candidates.parquet (+ CANDIDATES.md, funnel_all.csv, sidecars).
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from pocketevo import recon  # noqa: E402

SNAP = recon.cli_value("--snapshot") or recon.default_snapshot()   # --snapshot / $PE_SNAPSHOT / prereg


def md(df, cols=None):
    cols = cols or list(df.columns)
    out = ["| " + " | ".join(map(str, cols)) + " |", "|" + "---|" * len(cols)]
    for _, r in df.iterrows():
        out.append("| " + " | ".join("–" if (isinstance(r[c], float) and np.isnan(r[c])) else
                                     (f"{r[c]:.3f}" if isinstance(r[c], float) else str(r[c])) for c in cols) + " |")
    return "\n".join(out)


def key(df):
    return list(zip(df["og_id"], df["type"], df["clade_species"]))


def main():
    ap = argparse.ArgumentParser()
    recon.add_partb_args(ap, skip=("sog_run",))
    ap.add_argument("--release", "--sog-run", dest="release", default=None, help="e.g. O1_r3: sets --base/--nulls/--m3/--out defaults")
    ap.add_argument("--base", default=str(ROOT / f"results/partB/events/{SNAP}/O1/d11"))
    ap.add_argument("--thresholds", default=str(ROOT / "results/partB/sim/measured/d11_primary_known_gt10/thresholds.json"))
    ap.add_argument("--thresholds-sens", default=str(ROOT / "results/partB/sim/measured/d11_primary_redetect_gt10/thresholds.json"))
    ap.add_argument("--nulls", default=str(ROOT / f"results/partB/nulls/{SNAP}/O1/d11/tau0.2/primary/funnel_real_vs_null.csv"))
    ap.add_argument("--m3", default=str(ROOT / f"results/partB/m3/{SNAP}/O1/tau0.2/m3.json"))
    ap.add_argument("--out", default=str(ROOT / f"results/partB/{SNAP}/O1"))
    ap.add_argument("--results-root", default=str(ROOT / "results/partB"), help="Part B results root (tests: tmp)")
    a = ap.parse_args()
    RB = Path(a.results_root)
    if a.release:
        a.base = str(RB / f"events/{SNAP}/{a.release}")
        a.nulls = str(RB / f"nulls/{SNAP}/{a.release}/tau0.2/primary/funnel_real_vs_null.csv")
        a.m3 = str(RB / f"m3/{SNAP}/{a.release}/tau0.2/m3.json")
        a.out = str(RB / f"{SNAP}/{a.release}")
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    base = Path(a.base)
    th = json.loads(Path(a.thresholds).read_text())
    ths = json.loads(Path(a.thresholds_sens).read_text()) if Path(a.thresholds_sens).exists() else None
    fun_rows, cand_all = [], []
    rescue_s = {}
    for tau in ("0.2", "0.3"):
        pre = base / f"tau{tau}/primary"
        post = pre / "rescue/events"
        if not (pre / "funnel.csv").exists():
            continue
        f = pd.read_csv(pre / "funnel.csv"); f["tau"] = tau; f["phase"] = "pre-rescue"; fun_rows.append(f)
        man = json.loads((pre / "manifest.json").read_text())
        if (post / "funnel.csv").exists():
            f = pd.read_csv(post / "funnel.csv"); f["tau"] = tau; f["phase"] = "post-rescue"; fun_rows.append(f)
            ev = pd.read_parquet(post / "events.parquet")
            ev = ev[ev["called"]].copy()
            ev["tau"] = float(tau)
            cand_all.append(ev)
            rescue_s[tau] = json.loads((pre / "rescue/rescue_summary.json").read_text())
        rescue_s.setdefault(tau, None)
        rescue_s[f"{tau}_dropped"] = man.get("n_dropped_by_reason")
    fun = pd.concat(fun_rows, ignore_index=True)
    fun.to_csv(out / "funnel_all.csv", index=False)
    cand = pd.concat(cand_all, ignore_index=True) if cand_all else pd.DataFrame()
    if len(cand):
        # sensitivity flags: same (og, type, clade species) called under other costs / tau
        prim = cand[(cand["gain_cost"] == 2) & (cand["loss_cost"] == 1)]
        for (g, l) in ((1.0, 1.0), (3.0, 1.0)):
            for tau in (0.2, 0.3):
                s = set(key(cand[(cand["gain_cost"] == g) & (cand["loss_cost"] == l) & (cand["tau"] == tau)]))
                cand[f"also_{int(g)}to1_tau{tau}"] = [k in s for k in key(cand)]
        for tau in (0.2, 0.3):
            s = set(key(prim[prim["tau"] == tau]))
            cand[f"also_2to1_tau{tau}"] = [k in s for k in key(cand)]
        cand["primary"] = (cand["gain_cost"] == 2) & (cand["loss_cost"] == 1) & (cand["tau"] == 0.2)
        # D21: a type is FDR-supported only if its calibrated cut reaches the 10 % FDR with non-zero recall in the
        # measured-noise simulation; losses are audit-level candidates ("not FDR-calibrated")
        sup = {}
        for c, v in th["chosen"].items():
            for typ, x in v.items():
                sup[(c, typ)] = bool(x and x["recall"] >= 1e-4 and x["fdr"] <= th["event_fdr_target"])
        cand["fdr_supported"] = [sup.get((f"{int(g)}:{int(l)}", t), False) for g, l, t in
                                 zip(cand["gain_cost"], cand["loss_cost"], cand["type"])]
        cand["label"] = np.where(cand["fdr_supported"], "FDR-calibrated candidate", "not FDR-calibrated (audit)")
        keepcols = ["primary", "fdr_supported", "label", "tau", "gain_cost", "loss_cost", "og_id", "char", "type", "node", "node_name", "level",
                    "parent_level", "parent_kind", "clade_lca_species_tree", "n_species", "n_species_P",
                    "n_species_A", "n_P", "n_W", "n_A", "frac_present", "frac_absent", "out_A_members",
                    "out_A_lineages", "out_P_lineages", "margin", "n_clade_events_char", "f1_single_event",
                    "min_conf_clade", "score_S_partial", "deep", "footprint_size", "clade_species",
                    "out_A_species", "out_P_species"] + [c for c in cand.columns if c.startswith("also_")]
        cand = cand[[c for c in keepcols if c in cand.columns]].rename(columns={"char": "sog_id"})
        cand = cand.sort_values(["primary", "type", "score_S_partial"], ascending=[False, True, False])
    nul = pd.read_csv(a.nulls) if Path(a.nulls).exists() else None
    nullcmp = None
    if nul is not None and len(cand):
        # real vs null at the same cut (pre-rescue D11 states, tau 0.2, 2:1): calls per input
        rows = []
        for typ in ("gain", "loss"):
            r = nul[(nul["type"] == typ)]
            real = r[r["input"] == "real"]
            for inp in ("null_curveball", "null_pseudo_clade"):
                q = r[r["input"] == inp]
                rows.append(dict(type=typ, null=inp, real_called=int(real["called"].iloc[0]) if len(real) else None,
                                 null_called_mean=float(q["called"].mean()) if len(q) else None,
                                 null_called_max=int(q["called"].max()) if len(q) else None,
                                 real_f1=int(real["called_f1_single_event"].iloc[0]) if len(real) else None,
                                 null_f1_mean=float(q["called_f1_single_event"].mean()) if len(q) else None,
                                 n_null_seeds=int(len(q))))
                rows[-1]["empirical_null_FDR"] = (min(1.0, rows[-1]["null_called_mean"] / rows[-1]["real_called"])
                                                  if rows[-1]["real_called"] else None)
                sim = th["chosen"].get("2:1", {}).get(typ)
                # no simulated calls (recall ~0) -> simulated FDR undefined, not 0
                rows[-1]["simulated_FDR"] = sim["fdr"] if sim and sim["recall"] >= 1e-4 else None
        nullcmp = pd.DataFrame(rows)
        for _, r in nullcmp.iterrows():
            m = cand["type"] == r["type"]
            cand.loc[m, f"{r['null']}_called_mean_same_cut"] = r["null_called_mean"]
            cand.loc[m, "real_called_same_cut_prerescue"] = r["real_called"]
        nullcmp.to_csv(out / "real_vs_null_same_cut.csv", index=False)
    # false-loss negative controls: losses called pre-rescue (tau 0.2, 2:1) not called after rescue
    ctrl = pd.DataFrame()
    pre = base / "tau0.2/primary"
    if a.release and (out / "false_loss_controls.parquet").exists():
        ctrl = pd.read_parquet(out / "false_loss_controls.parquet")      # scripts/partB/false_loss_controls.py
    elif (pre / "rescue/events/events.parquet").exists():
        e0 = pd.read_parquet(pre / "events.parquet")
        e1 = pd.read_parquet(pre / "rescue/events/events.parquet")
        sel = lambda e: e[e["called"] & (e["type"] == "loss") & (e["gain_cost"] == 2)]  # noqa: E731
        post_keys = set(zip(sel(e1)["char"], sel(e1)["node"]))
        lost = sel(e0)[[k not in post_keys for k in zip(sel(e0)["char"], sel(e0)["node"])]].copy()
        rh = pd.read_parquet(pre / "rescue/rescue_hits.parquet")
        rs = pd.read_parquet(pre / "rescue/sog_states/bucket=all.parquet")
        trees = {t.tree_id: t for t in recon.iter_orthoxml_trees(
            ROOT / "data/orthology/ancestra_tier1_v0/accurate/hogs.orthoxml.gz", keep=set(lost["og_id"]))}
        st0 = {}
        rows = []
        post_state = dict(zip(zip(rs["sog_id"], rs["accession"]), rs["state"]))
        nh = dict(zip(zip(rh["sog_id"], rh["accession"]), rh["n_settings_hit"]))
        for _, r in lost.iterrows():
            t = trees[r["og_id"]]
            leaves = t.leaf_names()[t.leaf_lo[r["node"]]:t.leaf_hi[r["node"]]]
            after = [post_state.get((r["char"], x)) for x in leaves]
            hits = [nh.get((r["char"], x)) for x in leaves]
            rows.append(dict(sog_id=r["char"], og_id=r["og_id"], node=r["node"], level=r["level"],
                             clade_species=r["clade_species"], n_A_before=r["n_A"],
                             n_A_after=sum(x == "A" for x in after), n_W_after=sum(x == "W" for x in after),
                             max_settings_hit=max([h for h in hits if h is not None], default=0),
                             reason="rescue: clade A calls re-detected by alternative BURROW settings (A->W)"))
        ctrl = pd.DataFrame(rows)
        if len(ctrl):
            ctrl["frac_A_rescued"] = ctrl["n_W_after"] / ctrl["n_A_before"].clip(lower=1)
            ctrl = ctrl.sort_values(["frac_A_rescued", "n_A_before"], ascending=False)
            ctrl["negative_control_example"] = False
            ctrl.iloc[:5, ctrl.columns.get_loc("negative_control_example")] = True
        if not a.release:
            ctrl.to_parquet(out / "false_loss_controls.parquet", index=False)
    f = out / "candidates.parquet"
    cand.to_parquet(f, index=False)
    m3 = json.loads(Path(a.m3).read_text()) if Path(a.m3).exists() else None
    d11 = {t: json.loads((out / f"d11_t{t}.json").read_text()) for t in ("0.2", "0.3") if (out / f"d11_t{t}.json").exists()}
    L = [f"# Part B candidates — O1, snapshot {SNAP} (development subset, preliminary)\n",
         f"Generated {time.strftime('%Y-%m-%d %H:%M UTC', time.gmtime())} by `scripts/partB/write_candidates.py`; "
         f"prereg sha256 `{recon.prereg_sha256()}`; git `{recon.git_commit()}`. Producer: phylogenetics agent; "
         "**not yet verified**. Gain/loss only (specificity needs footprint residues + fingerprints per member, "
         "not yet exported by the SOG stage).\n",
         "## Pipeline\n",
         "1. States: SOG agent, τ = 0.2 (primary) and 0.3 (sensitivity), primary A (pLDDT ≥ 90). D8: SOGs with "
         "footprint ≤ 3 columns excluded. **D11**: A requires the member footprint to span > 10 residues, else U "
         + (f"(applied by the SOG stage in release `{a.release}`; small_footprint / empty_footprint SOGs excluded).\n"
            if a.release else
            f"(`scripts/partB/apply_d11.py`; A→U τ0.2: {d11.get('0.2', {}).get('A_to_U')}, τ0.3: "
            f"{d11.get('0.3', {}).get('A_to_U')}).\n") +
         "2. Event-eligible SOGs (≥ 2 P and ≥ 2 A members; necessary for any candidate) → Sankoff on the O1 HOG "
         "gene tree, costs 2:1 (primary), 1:1, 3:1 → prereg gain/loss rules + FDR-calibrated cuts.\n"
         "3. Rescue (k = 1, D9; hit = Jaccard ≥ 0.3, the rule k was calibrated with) of **every P and A member** of "
         "every SOG with a called candidate (A→W at ≥ 1/6 settings, P→W at 0/6), then events recomputed.\n",
         "## Thresholds (re-derived on measured noise; event FDR 0.10)\n",
         "Noise model (`results/partB/sim/measured/noise_primary_known_gt10.json`): per-species miss rate on known "
         "sites (> 10 residues, A_primary), rescue residual from error_rates.parquet (false-A after rescue k=1 "
         "= 0.054 for > 10 residues), **false-W on true absences 0.35** (control regions rescued), A-eligibility "
         "per species from the D11 states (pooled ≈ 1.4 % of non-P members are A → u_absent ≈ 0.98), false-P "
         "0.02 (placeholder, not measured).\n"]
    rows = []
    for nm, t_ in (("known-site miss (primary)", th), ("re-detection miss (sens.)", ths)):
        if not t_:
            continue
        for c, v in t_["chosen"].items():
            for typ, x in v.items():
                rows.append(dict(model=nm, cost=c, type=typ, min_species=x and x["min_species"],
                                 min_out_lineages=x and x["min_out_lineages"], min_margin=x and x["min_margin"],
                                 recall=x and x["recall"], FDR=x and x["fdr"]))
    L.append(md(pd.DataFrame(rows)) + "\n")
    L.append("Recall is per planted event on a random eligible branch of a real O1 tree; it is ≈ 0.1 % for gains "
             "and ≈ 0 for losses because true absences are almost never observable as A under the primary rule "
             "+ D11 (u_absent ≈ 0.98) and 35 % of the observable ones are rescued to W. The loss cuts that meet "
             "the FDR target call essentially nothing in simulation: **real loss calls are not FDR-supported** and "
             "are listed for audit only.\n")
    L.append("## Funnel (real)\n")
    cols = ["tau", "phase", "cost", "type", "n_sogs", "clade_size_ok", "unambiguous", "species", "frac", "outgroup",
            "candidate", "called", "called_f1_single_event"]
    L.append(md(fun[[c for c in cols if c in fun.columns]]) + "\n")
    L.append("Stage columns are cumulative: clade_size_ok = MPR event whose clade reaches the frozen species minimum; "
             "species/frac/outgroup = prereg gain/loss rules; candidate = all frozen rules; called = + calibrated "
             "cuts; called_f1 = + evidence filter 1 (single event spanning ≥ 2 species). Post-rescue rows cover "
             "only SOGs that had a called candidate before rescue.\n")
    for tau, r in rescue_s.items():
        if r and not str(tau).endswith("_dropped"):
            L.append(f"Rescue τ{tau}: {r['n_sogs']} SOGs, {r['n_proteins']} proteins, A {r['A_before']}→{r['A_after']} "
                     f"(A→W {r['A_to_W']}), P {r['P_before']}→{r['P_after']} (P→W {r['P_to_W']}), "
                     f"BURROW {r['burrow_s']:.0f} s.\n")
        elif str(tau).endswith("_dropped"):
            L.append(f"SOGs dropped τ{tau.split('_')[0]}: {r}\n")
    if nul is not None:
        L.append("## Real vs null inputs (pre-rescue D11 states, τ = 0.2, cost 2:1, same code and thresholds)\n")
        L.append(md(nul[[c for c in ["input", "seed", "type", "n_sogs", "clade_size_ok", "unambiguous", "candidate",
                                     "called", "called_f1_single_event", "called_false_by_construction", "note"]
                         if c in nul.columns]]) + "\n")
    if nullcmp is not None:
        L.append("### Real vs null at the same cut (D21, D24): empirical null FDR = mean null calls / real calls\n")
        L.append(md(nullcmp) + "\n")
        L.append("Each candidate row carries these counts (`null_*_called_mean_same_cut`, `real_called_same_cut_prerescue`)."
                 " Pre-rescue comparison: the null inputs are permutations of the same pre-rescue D11 states.\n")
    if len(ctrl) and "kind" in ctrl.columns:
        L.append("## False-loss negative controls (D24)\n")
        for kind, txt in (("bad_model", "losses called on states before the quality filters (O1_r2 without D11; "
                           "primary and pLDDT >= 70 strata) whose derived-clade A members are all poor models and "
                           "whose call disappears once A_primary + D11 apply"),
                          ("rescue_instability", "losses called before rescue on the release and not after rescue "
                           "(high-quality models; pocket-finder instability)")):
            c = ctrl[ctrl["kind"] == kind]
            L.append(f"**{kind}** ({len(c)} found; 5 marked): {txt}.\n")
            cols = [x for x in ["negative_control_example", "pre_stratum", "sog_id", "level", "n_A_clade", "frac_A_bad",
                                "min_A_footprint_plddt", "bad_reasons", "frac_A_rescued", "clade_species"]
                    if x in c.columns]
            L.append(md(c[c["negative_control_example"]], cols) + "\n")
    elif len(ctrl):
        L.append("## False-loss negative controls (losses overturned by rescue)\n")
        L.append(f"{len(ctrl)} losses called before rescue (τ 0.2, 2:1) are not called after rescue. The five with the "
                 "largest rescued fraction of clade A calls are the designated negative-control examples "
                 "(`false_loss_controls.parquet`, `negative_control_example`). The other artifact checks were not run.\n")
        L.append(md(ctrl.head(8), ["negative_control_example", "sog_id", "level", "n_A_before", "n_A_after",
                                   "n_W_after", "frac_A_rescued", "max_settings_hit", "clade_species"]) + "\n")
    if len(cand):
        p = cand[cand["primary"]]
        yp = int(len(p))
        L.append("## Yield preflight (pi_freeze.yield_preflight)\n")
        L.append(f"Primary stratum (pLDDT ≥ 90 A, D11, τ 0.2, 2:1, post-rescue): **{int((p['type'] == 'gain').sum())} "
                 f"gains (FDR-calibrated) + {int((p['type'] == 'loss').sum())} losses (not FDR-calibrated, D21) = {yp} "
                 f"called** "
                 f"({int(p['f1_single_event'].sum())} also pass f1). Threshold 5 → secondary stratum "
                 f"{'REQUIRED' if yp < 5 else 'not triggered'}.\n")
        L.append("## Primary candidates (top 25 by S_partial = persistence × min state confidence; evidence count "
                 "pending filters 2/3)\n")
        L.append(md(p.head(25), ["type", "label", "og_id", "sog_id", "level", "clade_lca_species_tree", "n_species_P",
                                 "n_species_A", "out_A_lineages", "out_P_lineages", "f1_single_event",
                                 "score_S_partial", "also_2to1_tau0.3", "also_1to1_tau0.2", "also_3to1_tau0.2"]) + "\n")
        st = {c: float(p[c].mean()) for c in p.columns if c.startswith("also_")}
        L.append(f"Robustness of primary calls (fraction also called): {json.dumps({k: round(v, 3) for k, v in st.items()})}\n")
    if m3:
        r = m3["results"]["multi_species"]
        L.append("## M3 tree consistency (D11 states, τ 0.2; `numbers.json` partA_descent.M3)\n")
        L.append(f"{r['n_sogs']} multi-species SOGs: M3 = {r['m3']:.3f}; checkerboard/curveball null "
                 f"{r['curveball_mean']:.3f} [{r['curveball_q025']:.3f}, {r['curveball_q975']:.3f}] "
                 f"({r['curveball_replicates']} reps), p = {r['p_null_ge_real']:.4f}; noisy-vertical "
                 f"{r['vertical_noisy_mean']:.3f} ± {r['vertical_noisy_sd']:.3f}. Swaps: {m3['results']['multi_species']['swaps_rule']}.\n")
    mix = []
    for mlt in (50, 100, 200):
        fm = Path(a.m3).with_name("m3.json" if mlt == 100 else f"m3_x{mlt}_multi_species.json")
        if fm.exists():
            r_ = json.loads(fm.read_text())["results"]["multi_species"]
            mix.append(dict(swap_mult=mlt, replicates=r_["curveball_replicates"], swaps=r_["swaps_target"],
                            null_mean=r_["curveball_mean"], null_q025=r_["curveball_q025"],
                            null_q975=r_["curveball_q975"], m3_real=r_["m3"]))
    if len(mix) > 1:
        L.append("### D22 swap-mixing check (null M3 at 50x / 100x / 200x min(#1, #0))\n")
        L.append(md(pd.DataFrame(mix)) + "\n")
        L.append("A stable null mean across multipliers means the swap chain has mixed.\n")
    L.append("## Caveats\n- Development subset; O1 only (46 species). Deep (Eukaryota-level) events are flagged `deep` "
             "(orthology-uncertain).\n- Artifact cascade (interface, TED/PAE, MSA-method flip, signal peptide "
             "already in states) not yet applied; required before any case is labelled.\n- false-P is a placeholder; "
             "loss calls carry no FDR support.\n")
    (out / "CANDIDATES.md").write_text("\n".join(L))
    for x in [f, out / "CANDIDATES.md", out / "funnel_all.csv"] + [y for y in (out / "real_vs_null_same_cut.csv",
                                                                              out / "false_loss_controls.parquet")
                                                                    if y.exists()]:
        recon.write_sidecar(x, SNAP, 0, {"script": "scripts/partB/write_candidates.py"})
    print(out / "CANDIDATES.md")


if __name__ == "__main__":
    main()
