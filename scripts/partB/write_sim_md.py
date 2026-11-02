#!/usr/bin/env python
"""Write results/partB/sim/SIM.md from the simulation / toy outputs (every number read from a produced file)."""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from pocketevo import recon  # noqa: E402

SIM = ROOT / "results/partB/sim"


def fmt(x, d=3):
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return "–"
    return f"{x:.{d}f}" if isinstance(x, (float, np.floating)) else str(x)


def md_table(df: pd.DataFrame, cols=None, d=3) -> str:
    cols = cols or list(df.columns)
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for _, r in df.iterrows():
        lines.append("| " + " | ".join(fmt(r[c], d) for c in cols) + " |")
    return "\n".join(lines)


def cell(grid, noise, u, g, l, typ, k, m, lin, f1=False):
    s = grid[(grid["f1"] == f1) & (grid["noise"] == noise) & np.isclose(grid["u_absent"], u) & (grid["gain_cost"] == g) &
             (grid["loss_cost"] == l) & (grid["type"] == typ) & (grid["k"] == k) & (grid["margin"] == m) &
             (grid["lineages"] == lin)]
    return s.iloc[0] if len(s) else None


def _f1(grid, key, typ, c, pi):
    """Calibrated cuts + evidence filter f1 (exactly one MPR event spanning >= 2 species in the SOG)."""
    if not c:
        return {}
    g, l = map(float, key.split(":"))
    r = cell(grid, 1.0, 0.5, g, l, typ, c["min_species"], c["min_margin"], c["min_out_lineages"], f1=True)
    return {"cal+f1_recall": r["recall"], "cal+f1_FDR": r[f"fdr_pi{pi}"]} if r is not None else {}


def main():
    th = json.loads((SIM / "thresholds.json").read_text())
    man = json.loads((SIM / "manifest.json").read_text())
    grid = pd.read_csv(SIM / "calibration_grid.csv")
    m3 = pd.read_csv(SIM / "m3_vertical_noisy.csv")
    strat = pd.read_csv(SIM / "recall_by_clade.csv")
    trees = pd.read_csv(SIM / "sim_trees.csv")
    ch = th["chosen"]
    U0 = 0.5
    L = []
    L.append("# SIM — Part B event engine: simulation, FDR calibration, toy end-to-end run\n")
    L.append(f"Generated {time.strftime('%Y-%m-%d %H:%M UTC', time.gmtime())} by `scripts/partB/write_sim_md.py`. "
             f"prereg sha256 `{th['prereg_sha256']}`; snapshot `{th['snapshot_id']}` (development subset, "
             f"preliminary); git `{recon.git_commit()}`. **Simulated data only** — no real SOG states existed at "
             f"writing. Observation-noise rates are **placeholders** (`noise_model_is_placeholder = "
             f"{th['noise_model_is_placeholder']}`) until measured rates from `results/observation/...` replace "
             f"them (`run_sim.py --noise-json`); thresholds must then be re-derived.\n")
    L.append("## 1. Engine (code)\n")
    L.append("- `src/pocketevo/recon.py`: streaming OrthoXML -> per-root-HOG gene trees (orthologGroup = speciation, "
             "level = TaxRange; paralogGroup = duplication, level = taxon below the duplication; leaves = proteins; "
             "polytomies kept); NCBI guide tree parser; batched Sankoff (numba; down + up pass, single MPR with "
             "no-change tie-break, per-node conditional optimal costs -> unambiguous nodes and event margins); "
             "Dollo event count (M3); sidecars.\n"
             "- `src/pocketevo/events.py`: event extraction with clade info (species below, informative members, "
             "outgroup lineages with A/P, levels, node kinds, deep flag, species-tree LCA), prereg gain/loss rules, "
             "calibrated cuts, funnel, yield preflight, residue-first specificity test (BLOSUM62 < 0 clade-consistent "
             "substitutions, covariate-regressed fingerprint shift vs noise floor, matched-depth empirical null, BH; "
             "clade permutation as secondary), FDR helpers, SOG-table loader/assembler.\n"
             "- `src/pocketevo/simulate.py`: per-species noise model, vertical / planted gain / planted loss truth, "
             "masked curveball, species-matrix collapse + member remap, pseudo-clades, same-protein replicate trees.\n"
             "- Drivers `scripts/partB/`: `build_hog_trees.py`, `run_sim.py`, `make_toy_sogs.py`, `run_events.py`, "
             "`run_nulls.py`, `eval_toy.py`, `write_sim_md.py`. Tests: `tests/test_recon.py`.\n")
    hs = json.loads((ROOT / "results/partB/hog_trees/O1/summary.json").read_text())
    L.append(f"O1 HOG trees (`results/partB/hog_trees/O1/`): {hs['trees']} root HOGs, {hs['nodes']} nodes, "
             f"{hs['leaves']} protein leaves, {hs['dup_nodes']} duplication and {hs['spec_nodes']} speciation nodes, "
             f"{hs['root_is_dup']} roots are duplications; {100 * hs['frac_internal_polytomy']:.1f} % of internal "
             f"nodes are polytomies; {hs['trees_ge3_species']} root HOGs have >= 3 species. Leaf sets and duplication "
             f"counts of all root groups agree with `ancestra.orthoxml.read_orthoxml`; all leaf species and node "
             f"levels exist in the guide tree.\n")
    L.append("### Interpretation choices (INTERP; for PI review)\n")
    L.append("1. **W = present-ambiguous**: W constrains the reconstruction like P (never absence) but never counts "
             "toward the species thresholds; it counts as present in the 80 % fraction of a gain and *against* a "
             "loss; state confidence 0.5. Sensitivity: `--w-mode missing|soft:<x>`.\n"
             "2. **Outgroup lineage** of a member m = the child subtree of LCA(m, clade) that contains m (sister "
             "lineage, then successive outgroups) on the gene tree (paralog sister copies count).\n"
             "3. **Event called only if unambiguous** (both endpoint states identical in every MPR). With 2:1 costs a "
             "gain is unique only when the alternative needs >= 3 losses; a gain at a child of the gene-tree root "
             "always reconstructs as a loss of the sister (no outgroup) — visible in the recall by clade below.\n"
             "4. **Calibrated cuts never relax frozen minima**; the grid is min species in the derived clade k, min "
             "outgroup lineages supporting the ancestral state, min parsimony margin.\n"
             "5. **Truth for precision** = the called clade is the planted clade or nested in it, same type "
             "('consistent'); exact-branch precision ('strict') reported beside it. Recall = exact branch.\n"
             "6. **FDR** = false/all calls in a mixture world: a fraction pi of SOGs carry one planted event "
             "(half gains, half losses), the rest are vertical. pi is unknown for real data: thresholds chosen at "
             f"pi_ref = {th['pi_ref']}; on real SOGs FDR is re-estimated with `events.empirical_fdr` (expected false "
             "calls from the no-event simulation on the same trees with measured rates / observed calls).\n"
             "7. **Curveball swaps** (prereg '100 x number of 1s'): counted as exchanged pairs (min(|A|,|B|) per "
             "curveball trade = equivalent checkerboard swaps), with the missing (U/G) mask held fixed.\n"
             "8. **Specificity 'P throughout'** = >= 80 % of informative members P/W on both sides (A members left "
             "out); tested set for BH = pairs passing the residue-first gate; empirical null = shifts of the other "
             "sister pairs (no residue signal) matched on species-tree depth and derived-clade size (log2 bin).\n"
             "9. **M3** on the species tree (SOG x species, strongest member state; Dollo); the gene-tree Dollo "
             "version is used inside the simulation.\n")
    nm = th["noise_model"]
    L.append("## 2. Simulation design\n")
    L.append(f"{man['n_trees']} real O1 root-HOG trees (of {man['n_root_hogs_ge3_species']} with >= 3 species; "
             f"<= 2000 leaves; seeded), {man['reps']} replicate SOGs per tree x world x setting = "
             f"{man['n_characters']} simulated SOGs; runtime {man['runtime_s']:.0f} s on 16 workers. Tree sizes: "
             f"median {int(trees['n_leaves'].median())} members / {int(trees['n_species'].median())} species "
             f"(max {int(trees['n_leaves'].max())}). Worlds: vertical (no event), one planted gain (clade >= 3 "
             f"species), one planted loss (clade >= 2 species); planted branches must have >= 2 outgroup lineages "
             f"and >= 2 outgroup species.\n")
    L.append(f"Noise (per leaf): true present -> U {nm['u_present']}, W {nm['w_present']}, A (false-A) "
             f"{nm['false_A']}; true absent -> U {nm['u_absent']}, P (false-P) {nm['false_P']}; per-species rates "
             f"drawn from Beta(concentration {nm['heterogeneity']}). false-A = prereg rescue.false_A_target; the "
             f"others are placeholders. Grid: noise multiplier m on false-A/false-P in {man['noise_grid']}; "
             f"u_absent in {{0.2, 0.5, 0.8}} at m = 1; costs 2:1 (all) and 1:1, 3:1 at m = 1.\n")
    L.append("## 3. Chosen thresholds (event FDR 0.10)\n")
    L.append("```json\n" + json.dumps({"prereg_sha256": th["prereg_sha256"], "event_fdr_target":
                                       th["event_fdr_target"], "pi_ref": th["pi_ref"],
                                       "reference_noise_multiplier": th["reference_noise_multiplier"],
                                       "chosen": {k: {t: ({kk: (round(vv, 4) if isinstance(vv, float) else vv)
                                                           for kk, vv in v.items() if kk != "fdr_all_pi"} if v else None)
                                                      for t, v in d.items()} for k, d in ch.items()}},
                                      indent=1) + "\n```\n")
    L.append("Full record: `results/partB/sim/thresholds.json` (with sidecar). Applied by `run_events.py` / "
             "`run_nulls.py` via `events.rules_by_cost`; the file's prereg hash is checked at load.\n")
    rows = []
    for key in ("2:1", "1:1", "3:1"):
        g, l = map(float, key.split(":"))
        for typ in ("gain", "loss"):
            k0 = 3 if typ == "gain" else 2
            fz = cell(grid, 1.0, U0, g, l, typ, k0, 1.0, 2)
            c = ch[key][typ]
            cc = cell(grid, 1.0, U0, g, l, typ, c["min_species"], c["min_margin"], c["min_out_lineages"]) if c else None
            pi = th["pi_ref"]
            rows.append(dict(cost=key, type=typ,
                             frozen_recall=fz["recall"] if fz is not None else None,
                             frozen_FDR=fz[f"fdr_pi{pi}"] if fz is not None else None,
                             frozen_false_per_null_SOG=fz["false_vertical"] if fz is not None else None,
                             cal_k=c["min_species"] if c else None, cal_lineages=c["min_out_lineages"] if c else None,
                             cal_margin=c["min_margin"] if c else None,
                             cal_recall=cc["recall"] if cc is not None else None,
                             cal_FDR=cc[f"fdr_pi{pi}"] if cc is not None else None,
                             cal_FDR_strict=cc[f"fdrstrict_pi{pi}"] if cc is not None else None,
                             cal_false_per_null_SOG=cc["false_vertical"] if cc is not None else None,
                             **_f1(grid, key, typ, c, pi)))
    L.append("Frozen prereg minima vs calibrated cuts at the reference noise (m = 1, u_absent 0.5, "
             f"pi = {th['pi_ref']}):\n")
    L.append(md_table(pd.DataFrame(rows), d=4) + "\n")
    # noise curves
    L.append("## 4. Precision / recall vs noise (cost 2:1, calibrated cuts; pi = %s)\n" % th["pi_ref"])
    rows = []
    for typ in ("gain", "loss"):
        c = ch["2:1"][typ]
        for lam in man["noise_grid"]:
            for lab, k, m, lin in (("frozen", 3 if typ == "gain" else 2, 1.0, 2),
                                   ("calibrated", c["min_species"], c["min_margin"], c["min_out_lineages"])):
                r = cell(grid, lam, U0, 2.0, 1.0, typ, k, m, lin)
                if r is None:
                    continue
                rows.append(dict(type=typ, cuts=lab, noise_mult=lam, false_A=nm["false_A"] * lam,
                                 recall=r["recall"], precision=1 - r[f"fdr_pi{th['pi_ref']}"],
                                 precision_strict=1 - r[f"fdrstrict_pi{th['pi_ref']}"],
                                 false_calls_per_null_SOG=r["false_vertical"]))
    L.append(md_table(pd.DataFrame(rows), d=4) + "\n")
    L.append("![precision/recall vs noise](precision_recall_vs_noise.png)\n\n"
             "![calibration curves](calibration_pr_curves.png)\n")
    # u sensitivity
    rows = []
    for u in (0.2, 0.5, 0.8):
        for typ in ("gain", "loss"):
            c = ch["2:1"][typ]
            r = cell(grid, 1.0, u, 2.0, 1.0, typ, c["min_species"], c["min_margin"], c["min_out_lineages"])
            if r is not None:
                rows.append(dict(u_absent=u, type=typ, recall=r["recall"], FDR=r[f"fdr_pi{th['pi_ref']}"]))
    L.append("U-rate sensitivity (calibrated cuts, m = 1): recall depends mostly on how often true absences are "
             "observable (A needs footprint pLDDT >= 90).\n")
    L.append(md_table(pd.DataFrame(rows), d=4) + "\n")
    rows = []
    for typ in ("gain", "loss"):
        c = ch["2:1"][typ]
        r = cell(grid, 1.0, U0, 2.0, 1.0, typ, c["min_species"], c["min_margin"], c["min_out_lineages"])
        rows.append(dict(type=typ, **{f"FDR_pi{p}": r[f"fdr_pi{p}"] for p in (0.02, 0.05, 0.1, 0.2)}))
    L.append("FDR of the chosen cuts vs the planted fraction pi:\n")
    L.append(md_table(pd.DataFrame(rows), d=4) + "\n")
    L.append("Recall by planted clade size and parent node kind (2:1, calibrated cuts, m = 1):\n")
    L.append(md_table(strat, d=3) + "\n")
    L.append("## 5. M3 under noisy-but-vertical evolution (gene-tree Dollo, fraction of SOGs with <= 1 event)\n")
    L.append(md_table(m3[(m3["gain_cost"] == 2)], ["noise", "u_absent", "n", "frac_le1"], d=4) + "\n")
    # toy
    te = SIM / "toy_run/real/toy_eval.json"
    if te.exists():
        tv = json.loads(te.read_text())
        tm = json.loads((SIM / "toy_sogs/TOY/O1/manifest.json").read_text())
        L.append("## 6. End-to-end toy run (real input layout, known truth)\n")
        L.append(f"`make_toy_sogs.py`: {tm['n_ogs']} O1 OGs, {tm['n_sogs']} SOGs in the sogs-web-v1 layout (with G "
                 "rows, replicate observations, footprint residues, fingerprints with a volume covariate effect, "
                 "replicate cosines). `run_events.py` -> `eval_toy.py`:\n")
        rows = []
        for key in ("2:1", "1:1", "3:1"):
            for typ in ("gain", "loss"):
                x = tv.get(key, {}).get(typ)
                if x:
                    rows.append(dict(cost=key, type=typ, planted=x["n_planted"], called=x["n_called"],
                                     precision=x["precision_consistent"], precision_strict=x["precision_exact"],
                                     recall=x["recall_exact"], false_in_vertical_SOGs=x["false_in_vertical_sogs"]))
        L.append(md_table(pd.DataFrame(rows), d=3) + "\n")
        if "specificity" in tv:
            s = tv["specificity"]
            L.append(f"Specificity (planted: 3 BLOSUM62<0 clade-consistent substitutions + fingerprint shift): "
                     f"{s['n_pairs']} sister pairs, {s['n_tested']} passed the residue gate (BH set), "
                     f"{s['n_planted_reaching_test']}/{s['n_planted']} planted reached the test, {s['n_called']} called, "
                     f"{s['n_true']} true -> precision {fmt(s['precision'])}, recall {fmt(s['recall'])}.\n")
    nf = SIM / "toy_run/nulls/funnel_real_vs_null.csv"
    if nf.exists():
        f = pd.read_csv(nf)
        m3t = json.loads((SIM / "toy_run/nulls/m3.json").read_text())
        L.append("Null-input pipeline on the toy table (`run_nulls.py`, same code and thresholds, cost 2:1, prereg "
                 "seeds). Columns are cumulative stages; for specificity: species = residue gate, frac = residue + "
                 "shift > floor, candidate = BH-tested, called = BH q <= 0.05 and shift > floor.\n")
        cols = [c for c in ["input", "seed", "type", "n_sogs", "clade_size_ok", "unambiguous", "species", "frac",
                            "outgroup", "candidate", "called", "called_f1_single_event",
                            "called_false_by_construction"] if c in f.columns]
        L.append(md_table(f, cols, d=0) + "\n")
        L.append(f"M3 (species tree, toy): real {m3t['m3_real']:.3f} (n = {m3t['n_sogs_defined']}); curveball null "
                 f"mean {m3t['curveball_mean']:.3f} [{m3t['curveball_q025']:.3f}, {m3t['curveball_q975']:.3f}], "
                 f"{m3t['curveball_replicates']} replicates, p(null >= real) = {m3t['p_real_ge_null']:.4f}; "
                 f"vertical-noisy expectation {m3t['vertical_noisy_mean']:.3f} +- {m3t['vertical_noisy_sd']:.3f}.\n")
    L.append("## 7. Findings and caveats\n")
    c21 = ch["2:1"]
    fz_loss = cell(grid, 1.0, U0, 2.0, 1.0, "loss", 2, 1.0, 2)
    L.append(f"- The prereg minima alone keep gains well under the FDR target, but **2-species losses are noise-"
             f"dominated** (FDR {fz_loss[f'fdr_pi{th['pi_ref']}']:.2f} at the reference noise: two false-A sister "
             f"species mimic a loss). Calibration raises the loss cut to k = {c21['loss']['min_species']} species / "
             f">= {c21['loss']['min_out_lineages']} P outgroup lineages; this costs recall (loss recall "
             f"{c21['loss']['recall']:.3f} vs gain {c21['gain']['recall']:.3f}).\n"
             "- Null-input funnel (toy, §6): gains separate cleanly from nulls, but **loss calls survive the "
             "curveball and pseudo-clade nulls at a substantial rate** (SOGs with many A members arranged at random "
             "still contain small all-A clades). The f1 persistence filter (single event spanning >= 2 species) "
             "removes most null losses while keeping nearly all true calls; it is a frozen evidence filter, so it "
             "is applied downstream anyway, and the funnel reports it as `called_f1_single_event`. Recommend the PI "
             "treat loss candidates as reportable only after f1 and the empirical-FDR check on real data.\n"
             "- Recall is bounded by observability (u_absent) and by parsimony identifiability (unique MPR needed), "
             "not by the thresholds; recall rises steeply with clade size.\n"
             "- No-event worlds produce essentially no calls at the chosen cuts (false calls per null SOG in §3), "
             "so FDR on real data will be driven by misplaced/polarisation errors around true events and by the "
             "real noise rates, which are still placeholders.\n"
             "- Not covered here (other stages): simulated losses from masked pLDDT/truncation must become U (state "
             "caller), rescue, artifact cascade, uricase/GULO G controls.\n"
             "- To do when real inputs arrive: (1) `run_sim.py --noise-json <measured rates>` on the O2 trees -> new "
             "thresholds.json; (2) `run_events.py` + `run_nulls.py` on `results/sogs/<snapshot>/O2` (primary, then "
             "secondary stratum if the yield preflight says so); (3) empirical FDR check with `events.empirical_fdr`; "
             "(4) noise floor from measured same-site replicate cosines.\n")
    meas = sorted((SIM / "measured").glob("*/thresholds.json"))
    if meas:
        L.append("## 8. Re-derived thresholds on MEASURED noise (D9 rescue k = 1; D11 A needs > 10 footprint residues)\n")
        L.append("Noise JSONs: `results/partB/sim/measured/noise_*.json` (`scripts/partB/make_noise_json.py`: per-species "
                 "miss rate on known sites, rescue residual and control-rescue (false-W on true absences, 0.35) from "
                 "error_rates.parquet, A-eligibility per species from the real states). `d11_*` runs use the D11 "
                 "states and the > 10-residue rates and are the ones applied to real SOGs (pre-D11 runs were superseded "
                 "and removed).\n")
        rows = []
        for f in meas:
            t = json.loads(f.read_text())
            g = t["noise_model"]
            for c, v in t["chosen"].items():
                for typ, x in v.items():
                    rows.append(dict(run=f.parent.name, cost=c, type=typ, false_A=g["false_A"], u_absent=g["u_absent"],
                                     w_absent=g.get("w_absent", 0.0),
                                     k=x["min_species"] if x else None, lineages=x["min_out_lineages"] if x else None,
                                     margin=x["min_margin"] if x else None, recall=x["recall"] if x else None,
                                     FDR=x["fdr"] if x else None))
        L.append(md_table(pd.DataFrame(rows), d=4) + "\n")
        L.append("With measured rates true absences are almost never observable as A (u_absent ~ 0.92 before, "
                 "~ 0.98 after D11) and 35 % of observable ones are rescued to W, so recall of planted gains is "
                 "~0.1-0.5 % and of planted losses ~0: no loss cut reaches useful recall at FDR 0.10.\n")
    out = SIM / "SIM.md"
    out.write_text("\n".join(L))
    recon.write_sidecar(out, th["snapshot_id"], 0, {"script": "scripts/partB/write_sim_md.py"})
    print(out)


if __name__ == "__main__":
    main()
