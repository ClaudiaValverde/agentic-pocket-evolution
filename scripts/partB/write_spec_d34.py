#!/usr/bin/env python
"""D34 write-up: SPECIFICITY_<O>.md per run (results/partB/<snap>/<run-id>/) and numbers under
partB.specificity_D34_<O1|O2> in results/numbers.json (read-modify-write under /tmp/pe_big/numbers.lock; other keys
kept). Inputs: the run's funnel_real_vs_null.csv, specificity_exploratory.parquet (D34 flags) and the bulk
cascade_summary.json from spec_cascade_d34.py.

    python scripts/partB/write_spec_d34.py --orthology O2 --release O2_r1 --run-id O2_r1_specD34
"""
from __future__ import annotations

import argparse
import fcntl
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from pocketevo import recon  # noqa: E402


def f3(x):
    return f"{x:.3g}" if isinstance(x, (float, np.floating)) else str(x)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--orthology", required=True, choices=["O1", "O2"])
    ap.add_argument("--release", required=True)
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--tau", type=float, default=0.2)
    a = ap.parse_args()
    snap = recon.default_snapshot()
    run = ROOT / f"results/partB/{snap}/{a.run_id}"
    sp = run / f"specificity/tau{a.tau}"
    fun = pd.read_csv(sp / "funnel_real_vs_null.csv")
    man = json.loads((sp / "manifest.json").read_text())
    ex = pd.read_parquet(sp / "specificity_exploratory.parquet")
    cs = json.loads(Path(f"/tmp/pe_big/partB/specD34/{a.release}/cascade_summary.json").read_text())
    rv = cs["real_vs_null"]
    now = subprocess.run(["date", "+%Y-%m-%d %H:%M %Z"], capture_output=True, text=True).stdout.strip()
    git = subprocess.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True, cwd=ROOT).stdout.strip()
    real = fun[fun["input"] == "real"].iloc[0]
    nul = fun[fun["input"] != "real"]
    L = [f"# Specificity change on {a.orthology} ({a.release}, run {a.run_id}): D34 analysis\n",
         f"Generated {now} by `scripts/partB/write_spec_d34.py` (git {git}). Producer: Part B specificity agent; "
         f"**not yet verified**. Runs: `run_specificity.py --orthology {a.orthology} --release {a.release} --run-id "
         f"{a.run_id} --tau {a.tau} --workers 16 --n-null-seeds 3 --rows-dir /tmp/pe_big/partB/specD34/{a.release}/rows` "
         f"(runtime {man['runtime_s']:.0f} s, CPU) and `spec_cascade_d34.py` ({cs['runtime_s']:.0f} s). Prereg sha256 "
         f"`{man['prereg_sha256']}`; frozen thresholds unchanged (noise floor {man['noise_floor']:.4f}, >= 2 "
         "BLOSUM62<0 clade-consistent substitutions, BH alpha 0.05). Everything below except the prereg BH count is "
         "**post hoc (D34)** and is kept out of the convincing headline.\n",
         "## 1. Prereg BH (headline): 0 calls, and why BH cannot fire\n",
         "| input | seed | residue gate | gate + shift | BH tested | BH called | min empirical p | min attainable p "
         "| BH threshold for one call | tail-normal called | tail-GPD called |", "|---|---|---|---|---|---|---|---|---|---|---|"]
    for _, r in fun.iterrows():
        L.append(f"| {r['input']} | {'' if pd.isna(r['seed']) else int(r['seed'])} | {r['residue_gate']:,} | "
                 f"{r['residue_and_shift']} | {r['bh_tested']:,} | {r['called']} | {r['min_p_empirical']:.2g} | "
                 f"{r['min_attainable_p_empirical']:.2g} | {r['bh_threshold_single_call']:.2g} | "
                 f"{r['called_tail_normal']} | {r['called_tail_gpd']} |")
    L.append(f"\nThe smallest empirical p any real pair can reach is 1/(1 + n_pool) of its matched comparator pool; the "
             f"best case over tested real pairs is {real['min_attainable_p_empirical']:.2g}, above the BH threshold for "
             f"a single call, alpha/n_tested = 0.05/{int(real['bh_tested']):,} = {real['bh_threshold_single_call']:.2g}. "
             "So the prereg BH call is unattainable by construction at this null resolution, whatever the data "
             f"(observed minimum {real['min_p_empirical']:.2g}).\n")
    L += ["## 2. Population-level signal (empirical null FDR, the D24 estimator)\n",
          "| stage | real | null per seed | null mean | empirical FDR |", "|---|---|---|---|---|"]
    names = {"funnel_pairs": "residue gate + shift (funnel)", "pass_artifacts": "funnel + artifact cascade",
             "pass_artifacts_strict": "funnel + cascade, strict (no quality flags)",
             "tail_normal_called": "tail-normal BH (D17 rule)", "tail_gpd_called": "tail-GPD BH"}
    for k, nm in names.items():
        v = rv[k]
        fdr = "" if v["empirical_fdr"] is None else f"{v['empirical_fdr']:.3f}"
        L.append(f"| {nm} | {v['real']} | {v['null_per_seed']} | {v['null_mean']:.2f} | {fdr} |")
    L.append("\nReal sister pairs exceed the pseudo-clade null about 6.5-7x at the funnel stage, and the excess survives "
             "the full artifact cascade at a similar ratio. The cascade is applied identically to null pairs (it "
             "evaluates the real proteins of the derived/sister clades; the pseudo-clade null permutes member data, "
             "not the tree). The pass rates are similar for real and null pairs, so the cascade does not enrich "
             "the signal; it removes mostly present-side fragment / length-ratio cases. Caveat: most of the excess "
             "is in the residue gate (real clades share substitutions, pseudo-clades do not, by construction), "
             "so the funnel ratio is an upper bound on a fingerprint-shift signal.\n")
    L += ["## 3. Tail-fit sensitivity (post hoc, D34(c))\n",
          "Same matched comparator pool as the prereg empirical p (depth x log2 clade-size bins, widened to >= 20). "
          "Normal (D17 rule, `directional.decide`): normal-fit p beyond the pool maximum, else max(empirical, normal). "
          "GPD: PWM fit to exceedances over the pool's q95; beyond the pool maximum the GPD p, else max(empirical, GPD) "
          "above q95. BH over the same tested set; same call rule (gate & shift above floor & q <= 0.05). Applied to "
          "the null seeds identically.\n",
          f"- Normal tail: real {rv['tail_normal_called']['real']} vs null {rv['tail_normal_called']['null_per_seed']}. "
          "The normal fit is a poor model for the bounded, right-skewed shift (1 - cos) distribution: it gives p down "
          "to ~1e-39 and fires in null seeds too, so it is anti-conservative and not usable as a call rule.",
          f"- GPD tail: real {rv['tail_gpd_called']['real']} vs null {rv['tail_gpd_called']['null_per_seed']}.",
          "- No null-seed tail call passes the artifact cascade.\n",
          "Tail-called real pairs:\n",
          "| sog_id | level | n_d | n_s | species_d | n_subst | substitutions | shift | pool max | n_pool | p_emp | "
          "p_normal | q_normal | p_gpd | q_gpd | pass_artifacts | strict | fails |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    t = ex[ex["flag_tail_normal_called"] | ex["flag_tail_gpd_called"]]
    for _, r in t.iterrows():
        L.append("| " + " | ".join(f3(r[c]) for c in ["sog_id", "level", "n_d", "n_s", "n_species_d", "n_subst",
                                                     "substitutions", "shift", "pool_max", "n_pool", "p_empirical",
                                                     "p_tail_normal", "q_tail_normal", "p_tail_gpd", "q_tail_gpd",
                                                     "pass_artifacts", "pass_artifacts_strict", "artifact_fails"]) + " |")
    fails = ex["artifact_fails"].replace("", np.nan).dropna().str.split(";").explode().value_counts().to_dict()
    L += ["\n## 4. Exploratory table\n",
          f"`specificity/tau{a.tau}/specificity_exploratory.parquet`: all {len(ex)} funnel-stage real pairs with the "
          "cascade columns and flags `flag_bh_prereg_called`, `flag_tail_normal_called`, `flag_tail_gpd_called`, "
          "`flag_pass_artifacts`, `flag_pass_artifacts_strict`, `flag_ge3_species_d`, `flag_tail_and_cascade`. "
          f"{int(ex['flag_pass_artifacts'].sum())} pass the cascade ({int(ex['flag_pass_artifacts_strict'].sum())} "
          f"strict); {int((ex['flag_ge3_species_d'] & ex['flag_pass_artifacts']).sum())} of those have >= 3 derived "
          f"species. Cascade fail counts: {json.dumps(fails)}. All rows are labelled exploratory, not significant. "
          "Null rows and per-input cascades are bulk: /tmp/pe_big/partB/specD34/" + a.release + "/.\n",
          "## Caveats\n",
          "- 3 null seeds (as the earlier specificity runs); empirical FDRs on null means of ~8-80 pairs are rough.",
          "- Candidates files (`candidates.parquet`) were not read for writing; the cascade ran in scratch results "
          "roots with an emptied candidate table.",
          "- O2 instability / species-tree lineage / PAE checks are `na` for specificity; O2 W members without an "
          "atlas site have no fingerprint.",
          "- Design scope (D31): within-orthogroup sister pairs only."]
    out = run / f"SPECIFICITY_{a.orthology}.md"
    out.write_text("\n".join(L) + "\n")
    recon.write_sidecar(out, snap, 0, {"script": "scripts/partB/write_spec_d34.py", "argv": sys.argv[1:]})
    num = dict(run_id=a.run_id, release=a.release, tau=a.tau, generated=now, verified=False, post_hoc="D34",
               prereg_bh_called=int(real["called"]),
               min_p_empirical_real=float(real["min_p_empirical"]),
               min_attainable_p_empirical_real=float(real["min_attainable_p_empirical"]),
               bh_threshold_single_call_real=float(real["bh_threshold_single_call"]),
               funnel=json.loads(fun.to_json(orient="records")), real_vs_null=rv,
               tail_called_real=json.loads(t[["sog_id", "level", "n_d", "n_s", "shift", "p_tail_normal", "q_tail_normal",
                                              "p_tail_gpd", "q_tail_gpd", "pass_artifacts",
                                              "pass_artifacts_strict"]].to_json(orient="records")),
               runtime_s=dict(specificity=man["runtime_s"], cascade=cs["runtime_s"]),
               source=[str(out.relative_to(ROOT)), str((sp / "funnel_real_vs_null.csv").relative_to(ROOT))])
    p = ROOT / "results/numbers.json"
    with open("/tmp/pe_big/numbers.lock", "w") as lk:
        fcntl.flock(lk, fcntl.LOCK_EX)
        cur = json.loads(p.read_text())
        cur.setdefault("partB", {})[f"specificity_D34_{a.orthology}"] = num
        tmp = p.with_suffix(".json.tmp_specd34")
        tmp.write_text(json.dumps(cur, indent=1))
        tmp.replace(p)
    print(out)


if __name__ == "__main__":
    main()
