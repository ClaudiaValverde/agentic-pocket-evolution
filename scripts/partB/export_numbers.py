#!/usr/bin/env python
"""Export Part B numbers to results/numbers.json under the key `partB` (other keys untouched). Every value is read
from a produced file; the source path is recorded next to each block."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from pocketevo import recon  # noqa: E402

SNAP = recon.cli_value("--snapshot") or recon.default_snapshot()


def rel(p):
    return str(Path(p).relative_to(ROOT))


def jl(p):
    return json.loads(Path(p).read_text())


def main():
    ap = argparse.ArgumentParser()
    recon.add_partb_args(ap, skip=("sog_run", "orthology"))
    ap.add_argument("--o1-run", default="O1_r3")
    ap.add_argument("--o2-run", default="O2_r1")
    a = ap.parse_args()
    R = ROOT / "results/partB"
    o1, o2 = a.o1_run, a.o2_run
    ev1 = R / f"events/{SNAP}/{o1}"
    nb = {"snapshot_id": SNAP, "prereg_sha256": recon.prereg_sha256(), "producer": "scripts/partB/export_numbers.py"}
    # ---- wall times per stage
    wt = {}
    for t in ("0.2", "0.3"):
        d = ev1 / f"tau{t}/primary"
        wt[f"O1_tau{t}"] = {"events_pre_s": jl(d / "manifest.json")["runtime_s"],
                            "rescue_s": jl(d / "rescue/rescue_summary.json")["runtime_s"],
                            "rescue_burrow_s": jl(d / "rescue/rescue_summary.json")["burrow_s"],
                            "events_post_s": jl(d / "rescue/events/manifest.json")["runtime_s"]}
    wt["O1_nulls_s"] = jl(R / f"nulls/{SNAP}/{o1}/tau0.2/primary/manifest.json")["runtime_s"]
    wt["O1_artifact_cascade_s"] = jl(R / f"{SNAP}/{o1}/artifact_summary.json")["runtime_s"]
    wt["O1_specificity_s"] = jl(R / f"{SNAP}/{o1}/specificity/tau0.2/manifest.json")["runtime_s"]
    wt["O2_events_s"] = jl(R / f"events/{SNAP}/{o2}/tau0.2/primary/manifest.json")["runtime_s"]
    wt["O2_nulls_s"] = jl(R / f"nulls/{SNAP}/{o2}/tau0.2/primary/manifest.json")["runtime_s"]
    wt["O2_artifact_cascade_s"] = jl(R / f"{SNAP}/{o2}/artifact_summary.json")["runtime_s"]
    wt["O2_rescue2_s"] = jl(ROOT / f"results/sogs/{SNAP}/{o2}_rescue2/rescue_summary.json").get("wall_s") \
        if (ROOT / f"results/sogs/{SNAP}/{o2}_rescue2/rescue_summary.json").exists() else None
    rb_side = jl(R / f"{SNAP}/{o1}/rbh_outgroups.parquet.json")
    wt["RBH_diamond_s"] = rb_side.get("runtime_s")
    nb["wall_times"] = wt
    # ---- RBH + Pfam architecture
    rbh = {}
    for run in (o1, o2):
        d = pd.read_parquet(R / f"{SNAP}/{run}/rbh_outgroups.parquet")
        d = d[d["testable"].astype(bool)]
        am = d["arch_match"].dropna().astype(bool)
        rbh[run] = {"n_outgroup_supporters": int(len(d)), "frac_rbh": float(d["rbh"].astype(bool).mean()),
                    "frac_rbh_by_type": {k: float(g["rbh"].astype(bool).mean()) for k, g in d.groupby("type")},
                    "n_arch_compared": int(len(am)), "frac_pfam_arch_match": float(am.mean()) if len(am) else None,
                    "source": rel(R / f"{SNAP}/{run}/rbh_outgroups.parquet")}
    nb["rbh_orthology"] = rbh
    # ---- duplication-associated gains
    dup = {}
    for run in (o1, o2):
        df = pd.read_parquet(R / f"{SNAP}/{run}/duplication_gains_exploratory.parquet")
        md = (R / f"{SNAP}/{run}/DUPLICATION_GAINS.md").read_text()
        import re
        m = re.search(r"orientations: (\d+) copies, (\d+) \(", md)
        dup[run] = {"n_tested": int(len(df)), "n_pass": int(df["duplication_associated_gain"].sum()),
                    "pattern_species_total": int(df["n_species_pattern"].sum()),
                    "swapped_species_total": int(df["n_species_swapped"].sum()),
                    "background_copies": int(m.group(1)) if m else None,
                    "background_ge3": int(m.group(2)) if m else None,
                    "label": "exploratory, not preregistered",
                    "source": rel(R / f"{SNAP}/{run}/duplication_gains_exploratory.parquet")}
    nb["duplication_gains"] = dup
    # ---- false-loss controls
    fl = pd.read_parquet(R / f"{SNAP}/{o1}/false_loss_controls.parquet")
    nb["false_loss_controls"] = {k: {"pool": int(len(g)), "examples": int(g["negative_control_example"].sum())}
                                 for k, g in fl.groupby("kind")}
    nb["false_loss_controls"]["source"] = rel(R / f"{SNAP}/{o1}/false_loss_controls.parquet")
    # ---- robustness of O1 primary gains
    c1 = pd.read_parquet(R / f"{SNAP}/{o1}/candidates.parquet")
    p = c1[c1["primary"].astype(bool) & (c1["type"] == "gain")]
    nb["O1_primary_gains"] = {"n": int(len(p)), "pass_artifacts": int(p["pass_artifacts"].sum()),
                              "robustness_frac_also_called": {c.replace("also_", ""): float(p[c].mean())
                                                              for c in p.columns if c.startswith("also_")},
                              "source": rel(R / f"{SNAP}/{o1}/candidates.parquet")}
    # ---- funnels + empirical null FDR (O1 and O2)
    th = jl(R / "sim/measured/d11_primary_known_gt10/thresholds.json")
    for key, run, evd in (("O1", o1, ev1 / "tau0.2/primary"), ("O2", o2, R / f"events/{SNAP}/{o2}/tau0.2/primary")):
        fun = pd.read_csv(evd / "funnel.csv")
        nul = pd.read_csv(R / f"nulls/{SNAP}/{run}/tau0.2/primary/funnel_real_vs_null.csv")
        blk = {"funnel_2to1": {r.type: {c: (None if pd.isna(getattr(r, c)) else int(getattr(r, c)))
                                        for c in ("n_sogs", "clade_size_ok", "unambiguous", "candidate", "called",
                                                  "called_f1_single_event")}
                               for r in fun[fun["cost"] == "2:1"].itertuples()},
               "empirical_null_FDR": {}, "simulated_FDR": {},
               "source": [rel(evd / "funnel.csv"), rel(R / f"nulls/{SNAP}/{run}/tau0.2/primary/funnel_real_vs_null.csv")]}
        for typ in ("gain", "loss"):
            real = float(nul[(nul["input"] == "real") & (nul["type"] == typ)]["called"].iloc[0])
            for inp in ("null_curveball", "null_pseudo_clade"):
                q = nul[(nul["input"] == inp) & (nul["type"] == typ)]
                blk["empirical_null_FDR"][f"{typ}_{inp}"] = {"null_mean": float(q["called"].mean()),
                                                             "n_seeds": int(len(q)), "real": real,
                                                             "fdr": min(1.0, float(q["called"].mean()) / real) if real else None}
            x = th["chosen"]["2:1"].get(typ)
            blk["simulated_FDR"][typ] = x["fdr"] if x and x["recall"] >= 1e-4 else None
        if key == "O2":
            c2 = pd.read_parquet(R / f"{SNAP}/{o2}/candidates.parquet")
            q2 = c2[c2["primary"].astype(bool)]
            blk["called_post_rescue"] = {t: int((q2["type"] == t).sum()) for t in ("gain", "loss")}
            blk["pass_artifacts"] = {t: int(q2.loc[q2["type"] == t, "pass_artifacts"].sum()) for t in ("gain", "loss")}
        nb[f"{key}_funnel_and_null"] = blk
    f = ROOT / "results/numbers.json"
    d = json.loads(f.read_text()) if f.exists() else {}
    d["partB"] = nb
    f.write_text(json.dumps(d, indent=1, default=float))
    print(json.dumps({k: (v if k != "wall_times" else "...") for k, v in nb.items()}, indent=1, default=float)[:3000])


if __name__ == "__main__":
    main()
