#!/usr/bin/env python3
"""Merge the co-folding numbers into results/numbers.json under the `cofold` key.

Read-modify-write immediately before the write, and touch ONLY `cofold`: other agents hold
uncommitted edits to the other top-level keys, so the file must never be rebuilt wholesale.
A stale mtime check guards against a write landing on top of a change made while this ran.

Usage: export_numbers.py
"""
from __future__ import annotations
import csv, json, time
from pathlib import Path

P = Path(__file__).resolve().parents[2]  # repo root
NUM = P / "results/numbers.json"


def tsv(p: Path) -> list[dict]:
    with open(p) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def f(x):
    try:
        return round(float(x), 4)
    except (TypeError, ValueError):
        return None


def build() -> dict:
    th = json.loads((P / "results/cofold/calibration/threshold.json").read_text())
    out = {
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "calibration": {
            "frozen_utc": th["frozen_utc"],
            "threshold_file": "results/cofold/calibration/threshold.json",
            "threshold_sha256": (P / "results/cofold/calibration/threshold.json.sha256"
                                 ).read_text().split()[0],
            "panel": th["panel"],
            "single_readout_auc": th["single_readout_auc"],
            "auc_by_ligand_class": th["single_readout_auc_by_ligand_class"],
            "auc_by_negative_kind": th["single_readout_auc_by_negative_kind"],
            "d23_composite": th["d23_composite"],
            "runtime": th["runtime"],
            "caveats": th.get("caveats"),
            "decoy_caveat": ("21 of 48 panel decoys exceed ECFP4 Tanimoto 0.4 to their cognate "
                             "(two at 1.0), so some negatives are real binders. This biases "
                             "sensitivity and AUC DOWN and makes specificity 0.96 a lower bound. "
                             "Not retuned; see docs/COFOLD_CONTROLS.md Sec.7."),
        },
        "tier1_cases": {}, "apo_pool": {}, "deep_dives": {},
    }
    for d in sorted((P / "results/cofold/cases").glob("*/decision.json")):
        cid = d.parent.name
        dec = json.loads(d.read_text())
        g = dec.get("positive_control_gate", {})
        c = dec.get("d23_composite", {})
        out["tier1_cases"][cid] = {
            "verdict": dec.get("verdict"),
            "gate_pass": g.get("pass"), "gate_lig_rmsd_A": g.get("lig_rmsd"), "gate_ref": g.get("ref"),
            "blinding": (dec.get("blinding") or {}).get("status"),
            "derived_members_binding": c.get("derived_members_binding"),
            "ancestral_members_binding": c.get("ancestral_members_binding"),
            "ligand_off_footprint": c.get("ligand_off_footprint"),
        }
    for s in tsv(P / "results/cofold/apo_pool/apo_pool_summary.tsv"):
        out["apo_pool"][s["case"]] = {
            "n_present": int(s["n_present"]), "n_absent": int(s["n_absent"]),
            "present_recall_median": f(s["present_recall_median"]),
            "absent_recall_median": f(s["absent_recall_median"]),
            "delta_recall": f(s["delta_recall"]),
        }
    out["apo_pool"]["_note"] = ("site recall is continuous (never thresholded at 0.8); on the four "
                                "co-folding positive controls a true positive reached only 0.44-1.00")
    for dd in sorted((P / "results/cofold/deepdives").glob("*/result.json")):
        r = json.loads(dd.read_text())
        g = r.get("positive_control_gate", {})
        out["deep_dives"][dd.parent.name] = {
            "verdict": r["verdict"], "blinded": True,
            "criterion_committed_in": "21b5c9e (before either dive finished folding)",
            "gate_pass": g.get("pass"), "gate_lig_rmsd_A": g.get("lig_rmsd"), "gate_ref": g.get("ref"),
            "interaction_lig_plddt_per_seed": r["interaction_lig_plddt_per_seed"],
            "interaction_neg_pae_per_seed": r["interaction_neg_pae_per_seed"],
            "expected_sign": "positive",
            "cells_lig_plddt": {k: v["lig_plddt"] for k, v in r["cells"].items()},
            "cells_binds": {k: v["binds"] for k, v in r["cells"].items()},
        }
    dd = P / "results/cofold/deep_dives/deep_dives_summary.json"
    if dd.exists():
        out["deep_dives"]["_independent_rerun"] = {
            "source": "results/cofold/deep_dives/ (separate run, not produced by this agent)",
            **json.loads(dd.read_text())}
    return out


def main() -> int:
    before = NUM.stat().st_mtime
    cof = build()
    cur = json.loads(NUM.read_text())          # re-read as late as possible
    if NUM.stat().st_mtime != before:
        raise SystemExit("numbers.json changed while building; re-run")
    cur["cofold"] = cof                        # only this key is touched
    tmp = NUM.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(cur, indent=1) + "\n")
    tmp.replace(NUM)
    print(f"merged cofold into {NUM}; top-level keys now: {sorted(cur)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
