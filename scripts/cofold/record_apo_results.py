#!/usr/bin/env python3
"""Write the ligand-free refold + BURROW results into each case.yaml cofold section.

Frozen interpretation (docs/COFOLD_CONTROLS.md Sec.4): site recall is CONTINUOUS, never
thresholded at 0.8. On the four co-folding positive controls a true positive reached only
0.44-1.00 against a 4.5 A ligand footprint, so a present-side recall near 0.5 is normal and
an absent-side recall of 0 is the informative contrast.

For the five false-loss controls the prediction was explicit: the pocket should REAPPEAR in
a fresh ligand-free model of the absent-side (bad-model) proteins, which would show the
original loss call was a prediction artifact. This script states per control whether it did,
using absent-side median recall, and does not soften a negative.

Usage: record_apo_results.py [--dry-run]
"""
from __future__ import annotations
import argparse, csv, json
from pathlib import Path
import yaml

P = Path(__file__).resolve().parents[2]  # repo root
SUMM = P / "results/cofold/apo_pool/apo_pool_summary.tsv"
PER = P / "results/cofold/apo_pool/apo_burrow_pool.tsv"
FALSE_LOSS_PREFIX = "control_falseloss_"
REAPPEAR_MIN = 0.3          # same floor as ligand_pocket_match.py's "detected"


def rows(path: Path) -> list[dict]:
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def f(x):
    try:
        return round(float(x), 3)
    except (TypeError, ValueError):
        return None


def main() -> int:
    ap = argparse.ArgumentParser(); ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    per = rows(PER)
    out = []
    for s in rows(SUMM):
        cid = s["case"]
        cf = P / "cases" / cid / "case.yaml"
        if not cf.exists():
            print(f"[apo] no case.yaml for {cid}"); continue
        pres, absn = f(s["present_recall_median"]), f(s["absent_recall_median"])
        n_abs = int(s["n_absent"])
        members = [{"acc": r["acc"], "side": r["side"], "recall": f(r.get("recall")),
                    "precision": f(r.get("precision")), "burrow_rank": r.get("burrow_rank") or None,
                    "n_pockets": r.get("n_pockets"), "construct": r.get("construct"),
                    "atlas_footprint_plddt": f(r.get("atlas_plddt")),
                    "site_source": r.get("site_source")}
                   for r in per if r["case"] == cid]
        entry = {"model_id": f"{cid}_apo", "role": "apo_refold",
                 "status": None, "ligands": [], "seeds": 3, "samples_per_seed": 5,
                 "templates": "none", "msa": "ColabFold MMseqs2",
                 "method": ("protein-only AF3 refold, then pocketcompare describe --af-mode; "
                            "site recall of the best-covering BURROW pocket against the case "
                            "footprint, reported CONTINUOUS (never thresholded at 0.8)"),
                 "present_recall_median": pres, "absent_recall_median": absn,
                 "delta_recall": f(s["delta_recall"]),
                 "n_present": int(s["n_present"]), "n_absent": n_abs,
                 "members": members,
                 "table": "results/cofold/apo_pool/apo_burrow_pool.tsv"}
        if cid.startswith(FALSE_LOSS_PREFIX):
            if n_abs == 0:
                entry["status"] = "control_not_testable"
                entry["prediction"] = ("the pocket should reappear in a fresh ligand-free model of "
                                       "the absent-side (bad-model) proteins")
                entry["outcome"] = ("NOT TESTED: no pre-D11 state-A member could be recovered for "
                                    "this SOG, so there is no absent side to refold")
            elif absn is not None and absn >= REAPPEAR_MIN:
                entry["status"] = "pocket_reappeared"
                entry["outcome"] = (f"CONFIRMED: absent-side median recall {absn} >= {REAPPEAR_MIN}; "
                                    f"the pocket is recovered in fresh models of the bad-model "
                                    f"proteins, so the original loss call was a prediction artifact")
            else:
                entry["status"] = "pocket_did_not_reappear"
                entry["outcome"] = (f"NOT CONFIRMED: absent-side median recall {absn}; the pocket does "
                                    f"NOT reappear in fresh ligand-free models of the bad-model "
                                    f"proteins. The false-loss prediction fails for this control, so "
                                    f"the absence is not explained as a modelling artifact by this test")
            entry["prediction"] = entry.get("prediction", "the pocket should reappear in a fresh "
                                            "ligand-free model of the absent-side (bad-model) proteins")
        else:
            if pres is not None and absn is not None and pres >= REAPPEAR_MIN and absn < REAPPEAR_MIN:
                entry["status"] = "apo_consistent_with_call"
                entry["outcome"] = (f"present-side median recall {pres}, absent-side {absn}: BURROW on "
                                    f"fresh ligand-free models finds the pocket on the present side and "
                                    f"not on the absent side, independently of the atlas")
            elif pres is not None and pres < REAPPEAR_MIN:
                entry["status"] = "apo_pocket_not_recovered"
                entry["outcome"] = (f"present-side median recall only {pres}: the pocket is not clearly "
                                    f"recovered even on the present side, so this test is uninformative")
            else:
                entry["status"] = "apo_inconclusive"
                entry["outcome"] = (f"present-side median recall {pres}, absent-side {absn}: the two "
                                    f"sides are not separated, so the refold neither supports nor "
                                    f"contradicts the call")
        doc = yaml.safe_load(cf.read_text())
        cof = [c for c in (doc.get("cofold") or []) if isinstance(c, dict)
               and c.get("role") != "apo_refold"]
        doc["cofold"] = cof + [entry]
        if not a.dry_run:
            tmp = cf.with_suffix(".yaml.tmp"); tmp.write_text(yaml.safe_dump(doc, sort_keys=False, width=100))
            tmp.replace(cf)
        out.append((cid, entry["status"], pres, absn))
        print(f"{cid:42s} {entry['status']:28s} present={pres} absent={absn}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
