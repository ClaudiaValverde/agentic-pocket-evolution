#!/usr/bin/env python3
"""Score the known-answer deep dives. Committed BEFORE any deep-dive result existed.

2 x 2 design (protein x ligand), 3 seeds x 5 samples, no templates, ColabFold MSAs,
HMAC-blinded job names with a persisted salt. Readouts are computed blind by collect();
labels are attached only in this script.

  steroid: MR LBD, GR LBD x aldosterone (AS4), cortisol (HCY).
           Known: MR binds both; GR prefers cortisol.  I = (GR.HCY-GR.AS4)-(MR.HCY-MR.AS4) > 0
  LDH/MDH: PfLDH, hMDH1 x pyruvate (PYR), oxaloacetate (OAA), + NADH.
           Known: LDH prefers PYR, MDH prefers OAA.   I = (LDH.PYR-LDH.OAA)-(MDH.PYR-MDH.OAA) > 0

PRE-DECLARED CRITERION: "resolved" only if I has the expected sign in ALL 3 seeds for
ligand pLDDT AND for min P-L PAE (sign-flipped, lower is better), and the frozen D23
composite does not call a dispreferred pairing "binds" while a preferred one "does not".
One readout or not all seeds -> "partial"; otherwise "not_resolved". The positive-control
gate (ligand RMSD <= 2 A) is applied first; failure -> "uninterpretable".
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import case_cofold as CC

P = Path(__file__).resolve().parents[2]  # repo root
DESIGN = {
 "deepdive_steroid_MR_vs_GR": dict(prot={"MR_LBD": "MR", "GR_LBD": "GR"}, ligs=["AS4", "HCY"],
     terms=[("GR","HCY",1),("GR","AS4",-1),("MR","HCY",-1),("MR","AS4",1)],
     preferred=[("GR","HCY"),("MR","AS4"),("MR","HCY")], dispreferred=[("GR","AS4")]),
 "deepdive_LDH_vs_MDH": dict(prot={"PfLDH": "LDH", "hMDH1": "MDH"}, ligs=["PYR", "OAA"],
     terms=[("LDH","PYR",1),("LDH","OAA",-1),("MDH","PYR",-1),("MDH","OAA",1)],
     preferred=[("LDH","PYR"),("MDH","OAA")], dispreferred=[("LDH","OAA"),("MDH","PYR")]),
}


def main() -> int:
    ap = argparse.ArgumentParser(); ap.add_argument("--case", required=True, choices=sorted(DESIGN))
    a = ap.parse_args(); D = DESIGN[a.case]
    wd = Path("/tmp/pe_big/cofold/cases") / a.case
    spec = CC.load_spec(P / "cases/cofold" / f"{a.case}.yaml")
    th = json.loads((P / "results/cofold/calibration/threshold.json").read_text())["d23_composite"]["thresholds"]
    key = json.loads((wd / "blind_key.json").read_text())["map"]
    ps = json.loads((wd / "blinded_per_seed.json").read_text())
    cell = {}
    for bid, k in key.items():
        pr = D["prot"].get(k["member"]); lg = next((x for x in D["ligs"] if x in k["ligands"].split("+")), None)
        if pr and lg and bid in ps:
            binds, n = CC.d23_binds(ps[bid], th)
            cell[(pr, lg)] = {"per_seed": ps[bid], "binds": binds, "n_seeds_binding": n}
    seeds = ["1", "2", "3"]
    def I(metric, sign):
        out = []
        for s in seeds:
            try:
                out.append(round(sign * sum(w * cell[(p, l)]["per_seed"][metric][s] for p, l, w in D["terms"]), 3))
            except KeyError:
                out.append(None)
        return out
    Ip, Ie = I("lig_plddt", +1), I("pae_min", -1)
    okp = all(x is not None and x > 0 for x in Ip); oke = all(x is not None and x > 0 for x in Ie)
    contra = any(cell.get(d, {}).get("binds") for d in D["dispreferred"]) and \
             any(not cell.get(p, {}).get("binds", False) for p in D["preferred"])
    gate = CC.gate(wd, spec, CC.load_prereg())
    if not gate.get("pass"):
        verdict = "uninterpretable"
    elif okp and oke and not contra:
        verdict = "resolved"
    elif okp or oke:
        verdict = "partial"
    else:
        verdict = "not_resolved"
    res = {"case": a.case, "verdict": verdict, "positive_control_gate": gate,
           "interaction_lig_plddt_per_seed": Ip, "interaction_neg_pae_per_seed": Ie,
           "expected_sign": "positive", "d23_contradiction": bool(contra),
           "cells": {f"{p}+{l}": {"binds": v["binds"], "n_seeds_binding": v["n_seeds_binding"],
                                  "lig_plddt": v["per_seed"]["lig_plddt"], "pae_min": v["per_seed"]["pae_min"],
                                  "centroid": v["per_seed"].get("centroid"), "recall": v["per_seed"].get("recall")}
                     for (p, l), v in sorted(cell.items())},
           "blinding": "BLINDED: readouts computed by collect() on hashed job names in one automated "
                       "pass; labels attached only in deepdive_score.py, whose criterion was committed "
                       "before the folds finished",
           "criterion": __doc__.split("PRE-DECLARED CRITERION:")[1].split('"""')[0].strip()}
    out = P / "results/cofold/deepdives" / a.case; out.mkdir(parents=True, exist_ok=True)
    tmp = out / "result.json.tmp"; tmp.write_text(json.dumps(res, indent=1)); tmp.replace(out / "result.json")
    print(json.dumps({k: res[k] for k in ("verdict", "interaction_lig_plddt_per_seed",
                                          "interaction_neg_pae_per_seed", "d23_contradiction")}))
    print("gate", gate.get("pass"), gate.get("lig_rmsd"))
    for c, v in res["cells"].items():
        print(f"  {c:9s} binds={v['binds']} pLDDT={v['lig_plddt']} PAE={v['pae_min']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
