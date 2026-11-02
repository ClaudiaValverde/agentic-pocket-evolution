#!/usr/bin/env python3
"""Apply the FROZEN D23 composite to a scored case. Thresholds are read, never fitted.

Per-member call (results/cofold/calibration/threshold.json, frozen 2026-10-30T16:13:27Z):

    binds  <=>  lig_plddt >= t_plddt  AND  (centroid <= t_centroid OR recall >= t_recall)
                holding in >= 2 of the 3 seeds

Case call, from the project's interpretation rules (which follow from the panel numbers):
  * "binds" is informative          -- specificity 0.96
  * "does not bind" is WEAK         -- sensitivity 0.42
  * so a gain/loss is co-fold-consistent only when the PRESENT side binds and the ABSENT
    side does not, and the absent-side half is still reported as weak;
  * metal-only ligands are uninterpretable (panel AUC 0.58, at chance).

A case whose positive-control gate failed is `uninterpretable` whatever the members show.

Usage: d23_call.py --case <case_id> [--all]
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path

P = Path(__file__).resolve().parents[2]  # repo root
THRESH = P / "results/cofold/calibration/threshold.json"
IONS = {"ZN", "MG", "CA", "MN", "FE", "FE2", "CU", "NI", "CO", "NA", "K"}


def load_thresholds() -> dict:
    d = json.loads(THRESH.read_text())
    t = d["d23_composite"]["thresholds"]
    return {"t_plddt": float(t["t_plddt"]), "t_centroid": float(t["t_centroid"]),
            "t_recall": float(t["t_recall"]),
            "specificity": d["d23_composite"]["specificity"],
            "sensitivity": d["d23_composite"]["sensitivity"],
            "frozen_utc": d["frozen_utc"]}


def binds(per_seed: dict, th: dict, min_seeds: int = 2):
    """Return (call, n_seeds_passing, per-seed detail)."""
    n, detail = 0, {}
    for sd, v in sorted(per_seed.items()):
        p, c, r = v.get("lig_plddt"), v.get("centroid"), v.get("site_recall")
        if r is None:
            r = v.get("recall")
        geom = ((c is not None and c <= th["t_centroid"])
                or (r is not None and r >= th["t_recall"]))
        ok = p is not None and p >= th["t_plddt"] and geom
        detail[sd] = {"lig_plddt": p, "centroid": c, "recall": r,
                      "plddt_ok": p is not None and p >= th["t_plddt"],
                      "geometry_ok": geom, "binds": ok}
        n += ok
    return (n >= min_seeds), n, detail


def call_case(cid: str) -> dict:
    th = load_thresholds()
    wd = Path("/tmp/pe_big/cofold/cases") / cid
    res = P / "results/cofold/cases" / cid
    key = json.loads((wd / "blind_key.json").read_text())["map"]
    ps = json.loads((wd / "blinded_per_seed.json").read_text())
    dec = json.loads((res / "decision.json").read_text()) if (res / "decision.json").exists() else {}
    members = []
    for bid, k in key.items():
        if bid not in ps:
            continue
        pseed = {s: {"lig_plddt": ps[bid]["lig_plddt"].get(s),
                     "centroid": ps[bid]["centroid"].get(s),
                     "site_recall": (ps[bid].get("site_recall") or ps[bid].get("recall") or {}).get(s)}
                 for s in ps[bid]["lig_plddt"]}
        b, n, detail = binds(pseed, th)
        members.append({"blind_id": bid, "member": k["member"], "state": k["state"],
                        "ligand_kind": k["ligand_kind"], "ligands": k["ligands"],
                        "binds": b, "n_seeds_binding": n, "per_seed": detail})
    cog = [m for m in members if m["ligand_kind"] == "cognate"]
    pres = [m for m in cog if m["state"] == "derived"]
    absn = [m for m in cog if m["state"] == "ancestral"]
    lig = {l for m in cog for l in m["ligands"].split("+")}
    metal_only = bool(lig) and lig <= IONS
    gate = dec.get("positive_control_gate", {})
    pres_b = sum(m["binds"] for m in pres)
    absn_b = sum(m["binds"] for m in absn)
    if metal_only:
        verdict, why = "uninterpretable_by_cofolding_metal", \
            "metal-only ligand; panel AUC 0.58 for ligand pLDDT, at chance"
    elif not gate.get("pass"):
        verdict, why = "uninterpretable", \
            f"positive-control gate failed (ligand RMSD {gate.get('lig_rmsd')} A > 2.0)"
    elif pres_b > 0 and absn_b == 0:
        verdict, why = "cofold_consistent", \
            (f"{pres_b}/{len(pres)} present-side members bind (informative, spec 0.96); "
             f"0/{len(absn)} absent-side bind (WEAK, sens 0.42)")
    elif pres_b == 0:
        verdict, why = "uninformative_no_binding", \
            (f"no present-side member binds ({pres_b}/{len(pres)}); at sensitivity 0.42 a "
             "non-binding call is weak, so this is absence of evidence, not evidence of "
             "absence")
    else:
        verdict, why = "not_cofold_consistent", \
            f"{pres_b}/{len(pres)} present-side and {absn_b}/{len(absn)} absent-side bind"
    out = {"case": cid, "d23_verdict": verdict, "reason": why,
           "thresholds": th, "threshold_file": str(THRESH.relative_to(P)),
           "positive_control_gate": gate,
           "n_present_binding": pres_b, "n_present": len(pres),
           "n_absent_binding": absn_b, "n_absent": len(absn),
           "members": sorted(members, key=lambda m: (m["state"], m["ligand_kind"])),
           "legacy_dvx_verdict": dec.get("verdict")}
    tmp = res / "d23_call.json.tmp"
    tmp.write_text(json.dumps(out, indent=1)); tmp.replace(res / "d23_call.json")
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", nargs="*")
    ap.add_argument("--all", action="store_true")
    a = ap.parse_args()
    cases = a.case or []
    if a.all:
        cases = sorted(p.name for p in (P / "results/cofold/cases").iterdir()
                       if (p / "decision.json").exists())
    for c in cases:
        try:
            o = call_case(c)
            print(f"{c}: {o['d23_verdict']}  present {o['n_present_binding']}/{o['n_present']} "
                  f"absent {o['n_absent_binding']}/{o['n_absent']}  -- {o['reason']}")
        except Exception as e:  # noqa: BLE001
            print(f"{c}: FAILED {e}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
