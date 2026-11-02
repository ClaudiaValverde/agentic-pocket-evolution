#!/usr/bin/env python3
"""Pick property-matched decoy ligands that are not near-cognates.

Decoy rule (prereg cofold.decoy_matching + a structural floor):
  1. candidates = CCD ligands of the SAME ligand class seen in BioLiP2;
  2. drop anything BioLiP2 records as binding any member of the case;
  3. drop anything with ECFP4 Tanimoto > MAX_DECOY_TANIMOTO to the cognate ligand;
  4. of what remains, take the nearest in z-scored (MW, logP, HBD, HBA, rotatable bonds).

Step 3 was added after the calibration panel: pure property matching returned decoys that
were the SAME molecule as the cognate (MET/MED are enantiomers, Tanimoto 1.0; APC/ACP are
both AMP-PCP, 0.79), which are not non-binders at all.

Usage: pick_decoys.py --ligand OAA --exclude-acc P0ABH7 Q8S4Y1 --n 1
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import calibration_panel as CP

P = Path(__file__).resolve().parents[2]  # repo root
KEYS = ["mw", "logp", "hbd", "hba", "rotb"]


def build_pool():
    args = CP.ccd_props()
    bl = pd.read_parquet(P / "data/annotations/biolip_sites.parquet",
                         columns=["uniprot_acc", "ligand", "ligand_class"])
    bl = bl[bl.ligand_class.isin(["regular", "metal"])]
    pr = {}
    for c in bl.ligand.unique():
        if isinstance(c, str) and 1 <= len(c) <= 3 and c.isalnum():
            v = CP.props_for(c, *args)
            if v:
                pr[c] = v
    cls = {c: CP.classify(c, pr[c]["n_heavy"]) for c in pr}
    binds = bl.groupby("uniprot_acc").ligand.apply(set).to_dict()
    return args, pr, cls, binds


def pick(ligand, exclude_accs, n, args, pr, cls, binds, max_tan=None):
    max_tan = CP.MAX_DECOY_TANIMOTO if max_tan is None else max_tan
    k = cls.get(ligand)
    pool = [c for c in pr if cls.get(c) == k]
    M = np.array([[pr[c][x] for x in KEYS] for c in pool], float)
    Z = (M - M.mean(0)) / np.where(M.std(0) == 0, 1, M.std(0))
    z = {c: Z[i] for i, c in enumerate(pool)}
    own = set().union(*[binds.get(a, set()) for a in exclude_accs]) if exclude_accs else set()
    out, rejected = [], []
    for dist, c in sorted((float(np.linalg.norm(z[c] - z[ligand])), c)
                          for c in pool if c != ligand and c not in own):
        tm = CP.tanimoto(ligand, c, *args)
        if tm is not None and tm > max_tan:
            rejected.append((c, round(tm, 2)))
            continue
        out.append({"decoy": c, "prop_dist": round(dist, 3), "tanimoto": round(tm, 3)
                    if tm is not None else None,
                    "props": {x: round(pr[c][x], 2) for x in KEYS}})
        if len(out) >= n:
            break
    return {"ligand": ligand, "ligand_class": k,
            "cognate_props": {x: round(pr[ligand][x], 2) for x in KEYS},
            "decoys": out, "rejected_as_near_cognate": rejected[:8],
            "max_tanimoto": max_tan}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ligand", required=True)
    ap.add_argument("--exclude-acc", nargs="*", default=[])
    ap.add_argument("--n", type=int, default=1)
    a = ap.parse_args()
    args, pr, cls, binds = build_pool()
    print(json.dumps(pick(a.ligand, a.exclude_acc, a.n, args, pr, cls, binds), indent=1))


if __name__ == "__main__":
    main()
