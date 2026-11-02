#!/usr/bin/env python
"""Pocket-finder stability: compare BURROW pocket lists of one structure across settings by residue Jaccard.

For every pocket of the reference run (first --run), report the best residue-Jaccard match in each other run.

    python scripts/burrow_settings_jaccard.py --stem AF-P04035-F1-model_v6 --top 5 \
        --run default=out_def --run noaf=out_noaf --run fast=out_fast
"""
import argparse

import pandas as pd


def load(d, stem):
    t = pd.read_csv(f"{d}/{stem}.csv")
    t["rs"] = t.residues.fillna("").str.split().apply(set)
    return t


def jac(a, b):
    return len(a & b) / max(1, len(a | b))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stem", required=True)
    ap.add_argument("--run", action="append", required=True, help="name=output_dir (first = reference)")
    ap.add_argument("--top", type=int, default=5)
    a = ap.parse_args()
    runs = [r.split("=", 1) for r in a.run]
    ref_name, ref = runs[0][0], load(runs[0][1], a.stem)
    rows = []
    for _, p in ref.head(a.top).iterrows():
        row = {"ref_rank": int(p["rank"]), "ref_score": round(p.score, 3), "n_res": len(p.rs)}
        for name, d in runs[1:]:
            t = load(d, a.stem)
            js = [jac(p.rs, q) for q in t.rs]
            k = max(range(len(js)), key=js.__getitem__) if js else None
            row[f"{name}_J"] = round(js[k], 2) if js else 0.0
            row[f"{name}_rank"] = int(t["rank"].iloc[k]) if js else None
        rows.append(row)
    print(f"reference = {ref_name}; n pockets per run: " +
          ", ".join(f"{n}={len(load(d, a.stem))}" for n, d in runs))
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    main()
