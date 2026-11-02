#!/usr/bin/env python
"""Score run_events.py output on the toy SOG table against its planted truth (make_toy_sogs.py).

    python scripts/partB/eval_toy.py --toy results/partB/sim/toy_sogs/TOY/O1 --run results/partB/sim/toy_run/real
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from pocketevo import recon  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--toy", default=str(ROOT / "results/partB/sim/toy_sogs/TOY/O1"))
    ap.add_argument("--run", default=str(ROOT / "results/partB/sim/toy_run/real"))
    ap.add_argument("--orthoxml", default=str(ROOT / "data/orthology/ancestra_tier1_v0/accurate/hogs.orthoxml.gz"))
    a = ap.parse_args()
    truth = pd.read_parquet(Path(a.toy) / "truth.parquet")
    ev = pd.read_parquet(Path(a.run) / "events.parquet")
    trees = {t.tree_id: t for t in recon.iter_orthoxml_trees(a.orthoxml, keep=set(truth["og_id"]))}
    tr = truth.set_index("sog_id")
    res = {}
    for (g, l), sub in ev[ev["called"]].groupby(["gain_cost", "loss_cost"]):
        key = f"{int(g)}:{int(l)}"
        res[key] = {}
        for typ in ("gain", "loss"):
            c = sub[sub["type"] == typ]
            ex = cons = 0
            for _, e in c.iterrows():
                t = tr.loc[e["char"]]
                if t["scenario"] != typ:
                    continue
                pl = int(t["planted_node"])
                tt = trees[t["og_id"]]
                if e["node"] == pl:
                    ex += 1
                if tt.leaf_lo[pl] <= tt.leaf_lo[e["node"]] and tt.leaf_hi[e["node"]] <= tt.leaf_hi[pl]:
                    cons += 1
            n_pl = int((truth["scenario"] == typ).sum())
            hit = len(set(c["char"]) & set(truth.loc[truth["scenario"] == typ, "sog_id"]))
            res[key][typ] = dict(n_called=len(c), n_exact=ex, n_consistent=cons, n_planted=n_pl,
                                 precision_consistent=cons / len(c) if len(c) else None,
                                 precision_exact=ex / len(c) if len(c) else None,
                                 recall_exact=ex / n_pl if n_pl else None, sogs_hit=hit,
                                 false_in_vertical_sogs=int(sum(tr.loc[x, "scenario"] == "vertical"
                                                                for x in c["char"])))
    sp = Path(a.run) / "specificity.parquet"
    if sp.exists():
        s = pd.read_parquet(sp)
        planted = truth[truth["specificity_node"] >= 0]
        key = set(zip(planted["sog_id"], planted["specificity_parent"]))
        called = s[s["called"]]
        tp = sum((r.sog_id, r.node) in key for r in called.itertuples())
        tested_planted = sum((r.sog_id, r.node) in key for r in s[s["tested"]].itertuples())
        res["specificity"] = dict(n_planted=len(planted), n_pairs=len(s), n_tested=int(s["tested"].sum()),
                                  n_planted_reaching_test=tested_planted, n_called=len(called), n_true=tp,
                                  precision=tp / len(called) if len(called) else None,
                                  recall=tp / len(planted) if len(planted) else None)
    out = Path(a.run) / "toy_eval.json"
    out.write_text(json.dumps(res, indent=1))
    recon.write_sidecar(out, "TOY", 0, {"script": "scripts/partB/eval_toy.py"})
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
