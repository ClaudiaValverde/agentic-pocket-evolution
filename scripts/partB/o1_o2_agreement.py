#!/usr/bin/env python
"""O1 vs O2 agreement for gain calls (misassigned-orthology check; PI 16:xx).

O1 gain = primary call on O1_r3 (tau 0.2, 2:1, post-rescue). Its derived-clade P members (accessions) are looked up
in the O2 events (tau 0.2, 2:1, O2_r1 + rescue overlay):
  recovered_called     an O2 *called* gain whose derived clade holds >= 50 % of the O1 derived P members
  recovered_candidate  same with an O2 frozen-rule candidate (before calibrated cuts)
  recovered_event      same with any unambiguous O2 MPR gain event (clade size ok)
  o2_same_sog          the O1 derived P members sit (>= 50 %) in one O2 SOG (any state)
Output: results/partB/<snap>/<sog-run>/o1_o2_gain_agreement.parquet + .json
"""
from __future__ import annotations

import argparse
import glob
import json
import sys
from collections import Counter
from pathlib import Path

import pandas as pd
import pyarrow.dataset as ds

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from pocketevo import recon  # noqa: E402

SNAP = recon.cli_value("--snapshot") or recon.default_snapshot()   # --snapshot / $PE_SNAPSHOT / prereg


def clade_P(ev, trees, states):
    out = []
    for r in ev.itertuples():
        t = trees[r.og_id]
        leaves = t.leaf_names()[t.leaf_lo[r.node]:t.leaf_hi[r.node]]
        out.append({x for x in leaves if states.get((r.char, x)) == "P"})
    return out


def main():
    ap = argparse.ArgumentParser()
    recon.add_partb_args(ap, sog_run="O1_r3")
    ap.add_argument("--o2-run", default="O2_r1")
    a = ap.parse_args()
    R1, R2 = a.sog_run, a.o2_run
    e1 = pd.read_parquet(ROOT / f"results/partB/events/{SNAP}/{R1}/tau0.2/primary/rescue/events/events.parquet")
    g1 = e1[e1["called"] & (e1["type"] == "gain") & (e1["gain_cost"] == 2)].copy()
    rs1 = pd.read_parquet(ROOT / f"results/partB/events/{SNAP}/{R1}/tau0.2/primary/rescue/sog_states/bucket=all.parquet")
    s1 = dict(zip(zip(rs1["sog_id"], rs1["accession"]), rs1["state"]))
    t1 = {t.tree_id: t for t in recon.iter_orthoxml_trees(
        ROOT / "data/orthology/ancestra_tier1_v0/accurate/hogs.orthoxml.gz", keep=set(g1["og_id"]))}
    g1["P_members"] = clade_P(g1, t1, s1)
    e2 = pd.read_parquet(ROOT / f"results/partB/events/{SNAP}/{R2}/tau0.2/primary/events.parquet")
    e2 = e2[(e2["type"] == "gain") & (e2["gain_cost"] == 2)]
    acc1 = set().union(*g1["P_members"])
    st2 = ds.dataset(sorted(glob.glob(f"/tmp/pe_big/partB/{SNAP}/{R2}/states_resc_t0.2/shard=*.parquet")),
                     format="parquet").to_table(columns=["og_id", "sog_id", "accession", "state"]).to_pandas()
    hits = st2[st2["accession"].isin(acc1)]
    sog_of = hits.groupby("accession")["sog_id"].apply(set).to_dict()
    st2 = st2[st2["sog_id"].isin(set(e2["char"]))]
    s2 = dict(zip(zip(st2["sog_id"], st2["accession"]), st2["state"]))
    t2 = {t.tree_id: t for t in recon.iter_orthoxml_trees(
        ROOT / "data/orthology/ancestra_O2/accurate/hogs.orthoxml.gz", keep=set(e2["og_id"]))}
    e2 = e2.copy()
    e2["P_members"] = clade_P(e2, t2, s2)
    lev = {"called": e2[e2["called"]], "candidate": e2[e2["candidate"]], "event": e2[e2["unambiguous"]]}
    rows = []
    for r in g1.itertuples():
        m = r.P_members
        rec = dict(sog_id=r.char, og_id=r.og_id, level=r.level, n_P=len(m), parent_kind=r.parent_kind)
        for k, df in lev.items():
            best = 0.0; best_sog = None
            for s, pm in zip(df["char"], df["P_members"]):
                if pm:
                    f = len(m & pm) / max(1, len(m))
                    if f > best:
                        best, best_sog = f, s
            rec[f"recovered_{k}"] = best >= 0.5
            rec[f"best_overlap_{k}"] = best
            rec[f"o2_sog_{k}"] = best_sog
        c = Counter(s for x in m for s in sog_of.get(x, ()))
        rec["o2_same_sog"] = bool(c and c.most_common(1)[0][1] / max(1, len(m)) >= 0.5)
        rec["n_P_in_O2"] = sum(1 for x in m if x in sog_of)
        rows.append(rec)
    df = pd.DataFrame(rows)
    out = ROOT / f"results/partB/{SNAP}/{R1}"
    f = out / "o1_o2_gain_agreement.parquet"
    df.to_parquet(f, index=False)
    summ = {"n_O1_gains": len(df), **{c: int(df[c].sum()) for c in df.columns if c.startswith("recovered_") or c == "o2_same_sog"},
            "by_parent_kind": {k: {c: int(g[c].sum()) for c in ("recovered_called", "recovered_event")} | {"n": len(g)}
                               for k, g in df.groupby("parent_kind")},
            "n_O2_gains_called": int(len(lev["called"])), "rule": "O2 clade holds >= 50 % of the O1 derived P members"}
    # diagnosis: where do the O1 members of these SOGs sit on O2 (modal O2 SOG of the derived P members)?
    rsx = rs1[rs1["sog_id"].isin(set(g1["char"])) & rs1["state"].isin(["P", "A", "W"])][["sog_id", "accession", "state"]]
    import pyarrow as pa
    import pyarrow.compute as pc
    raw = ds.dataset(sorted(glob.glob(f"/tmp/pe_big/sogs/{SNAP}/{R2}/states_t0.2/shard=*.parquet")), format="parquet")
    r2 = raw.to_table(filter=pc.is_in(pc.field("accession"), value_set=pa.array(sorted(set(rsx["accession"])))),
                      columns=["og_id", "sog_id", "accession", "state", "reason"]).to_pandas()
    p1 = rsx[rsx["state"] == "P"].merge(r2[r2["state"] == "P"][["sog_id", "accession"]], on="accession", suffixes=("1", "2"))
    best = p1.groupby(["sog_id1", "sog_id2"]).size().reset_index(name="n").sort_values("n").drop_duplicates("sog_id1", keep="last")
    y = rsx.rename(columns={"sog_id": "sog_id1", "state": "state1"}).merge(best[["sog_id1", "sog_id2"]], on="sog_id1") \
        .merge(r2.rename(columns={"sog_id": "sog_id2", "state": "state2"}), on=["sog_id2", "accession"], how="left")
    y["state2"] = y["state2"].fillna("not_in_O2_SOG_family")
    o2fam = dict(zip(r2["accession"], r2["og_id"]))
    ct = pd.crosstab(y["state1"], y["state2"])
    ua = y[(y["state1"] == "A") & (y["state2"] == "U")]["reason"].fillna("").str.split(":").str[0].value_counts()
    outA = y[y["state1"] == "A"]
    diag = {"O1_state_vs_raw_O2_state_in_modal_O2_SOG": {i: {c: int(ct.at[i, c]) for c in ct.columns} for i in ct.index},
            "O1_A_that_are_U_on_O2_by_reason": {k: int(v) for k, v in ua.items()},
            "O1_A_members_in_another_O2_root_HOG": int((outA["state2"] == "not_in_O2_SOG_family").sum()),
            "O1_A_members_total": int(len(outA)),
            "reading": "O1 gains rest on outgroup absences; the share of those absent proteins that O2 places in a "
                       "different root HOG than the derived clade measures orthology disagreement between O1 and O2."}
    summ["diagnosis"] = diag
    fj = out / "o1_o2_gain_agreement.json"
    fj.write_text(json.dumps(summ, indent=1))
    for x in (f, fj):
        recon.write_sidecar(x, SNAP, 0, {"script": "scripts/partB/o1_o2_agreement.py"})
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()
