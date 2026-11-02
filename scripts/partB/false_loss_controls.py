#!/usr/bin/env python
"""False-loss negative controls (DESIGN §5.4 control (ii); D24 item 4).

Two kinds, both from the same code and cuts (2:1, tau 0.2):
  bad_model           losses called on states BEFORE the quality filters (O1_r2, no D11; primary and pLDDT >= 70
                      strata) whose derived-clade A members are demonstrably poor models (footprint pLDDT < 90,
                      footprint <= 10 residues (D11), length ratio outside 0.7-1.3 or UniProt fragment) and whose
                      loss call disappears once the quality filters apply (not called in O1_r3 primary, pre-rescue).
  rescue_instability  losses called in O1_r3 before rescue and not after rescue (clade A calls re-detected by the
                      alternative BURROW settings); high-quality models, so these test pocket-finder instability.
The 5 worst of each kind are marked `negative_control_example`.

    python scripts/partB/false_loss_controls.py --release O1_r3
"""
from __future__ import annotations

import argparse
import glob
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.dataset as ds

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from pocketevo import recon  # noqa: E402

SNAP = recon.cli_value("--snapshot") or recon.default_snapshot()   # --snapshot / $PE_SNAPSHOT / prereg


def loss_calls(path):
    e = pd.read_parquet(path)
    return e[e["called"] & (e["type"] == "loss") & (e["gain_cost"] == 2)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--release", "--sog-run", dest="release", default="O1_r3")
    recon.add_partb_args(ap, skip=("sog_run",))
    ap.add_argument("--pre-d11-run", default="O1_r2", help="release before D11 (bad-model controls)")
    ap.add_argument("--results-root", default=str(ROOT / "results/partB"), help="Part B results root (tests: tmp)")
    a = ap.parse_args()
    RB = Path(a.results_root)
    ev = RB / f"events/{SNAP}"
    out = RB / f"{SNAP}/{a.release}"
    after = loss_calls(ev / f"{a.release}/tau0.2/primary/events.parquet")
    after_keys = set(zip(after["og_id"], after["clade_species"]))
    pre = []
    for sc in ("state", "state_s70"):
        e = loss_calls(ev / f"{a.pre_d11_run}_preD11/tau0.2/{sc}/events.parquet").copy()
        e["pre_stratum"] = "primary_pLDDT90" if sc == "state" else "secondary_pLDDT70"
        pre.append(e)
    pre = pd.concat(pre, ignore_index=True)
    pre = pre[[k not in after_keys for k in zip(pre["og_id"], pre["clade_species"])]]
    pre = pre.drop_duplicates(["og_id", "char", "node"])
    trees = {t.tree_id: t for t in recon.iter_orthoxml_trees(
        ROOT / "data/orthology/ancestra_tier1_v0/accurate/hogs.orthoxml.gz", keep=set(pre["og_id"]) | set(after["og_id"]))}
    r2 = ds.dataset(sorted(glob.glob(f"/tmp/pe_big/sogs/{SNAP}/{a.pre_d11_run}/states_t0.2/shard=*.parquet")), format="parquet")
    sogs = set(pre["char"])
    cols = ["sog_id", "accession", "state", "state_s70", "footprint_plddt", "frac_lt70", "len_ratio"]
    b = r2.to_table(columns=cols).to_pandas()
    b = b[b["sog_id"].isin(sogs)]
    r3 = ds.dataset(sorted(glob.glob(f"/tmp/pe_big/sogs/{SNAP}/{a.release}/states_t0.2/shard=*.parquet")), format="parquet")
    c3 = r3.to_table(columns=["sog_id", "accession", "reason"]).to_pandas()
    c3 = c3[c3["sog_id"].isin(sogs)]
    small = set(zip(c3.loc[c3["reason"].fillna("").str.startswith("footprint_residues"), "sog_id"],
                    c3.loc[c3["reason"].fillna("").str.startswith("footprint_residues"), "accession"]))
    info = {(s, x): r for s, x, *r in zip(b["sog_id"], b["accession"], b["state"], b["state_s70"],
                                          b["footprint_plddt"], b["frac_lt70"], b["len_ratio"])}
    u = pd.read_parquet(ROOT / "data/annotations/uniprot_proteins.parquet", columns=["accession", "fragment"])
    frag = dict(zip(u["accession"], u["fragment"]))
    rows = []
    for r in pre.itertuples():
        t = trees[r.og_id]
        sc = "state" if r.pre_stratum.startswith("primary") else "state_s70"
        leaves = t.leaf_names()[t.leaf_lo[r.node]:t.leaf_hi[r.node]]
        aa = [x for x in leaves if (r.char, x) in info and info[(r.char, x)][0 if sc == "state" else 1] == "A"]
        if not aa:
            continue
        bad = []
        pls = []
        for x in aa:
            st, st70, pl, fr, lr = info[(r.char, x)]
            pls.append(pl)
            reasons = []
            if pl is not None and pl < 90:
                reasons.append("footprint_pLDDT<90")
            if (r.char, x) in small:
                reasons.append("footprint<=10res")
            if lr is not None and (lr < 0.7 or lr > 1.3):
                reasons.append("length_ratio")
            if frag.get(x) is True:
                reasons.append("uniprot_fragment")
            bad.append(";".join(reasons))
        nbad = sum(1 for x in bad if x)
        rows.append(dict(kind="bad_model", pre_stratum=r.pre_stratum, sog_id=r.char, og_id=r.og_id, node=r.node,
                         level=r.level, clade_species=r.clade_species, n_A_clade=len(aa), n_A_bad=nbad,
                         frac_A_bad=nbad / len(aa), min_A_footprint_plddt=float(np.nanmin(np.array(pls, float))),
                         bad_reasons=";".join(sorted({y for x in bad for y in x.split(";") if y})),
                         overturned_by="quality filters (A_primary pLDDT>=90, D11 > 10 residues)"))
    bm = pd.DataFrame(rows).sort_values(["frac_A_bad", "min_A_footprint_plddt"], ascending=[False, True])
    bm["negative_control_example"] = False
    full = bm[bm["frac_A_bad"] >= 0.999]
    # diverse examples: 3 low-pLDDT absences (pLDDT>=70 stratum) + 2 from the primary stratum (D11 / length)
    sel = list(full[full["pre_stratum"].str.startswith("secondary")].head(3).index) + \
        list(full[full["pre_stratum"].str.startswith("primary")].head(2).index)
    sel = sel + [i for i in full.index if i not in sel][:5 - len(sel)]
    bm.loc[sel, "negative_control_example"] = True
    # rescue-instability controls on the release
    e0 = after
    e1 = loss_calls(ev / f"{a.release}/tau0.2/primary/rescue/events/events.parquet")
    k1 = set(zip(e1["char"], e1["node"]))
    lost = e0[[k not in k1 for k in zip(e0["char"], e0["node"])]]
    rs = pd.read_parquet(ev / f"{a.release}/tau0.2/primary/rescue/sog_states/bucket=all.parquet")
    post = dict(zip(zip(rs["sog_id"], rs["accession"]), rs["state"]))
    ri = []
    for r in lost.itertuples():
        t = trees[r.og_id]
        leaves = t.leaf_names()[t.leaf_lo[r.node]:t.leaf_hi[r.node]]
        aft = [post.get((r.char, x)) for x in leaves]
        ri.append(dict(kind="rescue_instability", pre_stratum="O1_r3 primary (D11)", sog_id=r.char, og_id=r.og_id,
                       node=r.node, level=r.level, clade_species=r.clade_species, n_A_clade=int(r.n_A),
                       n_W_after_rescue=sum(x == "W" for x in aft),
                       frac_A_rescued=sum(x == "W" for x in aft) / max(1, int(r.n_A)),
                       overturned_by="rescue (k=1): clade A re-detected by alternative BURROW settings"))
    ri = pd.DataFrame(ri).sort_values(["frac_A_rescued", "n_A_clade"], ascending=False)
    ri["negative_control_example"] = False
    ri.loc[ri.head(5).index, "negative_control_example"] = True
    allc = pd.concat([bm, ri], ignore_index=True)
    f = out / "false_loss_controls.parquet"
    allc.to_parquet(f, index=False)
    recon.write_sidecar(f, SNAP, 0, {"script": "scripts/partB/false_loss_controls.py", "argv": sys.argv[1:]})
    print(allc.groupby("kind").size().to_dict())
    print(allc[allc["negative_control_example"]][["kind", "pre_stratum", "sog_id", "level", "n_A_clade",
                                                  "frac_A_bad", "min_A_footprint_plddt", "bad_reasons",
                                                  "frac_A_rescued"]].to_string())


if __name__ == "__main__":
    main()
