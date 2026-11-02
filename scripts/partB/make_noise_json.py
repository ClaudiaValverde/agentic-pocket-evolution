#!/usr/bin/env python
"""Measured observation-noise models for run_sim.py --noise-json (D9: rescue k = 1).

Inputs: results/observation/<snap>/{false_a_sites.parquet, rescue_hits.parquet, error_rates.parquet} and the real
SOG state table (A-eligibility per species). One JSON per footprint size class (<= 10 / > 10 residues; the
rescue residual differs ~3x between them) and per miss model:
  known  (primary): miss = false-A rate before rescue on experimentally known sites, per species (shrunk to the
                    pooled rate, m = 50 pseudo-sites), A_primary criteria;
  redetect (sensitivity): miss = 1 - same-source re-detection (BURROW seeds s1/s2, tau = 0.2), per species.
Per species s, class c (q_s = A / (A + U) among non-P members of the real SOG states = A-eligibility):
  true present: U = miss_s (1 - q_s); among non-U: A = miss_s q_s stillA_c / (1 - U), W = miss_s q_s (1 - stillA_c) / (1 - U)
  true absent : U = 1 - q_s; among non-U: W = control_rescued_c (true negatives rescued at k = 1),
                P = false_P (placeholder 0.02, not measured)
stillA_c = fraction of false-A cases with 0/6 rescue hits (k = 1, jaccard rule); control_rescued_c = fraction of
true-negative control patches with >= 1 hit; both A_primary, all halves, by known-site size class.

    python scripts/partB/make_noise_json.py --stratum primary
"""
from __future__ import annotations

import argparse
import glob
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.dataset as ds

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from pocketevo import recon  # noqa: E402

SNAP = recon.cli_value("--snapshot") or recon.default_snapshot()   # --snapshot / $PE_SNAPSHOT / prereg


def shrink(k, n, p0, m=50.0):
    return (k + m * p0) / (n + m)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--obs", default=str(ROOT / f"results/observation/{SNAP}"))
    ap.add_argument("--states", default=f"/tmp/pe_big/sogs/{SNAP}/O1/states_t0.2")
    ap.add_argument("--guide", default=str(ROOT / "data/orthology/ancestra_tier1_v0/tier1_guide_tree.nwk"))
    ap.add_argument("--stratum", default="primary", choices=["primary", "secondary"])
    ap.add_argument("--false-p", type=float, default=0.02)
    ap.add_argument("--out", default=str(ROOT / "results/partB/sim/measured"))
    recon.add_partb_args(ap)
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    obs = Path(a.obs)
    crit = "A_primary" if a.stratum == "primary" else "A_sens"
    fcol = "false_A" if a.stratum == "primary" else "false_A_sens"
    st = recon.SpeciesTree.from_file(a.guide)
    sp_of_tax = {int(st.tree.taxid[v]): st.names[v] for v in st.tree.leaves}
    # ---- known-site false-A before rescue, by species and size class
    f = pd.read_parquet(obs / "false_a_sites.parquet")
    f["cls"] = np.where(f["size_bin"].isin(["3-5", "6-10"]), "le10", "gt10")
    el = f[f[crit] == True]  # noqa: E712
    pooled = {c: float(g[fcol].mean()) for c, g in el.groupby("cls")}
    pooled["all"] = float(el[fcol].mean())
    # ---- rescue residual / control rescue by size class (k = 1, jaccard)
    r = pd.read_parquet(obs / "rescue_hits.parquet")
    r = r[(r["rule"] == "jaccard")]
    size = f.drop_duplicates("known_id").set_index("known_id")["cls"]
    r["cls"] = r["known_id"].map(size)
    if r["cls"].isna().any():   # controls: known_id of the matched case (if not resolvable -> pooled only)
        r["cls"] = r["cls"].fillna("na")
    cases = r[(r["kind"] == "false_A_case") & (r["A_primary"] == (a.stratum == "primary"))] \
        if a.stratum == "primary" else r[r["kind"] == "false_A_case"]
    ctrl = r[r["kind"] == "control"]
    stillA = {c: float((g["n_hits"] == 0).mean()) for c, g in cases.groupby("cls") if c != "na"}
    stillA["all"] = float((cases["n_hits"] == 0).mean())
    crs = {c: float((g["n_hits"] >= 1).mean()) for c, g in ctrl.groupby("cls") if c != "na" and len(g) >= 50}
    crs["all"] = float((ctrl["n_hits"] >= 1).mean())
    nctl = {c: int(len(g)) for c, g in ctrl.groupby("cls")}
    # prefer the observation agent's published stratified table (error_rates.parquet, k = 1) when present:
    # stillA_c = false_A_after_rescue / false_A (same criterion, size class); false-W on true absences =
    # pooled control_rescued (PI D11 note: 35 % of true-negative regions rescued at k = 1)
    er = pd.read_parquet(obs / "error_rates.parquet")
    er = er[(er["kingdom"] == "all") & (er["ligand_class"].isin(["all"])) & (er["plddt_band"] == "all") &
            (er["criterion"] == crit)]
    cmap = {"<=10": "le10", ">10": "gt10", "all": "all"}
    src_rates = {}
    for c_er, c in cmap.items():
        b = er[(er["metric"] == "false_A") & (er["k"] == 0) & (er["site_size_class"] == c_er)]["rate"]
        af = er[(er["metric"] == "false_A_after_rescue") & (er["k"] == 1) & (er["site_size_class"] == c_er)]["rate"]
        if len(b) and len(af):
            stillA[c] = float(af.iloc[0] / b.iloc[0])
            pooled[c] = float(b.iloc[0])
            src_rates[c] = dict(false_A_before=float(b.iloc[0]), false_A_after_k1=float(af.iloc[0]))
    cr = er[(er["metric"] == "control_rescued") & (er["k"] == 1) & (er["site_size_class"] == "all")]["rate"]
    if len(cr):
        crs = {"all": float(cr.iloc[0])}
        src_rates["control_rescued_all_k1"] = float(cr.iloc[0])
    # ---- re-detection (seeds) per taxid
    e = pd.read_parquet(obs / "error_rates.parquet")
    rd = e[(e["metric"] == "redetect") & e["model_source"].isin(["seed_s1", "seed_s2"]) & (e["tau"] == 0.2) &
           (e["plddt_band"] == "all") & (e["site_size_class"] == "all")]
    rd_pool = rd[rd["kingdom"] == "all"]
    rd0 = float(rd_pool["n_events"].sum() / rd_pool["n"].sum())
    rdt = rd[(rd["kingdom"] == "by_taxid") & rd["taxid"].astype(str).str.isdigit()].copy()
    rdt["taxid"] = rdt["taxid"].astype(int)
    rdt = rdt.groupby("taxid")[["n", "n_events"]].sum()
    # ---- A-eligibility per species from real states
    scol = "state" if a.stratum == "primary" else "state_s70"
    tb = ds.dataset(sorted(glob.glob(f"{a.states}/shard=*.parquet")), format="parquet").to_table(
        columns=["taxid", scol, "reason"]).to_pandas()
    tb = tb[tb[scol].isin(["A", "U"])]
    q = tb.groupby("taxid")[scol].apply(lambda s: float((s == "A").mean()))
    q0 = float((tb[scol] == "A").mean())
    # missed present sites: the homologous region exists, so 'aligned<0.8' and 'partial_site_overlap' (a pocket
    # is there) do not apply -> eligibility among the remaining U reasons (conservative: higher false-A)
    tp = tb[~tb["reason"].fillna("").str.startswith(("aligned", "partial_site_overlap"))]
    qp = tp.groupby("taxid")[scol].apply(lambda s: float((s == "A").mean()))
    qp0 = float((tp[scol] == "A").mean())
    meta = dict(stratum=a.stratum, criterion=crit, pooled_false_A_before=pooled, stillA_k1=stillA,
                control_rescued_k1=crs, n_controls_by_class=nctl, error_rates_parquet=src_rates, redetect_seed_pooled=rd0,
                A_eligibility_pooled=q0, A_eligibility_missed_present_pooled=qp0, false_P_placeholder=a.false_p, rescue_k=1, decision="D9")
    files = []
    for model in ("known", "redetect"):
        for cls in ("le10", "gt10"):
            per = {}
            rows = []
            for tx, sp in sp_of_tax.items():
                g = el[(el["taxid"] == tx) & (el["cls"] == cls)]
                if model == "known":
                    miss = shrink(float(g[fcol].sum()), len(g), pooled[cls])
                else:
                    n, k = (rdt.loc[tx, "n"], rdt.loc[tx, "n_events"]) if tx in rdt.index else (0, 0)
                    miss = 1 - shrink(float(k), float(n), rd0)
                qs = float(q.get(tx, q0))
                qps = float(qp.get(tx, qp0))
                sa = stillA.get(cls, stillA["all"])
                u_p = miss * (1 - qps)
                fa = miss * qps * sa / (1 - u_p)
                w = miss * qps * (1 - sa) / (1 - u_p)
                per[sp] = dict(false_A=fa, w_present=w, u_present=u_p, u_absent=1 - qs,
                               w_absent=crs.get(cls, crs["all"]), false_P=a.false_p)
                rows.append(dict(species=sp, taxid=tx, miss=miss, q=qs, q_present=qps, **per[sp]))
            df = pd.DataFrame(rows)
            gl = {k: float(df[k].mean()) for k in ("false_A", "w_present", "u_present", "u_absent", "w_absent",
                                                   "false_P")}
            d = {"global": gl, "per_species": per,
                 "source": f"measured {SNAP} ({model} miss model, footprint {cls}, {a.stratum} A, rescue k=1 D9)",
                 "meta": meta, "residual_false_A_given_eligible": float((df['miss'] * stillA.get(cls)).mean())}
            fp = out / f"noise_{a.stratum}_{model}_{cls}.json"
            fp.write_text(json.dumps(d, indent=1))
            files.append(fp)
            df.to_csv(out / f"noise_{a.stratum}_{model}_{cls}.csv", index=False)
            files.append(out / f"noise_{a.stratum}_{model}_{cls}.csv")
    for fp in files:
        recon.write_sidecar(fp, SNAP, 0, {"script": "scripts/partB/make_noise_json.py", "argv": sys.argv[1:]})
    print(json.dumps(meta, indent=1))


if __name__ == "__main__":
    main()
