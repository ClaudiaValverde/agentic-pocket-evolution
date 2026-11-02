#!/usr/bin/env python
"""Part A stage 7: estimands (r_c, F_anc / F_offsite / F_H4 / F_analog) with expected FP, 3 weightings,
two-way bootstrap over root groups, E-threshold and tail-variant sensitivity, strata, raw top-k composition.

Expected false positives per pair (q0), declared rule (PARTA.md Deviations D3):
  nominal  q0 = lambda * E_thr / N_db, lambda = max(1, held-out observed/expected one-direction ratio at the deepest
           threshold with >= 100 expected pairs) -- the review's 'calibrated null' rate, never below nominal;
  bound    q0 = observed both-direction rate among the 10^9 random non-homologous pairs at the same E threshold
           (treats every random-pair tail hit as false: lower bound on F_analog; circular, see Deviations).
Outputs: results/partA/<snap>/estimands_<score>.csv, bootstrap_<score>.parquet (bulk), topk_composition.csv,
strata_<score>.csv, estimands_<score>.json
"""
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from _common import parser, resolve, logger, pairs
from pocketevo import estimands as est
from pocketevo import homology as hom

ap = parser(__doc__)
ap.add_argument("--score", choices=["cos", "csls"], default="cos")
ap.add_argument("--n-boot", type=int, default=200)
ap.add_argument("--no-strata", action="store_true")
ap.add_argument("--n-boot-strata", type=int, default=100)
ap.add_argument("--site-weight", choices=["none", "species_henikoff"], default="none",
                help="sensitivity: species_henikoff = 1 / (#universe proteins of the same taxid in the root group) per "
                     "site, pair weight w_i w_j, point estimates only (no bootstrap, no strata)")
ap.add_argument("--split-table", default=None, help="og_id -> split table (D20, e.g. O2_r1_split/og_split_o2.parquet)")
a = resolve(ap.parse_args())
log = logger(a.work / "partA.log")
T0 = time.time()
cd = a.work / f"cls_{a.score}{a.tag}"
sig = pd.read_parquet(cd / "pairs_cls.parquet")
smp = pd.read_parquet(cd / "sample_rest.parquet")
ws = pd.read_parquet(cd / "sample_within.parquet")
csumm = json.loads((cd / "summary.json").read_text())
c = np.load(a.work / "context.npz")
forest = hom.parse_orthoxml(str(Path(a.orthology) / "hogs.orthoxml.gz"))
ctx = hom.PairContext(c["site_prot"], c["site_lab"], c["prot_group"], c["prot_gene"], forest)
N = len(ctx.site_prot)
G = int(ctx.prot_group.max()) + 1
calt = pd.read_csv(a.results / f"null_{a.score}_calibration_heldout.csv")
ext_p = a.results / f"null_{a.score}_calibration_extended.csv"
ext = pd.read_csv(ext_p) if ext_p.exists() else None
PRIMARY_VARIANT = "pooled"
E_PRIMARY = float(a.prereg["pairs"]["null"]["significance"]["threshold"]) if "null" in a.prereg["pairs"] else \
    float(a.prereg["pairs"][None]["significance"]["threshold"])
E_GRID = [1e-4, E_PRIMARY, 1.0]
split = {k: float(v) for k, v in csumm["split_fraction"].items() if k in ("H1", "H2", "H12")}
for k in ("H1", "H2", "H12"):
    split.setdefault(k, 0.0)
log("split", split)


def lam(variant):
    d = calt[(calt.variant == variant) & (calt.expected_one_dir >= 100)]
    d = d.sort_values("p_threshold").iloc[0]
    return max(1.0, float(d.ratio_ij), float(d.ratio_ji)), float(d.p_threshold)


def q0_bound(variant, E):
    if ext is None:
        return np.nan
    d = ext[ext.variant == variant]
    d = d.iloc[(np.abs(np.log(d.E_threshold) - np.log(E))).argsort()[:1]]
    return float(d.obs_both.iloc[0] / d.n_pairs.iloc[0])


SW = None
if a.site_weight == "species_henikoff":
    uni_t = pd.read_parquet(a.work / "universe.parquet", columns=["taxid"]).taxid.to_numpy()
    prot_tax = pd.Series(uni_t).groupby(ctx.site_prot).first().reindex(np.arange(len(ctx.prot_group))).to_numpy()
    key = pd.Series(ctx.prot_group.astype(np.int64) * 10**8 + np.nan_to_num(prot_tax).astype(np.int64))
    n_same = key.map(key.value_counts()).to_numpy().astype(float)
    SW = (1.0 / n_same)[ctx.site_prot]
    sig = sig.assign(pair_w=SW[sig.i.values] * SW[sig.j.values])
    smp = smp.assign(pair_w=SW[smp.i.values] * SW[smp.j.values])
    a.n_boot, a.no_strata = 0, True
    a.tag = a.tag + "_henikoff"
subsets_all = hom.LabelSubsets.build(ctx, np.ones(N, bool), site_weight=SW)
totals_all = hom.exact_class_totals(ctx, subsets=subsets_all, site_weight=SW)
log("totals", json.dumps({k: f"{v:.4g}" for k, v in totals_all.items()}))
rows, payload = [], {"totals": totals_all, "split": split, "split_mode": csumm["split_mode"], "N_db": N}
SCHEMES = ("per_pair", "per_site", "per_family_pair")


def run(variant, E, q0kind, boot=False, sig_=None, smp_=None, totals=None, tag="all", subsets=None, mask=None,
        n_boot=None):
    sig_ = sig if sig_ is None else sig_
    smp_ = smp if smp_ is None else smp_
    totals = totals_all if totals is None else totals
    S = sig_[sig_[f"Emax_{variant}"] < E]
    S = S[S.cls != "H0"]
    lam_v, lam_p = lam(variant)
    q0 = lam_v * E / N if q0kind == "nominal" else q0_bound(variant, E)
    point = est.estimands(S, totals, split, smp_, q0, schemes=SCHEMES)
    rec = {"tag": tag, "variant": variant, "E": E, "q0kind": q0kind, "q0": q0, "lambda": lam_v,
           "lambda_p": lam_p, "n_sig": len(S)}
    flat = est.flatten(point)
    out = []
    reps = None
    if boot:
        sub = subsets if subsets is not None else subsets_all
        tf = (lambda w: hom.exact_class_totals(ctx, site_mask=mask, group_weight=w, subsets=sub))
        reps = est.bootstrap(S, tf, split, smp_, q0, G, n_boot or a.n_boot, a.seed, SCHEMES, log=log)
    for k, v in flat.items():
        r = dict(rec, key=k, value=v)
        if reps is not None and k in reps:
            col = pd.to_numeric(reps[k], errors="coerce").dropna()
            if len(col):
                r["ci_lo"], r["ci_hi"] = np.quantile(col, [0.025, 0.975])
        out.append(r)
    return out, reps


# primary: pooled tail, E < 0.01, nominal q0, bootstrap
prim, reps = run(PRIMARY_VARIANT, E_PRIMARY, "nominal", boot=a.n_boot > 0)
rows += prim
if reps is not None:
    reps.to_parquet(a.work / f"bootstrap_{a.score}{a.tag}.parquet", index=False)
log(f"primary done {time.time() - T0:.0f}s")
for variant in ("pooled", "perq"):
    for E in E_GRID:
        for kind in ("nominal", "bound"):
            if variant == PRIMARY_VARIANT and E == E_PRIMARY and kind == "nominal":
                continue
            r, _ = run(variant, E, kind)
            rows += r
df = pd.DataFrame(rows)
df.to_csv(a.results / f"estimands_{a.score}{a.tag}.csv", index=False)
log(f"sensitivity done {time.time() - T0:.0f}s")

# ---------------------------------------------------------------- strata (point estimates + bootstrap off)
if not a.no_strata:
    uni = pd.read_parquet(a.work / "universe.parquet", columns=["tier", "plddt_mean_rep", "rep_site_volume_A3"])
    vol = uni.rep_site_volume_A3.to_numpy()
    t1, t2 = np.quantile(vol, [1 / 3, 2 / 3])
    strata = {"both_tier1": uni.tier.to_numpy() == 1, "both_tier2": uni.tier.to_numpy() == 2,
              "both_plddt_ge90": uni.plddt_mean_rep.to_numpy() >= 90,
              "both_vol_small": vol <= t1, "both_vol_mid": (vol > t1) & (vol <= t2), "both_vol_large": vol > t2}
    # prereg tune/test split (root-HOG smallest-accession key; control families forced to test)
    prot_acc = pd.read_parquet(a.work / "proteins.parquet", columns=["accession"]).accession.to_numpy()
    forced = a.prereg["tune_test"]["forced_test_half"]["root_hogs_containing"]
    if a.split_table and Path(a.split_table).exists():
        half = hom.group_halves_from_table(int(ctx.prot_group.max()) + 1, forest, pd.read_parquet(a.split_table))
        payload["split_source"] = str(a.split_table)
    else:
        half = hom.group_halves(prot_acc, ctx.prot_group, forest, forced)
        payload["split_source"] = "computed: prereg smallest-accession hash (no crossover exclusion)"
    site_half = half[ctx.prot_group[ctx.site_prot]]
    np.save(a.work / "site_half.npy", site_half)
    strata = {"test_half": site_half == 1, "tune_half": site_half == 0, **strata}
    payload["half_sites"] = {"test": int((site_half == 1).sum()), "tune": int((site_half == 0).sum()),
                             "excluded_crossover_or_nonhog": int((site_half == -1).sum())}
    srows, sreps = [], []
    for name, mask in strata.items():
        sub = hom.LabelSubsets.build(ctx, mask)
        tot = hom.exact_class_totals(ctx, site_mask=mask, subsets=sub)
        sg = sig[mask[sig.i.values] & mask[sig.j.values]]
        sm = smp[mask[smp.i.values] & mask[smp.j.values]]
        r, rp = run(PRIMARY_VARIANT, E_PRIMARY, "nominal", sig_=sg, smp_=sm, totals=tot, tag=name,
                    boot=a.n_boot_strata > 0, subsets=sub, mask=mask, n_boot=a.n_boot_strata)
        if rp is not None:
            rp = rp.copy(); rp["tag"] = name; rp["rep"] = np.arange(len(rp)); sreps.append(rp)
        for x in r:
            x["n_sites_stratum"] = int(mask.sum()); x["n_rest_sample_stratum"] = len(sm)
        srows += r
        log(f"stratum {name} done")
    pd.DataFrame(srows).to_csv(a.results / f"strata_{a.score}{a.tag}.csv", index=False)
    if sreps:
        pd.concat(sreps, ignore_index=True).to_parquet(a.work / f"bootstrap_strata_{a.score}{a.tag}.parquet", index=False)
    payload["volume_tertiles_A3"] = [float(t1), float(t2)]

payload["wall_s"] = time.time() - T0
(a.results / f"estimands_{a.score}{a.tag}.json").write_text(json.dumps(payload, indent=1, default=float))
pairs.write_sidecar(a.results / f"estimands_{a.score}{a.tag}.csv", a.snapshot, a.seed,
                    {"stage": f"estimands_{a.score}{a.tag}", "n_boot": a.n_boot, "primary_variant": PRIMARY_VARIANT,
                     "E_primary": E_PRIMARY, "wall_s": payload["wall_s"]})
log("estimands done", f"{time.time() - T0:.0f}s")
