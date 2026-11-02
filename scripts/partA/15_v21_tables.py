#!/usr/bin/env python
"""Part A v2.1 stage 15: test half and tau = 0.3 sensitivity for every headline metric (prereg
sogs.tau_column_jaccard.reported "test half + full grid for every headline metric"; ledger D10).

Rows: run (tau0.2 primary O1 SOGs | tau0.3 O1_r2 sensitivity_tau0.3) x subset (all | test_half | tune_half).
Metrics: F_anc / F_analog (3 weightings, CI), r_adj for H1s H2s H12s H3j H3s H4 N, F_conv,func (per family pair /
per pair, full labels, CI), F_conv,dir (per pair / per family pair, CI), descent M1 ortholog / paralog and the
identity-standardised difference (descent_2d_summary.csv).
Output: headline_halves_tau.csv (long format: run, subset, metric, value, ci_lo, ci_hi)
"""
import numpy as np
import pandas as pd

from _common import parser, resolve, pairs

a = resolve(parser(__doc__).parse_args())
R = a.results
RUNS = {"tau0.2": "", "tau0.3": "_t0.3"}
rows = []


def add(run, subset, metric, v, lo=np.nan, hi=np.nan, src=""):
    rows.append({"run": run, "subset": subset, "metric": metric, "value": v, "ci_lo": lo, "ci_hi": hi, "source": src})


for run, tag in RUNS.items():
    est_p, st_p = R / f"estimands_cos{tag}.csv", R / f"strata_cos{tag}.csv"
    if not est_p.exists():
        continue
    e = pd.read_csv(est_p)
    e = e[(e.tag == "all") & (e.variant == "pooled") & (e.E == 0.01) & (e.q0kind == "nominal")]
    s = pd.read_csv(st_p) if st_p.exists() else pd.DataFrame(columns=e.columns)
    for subset, d in (("all", e), ("test_half", s[s.tag == "test_half"]), ("tune_half", s[s.tag == "tune_half"])):
        src = est_p.name if subset == "all" else st_p.name
        for sch in ("per_family_pair", "per_site", "per_pair"):
            for k in ("F_anc", "F_offsite", "F_H4", "F_analog"):
                r = d[d.key == f"{sch}.F_adj.{k}"]
                if len(r):
                    add(run, subset, f"{k}.{sch}", r.value.iloc[0], r.get("ci_lo", pd.Series([np.nan])).iloc[0],
                        r.get("ci_hi", pd.Series([np.nan])).iloc[0], src)
        for c in ("H1s", "H1x", "H2s", "H2x", "H12s", "H12x", "H3j", "H3s", "H4", "N"):
            r = d[d.key == f"per_pair.r_adj.{c}"]
            if len(r):
                add(run, subset, f"r_adj.{c}", r.value.iloc[0], r.get("ci_lo", pd.Series([np.nan])).iloc[0],
                    r.get("ci_hi", pd.Series([np.nan])).iloc[0], src)
    fc = R / f"functional_convergence{tag}.csv"
    if fc.exists():
        f = pd.read_csv(fc)
        f = f[f.exclude_train & (f.category == "any") & (f.label_variant == "full") & (f["set"] == "sig_N_vs_matched_random_N")]
        for subset in ("all", "test_half"):
            r = f[f.subset == subset]
            if len(r):
                r = r.iloc[0]
                for sch in ("per_family_pair", "per_pair"):
                    add(run, subset, f"F_conv_func.{sch}", r[f"F_conv_func.{sch}"], r.get(f"F_conv_func.{sch}.ci_lo"),
                        r.get(f"F_conv_func.{sch}.ci_hi"), fc.name)
                add(run, subset, "delta_func.family_pair", r["delta.family_pair"], r.get("delta.family_pair.ci_lo"),
                    r.get("delta.family_pair.ci_hi"), fc.name)
    dd = R / f"directional_by_class{tag}.csv"
    if dd.exists():
        d = pd.read_csv(dd)
        for subset in ("all", "test_half"):
            r = d[(d["class"] == "all") & (d.subset == subset)]
            if len(r):
                r = r.iloc[0]
                for sch in ("per_pair", "per_family_pair"):
                    add(run, subset, f"F_conv_dir.{sch}", r[f"F_conv_dir_{sch}"], r.get(f"F_conv_dir_{sch}.ci_lo"),
                        r.get(f"F_conv_dir_{sch}.ci_hi"), dd.name)
                add(run, subset, "F_conv_dir.passing_pairs", r.passing, src=dd.name)
ds = R / "descent_2d_summary.csv"
if ds.exists():
    d = pd.read_csv(ds)
    for dname, run in (("O1_tau0.2", "tau0.2"), ("O1_r2_tau0.2", "tau0.2_O1_r2"), ("O1_r2_tau0.3", "tau0.3")):
        for half in ("all", "test", "tune"):
            r = d[(d.dataset == dname) & (d.half == half)]
            if not len(r):
                continue
            r = r.iloc[0]
            sub = {"all": "all", "test": "test_half", "tune": "tune_half"}[half]
            for k in ("M1_ortholog", "M1_paralog", "diff_identity", "diff_identity_x_time", "diff_time"):
                add(run, sub, f"descent.{k}", r[k], r[k + ".ci_lo"], r[k + ".ci_hi"], ds.name)
out = pd.DataFrame(rows)
out.to_csv(R / "headline_halves_tau.csv", index=False)
pairs.write_sidecar(R / "headline_halves_tau.csv", a.snapshot, a.seed, {"stage": "v21_tables", "n_rows": len(out)})
print(out.pivot_table(index="metric", columns=["run", "subset"], values="value").to_string())
