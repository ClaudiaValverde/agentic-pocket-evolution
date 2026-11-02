#!/usr/bin/env python
"""Part A v3 stage 18: O2 (v3) vs O1 (v2.1) comparison and class moves caused by real placement.

  * class moves: crosstab of final classes (O1 run vs O2 run) over the union of pairs significant (pooled tail,
    E < 0.01 both directions) in either run; separately for pairs involving >= 1 tier-2 protein (where the O1 run
    used the PANTHER-subfamily proxy) and tier-1-only pairs; group-kind of tier-2 proteins in both runs.
  * one comparison table: per-family-pair F_anc / F_analog (all + test half), r_c, F_conv,func (headline and the
    specific-ligand sensitivity 'tokens_no_metal'), descent M1 ortholog / paralog.
Outputs (results/partA/<snap>/v3/): class_moves.csv, class_moves_summary.json, o2_vs_o1.csv
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

from _common import parser, resolve, pairs

ap = parser(__doc__)
ap.add_argument("--work-o1", default="/tmp/pe_big/partA/20261030T085136Z")
ap.add_argument("--results-o1", default=str(pairs.REPO / "results/partA/20261030T085136Z/v2_1"))
a = resolve(ap.parse_args())
R, R1 = a.results, Path(a.results_o1)
W1 = Path(a.work_o1)
cols = ["i", "j", "cls", "Emax_pooled"]
p1 = pd.read_parquet(W1 / "cls_cos" / "pairs_cls.parquet", columns=cols)
p2 = pd.read_parquet(a.work / "cls_cos" / "pairs_cls.parquet", columns=cols)
p1, p2 = p1[p1.Emax_pooled < 0.01], p2[p2.Emax_pooled < 0.01]
m = p1.merge(p2, on=["i", "j"], how="outer", suffixes=("_O1", "_O2"))
m["cls_O1"] = m.cls_O1.fillna("not_sig")
m["cls_O2"] = m.cls_O2.fillna("not_sig")
uni = pd.read_parquet(a.work / "universe.parquet", columns=["tier"])
tier = uni.tier.to_numpy()
m["has_tier2"] = (tier[m.i.values] == 2) | (tier[m.j.values] == 2)
rows = []
for sub, d in (("all", m), ("pairs_with_tier2", m[m.has_tier2]), ("tier1_only", m[~m.has_tier2])):
    ct = pd.crosstab(d.cls_O1, d.cls_O2)
    for c1 in ct.index:
        for c2 in ct.columns:
            if ct.loc[c1, c2]:
                rows.append({"subset": sub, "cls_O1": c1, "cls_O2": c2, "n": int(ct.loc[c1, c2])})
cm = pd.DataFrame(rows)
cm.to_csv(R / "class_moves.csv", index=False)
both = m[(m.cls_O1 != "not_sig") & (m.cls_O2 != "not_sig")]
summ = {
    "sig_pairs_O1": int(len(p1)), "sig_pairs_O2": int(len(p2)), "sig_both": int(len(both)),
    "moved_class_among_sig_both": int((both.cls_O1 != both.cls_O2).sum()),
    "moved_class_among_sig_both_with_tier2": int(((both.cls_O1 != both.cls_O2) & both.has_tier2).sum()),
    "O1_H12_with_tier2_sig_both": int(((both.cls_O1.str.startswith("H12")) & both.has_tier2).sum()),
    "O1_H12_with_tier2_to": both[both.cls_O1.str.startswith("H12") & both.has_tier2].cls_O2.value_counts().to_dict(),
    "into_O2_H12_from": both[both.cls_O2.str.startswith("H12")].cls_O1.value_counts().to_dict(),
}
pr1 = pd.read_parquet(W1 / "proteins.parquet", columns=["accession", "tier", "group_kind"])
pr2 = pd.read_parquet(a.work / "proteins.parquet", columns=["accession", "tier", "group_kind"])
summ["group_kind_O1"] = pr1.groupby(["tier", "group_kind"]).size().astype(int).to_dict().__repr__()
summ["group_kind_O2"] = pr2.groupby(["tier", "group_kind"]).size().astype(int).to_dict().__repr__()
(R / "class_moves_summary.json").write_text(json.dumps(summ, indent=1, default=str))

# comparison table
out = []


def est_rows(Rdir, run, est_name="estimands_cos.csv", strata_name="strata_cos.csv"):
    e = pd.read_csv(Rdir / est_name)
    e = e[(e.tag == "all") & (e.variant == "pooled") & (e.E == 0.01) & (e.q0kind == "nominal")]
    s = pd.read_csv(Rdir / strata_name)
    for subset, d in (("all", e), ("test_half", s[s.tag == "test_half"]), ("plddt_ge90", s[s.tag == "both_plddt_ge90"])):
        for key in ["per_family_pair.F_adj.F_anc", "per_family_pair.F_adj.F_analog", "per_pair.F_adj.F_anc",
                    "per_pair.F_adj.F_analog"] + [f"per_pair.r_adj.{c}" for c in
                                                  ("H1s", "H1x", "H2s", "H2x", "H12s", "H12x", "H3j", "H3s", "H4", "N")] \
                + [f"n.{c}" for c in ("H1s", "H2s", "H12s", "H12x", "H3j", "N")] + \
                [f"per_pair.K.{c}" for c in ("H1s", "H2s", "H12s", "H12x", "H3j", "H4", "N")]:
            r = d[d.key == key]
            if len(r):
                out.append({"run": run, "subset": subset, "metric": key, "value": r.value.iloc[0],
                            "ci_lo": r.ci_lo.iloc[0] if "ci_lo" in r else np.nan,
                            "ci_hi": r.ci_hi.iloc[0] if "ci_hi" in r else np.nan})


def func_rows(Rdir, run):
    f = pd.read_csv(Rdir / "functional_convergence.csv")
    f = f[f.exclude_train & (f["set"] == "sig_N_vs_matched_random_N") & (f.category == "any")]
    for subset in ("all", "test_half", "plddt_ge90"):
        for var in ("full", "tokens_no_metal"):
            r = f[(f.subset == subset) & (f.label_variant == var)]
            if len(r):
                r = r.iloc[0]
                out.append({"run": run, "subset": subset, "metric": f"F_conv_func.per_family_pair[{var}]",
                            "value": r["F_conv_func.per_family_pair"], "ci_lo": r.get("F_conv_func.per_family_pair.ci_lo"),
                            "ci_hi": r.get("F_conv_func.per_family_pair.ci_hi")})
                out.append({"run": run, "subset": subset, "metric": f"functional_n_test[{var}]", "value": r.n_test})


def descent_rows(path, run, ds):
    if not Path(path).exists():
        return
    d = pd.read_csv(path)
    d = d[d.dataset == ds]
    for r in d.itertuples(index=False):
        rr = r._asdict()
        sub = {"all": "all", "test": "test_half", "tune": "tune_half"}[rr["half"]]
        for k in ("M1_ortholog", "M1_paralog", "diff_identity", "diff_identity_x_time"):
            out.append({"run": run, "subset": sub, "metric": f"descent.{k}", "value": rr[k]})


est_rows(R1, "O1_v2.1")
est_rows(R, "O2_v3")
func_rows(R1, "O1_v2.1")
func_rows(R, "O2_v3")
descent_rows(R1 / "descent_2d_summary.csv", "O1_v2.1", "O1_tau0.2")
descent_rows(R / "descent_O2_summary.csv", "O2_v3", "O2_tau0.2")
o = pd.DataFrame(out)
o.to_csv(R / "o2_vs_o1.csv", index=False)
pairs.write_sidecar(R / "o2_vs_o1.csv", a.snapshot, a.seed, {"stage": "o2_vs_o1", **{k: v for k, v in summ.items()
                                                                                    if isinstance(v, int)}})
pd.set_option("display.width", 250)
print(json.dumps(summ, indent=1, default=str))
print(o.pivot_table(index="metric", columns=["subset", "run"], values="value").to_string())
