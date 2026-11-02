#!/usr/bin/env python
"""Part A: collect headline numbers into results/numbers.json under key "partA" (merge; other keys untouched)
and write results/partA/<snap>/headline.csv. Every value is read from a stage output of this snapshot."""
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from _common import parser, resolve, pairs

_ap = parser(__doc__)
_ap.add_argument("--numbers-key", default="partA")
_ap.add_argument("--numbers-json", default=None,
                 help="merged numbers file to update (default results/numbers.json); 'none' = write only "
                      "<results>/numbers_partA.json. results/numbers.json[partA] is only ever written for the "
                      "calibration snapshot (prereg meta.snapshot_id): other snapshots must use another key/file")
a = resolve(_ap.parse_args())
R = a.results
est = pd.read_csv(R / "estimands_cos.csv")
prim = est[(est.tag == "all") & (est.variant == "pooled") & (est.E == 0.01) & (est.q0kind == "nominal")]


def g(key, df=prim):
    r = df[df.key == key]
    if not len(r):
        return None
    r = r.iloc[0]
    out = {"value": float(r.value)}
    if "ci_lo" in r and pd.notna(r.ci_lo):
        out["ci95"] = [float(r.ci_lo), float(r.ci_hi)]
    return out


num = {"snapshot_id": a.snapshot, "prereg_sha256": a.prereg_sha, "git_commit": pairs.git_commit(), "seed": a.seed,
       "status": "development subset, preliminary",
       "written_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
num["n_sites"] = int(json.loads((R / "sig_cos_summary.json").read_text())["N"])
sig = json.loads((R / "sig_cos_summary.json").read_text())
num["n_sig_pairs_E0.01_pooled"] = sig["n_sig_pairs_E0.01"]["pooled"]
num["knn_saturated_frac"] = sig["saturated_frac"]
r0 = prim.iloc[0]
num["q0_per_pair"] = float(r0.q0)
num["lambda_calibration"] = float(r0["lambda"])
heads = {}
for sch in ("per_family_pair", "per_pair", "per_site"):
    heads[sch] = {k: g(f"{sch}.F_adj.{k}") for k in ("F_anc", "F_offsite", "F_H4", "F_analog")}
    heads[sch]["raw"] = {k: g(f"{sch}.F_raw.{k}") for k in ("F_anc", "F_offsite", "F_H4", "F_analog")}
num["F_headline_per_family_pair"] = heads["per_family_pair"]
num["F_per_pair"] = heads["per_pair"]
num["F_per_site"] = heads["per_site"]
classes = ["H1s", "H1x", "H2s", "H2x", "H12s", "H12x", "H3j", "H3s", "H4", "N", "unknown"]
num["per_class"] = {c: {"K": g(f"per_pair.K.{c}"), "expected_FP": g(f"per_pair.E0.{c}"), "n": g(f"n.{c}"),
                        "r_adj": g(f"per_pair.r_adj.{c}"), "r_raw": g(f"per_pair.r_raw.{c}")} for c in classes}
# sensitivity grid
sens = []
for (tag, v, E, kind), d in est.groupby(["tag", "variant", "E", "q0kind"]):
    row = {"variant": v, "E": E, "q0kind": kind, "n_sig": int(d.n_sig.iloc[0]), "q0": float(d.q0.iloc[0])}
    for sch in ("per_family_pair", "per_pair"):
        for k in ("F_anc", "F_offsite", "F_H4", "F_analog"):
            x = d[d.key == f"{sch}.F_adj.{k}"]
            row[f"{sch}.{k}"] = float(x.value.iloc[0]) if len(x) else None
    sens.append(row)
pd.DataFrame(sens).to_csv(R / "sensitivity_grid.csv", index=False)
num["sensitivity_grid_file"] = str((R / "sensitivity_grid.csv").relative_to(pairs.REPO))
cs = R / "estimands_csls.csv"
if cs.exists():
    e2 = pd.read_csv(cs)
    p2 = e2[(e2.variant == "pooled") & (e2.E == 0.01) & (e2.q0kind == "nominal")]
    num["csls_sensitivity"] = {sch: {k: g(f"{sch}.F_adj.{k}", p2) for k in ("F_anc", "F_offsite", "F_H4", "F_analog")}
                               for sch in ("per_family_pair", "per_pair")}
st = R / "strata_cos.csv"
if st.exists():
    s = pd.read_csv(st)
    num["strata"] = {}
    for tag, d in s.groupby("tag"):
        num["strata"][tag] = {f"{sch}.{k}": g(f"{sch}.F_adj.{k}", d) for sch in ("per_family_pair", "per_pair")
                              for k in ("F_anc", "F_H4", "F_analog")}
        num["strata"][tag]["n_sig"] = int(d.n_sig.iloc[0])
ho = pd.read_csv(R / "null_cos_calibration_heldout.csv")
ex = pd.read_csv(R / "null_cos_calibration_extended.csv")
num["calibration"] = {
    "heldout_ratio_pooled": {f"{p:g}": [float(x), float(y)] for p, x, y in
                             ho[ho.variant == "pooled"][["p_threshold", "ratio_ij", "ratio_ji"]].values},
    "extended_n_pairs": int(ex.n_pairs.iloc[0]),
    "extended_both_dir_E0.01_rate": {v: float(d.iloc[(np.abs(d.E_threshold - 0.01)).argsort()[:1]].obs_both.iloc[0]
                                          / d.n_pairs.iloc[0]) for v, d in ex.groupby("variant")},
}
covm = json.loads((a.work / "null_cos" / "covmodel.json").read_text())["fit"]
num["covariate_model"] = {"r2_covariates": covm["r2_cov"], "r2_size_only": covm["r2_size_only"]}
if (R / "ligand_enrichment.csv").exists():   # v1 stage 9 (superseded by 12_functional in v2)
    lig = pd.read_csv(R / "ligand_enrichment.csv")
    l0 = lig[(lig.set == "sig_N_vs_matched_random_N") & (lig.exclude_pocketcompare_train)].iloc[0]
    num["ligand_enrichment_N_v1"] = {k: (float(l0[k]) if pd.notna(l0[k]) else None) for k in
                                     ("n_test", "shared_test", "shared_control_matched", "enrichment")}
def _clean(records):
    return [{k: (None if (isinstance(v, float) and not np.isfinite(v)) else v) for k, v in r.items()} for r in records]


fc = R / "functional_convergence.csv"
if fc.exists():   # named keys (v2.1 fix: itertuples()._asdict() renamed dotted columns to _13/_14/...)
    num["functional"] = _clean(pd.read_csv(fc).to_dict("records"))
for extra in ("headline_halves_tau.csv", "directional_toy_seeds.csv", "descent_2d_strata.csv",
              "descent_2d_summary.csv"):
    if (R / extra).exists():
        num[extra.replace(".csv", "")] = _clean(pd.read_csv(R / extra).to_dict("records"))
hh = R / "hidden_homology.csv"
if hh.exists():
    h = pd.read_csv(hh)
    num["hidden_homology"] = h.to_dict(orient="records")
    eps = float(h[h.set == "sig_N"].eps_lenient.iloc[0])
    num["F_analog_hidden_homology_corrected"] = {
        sch: {"value": heads[sch]["F_analog"]["value"] * (1 - eps),
              "ci95": [x * (1 - eps) for x in heads[sch]["F_analog"].get("ci95", [np.nan, np.nan])],
              "eps_lenient": eps, "note": "CI scales the bootstrap CI; epsilon uncertainty in hidden_homology"}
        for sch in ("per_family_pair", "per_pair", "per_site")}
dd = R / "directional_by_class.csv"
if dd.exists():
    num["directional"] = _clean(pd.read_csv(dd).to_dict("records"))
    num["directional_summary"] = json.loads((R / "directional_summary.json").read_text())
tk = pd.read_csv(R / "topk_composition.csv")
num["topk_composition"] = tk.to_dict(orient="records")
cl = json.loads((R / "cls_cos_summary.json").read_text())
num["split_mode"] = cl["split_mode"]
num["split_fraction"] = cl["split_fraction"]
num["rest_sample_classes"] = cl["rest_sample_classes"]
ej = json.loads((R / "estimands_cos.json").read_text())
num["exact_totals"] = ej["totals"]

pairs.write_sidecar(R / "numbers_partA.prov", a.snapshot, a.seed, {"stage": "numbers_partA"})
(R / "numbers_partA.json").write_text(json.dumps(num, indent=1, default=float))   # always: per-snapshot copy
path = None if (a.numbers_json or "").lower() == "none" else Path(a.numbers_json or pairs.REPO / "results" / "numbers.json")
from pocketevo import frozen  # noqa: E402
if (path is not None and path.resolve() == (pairs.REPO / "results" / "numbers.json").resolve()
        and a.numbers_key == "partA" and a.snapshot != frozen.calibration_snapshot()):
    raise SystemExit(f"refusing to overwrite results/numbers.json['partA'] with snapshot {a.snapshot} values (the "
                     f"reported numbers are snapshot {frozen.calibration_snapshot()}); {R / 'numbers_partA.json'} "
                     "was written; pass --numbers-json none or another --numbers-key")
if path is not None:
    allnum = json.loads(path.read_text()) if path.exists() else {}
    allnum[a.numbers_key] = num
    path.write_text(json.dumps(allnum, indent=1, default=float))
print(json.dumps({k: num[k] for k in ("F_headline_per_family_pair", "q0_per_pair")}, indent=1, default=float))
