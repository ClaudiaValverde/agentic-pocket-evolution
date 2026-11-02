#!/usr/bin/env python
"""Machine-readable observation error rates for the simulation agent + rescue breakdowns (PI request).

Inputs (results/observation/<snapshot>/): redetection_sites.parquet, false_a_sites.parquet, rescue_hits.parquet,
rescue_k.json. Outputs:
  error_rates.parquet          long table, one row per (metric, kingdom, taxid, plddt_band, model_source,
                               site_size_class, ligand_class, tau, k): n, rate, Wilson 95% CI
     metric = redetect          P(v6 atlas site re-detected | other source/seed), model_source in
                                {foldcomp_v4_kept, foldcomp_v4_all, seed_s1, seed_s2, seed_s0}
              false_A           P(no atlas site covers >= 50 % | known-present site, A-eligible), model_source afdb_v6
              false_A_after_rescue  false_A x P(< k rescue settings hit | false-A case), per k
              control_rescued   P(>= k settings hit | random non-site surface patch of equal size passing A)
     Rows with taxid = -1 / kingdom = 'all' / plddt_band = 'all' / ... are marginal (pooled) rows.
  rescue_breakdown.csv         false-A after rescue at each k by site size class (<=10 / >10) and ligand class,
                               tune / test / all halves, and control rescue rate at each k (same strata)

  python scripts/observation/error_rates.py
"""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd

P = Path(__file__).resolve().parents[2]
SNAP = "20261030T085136Z"


def wilson(k, n, z=1.96):
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def size_class(n):
    return "<=10" if n <= 10 else ">10"


def grouped(df, flag, dims, fixed):
    """Rates of boolean `flag` for every combination of dims marginalised in/out (value 'all')."""
    rows = []
    for keep in itertools.product([True, False], repeat=len(dims)):
        by = [d for d, k in zip(dims, keep) if k]
        groups = df.groupby(by, observed=True) if by else [((), df)]
        for key, g in groups:
            key = key if isinstance(key, tuple) else (key,)
            rec = {d: "all" for d in dims}
            rec.update(dict(zip(by, key)))
            k, n = int(g[flag].sum()), int(len(g))
            lo, hi = wilson(k, n)
            rows.append({**fixed, **rec, "n": n, "n_events": k, "rate": k / n if n else np.nan, "ci_lo": lo, "ci_hi": hi})
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default=str(P / f"results/observation/{SNAP}"))
    ap.add_argument("--taus", type=float, nargs="+", default=[0.2, 0.3, 0.4, 0.5])
    a = ap.parse_args()
    R = Path(a.dir)
    rows = []
    # --- re-detection
    S = pd.read_parquet(R / "redetection_sites.parquet")
    M = pd.read_parquet("/tmp/pe_big/obs/redetect_models.parquet")[["accession", "taxid"]] \
        if Path("/tmp/pe_big/obs/redetect_models.parquet").exists() else None
    if M is not None:
        S = S.merge(M, on="accession", how="left")
    S["site_size_class"] = S.n_res.map(size_class)
    S["plddt_band"] = S.band
    S["model_source"] = S.comparison
    S["ligand_class"] = "any"
    for tau in a.taus:
        S["ev"] = S.best_jaccard >= tau
        for src, g in S.groupby("model_source"):
            for r in grouped(g, "ev", ["kingdom", "plddt_band", "site_size_class"],
                             dict(metric="redetect", model_source=src, tau=tau, k=-1, ligand_class="any")):
                rows.append(r)
            for r in grouped(g, "ev", ["taxid", "plddt_band"],
                             dict(metric="redetect", model_source=src, tau=tau, k=-1, ligand_class="any",
                                  kingdom="by_taxid", site_size_class="all")):
                rows.append(r)
    # --- false-A (known-present sites)
    K = pd.read_parquet(R / "false_a_sites.parquet")
    K["site_size_class"] = K.residues.map(len).map(size_class)
    K["plddt_band"] = K.band
    kd = {}
    try:
        import sys
        sys.path.insert(0, str(P / "scripts/observation"))
        from sample_redetection import kingdoms
        kd = kingdoms(K.taxid.unique())
    except Exception as e:  # noqa: BLE001
        print("kingdom lookup failed:", e)
    K["kingdom"] = K.taxid.map(kd).fillna("unknown")
    for crit, flag in (("A_primary", "false_A"), ("A_sens", "false_A_sens")):
        E = K[K[crit]].copy()
        E["ev"] = ~E.detected
        for r in grouped(E, "ev", ["kingdom", "plddt_band", "site_size_class", "ligand_class"],
                         dict(metric="false_A", criterion=crit, model_source="afdb_v6", tau=-1, k=0)):
            rows.append(r)
        for r in grouped(E, "ev", ["taxid", "site_size_class"],
                         dict(metric="false_A", criterion=crit, model_source="afdb_v6", tau=-1, k=0,
                              kingdom="by_taxid", plddt_band="all", ligand_class="all")):
            rows.append(r)
    # --- rescue
    brk = []
    hp = R / "rescue_hits.parquet"
    if hp.exists():
        H = pd.read_parquet(hp)
        H = H[H.rule == "jaccard"]
        info = K.set_index("known_id")[["site_size_class", "ligand_class", "kingdom", "taxid", "residues", "half"]]
        H["base_id"] = np.where(H.kind == "control", H.known_id.str[5:], H.known_id)
        H = H.join(info[["site_size_class", "kingdom", "taxid"]], on="base_id")
        H["case_ligand_class"] = H.base_id.map(info.ligand_class)
        kgrid = json.loads((R / "rescue_k.json").read_text())["grid"]
        for crit, flag in (("A_primary", "false_A"), ("A_sens", "false_A_sens")):
            for half in ("tune", "test", "all"):
                Kh = K if half == "all" else K[K.half == half]
                for dims in ([], ["site_size_class"], ["ligand_class"], ["site_size_class", "ligand_class"],
                             ["plddt_band"], ["kingdom"]):
                    groups = Kh.groupby(dims) if dims else [((), Kh)]
                    for key, g in groups:
                        key = key if isinstance(key, tuple) else (key,)
                        elig = g[g[crit]]
                        ids = set(g.loc[g[flag], "known_id"])
                        fa0 = len(ids) / len(elig) if len(elig) else np.nan
                        hc = H[(H.kind == "false_A_case") & H.known_id.isin(ids)]
                        hn = H[(H.kind == "control") & H.base_id.isin(ids)]
                        for k in kgrid:
                            still = float((hc.n_hits < k).mean()) if len(hc) else np.nan
                            rec = dict(criterion=crit, half=half, strata="+".join(dims) or "all",
                                       **{d: v for d, v in zip(dims, key)}, k=k, n_A_eligible=len(elig),
                                       n_false_A=len(ids), false_A_before=fa0, n_cases_rescue_run=len(hc),
                                       frac_cases_still_A=still, false_A_after=fa0 * still,
                                       n_controls=len(hn),
                                       control_rescued=float((hn.n_hits >= k).mean()) if len(hn) else np.nan,
                                       control_survive_as_A=float((hn.n_hits < k).mean()) if len(hn) else np.nan)
                            brk.append(rec)
        B = pd.DataFrame(brk)
        B.to_csv(R / "rescue_breakdown.csv", index=False)
        # long-table rows for the simulation agent (all halves pooled)
        for r in B[(B.half == "all")].itertuples():
            base = dict(criterion=r.criterion, model_source="afdb_v6", tau=-1, k=r.k,
                        kingdom=getattr(r, "kingdom", "all") if isinstance(getattr(r, "kingdom", None), str) else "all",
                        plddt_band=r.plddt_band if isinstance(getattr(r, "plddt_band", None), str) else "all",
                        site_size_class=r.site_size_class if isinstance(getattr(r, "site_size_class", None), str) else "all",
                        ligand_class=r.ligand_class if isinstance(getattr(r, "ligand_class", None), str) else "all",
                        taxid="all")
            rows.append({**base, "metric": "false_A_after_rescue", "n": r.n_A_eligible,
                         "n_events": np.nan, "rate": r.false_A_after, "ci_lo": np.nan, "ci_hi": np.nan,
                         "n_rescue_sample": r.n_cases_rescue_run})
            rows.append({**base, "metric": "control_rescued", "n": r.n_controls, "n_events": np.nan,
                         "rate": r.control_rescued, "ci_lo": np.nan, "ci_hi": np.nan})
    T = pd.DataFrame(rows)
    for c in ("kingdom", "plddt_band", "site_size_class", "ligand_class", "model_source", "criterion", "taxid"):
        T[c] = T[c].astype(str)
    T["snapshot_id"] = SNAP
    T.to_parquet(R / "error_rates.parquet", index=False)
    print(T.groupby("metric").size())
    if brk:
        B = pd.DataFrame(brk)
        print(B[(B.criterion == "A_primary") & (B.half == "all") & B.strata.isin(["all", "site_size_class", "ligand_class"])]
              [["strata", "site_size_class", "ligand_class", "k", "false_A_before", "false_A_after", "control_rescued"]]
              .to_string())


if __name__ == "__main__":
    main()
