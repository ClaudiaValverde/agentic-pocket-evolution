#!/usr/bin/env python
"""D28 co-fold pool: combine the per-orthology case lists (build_cases.py --pool eligible), deduplicate events that
involve the same proteins, then apply prereg case_selection (top-k by S per category + n seeded audit) to the pool.

Pool = filter_label in {convincing, supported} AND full artifact cascade pass (pass_artifacts incl. orthology_split,
plus the M4 orthology_side_family check) AND not exploratory. Specificity rows that are exploratory (not BH-significant)
are kept as audit rows only. Negative controls are carried over from each list. No backfill with failing cases.

Dedup rule: two rows of the same category are duplicates if their present-side protein sets overlap with
Jaccard >= 0.5. The row with the higher S is kept; on ties the later orthology in --lists (O2) is kept.

  .venv/bin/python scripts/cases/combine_pool.py \
      --lists results/cases/<snap>/stage_O1.parquet results/cases/<snap>/stage_O2.parquet \
      --stage-dirs /tmp/pe_big/cases/stage/O1 /tmp/pe_big/cases/stage/O2 --out results/cases/<snap>/case_list_v2.parquet
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from pocketevo import evidence as E  # noqa: E402
from pocketevo import recon  # noqa: E402

DEDUP_JACCARD = 0.5


def accset(x) -> frozenset:
    return frozenset(a for a in str(x or "").split(";") if a and a != "None" and a != "nan")


def dedup(pool: pd.DataFrame) -> pd.DataFrame:
    """Greedy: sort by S desc, orthology rank desc (later list wins ties), case_id; drop rows whose present set
    overlaps a kept row of the same category with Jaccard >= DEDUP_JACCARD."""
    pool = pool.sort_values(["S", "orth_rank", "case_id"], ascending=[False, False, True])
    kept, dup_of = [], {}
    for i, r in pool.iterrows():
        a = accset(r.present_accs)
        hit = None
        for j in kept:
            k = pool.loc[j]
            if k.category != r.category:
                continue
            b = accset(k.present_accs)
            if a and b and len(a & b) / len(a | b) >= DEDUP_JACCARD:
                hit = k.case_id
                break
        if hit is None:
            kept.append(i)
        else:
            dup_of[r.case_id] = hit
    out = pool.copy()
    out["duplicate_of"] = out.case_id.map(dup_of)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--lists", nargs="+", required=True)
    ap.add_argument("--stage-dirs", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--cases-dir", default="cases")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--snapshot", default="20261030T085136Z")
    ap.add_argument("--candidates", nargs="*", default=[], help="orth=path of Part B candidate tables (null FDR)")
    a = ap.parse_args()
    t0 = time.time()
    rules = E.rules_from_prereg(recon.load_prereg())
    parts = []
    for k, (f, sd) in enumerate(zip(a.lists, a.stage_dirs)):
        d = pd.read_parquet(f)
        d["orth_rank"], d["source_list"], d["stage_dir"] = k, str(f), str(sd)
        parts.append(d)
    df = pd.concat(parts, ignore_index=True)
    elig = df.selection.isin(["eligible"])
    pool = dedup(df[elig])
    uniq = pool[pool.duplicate_of.isna()]
    sel = E.select_cases(uniq, rules, a.seed, bar=("convincing", "supported"), label_col="filter_label")
    df = df.merge(pool[["case_id", "duplicate_of"]], on="case_id", how="left")
    df = df.merge(sel[["case_id", "selected", "selection"]].rename(columns={"selected": "pool_selected",
                                                                          "selection": "pool_selection"}),
                  on="case_id", how="left")
    df["pool_selected"] = df.pool_selected.fillna(False).astype(bool)
    df["pool_selection"] = df.pool_selection.fillna("")
    keep_other = df.selection.isin(["negative_control", "audit_exploratory", "exploratory_duplication"])
    df.loc[keep_other, "pool_selected"] = True
    df.loc[keep_other, "pool_selection"] = df.loc[keep_other, "selection"]
    df["selected"], df["selection"] = df.pool_selected, df.pool_selection
    df = df.drop(columns=["pool_selected", "pool_selection"])
    out = ROOT / a.out if not Path(a.out).is_absolute() else Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    for c in df.columns:
        if df[c].dtype == object:
            df[c] = df[c].map(lambda v: v if v is None or isinstance(v, str) else str(v))
    df.to_parquet(out, index=False)
    # copy selected case dirs from staging, stamping the final selection
    n = 0
    for r in df[df.selected].itertuples():
        src = Path(r.stage_dir) / r.case_id
        if not (src / "case.yaml").exists():
            continue
        dst = (ROOT / a.cases_dir) / r.case_id
        dst.mkdir(parents=True, exist_ok=True)
        c = yaml.safe_load((src / "case.yaml").read_text())
        c["notes"] = (c.get("notes") or "") + f" Final co-fold pool selection (D28, {out.name}): {r.selection}."
        c.setdefault("provenance", {})["case_list"] = str(out.relative_to(ROOT)) if out.is_relative_to(ROOT) else str(out)
        (dst / "case.yaml").write_text(yaml.safe_dump(c, sort_keys=False, allow_unicode=True, width=120))
        (dst / "evidence").mkdir(exist_ok=True)
        if (src / "evidence/evidence_summary.json").exists():
            shutil.copy(src / "evidence/evidence_summary.json", dst / "evidence/evidence_summary.json")
        n += 1
    counts = {}
    for (orth, cat), g in df.groupby(["orthology", "category"]):
        counts[f"{orth}:{cat}"] = dict(
            n=int(len(g)), f1=int((g.f1 == "pass").sum()), f2=int((g.f2 == "pass").sum()),
            f3=int((g.f3 == "pass").sum()), artifacts_pass=int((g.artifacts == "pass").sum()),
            supported_or_convincing=int(g.filter_label.isin(["supported", "convincing"]).sum()),
            convincing_headline=int(((g.filter_label == "convincing") & (g.artifacts == "pass")).sum()),
            convincing_with_literature_post_hoc=int(((g.get("filter_label_with_literature") == "convincing")
                                                     & (g.artifacts == "pass")).sum()),
            eligible=int((g.selection != "").sum() if False else g.case_id.isin(pool.case_id).sum()),
            duplicates=int(g.duplicate_of.notna().sum()), selected=int(g.selected.sum()))
    fdr = {}
    for spec_ in a.candidates:
        orth, path = spec_.split("=", 1)
        c = pd.read_parquet(path)
        c = c[c["primary"]] if "primary" in c else c
        for typ, g in c.groupby("type"):
            real = (float(g["real_called_same_cut_prerescue"].iloc[0]) if "real_called_same_cut_prerescue" in g
                    else float(len(g)))
            fdr[f"{orth}:{typ}"] = {k: round(min(1.0, float(g[f"null_{k}_called_mean_same_cut"].iloc[0]) / real), 3)
                                   for k in ("curveball", "pseudo_clade") if f"null_{k}_called_mean_same_cut" in g}
            fdr[f"{orth}:{typ}"]["real_called"] = real
    md = [f"# case_list_v2 (D25c/D28 co-fold pool), snapshot {a.snapshot}", "",
          "Headline: convincing_only (all three filters AND artifact cascade pass). Losses are not FDR-calibrated (D21).",
          "", "## Empirical null FDR per category (mean null calls / real calls at the same cut, Part B)", "",
          "| orthology:type | curveball | pseudo-clade | real called |", "|---|---|---|---|"]
    md += [f"| {k} | {v.get('curveball')} | {v.get('pseudo_clade')} | {v['real_called']:.0f} |" for k, v in fdr.items()]
    md += ["", "Category `duplication_gain_exploratory`: not preregistered; null is selection-conditioned; never in "
           "headline numbers.", "", "## Selected", "", "| case_id | orthology | category | selection | S | filter_label | f1 | f2 | f3 | f2+lit |",
           "|---|---|---|---|---|---|---|---|---|---|"]
    for r in df[df.selected].sort_values(["category", "selection", "S"], ascending=[True, True, False]).itertuples():
        md.append(f"| {r.case_id} | {r.orthology} | {r.category} | {r.selection} | {r.S:g} | {r.filter_label} | {r.f1} | "
                  f"{r.f2} | {r.f3} | {getattr(r, 'f2_with_literature', '')} |")
    md += ["", "## Counts", "", "```", json.dumps(counts, indent=1), "```"]
    (out.with_suffix(".md")).write_text("\n".join(md) + "\n")
    recon.write_sidecar(out, a.snapshot, a.seed, {"empirical_null_fdr": fdr,
        "script": "scripts/cases/combine_pool.py", "argv": sys.argv[1:],
        "inputs": {f: recon.file_sha256(f) for f in a.lists},
        "rule": "D25c/D28: pool = filter_label in {convincing, supported} AND artifacts pass AND not exploratory; "
                f"dedup by present-side protein Jaccard >= {DEDUP_JACCARD} (higher S kept, tie -> later list); "
                f"prereg top_k={rules.top_k} + audit={rules.n_audit} per category, seed {a.seed}; no backfill",
        "headline": "convincing_only (convincing and artifacts pass)", "counts": counts,
        "n_case_dirs_written": n, "code_sha256": recon.code_sha256(__file__, E.__file__),
        "runtime_s": round(time.time() - t0, 1)})
    print(json.dumps(counts, indent=1))
    print(df[df.selected][["case_id", "orthology", "category", "selection", "S", "filter_label"]].to_string())


if __name__ == "__main__":
    main()
