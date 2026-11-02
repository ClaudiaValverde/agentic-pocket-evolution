#!/usr/bin/env python
"""PI decision D20: leakage-safe tune/test split for O2 families, and the §8.1 benchmark re-reported on it.

An O2 root OG is
  test       if its O2 half (prereg hash rule on O2 membership) is test AND none of its members was in an O1 TUNE family;
  tune       if its O2 half is tune AND none of its members was in an O1 TEST family;
  crossover  otherwise (members from the other O1 half): excluded from both tune and test reporting.
"O1 family" = an O1 root OG that was run (>= 3 species), i.e. one whose annotations could enter the tau tuning.
Members not in any run O1 OG impose no constraint.

    .venv/bin/python scripts/sog_split_o2.py --o1 /tmp/pe_big/sogs/20261030T085136Z/O1_r2 \\
        --o2 /tmp/pe_big/sogs/20261030T085136Z/O2_r1 --o2-results results/sogs/20261030T085136Z/O2_r1 \\
        --out results/sogs/20261030T085136Z/O2_r1_split
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

P = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(P / "scripts"))
sys.path.insert(0, str(P / "src"))
from build_sogs import git_commit  # noqa: E402
from pocketevo import sog  # noqa: E402
from sog_tau_benchmark import boot_ci, metrics  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--o1", required=True)
    ap.add_argument("--o2", required=True)
    ap.add_argument("--o2-results", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--prereg", default=str(P / "config/prereg.yaml"))
    ap.add_argument("--n-boot", type=int, default=200)
    a = ap.parse_args()
    prm = sog.prereg_params(yaml.safe_load(open(a.prereg)))
    seed = int(prm["seeds"][0])
    o1, o2 = Path(a.o1), Path(a.o2)
    og1 = pd.read_parquet(o1 / "inputs/ogs.parquet")
    m1 = pd.read_parquet(o1 / "inputs/members.parquet", columns=["og_id", "accession"])
    og1 = og1[og1.eligible]
    m1 = m1.merge(og1[["og_id", "half"]], on="og_id")
    half1 = m1.groupby("accession").half.first()
    og2 = pd.read_parquet(o2 / "inputs/ogs.parquet")
    m2 = pd.read_parquet(o2 / "inputs/members.parquet", columns=["og_id", "accession"])
    m2["o1_half"] = m2.accession.map(half1)
    agg = m2.groupby("og_id").o1_half.agg(n_o1_tune=lambda s: int((s == "tune").sum()),
                                          n_o1_test=lambda s: int((s == "test").sum()))
    t = og2[["og_id", "half", "forced_test", "eligible", "n_members", "n_species"]].merge(agg, on="og_id", how="left")
    t[["n_o1_tune", "n_o1_test"]] = t[["n_o1_tune", "n_o1_test"]].fillna(0).astype(int)
    t["split_o2"] = np.where((t.half == "test") & (t.n_o1_tune == 0), "test",
                             np.where((t.half == "tune") & (t.n_o1_test == 0), "tune", "crossover"))
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    t.rename(columns={"half": "half_o2_hash"}).to_parquet(out / "og_split_o2.parquet", index=False)
    el = t[t.eligible]
    counts = {"all_ogs": t.split_o2.value_counts().to_dict(), "eligible_ogs": el.split_o2.value_counts().to_dict(),
              "eligible_by_hash_half": {f"{k[0]}->{k[1]}": int(v) for k, v in el.groupby(["half", "split_o2"]).size().items()},
              "forced_test_crossover": int((el.forced_test & (el.split_o2 == "crossover")).sum())}
    # benchmark re-report on the clean halves (pairs from O2_r1, unchanged SOGs)
    pairs = pd.read_parquet(Path(a.o2_results) / "tau_benchmark_pairs.parquet")
    pairs = pairs.merge(t[["og_id", "split_o2"]], on="og_id")
    rows = []
    for split in ["tune", "test", "crossover"]:
        for subset, msk in [("primary", pairs.primary), ("with_metals", pd.Series(True, index=pairs.index))]:
            d = pairs[msk & (pairs.split_o2 == split)]
            for tau in prm["tau_grid"]:
                r = metrics(d, tau)
                r.update({"split_o2": split, "subset": subset})
                if subset == "primary":
                    r["f1_ci_lo"], r["f1_ci_hi"] = boot_ci(d, tau, a.n_boot, seed)
                rows.append(r)
    R = pd.DataFrame(rows)
    R.to_csv(out / "tau_benchmark_o2_split.csv", index=False)
    counts["benchmark_pairs"] = {f"{k[0]}|{k[1]}": int(v) for k, v in pairs.groupby(["split_o2", "kind"]).size().items()}
    side = {"snapshot_id": o2.parent.name, "prereg_sha256": sog.prereg_sha256(a.prereg), "git_commit": git_commit(),
            "seed": seed, "code_version": git_commit(), "decision": "D20", "o1_run": str(o1), "o2_run": str(o2)}
    json.dump({**side, **counts}, open(out / "split_o2_summary.json", "w"), indent=1, default=int)
    for f in ["og_split_o2.parquet", "tau_benchmark_o2_split.csv"]:
        json.dump(side, open(out / (f + ".json"), "w"), indent=1)
    print(json.dumps(counts, indent=1, default=int))
    print(R[R.subset == "primary"][["split_o2", "tau", "n_pos_detected", "n_neg", "recall", "precision", "f1",
                                    "f1_ci_lo", "f1_ci_hi", "n_ogs"]].round(3).to_string())


if __name__ == "__main__":
    main()
