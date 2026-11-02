#!/usr/bin/env python
"""G2 functional slice (5 control families): trace literature binding-site residues of the human reference protein
through the SOG pipeline and summarise states. Functional check only (no tuning).

Known-site residues below are literature reference residues used ONLY for this functional check (tests/slice
input, not pipeline logic): CDK2 ATP site, CA2 zinc/active site, ESR1/AR ligand-binding pocket, LDHA NADH/substrate.
"""
import argparse
import glob
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

P = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(P / "scripts"))
from build_sogs import git_commit  # noqa: E402

SLICE = {  # accession: (name, residues)  -- UniProt numbering
    "P24941": ("CDK2 ATP site", [10, 18, 31, 33, 64, 80, 81, 82, 83, 84, 86, 89, 131, 132, 134, 144, 145]),
    "P00918": ("CA2 active site", [64, 92, 94, 96, 119, 121, 131, 141, 143, 198, 199, 200, 209]),
    "P03372": ("ESR1 ligand pocket", [343, 346, 347, 349, 350, 353, 383, 384, 387, 388, 391, 394, 404, 421, 424, 425, 428, 521, 524, 525, 528]),
    "P10275": ("AR ligand pocket", [701, 704, 705, 707, 711, 742, 745, 746, 749, 752, 764, 780, 787, 873, 877, 891, 895]),
    "P00338": ("LDHA NADH/substrate", [28, 29, 30, 52, 53, 97, 98, 99, 106, 138, 139, 140, 163, 169, 193, 194, 238, 246, 248]),
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sogs", required=True)
    ap.add_argument("--tau", type=float, required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    root = Path(a.sogs)
    M = pd.read_parquet(root / "inputs/members.parquet")
    og = pd.read_parquet(root / "inputs/ogs.parquet")
    sites = pd.read_parquet(root / "inputs/shard_sites", columns=["site_id", "accession", "residues", "og_id"])
    ss = pd.concat([pd.read_parquet(f) for f in glob.glob(str(root / "sogs/sites_shard=*.parquet"))])
    sg = pd.concat([pd.read_parquet(f) for f in glob.glob(str(root / "sogs/shard=*.parquet"))])
    sg = sg[np.isclose(sg.tau, a.tau)]
    stt = pd.concat([pd.read_parquet(f) for f in glob.glob(str(root / f"states_t{a.tau:g}/shard=*.parquet"))])
    stats = pd.concat([pd.read_parquet(f) for f in glob.glob(str(root / "sogs/ogstats_shard=*.parquet"))])
    rows, og_rows = [], []
    for acc, (name, res) in SLICE.items():
        og_id = M[M.accession == acc].og_id.iloc[0]
        o = og[og.og_id == og_id].iloc[0]
        s = sites[sites.accession == acc].copy()
        s["coverage"] = s.residues.map(lambda r: len(set(res) & set(r)) / len(res))
        b = s.sort_values(["coverage", "site_id"], ascending=[False, True]).iloc[0]
        lab = ss.set_index("site_id").loc[b.site_id, f"sog_t{a.tau:g}"]
        sog_id = f"{og_id}:g{lab}"
        g = sg[sg.sog_id == sog_id].iloc[0]
        st = stt[stt.sog_id == sog_id]
        vc = st.state.value_counts().to_dict()
        vc70 = st.state_s70.value_counts().to_dict()
        pspec = st[st.state == "P"].taxid.nunique()
        rows.append({"accession": acc, "site": name, "og_id": og_id, "half": o.half, "forced_test": bool(o.forced_test),
                     "best_site": b.site_id, "site_n_res": len(b.residues), "known_coverage": round(b.coverage, 3),
                     "sog_id": sog_id, "sog_n_sites": int(g.n_sites), "sog_n_species": int(g.n_species),
                     "footprint_size": int(g.footprint_size), "n_P": vc.get("P", 0), "n_P_species": int(pspec),
                     "n_A": vc.get("A", 0), "n_U": vc.get("U", 0), "n_G": vc.get("G", 0), "n_A_s70": vc70.get("A", 0),
                     "U_reasons": "; ".join(f"{k}:{v}" for k, v in st[st.state == "U"].reason.value_counts().items())})
    for og_id in sorted(set(r["og_id"] for r in rows)):
        o = og[og.og_id == og_id].iloc[0]
        t = stats[stats.og_id == og_id].iloc[0]
        st = stt[stt.og_id == og_id]
        sgo = sg[sg.og_id == og_id]
        og_rows.append({"og_id": og_id, "level": o.level, "half": o.half, "n_members": int(o.n_members),
                        "n_species": int(o.n_species), "n_with_model": int(t.n_with_model), "n_seq_mismatch": int(t.n_seq_mismatch),
                        "n_cols": int(t.n_cols), "identity": round(float(t.identity), 3), "n_sites": int(t.n_sites),
                        "n_sogs": int(len(sgo)), "n_sogs_multispecies": int((sgo.n_species >= 2).sum()),
                        "msa_s": round(float(t.msa_s), 2), "graph_s": round(float(t.graph_s), 2),
                        **{f"rows_{k}": int(v) for k, v in st.state.value_counts().items()}})
    R, O = pd.DataFrame(rows), pd.DataFrame(og_rows)
    out = Path(a.out)
    R.to_csv(out / "slice_G2_sites.csv", index=False)
    O.to_csv(out / "slice_G2_ogs.csv", index=False)
    side = {"snapshot_id": root.parent.name, "prereg_sha256": json.load(open(root / "inputs/members.parquet.json"))["prereg_sha256"],
            "git_commit": git_commit(), "seed": 0, "tau": a.tau}
    for f in ["slice_G2_sites.csv", "slice_G2_ogs.csv"]:
        json.dump(side, open(out / (f + ".json"), "w"), indent=1)
    print(O.to_string())
    print(R.drop(columns=["U_reasons"]).to_string())
    print(R[["accession", "U_reasons"]].to_string())


if __name__ == "__main__":
    main()
