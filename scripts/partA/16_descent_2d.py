#!/usr/bin/env python
"""Part A v2.1 stage 16: descent M1 (prereg part_a.M1) — ortholog vs paralog at matched identity AND divergence
time (2-D strata), test/tune halves, and the tau = 0.3 sensitivity (PI v2.1 items 1 and 5).

Reads the descent pair tables written by scripts/partA_descent/build_descent_pairs.py (not modified here):
  O1      tau 0.2  /tmp/pe_big/partA_descent/<snap>/O1/pairs.parquet          (the verified DESCENT.md run)
  O1_r2   tau 0.2  /tmp/pe_big/partA_descent_r2t02/<snap>/O1/pairs.parquet   (O1_r2 states; like-for-like with 0.3)
  O1_r2   tau 0.3  /tmp/pe_big/partA_descent_t03/<snap>/O1/pairs.parquet
Primary variant qm90, capped sampling (in_cap), as DESCENT.md. Divergence time per species pair from
analyze_descent.divergence_table (TimeTree; same-species = 0). Strata: identity deciles x time bins
{0 (same species), 0.5-100, 100-300, 300-700, 700-1200, 1200-2000, >2000 Mya}.
Standardised difference (ortholog - paralog) over cells where both classes have >= 20 informative observations,
cell weight = total informative observations of the cell (both classes); 200 OG-bootstrap reps
(descent.bootstrap_weights, seed*1000+rep). Reported side by side: identity-only, time-only, identity x time.
Halves: prereg tune/test key (smallest accession among root-HOG members, sha1 % 2), control families forced test.
Wording (PI): "ortholog conjecture not supported for presence", not "orthologs less conserved".
Outputs: descent_2d_strata.csv (cells, primary dataset/all), descent_2d_summary.csv
"""
import importlib.util
import time
from pathlib import Path

import numpy as np
import pandas as pd

from _common import parser, resolve, logger, pairs
from pocketevo import descent, frozen
from pocketevo import homology as hom

ap = parser(__doc__)
ap.add_argument("--n-boot", type=int, default=200)
ap.add_argument("--min-cell", type=int, default=20)
ap.add_argument("--datasets", default=None, help="name=path[,name=path] (default: the O1 / O1_r2 tables)")
ap.add_argument("--split-table", default=None, help="og_id -> split (D20); default: prereg hash on --orthology")
ap.add_argument("--out-prefix", default="descent_2d")
ap.add_argument("--timetree", default=None, help="TimeTree pairwise CSV (default: the frozen species-level table, "
                "pocketevo.frozen 'timetree_O1'; snapshot-invariant)")
ap.add_argument("--species-tree", default=None, help="species tree (default results/sogs/<snap>/O1/web/species_tree.nwk)")
ap.add_argument("--variant", default="qm90", help="descent informative variant (qm90 primary; qm90d11 = D11 symmetric)")
a = resolve(ap.parse_args())
log = logger(a.work / "partA.log")
T0 = time.time()
spec = importlib.util.spec_from_file_location("analyze_descent", pairs.REPO / "scripts/partA_descent/analyze_descent.py")
AD = importlib.util.module_from_spec(spec)
spec.loader.exec_module(AD)
snap = a.snapshot
TT = a.timetree or str(frozen.path("timetree_O1"))      # species-level: not keyed on the atlas snapshot
SPT = a.species_tree or str(pairs.REPO / f"results/sogs/{snap}/O1/web/species_tree.nwk")
div = AD.divergence_table(TT, SPT)
DATASETS = {"O1_tau0.2": f"/tmp/pe_big/partA_descent/{snap}/O1/pairs.parquet",
            "O1_r2_tau0.2": f"/tmp/pe_big/partA_descent_r2t02/{snap}/O1/pairs.parquet",
            "O1_r2_tau0.3": f"/tmp/pe_big/partA_descent_t03/{snap}/O1/pairs.parquet"}
if a.datasets:
    DATASETS = dict(x.split("=", 1) for x in a.datasets.split(","))
V = a.variant
T_EDGES = [0.5, 100, 300, 700, 1200, 2000]
T_LAB = ["0", "0.5-100", "100-300", "300-700", "700-1200", "1200-2000", ">2000"]

forest = hom.parse_orthoxml(str(Path(a.orthology) / "hogs.orthoxml.gz"))
has = forest.gene_parent >= 0
roots = forest.node_root_rank[forest.gene_parent[has]]
mn = pd.DataFrame({"r": roots, "acc": forest.gene_acc[has].astype(str)}).groupby("r").acc.min()
import hashlib
forced = set(a.prereg["tune_test"]["forced_test_half"]["root_hogs_containing"])
forced_roots = set(forest.node_root_rank[forest.gene_parent[has & np.isin(forest.gene_acc.astype(str), list(forced))]])
og_half = {forest.root_names[r]: (1 if (r in forced_roots or int(hashlib.sha1(k.encode()).hexdigest(), 16) % 2) else 0)
           for r, k in mn.items()}
if a.split_table:
    st = pd.read_parquet(a.split_table)
    og_half = {o: (1 if v == "test" else 0 if v == "tune" else -1) for o, v in zip(st.og_id, st.split_o2)}


def load(path):
    cols = ["og_id", "rel", "taxid_a", "taxid_b", "identity", "in_cap", f"n_inf_{V}", f"n_match_{V}"]
    P = pd.read_parquet(path, columns=cols)
    P = P[P.in_cap & (P[f"n_inf_{V}"] > 0) & P.identity.notna()]
    P = P.merge(div[["taxid_a", "taxid_b", "age_mya"]], on=["taxid_a", "taxid_b"], how="left")
    P["ib"] = np.minimum((P.identity * 10).astype(int), 9)
    P["tb"] = np.digitize(P.age_mya.fillna(-1).to_numpy(), T_EDGES)
    P.loc[P.age_mya.isna(), "tb"] = -1
    P["half"] = P.og_id.map(og_half).fillna(-1).astype(int)
    return P


def agg(P):
    """aggregate pairs to (og, identity bin, time bin, cross-species, relation) cells -> s (matches), t (informative)."""
    A = P.assign(cross=(P.taxid_a != P.taxid_b)).groupby(["og_id", "ib", "tb", "cross", "rel"], sort=False)[
        [f"n_match_{V}", f"n_inf_{V}"]].sum().reset_index()
    return A.rename(columns={f"n_match_{V}": "s", f"n_inf_{V}": "t"})


def std_diff(A, keys, w):
    """ortholog - paralog M1 standardised over cells (keys) where both classes have >= min_cell informative."""
    q = A.assign(s=A.s * w, t=A.t * w)
    g = q.groupby(keys + ["rel"])[["s", "t"]].sum().unstack("rel")
    if 1 not in g["t"] or 2 not in g["t"]:
        return np.nan
    ok = (g["t"][1] >= a.min_cell) & (g["t"][2] >= a.min_cell)
    g = g[ok]
    m1_ = g["s"][1] / g["t"][1]; m2_ = g["s"][2] / g["t"][2]
    wt = g["t"][1] + g["t"][2]
    return float(((m1_ - m2_) * wt).sum() / wt.sum()) if wt.sum() > 0 else np.nan


def m1(A, rel, w):
    d = A.rel.to_numpy() == rel
    return float((A.s.to_numpy() * w)[d].sum() / (A.t.to_numpy() * w)[d].sum())


rows = []
for name, path in DATASETS.items():
    if not Path(path).exists():
        log(f"descent 2d: {name} missing ({path})")
        continue
    P0 = load(path)
    for half in ("all", "test", "tune"):
        P = P0 if half == "all" else P0[P0.half == (1 if half == "test" else 0)]
        A = agg(P)
        codes, ogs = pd.factorize(A.og_id)
        tm = (A.tb >= 0).to_numpy()
        cr = tm & A.cross.to_numpy()
        stats = {
            "M1_ortholog": lambda w: m1(A, 1, w), "M1_paralog": lambda w: m1(A, 2, w),
            "diff_identity": lambda w: std_diff(A, ["ib"], w),
            "diff_time": lambda w: std_diff(A[tm], ["tb"], w[tm]),
            "diff_identity_x_time": lambda w: std_diff(A[tm], ["ib", "tb"], w[tm]),
            "diff_identity_x_time_cross_species": lambda w: std_diff(A[cr], ["ib", "tb"], w[cr]),
        }
        rec = {"dataset": name, "half": half, "n_pairs_ortholog": int((P.rel == 1).sum()),
               "n_pairs_paralog": int((P.rel == 2).sum()), "n_ogs": len(ogs)}
        reps = {k: [] for k in stats}
        for r in range(a.n_boot):
            w = descent.bootstrap_weights(codes, len(ogs), r, a.seed)
            for k, f in stats.items():
                reps[k].append(f(w))
        w1 = np.ones(len(A))
        for k, f in stats.items():
            rec[k] = f(w1)
            rec[k + ".ci_lo"], rec[k + ".ci_hi"] = np.nanquantile(np.array(reps[k], float), [0.025, 0.975])
        rows.append(rec)
        log(f"descent 2d {name} {half}: {rec['M1_ortholog']:.4f} / {rec['M1_paralog']:.4f} "
            f"id {rec['diff_identity']:+.4f} time {rec['diff_time']:+.4f} id x time {rec['diff_identity_x_time']:+.4f}")
        if half == "all" and name == list(DATASETS)[0]:
            g = A[tm].groupby(["ib", "tb", "rel"])[["s", "t"]].sum().reset_index()
            g["M1"] = g.s / g.t
            g["identity_bin"] = g.ib.map(lambda x: f"{x / 10:.1f}-{(x + 1) / 10:.1f}")
            g["time_bin_mya"] = g.tb.map(lambda x: T_LAB[x])
            g["relation"] = g.rel.map({1: "ortholog", 2: "paralog"})
            g.drop(columns=["rel"]).rename(columns={"s": "n_match", "t": "n_informative"}).to_csv(
                a.results / f"{a.out_prefix}_strata.csv", index=False)
if not rows:
    raise SystemExit(f"descent 2d: none of the descent pair tables exists ({DATASETS}); build them with "
                     "scripts/partA_descent/build_descent_pairs.py or pass --datasets")
out = pd.DataFrame(rows)
out.to_csv(a.results / f"{a.out_prefix}_summary.csv", index=False)
pairs.write_sidecar(a.results / f"{a.out_prefix}_summary.csv", a.snapshot, a.seed,
                    {"stage": "descent_2d", "datasets": DATASETS, "timetree": TT, "species_tree": SPT, "variant": V, "min_cell": a.min_cell,
                     "n_boot": a.n_boot, "wall_s": time.time() - T0})
log(f"descent 2d done {time.time() - T0:.0f}s")
