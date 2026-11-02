#!/usr/bin/env python
"""Part A v2 stage 13: cross-fitted directional convergence test (DESIGN v2 §4.2; src/pocketevo/directional.py).

Tested: significant pairs (pooled, E < 0.01) of classes N, H3j, H3s, H4; each orientation (a, b) in which b's site
belongs to an O1 SOG (tau from the SOG run) with >= 3 species is tested (a pair passes if any orientation passes;
BH over all tested orientations). Similarity is the covariate-standardised z of the primary null.
Nulls: (i) significant same-SOG ortholog pairs (H1s; a's protein excluded from the member sets), (ii) random
non-homologous partners a' of tested b sites. Toy calibration: tests/test_partA.py::test_directional_* (no-convergence
world reads 0, planted convergence detected).
F_conv,dir = passing / testable pairs (per pair; per family pair also), two-way bootstrap over root groups.
Outputs: results/partA/<snap>/<run>/directional_summary.json, directional_by_class.csv,
         convergence_candidates.parquet (passing pairs with proteins, sites, SOGs, similarity, delta, ligand evidence)
"""
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from _common import parser, resolve, logger, pairs
from pocketevo import directional as dr
from pocketevo import homology as hom
from pocketevo import nulls
from pocketevo.estimands import pair_group_weight

ap = parser(__doc__)
ap.add_argument("--sogs-dir", default=str(pairs.REPO / "results/sogs/20261030T085136Z/O1"))
ap.add_argument("--tree", default=str(pairs.REPO / "data/orthology/ancestra_tier1_v0/tier1_guide_tree.nwk"))
ap.add_argument("--tau", type=float, default=None, help="default: <sogs-dir>/tau_choice.json")
ap.add_argument("--n-null", type=int, default=300_000)
ap.add_argument("--n-boot", type=int, default=200)
ap.add_argument("--batch", type=int, default=8192)
a = resolve(ap.parse_args())
log = logger(a.work / "partA.log")
dev = torch.device("cuda:0")
T0 = time.time()
summ = {}
tau = a.tau if a.tau is not None else float(json.loads((Path(a.sogs_dir) / "tau_choice.json").read_text())["tau"])
summ["tau"] = tau
summ["sogs_dir"] = a.sogs_dir
uni = pd.read_parquet(a.work / "universe.parquet", columns=["site_id", "accession", "taxid", "tier"])
N = len(uni)
ss = pd.read_parquet(Path(a.sogs_dir) / "site_sog_all_tau.parquet", columns=["site_id", "og_id", f"sog_t{tau:g}"])
sg = pd.read_parquet(Path(a.sogs_dir) / "sogs_all_tau.parquet", columns=["og_id", "tau", "sog_index", "n_species"])
sg = sg[np.isclose(sg.tau, tau)]
sg_fp = pd.read_parquet(Path(a.sogs_dir) / "sogs_all_tau.parquet", columns=["og_id", "tau", "sog_index", "footprint_size"])
sg_fp = sg_fp[np.isclose(sg_fp.tau, tau)]
sg_fp["sog"] = sg_fp.og_id + ":g" + sg_fp.sog_index.astype(str)
ss["sog"] = ss.og_id + ":g" + ss[f"sog_t{tau:g}"].astype(str)
nsp = dict(zip(sg.og_id + ":g" + sg.sog_index.astype(str), sg.n_species))
si_of = pd.Series(np.arange(N), index=uni.site_id.values)
ss["si"] = si_of.reindex(ss.site_id.values).to_numpy()
ss = ss.dropna(subset=["si"])
ss["si"] = ss.si.astype(np.int64)
ss = ss[ss.sog.map(nsp).fillna(0) >= 3]
site_sog = np.full(N, -1, np.int64)
codes, sog_names = pd.factorize(ss.sog)
site_sog[ss.si.values] = codes
members = pd.Series(ss.si.values).groupby(codes).apply(np.array).to_dict()
tree = dr.SpeciesTree(Path(a.tree).read_text(), leaf_key=lambda nm: int(nm.split("_")[-1]))
taxid = uni.taxid.to_numpy()
leaf = np.array([tree.leaf_of.get(int(t), -1) for t in taxid])
summ["sites_in_multispecies_sogs"] = int((site_sog >= 0).sum())
summ["sog_sites_without_tree_leaf"] = int(((site_sog >= 0) & (leaf < 0)).sum())

c = np.load(a.work / "context.npz")
site_prot = c["site_prot"]
site_group = c["prot_group"][site_prot]
model = nulls.CovModel.from_dict(json.loads((a.work / "null_cos" / "covmodel.json").read_text()))
cov = np.load(a.work / "covariates.npy")
Fn = torch.from_numpy(model.norm(np.where(np.isfinite(cov), cov, np.nanmean(cov, 0))).astype(np.float32)).to(dev)
X = np.load(a.work / "fp.f16.npy", mmap_mode="r")
Xg = torch.empty((N, X.shape[1]), dtype=torch.float16, device=dev)
for s in range(0, N, 500_000):
    Xg[s:s + 500_000] = torch.from_numpy(np.ascontiguousarray(X[s:s + 500_000]))

set_cache: dict = {}


def sets_for(b, excl_prot=-1):
    key = (b, excl_prot)
    if key in set_cache:
        return set_cache[key]
    res = None
    if site_sog[b] >= 0 and leaf[b] >= 0:
        mem = members[site_sog[b]]
        mem = mem[(mem != b) & (leaf[mem] >= 0) & (site_prot[mem] != site_prot[b]) & (site_prot[mem] != excl_prot)]
        rng = np.random.default_rng([a.seed, int(b), int(excl_prot) + 1])
        sp = dr.derived_split(tree, leaf[b], leaf[mem], rng) if len(mem) else None
        if sp is not None:
            der, out = mem[sp[0]], mem[sp[1]]
            sel, tst = dr.halves(der, rng)
            res = (dr.cap(np.r_[sel, b], rng), dr.cap(tst, rng), dr.cap(out, rng))
            if len(res[1]) < 2:
                res = None
    if excl_prot < 0:
        set_cache[key] = res
    return res


def pad(lists, m=dr.MAX_SET + 1):
    arr = np.full((len(lists), m), -1, np.int64)
    for k, l in enumerate(lists):
        arr[k, :len(l)] = l
    return arr


def mean_z(av, idx):
    """av (B,) site idx of a; idx (B, M) member idx (-1 pad) -> mean z over valid."""
    va = idx >= 0
    ii = torch.where(va, idx, torch.zeros_like(idx))
    s = (Xg[ii].float() * Xg[av][:, None, :].float()).sum(-1)
    fa = Fn[av][:, None, :].expand(-1, idx.shape[1], -1)
    mu, ls = model.mu_logsd_pairs(fa.reshape(-1, fa.shape[-1]), Fn[ii].reshape(-1, Fn.shape[1]))
    z = model.z(s.reshape(-1), mu, ls).reshape(s.shape)
    z = torch.where(va, z, torch.zeros_like(z))
    return z.sum(1) / va.sum(1).clamp(min=1)


def deltas(A, B, excl=None):
    """orientations (A[k] -> B[k]) -> (delta_test, delta_sel, testable mask)."""
    out_t = np.full(len(A), np.nan); out_s = np.full(len(A), np.nan)
    sets = [sets_for(int(b), int(excl[k]) if excl is not None else -1) for k, b in enumerate(B)]
    ok = np.array([x is not None for x in sets])
    idx = np.flatnonzero(ok)
    for s0 in range(0, len(idx), a.batch):
        kk = idx[s0:s0 + a.batch]
        av = torch.from_numpy(A[kk].astype(np.int64)).to(dev)
        S_, T_, O_ = (torch.from_numpy(pad([sets[k][q] for k in kk])).to(dev) for q in (0, 1, 2))
        zo = mean_z(av, O_)
        out_t[kk] = (mean_z(av, T_) - zo).cpu().numpy()
        out_s[kk] = (mean_z(av, S_) - zo).cpu().numpy()
    return out_t, out_s, ok


sig = pd.read_parquet(a.work / f"cls_cos{a.tag}" / "pairs_cls.parquet", columns=["i", "j", "s", "z", "cls", "Emax_pooled",
                                                                         "gi", "gj"])
sig = sig[sig.Emax_pooled < 0.01].reset_index(drop=True)
TESTED = ["N", "H3j", "H3s", "H4"]
cand = sig[sig.cls.isin(TESTED)].reset_index(drop=True)
oa = np.r_[cand.i.values, cand.j.values].astype(np.int64)
ob = np.r_[cand.j.values, cand.i.values].astype(np.int64)
opair = np.r_[np.arange(len(cand)), np.arange(len(cand))]
has = (site_sog[ob] >= 0) & (leaf[ob] >= 0)
oa, ob, opair = oa[has], ob[has], opair[has]
log(f"directional: {len(cand)} candidate pairs, {len(oa)} orientations with a multi-species SOG on the b side")
dt, ds, ok = deltas(oa, ob)
summ["orientations_with_sog"] = int(len(oa))
summ["orientations_testable"] = int(ok.sum())
log(f"testable orientations {int(ok.sum())} ({time.time() - T0:.0f}s)")

# nulls
rng = np.random.default_rng(a.seed + 40)
h1 = sig[sig.cls == "H1s"]
h1 = h1.iloc[rng.permutation(len(h1))[:a.n_null]]
na = np.r_[h1.i.values, h1.j.values].astype(np.int64)
nbb = np.r_[h1.j.values, h1.i.values].astype(np.int64)
n1, _, ok1 = deltas(na, nbb, excl=site_prot[na])
null1 = n1[ok1]
ctx = hom.PairContext(site_prot, c["site_lab"], c["prot_group"], np.full(len(c["prot_group"]), -1), None)
tb = ob[ok]
pick = tb[rng.integers(0, len(tb), a.n_null)]
ra = rng.integers(0, N, a.n_null)
keep = hom.classify_pairs(ctx, ra, pick) == hom.REST
n2, _, ok2 = deltas(ra[keep], pick[keep])
null2 = n2[ok2]
summ["null_ortholog_n"] = int(len(null1)); summ["null_random_n"] = int(len(null2))
summ["null_ortholog_mean_sd_q95"] = [float(np.mean(null1)), float(np.std(null1)), float(np.quantile(null1, 0.95))]
summ["null_random_mean_sd_q95"] = [float(np.mean(null2)), float(np.std(null2)), float(np.quantile(null2, 0.95))]
log("nulls", json.dumps({k: summ[k] for k in summ if k.startswith("null")}))

tested = np.flatnonzero(ok)
passed_t, p_t, q95 = dr.decide(dt[tested], ds[tested], [null1, null2])
summ["q95_threshold"] = float(q95)
summ["selected_orientations"] = int((ds[tested] > 0).sum())
po = np.zeros(len(oa), bool); po[tested] = passed_t
dtt = dt[tested]
summ["delta_tested_mean_sd_q50_q95_q99"] = [float(np.mean(dtt)), float(np.std(dtt))] + \
    [float(x) for x in np.quantile(dtt, [0.5, 0.95, 0.99])]
summ["frac_selected_and_delta_gt_q95"] = float(((ds[tested] > 0) & (dtt > q95)).mean())
pr, _, q95r = dr.decide(dtt, ds[tested], [null2])
summ["sensitivity_random_null_only"] = {"q95": float(q95r), "passing_orientations": int(pr.sum())}
po_r = np.zeros(len(oa), bool); po_r[tested] = pr
pv0 = np.full(len(oa), np.nan); pv0[tested] = p_t
top = np.zeros(len(oa), bool)
sel_idx = tested[ds[tested] > 0]
top[sel_idx[np.argsort(pv0[sel_idx])[:500]]] = True       # top 500 selected orientations by p (flagged)
report = po | po_r | top
pv = np.full(len(oa), np.nan); pv[tested] = p_t
# pair level
testable_pair = np.zeros(len(cand), bool); testable_pair[opair[ok]] = True
pass_pair = np.zeros(len(cand), bool); np.logical_or.at(pass_pair, opair[po], True)
cand["testable"], cand["passed"] = testable_pair, pass_pair
summ["testable_pairs"] = int(testable_pair.sum())
summ["testable_frac_of_tested_classes"] = float(testable_pair.mean())
summ["passing_pairs"] = int(pass_pair.sum())

# F_conv,dir with two-way bootstrap (same Poisson seeds as 07)
G = int(c["prot_group"].max()) + 1
T = cand[cand.testable]


def f_dir(df, w=None, fam=False):
    pw = np.ones(len(df)) if w is None else pair_group_weight(w, df.gi.values, df.gj.values)
    if fam:   # each root-group cell weight 1
        cell = np.minimum(df.gi.values, df.gj.values).astype(np.int64) * (1 << 31) + np.maximum(df.gi.values, df.gj.values)
        cnt = pd.Series(cell).map(pd.Series(cell).value_counts()).to_numpy()
        pw = pw / cnt
    return float(np.sum(pw * df.passed.values) / np.sum(pw)) if np.sum(pw) > 0 else np.nan


site_half = np.load(a.work / "site_half.npy") if (a.work / "site_half.npy").exists() else None
rows = []
specs = [(cl, "all") for cl in ["all"] + TESTED]
if site_half is not None:
    specs.append(("all", "test_half"))
for cl, subset in specs:
    d = T if cl == "all" else T[T.cls == cl]
    if subset == "test_half":
        d = d[(site_half[d.i.values] == 1) & (site_half[d.j.values] == 1)]
    r = {"class": cl, "subset": subset, "testable_pairs": len(d), "passing": int(d.passed.sum()),
         "F_conv_dir_per_pair": f_dir(d), "F_conv_dir_per_family_pair": f_dir(d, fam=True)}
    bp, bf = [], []
    for k in range(a.n_boot):
        w = np.random.default_rng(a.seed * 1000 + k).poisson(1.0, G).astype(float)
        bp.append(f_dir(d, w)); bf.append(f_dir(d, w, fam=True))
    if len(d):
        r["F_conv_dir_per_pair.ci_lo"], r["F_conv_dir_per_pair.ci_hi"] = np.nanquantile(bp, [0.025, 0.975])
        r["F_conv_dir_per_family_pair.ci_lo"], r["F_conv_dir_per_family_pair.ci_hi"] = np.nanquantile(bf, [0.025, 0.975])
    r["testable_of_class"] = int(len(cand[cand.cls == cl])) if cl != "all" else len(cand)
    rows.append(r)
byc = pd.DataFrame(rows)
byc.to_csv(a.results / f"directional_by_class{a.tag}.csv", index=False)
log("directional\n" + byc.to_string())

# candidates
P = pd.DataFrame({"pair": opair[report], "a": oa[report], "b": ob[report], "delta": dt[report],
                  "delta_sel": ds[report], "p": pv[report], "passed_primary": po[report],
                  "passed_random_null_only": po_r[report]})
if len(P):
    from pocketevo.directional import bh  # noqa: F401
    P = P.merge(cand[["i", "j", "s", "z", "cls", "Emax_pooled"]].reset_index().rename(columns={"index": "pair"}), on="pair")
    P["site_a"], P["site_b"] = uni.site_id.values[P.a.values], uni.site_id.values[P.b.values]
    P["acc_a"], P["acc_b"] = uni.accession.values[P.a.values], uni.accession.values[P.b.values]
    P["sog_b"] = np.array(sog_names)[site_sog[P.b.values]]
    P["sog_a"] = [sog_names[x] if x >= 0 else None for x in site_sog[P.a.values]]
    fc = a.work / "functional_site_labels.parquet"
    if fc.exists():
        lab = pd.read_parquet(fc).set_index("si").tokens
        P["ligands_a"] = lab.reindex(P.a.values).values
        P["ligands_b"] = lab.reindex(P.b.values).values
        isl = lambda x: isinstance(x, (list, np.ndarray))
        P["shared_ligand_token"] = [bool(set(x) & set(y)) if isl(x) and isl(y) else None
                                    for x, y in zip(P.ligands_a, P.ligands_b)]
    # artifact-like flags (rule declared in v2.1 before inspecting new runs; thresholds from prereg A_primary pLDDT
    # and the atlas filters where they exist, else stated here): any of
    #   min site volume < 100 A^3 | volume ratio > 5 | pLDDT_rep < 80 or frac<70 > 0.2 on either side |
    #   SOG footprint <= 6 columns (b side) | either site < 8 residues | burrow score < 0.1 on either side
    U = pd.read_parquet(a.work / "universe.parquet", columns=["rep_site_volume_A3", "plddt_mean_rep",
                                                              "plddt_frac_lt70_max", "n_res", "max_burrow_score"])
    fp_size = dict(zip(sg_fp.sog, sg_fp.footprint_size))
    for side, col in (("a", "a"), ("b", "b")):
        idx = P[col].values
        for k in U.columns:
            P[f"{k}_{side}"] = U[k].values[idx]
    P["footprint_size_b"] = P.sog_b.map(fp_size).values
    vmin = np.minimum(P.rep_site_volume_A3_a, P.rep_site_volume_A3_b)
    vrat = np.maximum(P.rep_site_volume_A3_a, P.rep_site_volume_A3_b) / np.maximum(vmin, 1)
    f = {"flag_small_volume": vmin < 100, "flag_volume_ratio": vrat > 5,
         "flag_low_plddt": (np.minimum(P.plddt_mean_rep_a, P.plddt_mean_rep_b) < 80) |
                           (np.maximum(P.plddt_frac_lt70_max_a, P.plddt_frac_lt70_max_b) > 0.2),
         "flag_small_footprint": P.footprint_size_b <= 6,
         "flag_few_residues": np.minimum(P.n_res_a, P.n_res_b) < 8,
         "flag_low_burrow": np.minimum(P.max_burrow_score_a, P.max_burrow_score_b) < 0.1}
    for k, v in f.items():
        P[k] = np.asarray(v, bool)
    P["artifact_like"] = np.any(np.stack([P[k].values for k in f]), 0)
    prim_p = P[P.passed_primary]
    summ["primary_passes_artifact_like"] = int(prim_p.artifact_like.sum())
    summ["primary_passes_flags"] = {k: int(prim_p[k].sum()) for k in f}
    summ["all_listed_artifact_like_frac"] = float(P.artifact_like.mean())
    P = P.sort_values("p").drop(columns=["a", "b"])
P.to_parquet(a.results / f"convergence_candidates{a.tag}.parquet", index=False)
summ["wall_s"] = time.time() - T0
(a.results / f"directional_summary{a.tag}.json").write_text(json.dumps(summ, indent=1, default=float))
pairs.write_sidecar(a.results / f"convergence_candidates{a.tag}.parquet", a.snapshot, [a.seed, a.seed + 40],
                    {"stage": "directional", **summ})
log("directional done", json.dumps(summ, default=float))
