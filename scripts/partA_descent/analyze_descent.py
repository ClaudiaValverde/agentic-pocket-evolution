#!/usr/bin/env python
"""Part A §4.1 M1/M2 analysis: rates vs divergence time and identity, ortholog-conjecture test, M2 cosine vs noise
floor and non-homologous background, positive-control families, figures, DESCENT numbers (results/numbers.json).

Reads the pair table from build_descent_pairs.py. Everything is parameterised by snapshot / orthology / prereg.
"""
import argparse
import fcntl
import json
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pyarrow.parquet as pq  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))
from pocketevo import descent  # noqa: E402
from pocketevo.homology import lca_type, parse_orthoxml, share_any  # noqa: E402
from pocketevo.pairs import check_prereg, load_prereg, sha256_file, write_sidecar  # noqa: E402

C_ORTHO, C_PARA, C_BG, C_FLOOR, C_FLOOR2 = "#2a78d6", "#eb6834", "#8a8984", "#1baf7a", "#6250d6"
ID_EDGES = np.round(np.arange(0.0, 1.0001, 0.1), 2)
T_EDGES = np.array([-1, 0.5, 100, 300, 700, 1200, 2000, 5000])   # Mya; first bin = same species (0)
T_LABELS = ["0 (same sp.)", "0.5-100", "100-300", "300-700", "700-1200", "1200-2000", ">2000"]
FIG_SIDECAR = {}
# Rank proxy (named fallback timetree_to_rank_proxy, flagged per pair): TimeTree 5 has no study for the root of
# cellular life (API: "No TimeTree study info available for this MRCA", taxon 131567). LUCA age ~4,200 Mya
# (Moody et al. 2024, Nat Ecol Evol 8:1654, doi:10.1038/s41559-024-02461-1). Only used where TimeTree has no age.
RANK_PROXY_MYA = {"cellular organisms": 4200.0}


# ----------------------------------------------------------------------------------------------- divergence time
def parse_newick(s: str):
    """Minimal newick parser -> (parent dict, leaf labels, internal labels). Leaves are taxids."""
    s = s.strip().rstrip(";")
    parent, label, stack, cur, k, n = {}, {}, [], None, 0, 0
    tok = ""

    def new():
        nonlocal n
        n += 1
        return n - 1
    root = None
    i = 0
    while i < len(s):
        c = s[i]
        if c == "(":
            node = new()
            if stack:
                parent[node] = stack[-1]
            else:
                root = node
            stack.append(node)
            i += 1
            last = None
            continue
        if c in ",)":
            i += 1
            if c == ")":
                last = stack.pop()
            continue
        j = i
        if c == "'":
            j = s.index("'", i + 1) + 1
            name = s[i + 1:j - 1]
        else:
            while j < len(s) and s[j] not in ",()":
                j += 1
            name = s[i:j]
        if s[i - 1] == ")":           # label of the node just closed
            label[last] = name
        else:
            leaf = new()
            parent[leaf] = stack[-1]
            label[leaf] = name
        i = j
    return parent, label, root


def species_lca_table(nwk_path: str) -> dict:
    parent, label, root = parse_newick(Path(nwk_path).read_text())
    children = {}
    for c, p in parent.items():
        children.setdefault(p, []).append(c)
    leaves = {int(label[x]): x for x in label if x not in children and label[x].isdigit()}

    def anc(x):
        out = [x]
        while x in parent:
            x = parent[x]
            out.append(x)
        return out
    A = {t: anc(x) for t, x in leaves.items()}
    lca = {}
    for a in leaves:
        sa = A[a]
        for b in leaves:
            sb = set(A[b])
            node = next(x for x in sa if x in sb)
            lca[(a, b)] = (node, label.get(node, f"node{node}"))
    return lca


def divergence_table(tt_csv: str, nwk: str) -> pd.DataFrame:
    """Per ordered species pair: age (Mya) and source (timetree / timetree_node_median / rank_proxy)."""
    tt = pd.read_csv(tt_csv)
    tt["age"] = pd.to_numeric(tt.age_mya, errors="coerce")
    lca = species_lca_table(nwk)
    ages = {}
    for r in tt.itertuples():
        if r.age > 0:
            ages[(r.taxid_a, r.taxid_b)] = ages[(r.taxid_b, r.taxid_a)] = r.age
    node_ages = {}
    for (a, b), (node, _) in lca.items():
        if (a, b) in ages:
            node_ages.setdefault(node, []).append(ages[(a, b)])
    # rank proxy (flagged): fixed coarse ages for named NCBI ranks, used only when TimeTree has nothing at the node
    rows = []
    for (a, b), (node, name) in lca.items():
        if a == b:
            rows.append((a, b, 0.0, "same_species", name))
        elif (a, b) in ages:
            rows.append((a, b, ages[(a, b)], "timetree", name))
        elif node in node_ages:
            rows.append((a, b, float(np.median(node_ages[node])), "timetree_node_median", name))
        elif name in RANK_PROXY_MYA:
            rows.append((a, b, RANK_PROXY_MYA[name], "rank_proxy", name))
        else:
            rows.append((a, b, np.nan, "missing", name))
    return pd.DataFrame(rows, columns=["taxid_a", "taxid_b", "age_mya", "age_source", "lca_name"])


# ----------------------------------------------------------------------------------------------- helpers
def weights(P: pd.DataFrame, scheme: str) -> np.ndarray:
    if scheme == "cap":
        return P.in_cap.to_numpy(float)
    if scheme == "unc":
        return P.w_unc.to_numpy(float)
    raise ValueError(scheme)


def boot_ci(stat_fn, og_codes, n_og, n_rep, seed):
    vals = []
    for r in range(n_rep):
        w = descent.bootstrap_weights(og_codes, n_og, r, seed)
        vals.append(stat_fn(w))
    vals = np.array(vals, dtype=float)
    return np.nanpercentile(vals, [2.5, 97.5], axis=0), vals


def savefig(fig, path, inputs, script):
    fig.savefig(path, dpi=110, bbox_inches="tight")
    plt.close(fig)
    FIG_SIDECAR[str(path)] = inputs
    return path


def style(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#e6e5e0", lw=0.8)
    ax.set_axisbelow(True)


# ----------------------------------------------------------------------------------------------- M1
def m1_tables(P, v, scheme, codes, n_og, n_rep, seed):
    """Overall + by identity + by time, ortholog vs paralog, raw and ceiling-normalised, with OG-bootstrap CIs."""
    w0 = weights(P, scheme)
    rel = P.rel.to_numpy()
    num = P[f"n_match_{v}"].to_numpy(float)
    den = P[f"n_inf_{v}"].to_numpy(float)
    cei = P[f"ceil_sum_{v}"].to_numpy(float)
    crm = P[f"ceilrm_sum_{v}"].to_numpy(float)
    idb = np.clip(np.digitize(P.identity.fillna(-1).to_numpy(), ID_EDGES) - 1, -1, len(ID_EDGES) - 2)
    tb = np.digitize(P.age_mya.fillna(-9).to_numpy(), T_EDGES) - 1
    groups = {"overall": np.zeros(len(P), int)}
    groups["identity"] = idb
    groups["time"] = tb
    out = []
    for gname, g in groups.items():
        nb = 1 if gname == "overall" else (len(ID_EDGES) - 1 if gname == "identity" else len(T_LABELS))
        for r, rname in ((descent.REL_ORTHO, "ortholog"), (descent.REL_PARA, "paralog")):
            m = (rel == r) & (g >= 0) & (g < nb)
            gi = g[m]

            def stat(wb, m=m, gi=gi):
                ww = w0[m] * (1 if wb is None else wb[m])
                N = np.bincount(gi, ww * num[m], nb)
                D = np.bincount(gi, ww * den[m], nb)
                Cc = np.bincount(gi, ww * cei[m], nb)
                Cr = np.bincount(gi, ww * crm[m], nb)
                with np.errstate(invalid="ignore", divide="ignore"):
                    return np.r_[N / D, N / Cc, N / Cr]
            pt = stat(None)
            ci, _ = boot_ci(stat, codes, n_og, n_rep, seed)
            D = np.bincount(gi, w0[m] * den[m], nb)
            npairs = np.bincount(gi, (w0[m] > 0).astype(float), nb)
            nogs = [len(set(P.og_id.to_numpy()[m][(gi == b) & (w0[m] > 0)])) for b in range(nb)]
            for b in range(nb):
                lab = "all" if gname == "overall" else (f"{ID_EDGES[b]:.1f}-{ID_EDGES[b+1]:.1f}" if gname == "identity"
                                                        else T_LABELS[b])
                out.append({"variant": v, "weighting": scheme, "by": gname, "bin": lab, "bin_index": b, "rel": rname,
                            "m1": pt[b], "m1_lo": ci[0][b], "m1_hi": ci[1][b], "m1_norm": pt[nb + b],
                            "m1_norm_lo": ci[0][nb + b], "m1_norm_hi": ci[1][nb + b],
                            "m1_normrm": pt[2 * nb + b], "m1_normrm_lo": ci[0][2 * nb + b],
                            "m1_normrm_hi": ci[1][2 * nb + b], "n_obs": D[b],
                            "n_pairs": npairs[b], "n_ogs": nogs[b]})
    return pd.DataFrame(out)


def conjecture_cells(P, v, scheme, width=0.02, cross_species_only=False, by_time=False):
    w0 = weights(P, scheme)
    m = (w0 > 0) & P.identity.notna().to_numpy() & (P[f"n_inf_{v}"].to_numpy() > 0)
    if cross_species_only:
        m &= (P.taxid_a != P.taxid_b).to_numpy()
    q = P.loc[m, ["og_id", "rel", "identity", f"n_match_{v}", f"n_inf_{v}", f"ceil_sum_{v}", "age_mya"]].copy()
    ww = w0[m]
    q["ib"] = np.minimum((q.identity / width).astype(int), int(1 / width) - 1)
    q["tb"] = np.digitize(q.age_mya.fillna(0).to_numpy(), [0.5, 100, 300, 700, 1200, 2000]) if by_time else 0
    q["succ"] = q[f"n_match_{v}"] * ww
    q["trials"] = q[f"n_inf_{v}"] * ww
    q["ceil"] = q[f"ceil_sum_{v}"] * ww
    c = q.groupby(["og_id", "rel", "ib", "tb"], sort=False)[["succ", "trials", "ceil"]].sum().reset_index()
    c["identity"] = (c.ib + 0.5) * width
    return c


def conjecture_test(cells, codes_of, n_rep, seed, knots, with_time=False):
    og_codes = codes_of.reindex(cells.og_id).to_numpy()
    n_og = int(codes_of.max()) + 1
    edges = ID_EDGES

    def fit(w):
        if with_time:
            # identity spline + time-bin indicators + class
            x = cells.identity.to_numpy(float)
            B = descent.bspline_basis(x, knots, 0, 1)
            T = pd.get_dummies(cells.tb).to_numpy(float)[:, 1:]
            orth = (cells.rel.to_numpy() == descent.REL_ORTHO).astype(float)
            X = np.column_stack([np.ones(len(x)), B, T, orth])
            beta = descent.fit_logistic(X, cells.succ.to_numpy(float), cells.trials.to_numpy(float), w)
            X1, X0 = X.copy(), X.copy()
            X1[:, -1], X0[:, -1] = 1, 0
            ww = cells.trials.to_numpy(float) * (1 if w is None else w) * (orth == 0)
            d = np.sum(ww * (descent.expit(X1 @ beta) - descent.expit(X0 @ beta))) / ww.sum()
            return np.array([beta[-1], d, np.nan])
        e = descent.ortholog_effect(cells, knots, w)
        b = descent.binned_effect(cells, edges, w)
        return np.array([e["log_or"], e["marginal_diff"], b["std_diff"]])
    pt = fit(None)
    ci, vals = boot_ci(fit, og_codes, n_og, n_rep, seed)
    return {"log_or": pt[0], "log_or_ci": ci[:, 0].tolist(), "or": float(np.exp(pt[0])),
            "or_ci": np.exp(ci[:, 0]).tolist(), "marginal_diff": pt[1], "marginal_diff_ci": ci[:, 1].tolist(),
            "binned_std_diff": pt[2], "binned_std_diff_ci": ci[:, 2].tolist(),
            "p_boot_two_sided_marginal": float(min(1.0, 2 * min(np.mean(vals[:, 1] <= 0), np.mean(vals[:, 1] >= 0)))),
            "n_cells": int(len(cells)), "n_ogs": int(cells.og_id.nunique()),
            "n_obs_ortho": float(cells.trials[cells.rel == 1].sum()), "n_obs_para": float(cells.trials[cells.rel == 2].sum())}


# ----------------------------------------------------------------------------------------------- M2
def m2_background(M2, n, seed, universe, ctx_path, fp):
    """Random non-homologous size-matched tier-1 site pairs: one site drawn from the matched sites, partner a
    random tier-1 site with the same n_res bin, different root group and no shared protein-level label."""
    rng = np.random.default_rng(seed)
    z = np.load(ctx_path, allow_pickle=True)
    site_prot, prot_group, prot_lab = z["site_prot"], z["prot_group"], z["prot_lab"]
    nres = universe.n_res.to_numpy()
    tier1 = np.flatnonzero(universe.tier.to_numpy() == 1)
    bins = np.r_[np.arange(3, 31), np.arange(35, 101, 5), 1000, 100000]
    b_all = np.digitize(nres, bins)
    by_bin = {b: tier1[b_all[tier1] == b] for b in np.unique(b_all[tier1])}
    src = rng.choice(np.r_[M2.si_a.to_numpy(), M2.si_b.to_numpy()], n * 3)
    out_i, out_j = [], []
    for k in range(0, len(src), 20000):
        i = src[k:k + 20000]
        bi = b_all[i]
        j = np.array([by_bin[b][rng.integers(len(by_bin[b]))] for b in bi])
        pi, pj = site_prot[i], site_prot[j]
        ok = (prot_group[pi] != prot_group[pj]) & ~share_any(prot_lab[pi], prot_lab[pj])
        out_i.append(i[ok])
        out_j.append(j[ok])
        if sum(len(x) for x in out_i) >= n:
            break
    i = np.concatenate(out_i)[:n]
    j = np.concatenate(out_j)[:n]
    fa = fp[i].astype(np.float32)
    fb = fp[j].astype(np.float32)
    cos = (fa * fb).sum(1) / (np.linalg.norm(fa, axis=1) * np.linalg.norm(fb, axis=1) + 1e-12)
    return pd.DataFrame({"si_a": i, "si_b": j, "cosine": cos, "n_res_a": nres[i], "n_res_b": nres[j]})


def noise_floor(redet_path, tau):
    r = pd.read_parquet(redet_path)
    r = r[r.best_jaccard >= tau - 1e-9]
    out = {}
    for comp, g in r.groupby("comparison"):
        c = g.cos_best.dropna()
        out[comp] = {"n": int(len(c)), "median": float(c.median()), "p05": float(c.quantile(0.05)),
                     "p25": float(c.quantile(0.25)), "p75": float(c.quantile(0.75))}
    return out


TASK = {}


def _run_task(t):
    kind, v, scheme, extra = t
    Pm, codes_m, codes_of, n_og, n_boot, seed, knots = (TASK[k] for k in
                                                        ("Pm", "codes_m", "codes_of", "n_og", "n_boot", "seed", "knots"))
    if kind == "m1":
        return t, m1_tables(Pm, v, scheme, codes_m, n_og, n_boot, seed)
    cs = extra in ("cross_species", "cross_species|+time")
    cells = conjecture_cells(Pm, v, scheme, cross_species_only=cs, by_time=extra == "cross_species|+time")
    return t, conjecture_test(cells, codes_of, n_boot, seed, knots, with_time=extra == "cross_species|+time")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--snapshot", default="20261030T085136Z")
    ap.add_argument("--orthology-version", default="O1")
    ap.add_argument("--bulk", default="/tmp/pe_big/partA_descent")
    ap.add_argument("--timetree", default=None)
    ap.add_argument("--species-tree", default=None)
    ap.add_argument("--universe", default="/tmp/pe_big/partA/20261030T085136Z/universe.parquet")
    ap.add_argument("--context", default="/tmp/pe_big/partA/20261030T085136Z/context.npz")
    ap.add_argument("--fingerprints", default="/tmp/pe_big/partA/20261030T085136Z/fp.f16.npy")
    ap.add_argument("--redetection", default=None)
    ap.add_argument("--slice", default=None)
    ap.add_argument("--n-boot", type=int, default=200)
    ap.add_argument("--workers", type=int, default=24)
    ap.add_argument("--n-background", type=int, default=100000)
    ap.add_argument("--results", default="results/partA_descent")
    ap.add_argument("--numbers", default="results/numbers.json")
    ap.add_argument("--lock", default="/tmp/pe_big/numbers.lock")
    a = ap.parse_args()
    t0 = time.time()
    check_prereg()
    pre = load_prereg()
    seed = int(pre["meta"]["seeds"][0])
    snap, ov = a.snapshot, a.orthology_version
    a.timetree = a.timetree or f"data/timetree/{snap}_{ov}_pairwise.csv"
    a.species_tree = a.species_tree or f"results/sogs/{snap}/{ov}/web/species_tree.nwk"
    a.redetection = a.redetection or f"results/observation/{snap}/redetection_sites.parquet"
    a.slice = a.slice or f"results/sogs/{snap}/{ov}/slice_G2_sites.csv"
    bulk = Path(a.bulk) / snap / ov
    outdir = REPO / a.results / snap / ov
    figdir = outdir / "figures"
    figdir.mkdir(parents=True, exist_ok=True)
    side = json.loads((bulk / "pairs.parquet.json").read_text())
    tau = side["tau"]

    P = pd.read_parquet(bulk / "pairs.parquet")
    M2 = pd.read_parquet(bulk / "m2_matches.parquet")
    div = divergence_table(a.timetree, a.species_tree)
    P = P.merge(div[["taxid_a", "taxid_b", "age_mya", "age_source", "lca_name"]], on=["taxid_a", "taxid_b"], how="left")
    P = P[P.rel.isin([1, 2])].reset_index(drop=True)
    codes_of = pd.Series(np.arange(P.og_id.nunique()), index=sorted(P.og_id.unique()))
    codes = codes_of[P.og_id].to_numpy()
    n_og = len(codes_of)
    nonctrl = ~P.control.to_numpy() | P.in_cap.to_numpy() | P.in_unc.to_numpy()
    Pm = P[nonctrl].reset_index(drop=True)
    codes_m = codes_of[Pm.og_id].to_numpy()

    numbers = {"snapshot_id": snap, "orthology": ov, "tau": tau, "prereg_sha256": side["prereg_sha256"],
               "ceiling_source": side["ceiling_source"], "ceiling": side["ceiling"], "cap": side["cap"],
               "unc_max": side["unc_max"], "n_boot": a.n_boot, "seed": seed}
    fb = div[div.taxid_a < div.taxid_b].age_source.value_counts().to_dict()
    numbers["divergence_time_sources_species_pairs"] = fb
    numbers["named_fallback_counters"] = {"timetree_to_rank_proxy": int(fb.get("rank_proxy", 0)),
                                          "timetree_missing": int(fb.get("missing", 0)),
                                          "timetree_to_node_median": int(fb.get("timetree_node_median", 0))}
    numbers["pairs"] = {
        "n_rows": int(len(P)), "n_ogs": int(n_og),
        "cap_pairs": {"ortholog": int((Pm.in_cap & (Pm.rel == 1)).sum()), "paralog": int((Pm.in_cap & (Pm.rel == 2)).sum())},
        "unc_pairs_sampled": int(Pm.in_unc.sum()),
        "candidate_pairs_total": int(Pm.groupby("og_id").n_cand.first().sum())}

    # ---------------- M1 tables + ortholog conjecture (parallel over variant x weighting; fork shares the table)
    rmc = descent.rule_matched_ceiling(pd.read_parquet(a.redetection), side["ceiling_source"], tau)
    numbers["ceiling_rule_matched"] = {f"{k[0]}|{k[1]}": v for k, v in rmc.table.items()}
    for v in descent.VARIANTS:
        for D_ in (Pm, P):
            D_[f"ceilrm_sum_{v}"] = sum(rmc.value(b, s_) * D_[f"n_inf_{v}_s{k}"].to_numpy(float)
                                        for k, (b, s_) in enumerate(descent.STRATA))
    knots = np.array([0.2, 0.3, 0.4, 0.55, 0.75])
    TASK.update({"Pm": Pm, "codes_m": codes_m, "codes_of": codes_of, "n_og": n_og, "n_boot": a.n_boot, "seed": seed,
                 "knots": knots})
    tasks = [("m1", v, sc, "") for v in descent.VARIANTS for sc in ("cap", "unc")]
    tasks += [("conj", v, sc, "") for v in descent.VARIANTS for sc in ("cap", "unc")]
    tasks += [("conj", "qm90", "cap", "cross_species"), ("conj", "qm90", "cap", "cross_species|+time")]
    import multiprocessing as mp
    with mp.get_context("fork").Pool(min(a.workers, len(tasks))) as pool:
        res = dict(pool.map(_run_task, tasks, chunksize=1))
    T = pd.concat([res[t] for t in tasks if t[0] == "m1"], ignore_index=True)
    conj = {f"{t[1]}|{t[2]}" + (f"|{t[3]}" if t[3] else ""): res[t] for t in tasks if t[0] == "conj"}
    T.to_csv(outdir / "m1_rates.csv", index=False)
    write_sidecar(outdir / "m1_rates.csv", snap, seed, {"inputs": {"pairs": str(bulk / "pairs.parquet")}})
    ov_ = T[T.by == "overall"]
    numbers["M1_overall"] = {f"{r.variant}|{r.weighting}|{r.rel}": {
        "m1": r.m1, "ci": [r.m1_lo, r.m1_hi], "m1_norm": r.m1_norm, "norm_ci": [r.m1_norm_lo, r.m1_norm_hi],
        "m1_norm_rule_matched": r.m1_normrm, "norm_rule_matched_ci": [r.m1_normrm_lo, r.m1_normrm_hi],
        "n_obs": r.n_obs, "n_pairs": r.n_pairs, "n_ogs": r.n_ogs} for r in ov_.itertuples()}
    for by in ("time", "identity"):
        tt = T[(T.by == by) & (T.variant == "qm90") & (T.weighting == "cap")]
        numbers[f"M1_qm90_cap_by_{by}"] = {f"{r.rel}|{r.bin}": {
            "m1": r.m1, "ci": [r.m1_lo, r.m1_hi], "m1_norm": r.m1_norm, "m1_norm_rule_matched": r.m1_normrm,
            "norm_rule_matched_ci": [r.m1_normrm_lo, r.m1_normrm_hi], "n_pairs": r.n_pairs, "n_ogs": r.n_ogs}
            for r in tt.itertuples()}
    numbers["ortholog_conjecture"] = conj

    print("stage", ' M2', f"{time.time()-t0:.0f}s", flush=True)
    # ---------------- M2
    M2 = M2.merge(P[["pair_id", "rel", "identity", "age_mya", "in_cap", "in_all", "control"]], on="pair_id", how="inner")
    M2c = M2[M2.in_cap].reset_index(drop=True)
    universe = pq.read_table(a.universe, columns=["si", "site_id", "n_res", "tier"]).to_pandas()
    fp = np.load(a.fingerprints, mmap_mode="r")
    bg = m2_background(M2c, a.n_background, seed, universe, a.context, fp)
    bg.to_parquet(bulk / "m2_background.parquet", index=False)
    floor = noise_floor(a.redetection, tau)
    bg_med, bg_p95, bg_p99 = float(bg.cosine.median()), float(bg.cosine.quantile(0.95)), float(bg.cosine.quantile(0.99))
    fl = floor.get("foldcomp_v4_kept", {})
    numbers["M2"] = {"noise_floor": floor, "background": {"n": int(len(bg)), "median": bg_med, "p95": bg_p95,
                                                          "p99": bg_p99}}
    rows = []
    idb = np.clip(np.digitize(M2c.identity.fillna(-1).to_numpy(), ID_EDGES) - 1, -1, 9)
    tb = np.digitize(M2c.age_mya.fillna(-9).to_numpy(), T_EDGES) - 1
    codes2 = codes_of[M2c.og_id].to_numpy()
    for by, g, labs in (("overall", np.zeros(len(M2c), int), ["all"]), ("identity", idb,
                        [f"{ID_EDGES[b]:.1f}-{ID_EDGES[b+1]:.1f}" for b in range(10)]), ("time", tb, T_LABELS)):
        for r, rn in ((1, "ortholog"), (2, "paralog")):
            for b, lab in enumerate(labs):
                m = (M2c.rel.to_numpy() == r) & (g == b)
                if m.sum() == 0:
                    continue
                c = M2c.cosine.to_numpy()[m]

                def st(w, m=m):
                    ww = np.ones(len(M2c)) if w is None else w
                    x = M2c.cosine.to_numpy()[m]
                    wm = ww[m]
                    if wm.sum() == 0:
                        return np.array([np.nan, np.nan])
                    o = np.argsort(x)
                    cw = np.cumsum(wm[o]) / wm.sum()
                    med = x[o][np.searchsorted(cw, 0.5)]
                    return np.array([med, np.sum(wm * (x > bg_p95)) / wm.sum()])
                ci, _ = boot_ci(st, codes2, n_og, a.n_boot, seed)
                rows.append({"by": by, "bin": lab, "rel": rn, "n": int(m.sum()), "median": float(np.median(c)),
                             "median_lo": ci[0][0], "median_hi": ci[1][0], "q25": float(np.quantile(c, 0.25)),
                             "q75": float(np.quantile(c, 0.75)), "frac_above_bg_p95": float(np.mean(c > bg_p95)),
                             "frac_above_bg_p95_lo": ci[0][1], "frac_above_bg_p95_hi": ci[1][1],
                             "frac_above_floor_p05": float(np.mean(c > fl.get("p05", np.nan))),
                             "retained": float((np.median(c) - bg_med) / (fl.get("median", np.nan) - bg_med))})
    M2T = pd.DataFrame(rows)
    M2T.to_csv(outdir / "m2_cosine.csv", index=False)
    write_sidecar(outdir / "m2_cosine.csv", snap, seed, {"inputs": {"m2": str(bulk / "m2_matches.parquet")}})
    numbers["M2"]["by"] = {f"{r.by}|{r.bin}|{r.rel}": {"n": r.n, "median": r.median, "ci": [r.median_lo, r.median_hi],
                                                       "frac_above_bg_p95": r.frac_above_bg_p95,
                                                       "retained_fraction": r.retained} for r in M2T.itertuples()}
    # matched-identity M2 ortholog-paralog difference (binned standardisation, paralog weights)
    def m2diff(w):
        ww = np.ones(len(M2c)) if w is None else w
        d, tot = 0.0, 0.0
        for b in range(10):
            mo = (M2c.rel.to_numpy() == 1) & (idb == b)
            mp = (M2c.rel.to_numpy() == 2) & (idb == b)
            if ww[mo].sum() > 0 and ww[mp].sum() > 0:
                mo_ = np.sum(ww[mo] * M2c.cosine.to_numpy()[mo]) / ww[mo].sum()
                mp_ = np.sum(ww[mp] * M2c.cosine.to_numpy()[mp]) / ww[mp].sum()
                d += ww[mp].sum() * (mo_ - mp_)
                tot += ww[mp].sum()
        return np.array([d / tot if tot else np.nan])
    ci, _ = boot_ci(m2diff, codes2, n_og, a.n_boot, seed)
    numbers["M2"]["ortholog_minus_paralog_mean_cosine_at_matched_identity"] = {"diff": float(m2diff(None)[0]),
                                                                               "ci": [ci[0][0], ci[1][0]]}

    print("stage", ' positive controls', f"{time.time()-t0:.0f}s", flush=True)
    # ---------------- positive controls (all pairs of the control OGs)
    sl = pd.read_csv(a.slice)
    ctrl_rows = []
    for og_id, g in P[P.in_all].groupby("og_id"):
        fams = "; ".join(sorted(set(sl[sl.og_id == og_id].site)))
        for r, rn in ((1, "ortholog"), (2, "paralog")):
            q = g[g.rel == r]
            if len(q) == 0:
                continue
            m2q = M2[(M2.og_id == og_id) & (M2.rel == r) & M2.in_all]
            for v in ("lit", "qm90", "s70"):
                ctrl_rows.append({"og_id": og_id, "family": fams, "rel": rn, "variant": v, "n_pairs": len(q),
                                  "m1": q[f"n_match_{v}"].sum() / max(q[f"n_inf_{v}"].sum(), 1),
                                  "m1_norm": q[f"n_match_{v}"].sum() / max(q[f"ceil_sum_{v}"].sum(), 1e-9),
                                  "m1_norm_rule_matched": q[f"n_match_{v}"].sum() / max(q[f"ceilrm_sum_{v}"].sum(), 1e-9),
                                  "n_obs": float(q[f"n_inf_{v}"].sum()),
                                  "m2_median_cos": float(m2q.cosine.median()) if len(m2q) else np.nan,
                                  "m2_n": int(len(m2q)),
                                  "median_identity": float(q.identity.median())})
        # known-pocket SOG only (M2: matched sites in that SOG)
        for srow in sl[sl.og_id == og_id].itertuples():
            sidx = int(str(srow.sog_id).split(":g")[1])
            m2q = M2[(M2.og_id == og_id) & (M2.sog_index == sidx) & M2.in_all]
            for r, rn in ((1, "ortholog"), (2, "paralog")):
                mq = m2q[m2q.rel == r]
                ctrl_rows.append({"og_id": og_id, "family": srow.site, "rel": rn, "variant": f"SOG {srow.sog_id}",
                                  "n_pairs": np.nan, "m1": np.nan, "m1_norm": np.nan, "n_obs": np.nan,
                                  "m2_median_cos": float(mq.cosine.median()) if len(mq) else np.nan,
                                  "m2_n": int(len(mq)), "median_identity": float(mq.identity.median()) if len(mq) else np.nan,
                                  "sog_n_P": int(srow.n_P), "sog_n_A": int(srow.n_A)})
    CT = pd.DataFrame(ctrl_rows)
    CT.to_csv(outdir / "controls.csv", index=False)
    write_sidecar(outdir / "controls.csv", snap, seed, {})
    numbers["controls"] = CT.replace({np.nan: None}).to_dict(orient="records")

    print("stage", ' figures', f"{time.time()-t0:.0f}s", flush=True)
    # ---------------- figures
    figs(T, M2T, M2c, bg, floor, conj, CT, figdir, bulk, Pm, knots)
    for pth, inp in FIG_SIDECAR.items():
        write_sidecar(Path(pth), snap, seed, {"script": "scripts/partA_descent/analyze_descent.py", "inputs": inp})
    numbers["runtime_s_analysis"] = time.time() - t0
    numbers["runtime_s_build"] = side["runtime_s"]
    (outdir / "descent_numbers.json").write_text(json.dumps(numbers, indent=1, default=_jd))
    write_sidecar(outdir / "descent_numbers.json.sidecar", snap, seed, {})
    merge_numbers(REPO / a.numbers, a.lock, "partA_descent", json.loads(json.dumps(numbers, default=_jd)))
    print(json.dumps({k: numbers[k] for k in ("pairs", "M1_overall")}, indent=1, default=_jd)[:4000])
    print("conjecture", json.dumps(conj, indent=1, default=_jd)[:3000])
    print(f"done {time.time()-t0:.0f}s")


def _jd(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return None if np.isnan(o) else float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, float) and np.isnan(o):
        return None
    return str(o)


def merge_numbers(path, lock, key, value):
    """Read-modify-write under an exclusive lock (file may be written concurrently by other agents)."""
    Path(lock).parent.mkdir(parents=True, exist_ok=True)
    with open(lock, "w") as lk:
        fcntl.flock(lk, fcntl.LOCK_EX)
        cur = json.loads(Path(path).read_text()) if Path(path).exists() else {}
        cur[key] = value
        tmp = Path(str(path) + ".tmp_descent")
        tmp.write_text(json.dumps(cur, indent=1, default=_jd, allow_nan=False) if False else
                       json.dumps(cur, indent=1, default=_jd))
        tmp.replace(path)
        fcntl.flock(lk, fcntl.LOCK_UN)


def figs(T, M2T, M2c, bg, floor, conj, CT, figdir, bulk, Pm, knots):
    inp = {"pairs": str(bulk / "pairs.parquet"), "m2": str(bulk / "m2_matches.parquet")}
    # F1 M1 vs time (orthologs + paralogs), qm90 and lit, cap
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.8), sharey=True)
    for ax, col, ttl in ((axes[0], "m1", "M1 (raw)"), (axes[1], "m1_norm", "M1 / prereg ceiling (error_rates)"),
                         (axes[2], "m1_normrm", "M1 / rule-matched ceiling")):
        for v, ls in (("qm90", "-"), ("lit", ":")):
            for rel, c in (("ortholog", C_ORTHO), ("paralog", C_PARA)):
                q = T[(T.by == "time") & (T.variant == v) & (T.weighting == "cap") & (T.rel == rel) & (T.n_pairs >= 30)]
                x = q.bin_index.to_numpy()
                ax.plot(x, q[col], ls, color=c, lw=2, marker="o", ms=5, label=f"{rel}, {v}")
                if v == "qm90":
                    ax.fill_between(x, q[f"{col}_lo"], q[f"{col}_hi"], color=c, alpha=0.18, lw=0)
        if col != "m1":
            ax.axhline(1.0, color=C_BG, lw=1, ls="--")
        ax.set_xticks(range(len(T_LABELS)), T_LABELS, rotation=35, ha="right", fontsize=8)
        ax.set_xlabel("species divergence (Mya, TimeTree)")
        ax.set_title(ttl, fontsize=10)
        style(ax)
    axes[0].set_ylabel("P(match | informative)")
    axes[0].legend(fontsize=7, frameon=False)
    fig.suptitle("M1 vs divergence time (cap ≤200 pairs/OG; band = 95% OG bootstrap, qm90)", fontsize=10)
    savefig(fig, figdir / "m1_vs_time.png", inp, None)
    # F2 M1 vs identity + logistic fit
    fig, ax = plt.subplots(figsize=(5.6, 3.8))
    for rel, c in (("ortholog", C_ORTHO), ("paralog", C_PARA)):
        q = T[(T.by == "identity") & (T.variant == "qm90") & (T.weighting == "cap") & (T.rel == rel) & (T.n_pairs >= 30)]
        x = np.array([float(b.split("-")[0]) + 0.05 for b in q.bin])
        ax.errorbar(x, q.m1, yerr=[np.clip(q.m1 - q.m1_lo, 0, None), np.clip(q.m1_hi - q.m1, 0, None)], color=c, marker="o", ms=5, lw=2, capsize=2,
                    label=rel)
    ax.set_xlabel("global sequence identity (OG MSA)")
    ax.set_ylabel("M1 = P(match | informative)")
    k = conj["qm90|cap"]
    ax.set_title(f"Pocket ortholog conjecture (qm90, cap): Δ = {k['marginal_diff']:+.3f} "
                 f"[{k['marginal_diff_ci'][0]:+.3f}, {k['marginal_diff_ci'][1]:+.3f}]", fontsize=9)
    ax.legend(frameon=False, fontsize=8)
    style(ax)
    savefig(fig, figdir / "m1_vs_identity.png", inp, None)
    # F3 M2 cosine vs time / identity with floor and background
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.6), sharey=True)
    fl = floor.get("foldcomp_v4_kept", {})
    fs = floor.get("seed_s1", {})
    for ax, by in ((axes[0], "time"), (axes[1], "identity")):
        for rel, c in (("ortholog", C_ORTHO), ("paralog", C_PARA)):
            q = M2T[(M2T.by == by) & (M2T.rel == rel) & (M2T.n >= 30)]
            if by == "time":
                x = np.array([T_LABELS.index(b) for b in q.bin])
            else:
                x = np.array([float(b.split("-")[0]) + 0.05 for b in q.bin])
            ax.plot(x, q["median"], color=c, lw=2, marker="o", ms=5, label=f"{rel} median")
            ax.fill_between(x, q.q25, q.q75, color=c, alpha=0.15, lw=0)
        ax.axhline(fl.get("median", np.nan), color=C_FLOOR, lw=1.5, ls="--", label="same protein, v4 vs v6 (median)")
        ax.axhline(fs.get("median", np.nan), color=C_FLOOR2, lw=1.5, ls="--", label="same model, BURROW seed (median)")
        ax.axhline(float(bg.cosine.median()), color=C_BG, lw=1.5, ls="-.", label="non-homologous, size-matched (median)")
        ax.axhline(float(bg.cosine.quantile(0.95)), color=C_BG, lw=1, ls=":", label="background p95")
        if by == "time":
            ax.set_xticks(range(len(T_LABELS)), T_LABELS, rotation=35, ha="right", fontsize=8)
            ax.set_xlabel("species divergence (Mya)")
        else:
            ax.set_xlabel("global sequence identity")
        style(ax)
    axes[0].set_ylabel("fingerprint cosine, matched sites")
    axes[1].legend(fontsize=6.5, frameon=False, loc="lower right")
    fig.suptitle("M2: fingerprint similarity of positionally matched sites (band = IQR)", fontsize=10)
    savefig(fig, figdir / "m2_cosine.png", inp, None)
    # F4 controls
    q = CT[CT.variant == "qm90"]
    if len(q):
        fig, ax = plt.subplots(figsize=(6, 3.2))
        labs = []
        for k_, (og, g) in enumerate(q.groupby("og_id")):
            for r, c, dx in (("ortholog", C_ORTHO, -0.18), ("paralog", C_PARA, 0.18)):
                gg = g[g.rel == r]
                if len(gg):
                    ax.bar(k_ + dx, gg.m1.iloc[0], 0.34, color=c, label=r if k_ == 0 else None)
            labs.append(f"{og}\n{g.family.iloc[0][:28]}")
        ax.set_xticks(range(len(labs)), labs, fontsize=7)
        ax.set_ylabel("M1 (qm90, all pairs)")
        ax.set_ylim(0, 1.05)
        ax.legend(frameon=False, fontsize=8)
        ax.set_title("Positive-control families", fontsize=10)
        style(ax)
        savefig(fig, figdir / "controls_m1.png", inp, None)


if __name__ == "__main__":
    main()
