"""Part A §4.1 M1/M2: do pockets descend by orthology? (DESIGN v2 §4.1; prereg part_a.M1 / part_a.M2).

M1  For a protein pair (X, Y) of the same root OG and an SOG s in which X is P (has a site at the footprint), the
    observation is *informative* if Y's footprint is informative (Y in {P, A}); it is a *match* if Y is P.
    M1 = sum(match) / sum(informative), by relation (ortholog = LCA is a speciation node of the HOG, paralog = a
    duplication node), divergence time and global sequence identity. Every observation is counted in both directions
    (X->Y and Y->X), so a P/P observation counts as 2 matches out of 2.
    Normalisation: each informative observation carries the same-protein re-detection ceiling c(band_Y, size_X)
    (prereg ceiling_stratified_by [plddt_band, model_source]); M1_norm = sum(match) / sum(c).

    Variants of "informative" (see DESCENT.md):
      lit    prereg-literal: Y state in {P, A} (A primary, footprint pLDDT >= 90). P has no quality requirement.
      qm90   quality-matched: Y = A, or Y = P whose footprint passes the same A_primary quality criteria that are
             available for P rows (mean pLDDT >= 90, frac < 70 <= 0.2, aligned >= 0.8, length ratio in range).
             Removes the selection bias that a low-confidence Y counts if P but never if not P.
      s70    sensitivity stratum: as qm90 with state_s70 / pLDDT >= 70.
      qm90nf qm90 restricted to Y models without atlas-filtered pockets (n_filtered_pockets == 0), symmetric for
             P and A (A calls on such models are 'rescue:pending', SOGS.md).
      qm90fp qm90 restricted to SOGs with footprint > 3 columns (decision D8 small_footprint flag).
      qm90d11 (added 2026-10-30 by the Part A lead, PI request after ledger D29; the descent agent had finished)
             qm90 with ledger D11 applied SYMMETRICALLY: the Y member's footprint must span > 10 residues for P as well
             as for A. D11 (O1_r2 / O2_r1 states) already makes A with <= 10 footprint residues U, while P had no such
             requirement, which removes informative non-matches on small footprints and inflates M1. Footprint
             residues of a state row = footprint_size (SOG columns) x aligned_frac (non-gap fraction of the member at
             the footprint), rounded; rows without that information keep qm90 behaviour (fp_res = inf).

M2  For each P/P observation the cosine of the two sites' representative fingerprints (computed by the caller).

Everything here is pure numpy/scipy so that tests/test_descent.py can check it on synthetic data with known answers.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.interpolate import BSpline

VARIANTS = ("lit", "qm90", "s70", "qm90nf", "qm90fp", "qm90d11")
D11_MIN_FOOTPRINT_RESIDUES = 10   # ledger D11: A needs footprint residues > 10; qm90d11 applies it to P too
GAP = ord("-")
REL_ORTHO, REL_PARA = 1, 2


# ----------------------------------------------------------------------------------------------- ceiling
@dataclass
class Ceiling:
    """c[band][size] with band in {'>=90', '70-90'} and size in {'<=10', '>10'}."""
    table: dict
    source: str
    tau: float

    def value(self, band: str, size: str) -> float:
        return float(self.table[(band, size)])


def load_ceiling(error_rates: pd.DataFrame, source: str = "foldcomp_v4_kept", tau: float = 0.2) -> Ceiling:
    e = error_rates[(error_rates.metric == "redetect") & (error_rates.model_source == source)
                    & (np.isclose(error_rates.tau, tau)) & (error_rates.kingdom == "all")]
    tab = {}
    for band in (">=90", "70-90"):
        for size in ("<=10", ">10"):
            r = e[(e.plddt_band == band) & (e.site_size_class == size)]
            if len(r) != 1:
                raise ValueError(f"ceiling stratum missing/ambiguous: {source} {tau} {band} {size} ({len(r)})")
            tab[(band, size)] = float(r.rate.iloc[0])
    return Ceiling(tab, source, tau)


def rule_matched_ceiling(redet_sites: pd.DataFrame, source: str = "foldcomp_v4_kept", tau: float = 0.2) -> Ceiling:
    """Same-protein re-detection ceiling computed with the M1 observation rules (not a prereg number; reported next
    to the literal ceiling): a re-described site with best Jaccard >= tau is a match (P), one with no overlapping
    site at all (Jaccard 0) a non-match (A-like); partial overlaps 0 < J < tau are U and leave the denominator, exactly
    as `partial_site_overlap` members are U in the SOG states."""
    r = redet_sites[redet_sites.comparison == source]
    size = np.where(r.n_res.to_numpy() <= 10, "<=10", ">10")
    tab = {}
    for band in (">=90", "70-90"):
        for sz in ("<=10", ">10"):
            g = r[(r.band.to_numpy() == band) & (size == sz)]
            P = int((g.best_jaccard >= tau - 1e-9).sum())
            A = int((g.best_jaccard <= 0).sum())
            if P + A == 0:
                raise ValueError(f"empty rule-matched ceiling stratum {band} {sz}")
            tab[(band, sz)] = P / (P + A)
    return Ceiling(tab, source + ":rule_matched", tau)


STRATA = ((">=90", "<=10"), (">=90", ">10"), ("70-90", "<=10"), ("70-90", ">10"))


def ceil_from_strata(strata_counts: dict, ceiling: Ceiling):
    """sum_k c_k * n_k over the four (band_Y, size_X) strata."""
    return sum(ceiling.value(b, s) * strata_counts[(b, s)] for b, s in STRATA)


# ----------------------------------------------------------------------------------------------- informative masks
def informative_masks(rows: pd.DataFrame, a_crit: dict) -> dict:
    """Per state row: boolean 'informative' and 'P' per variant, plus Y band (>=90 or 70-90).

    rows needs: state, state_s70, footprint_plddt, frac_lt70, aligned_frac, len_ratio, n_filtered_pockets.
    a_crit: {'plddt': 90, 'plddt_s70': 70, 'frac_lt70': 0.2, 'aligned': 0.8, 'len_lo': 0.7, 'len_hi': 1.3}
    """
    st = rows["state"].to_numpy()
    st70 = rows["state_s70"].to_numpy()
    isP = st == "P"
    pl = rows["footprint_plddt"].to_numpy(float)
    qual_base = ((rows["frac_lt70"].to_numpy(float) <= a_crit["frac_lt70"] + 1e-9)
                 & (rows["aligned_frac"].to_numpy(float) >= a_crit["aligned"] - 1e-9)
                 & (rows["len_ratio"].to_numpy(float) >= a_crit["len_lo"] - 1e-9)
                 & (rows["len_ratio"].to_numpy(float) <= a_crit["len_hi"] + 1e-9))
    q90 = qual_base & (pl >= a_crit["plddt"] - 1e-9)
    q70 = qual_base & (pl >= a_crit["plddt_s70"] - 1e-9)
    nf = np.nan_to_num(rows["n_filtered_pockets"].to_numpy(float), nan=0.0) == 0
    inf = {
        "lit": isP | (st == "A"),
        "qm90": (isP & q90) | (st == "A"),
        "s70": (isP & q70) | (st70 == "A"),
        "qm90nf": ((isP & q90) | (st == "A")) & nf,
    }
    if "big_fp" in rows:   # D8 sensitivity: drop SOGs with footprint <= 3 columns
        inf["qm90fp"] = inf["qm90"] & rows["big_fp"].to_numpy(bool)
    fp_res = rows["fp_res"].to_numpy(float) if "fp_res" in rows else np.full(len(rows), np.inf)
    inf["qm90d11"] = inf["qm90"] & (np.nan_to_num(fp_res, nan=np.inf) > D11_MIN_FOOTPRINT_RESIDUES)
    band90 = pl >= a_crit["plddt"] - 1e-9
    return {"inf": inf, "isP": isP, "band90": band90}


# ----------------------------------------------------------------------------------------------- per-OG pair counts
def og_matrices(member_idx: np.ndarray, sog_idx: np.ndarray, n_mem: int, n_sog: int, masks: dict,
                site_small: np.ndarray, site_si: np.ndarray) -> dict:
    """Sparse member x SOG indicator matrices for one OG."""
    def m(sel, data=None):
        d = np.ones(sel.sum(), np.float64) if data is None else data[sel]
        return sp.csr_matrix((d, (member_idx[sel], sog_idx[sel])), shape=(n_mem, n_sog))
    isP = masks["isP"]
    out = {"P": m(isP), "P_small": m(isP & site_small), "P_large": m(isP & ~site_small),
           "SI": m(isP, site_si.astype(np.float64) + 1.0)}
    for v, inf in masks["inf"].items():
        out[f"I_{v}"] = m(inf)
        out[f"IP_{v}"] = m(inf & isP)
        out[f"I90_{v}"] = m(inf & masks["band90"])
        out[f"I70_{v}"] = m(inf & ~masks["band90"])
    return out


def pair_counts(M: dict, ceiling: Ceiling, variants=VARIANTS) -> dict:
    """Dense (n x n) ordered counts. inf[x, y] = #SOGs with X P and Y informative (x != y); match likewise with Y P.

    Returns per variant symmetric totals (both directions): n_inf, n_match, ceil_sum.
    """
    P = M["P"]
    out = {}
    for v in variants:
        Ci = (P @ M[f"I_{v}"].T).toarray().astype(np.float32)
        Cm = (P @ M[f"IP_{v}"].T).toarray().astype(np.float32)
        Cc = np.zeros_like(Ci)
        strata = {}
        for xs, key_s in (("P_small", "<=10"), ("P_large", ">10")):
            for yb, key_b in (("I90", ">=90"), ("I70", "70-90")):
                c = ceiling.value(key_b, key_s)
                n = (M[xs] @ M[f"{yb}_{v}"].T).toarray().astype(np.float32)
                np.fill_diagonal(n, 0)
                strata[(key_b, key_s)] = n + n.T
                Cc += np.float32(c) * n
        np.fill_diagonal(Ci, 0)
        np.fill_diagonal(Cm, 0)
        np.fill_diagonal(Cc, 0)
        out[v] = {"n_inf": Ci + Ci.T, "n_match": Cm + Cm.T, "ceil_sum": Cc + Cc.T, "strata": strata}
    return out


def sample_pairs(cand_i: np.ndarray, cand_j: np.ndarray, cap: int, rng: np.random.Generator) -> np.ndarray:
    """Indices into the candidate list: all if <= cap, else a seeded uniform sample of size cap (sorted)."""
    n = len(cand_i)
    if n <= cap:
        return np.arange(n)
    return np.sort(rng.choice(n, cap, replace=False))


def pair_identity(mat: np.ndarray, i: np.ndarray, j: np.ndarray, chunk: int = 2048) -> tuple[np.ndarray, np.ndarray]:
    """Global identity = identical / columns where both non-gap; also the number of such columns."""
    ident = np.empty(len(i))
    nb = np.empty(len(i), np.int64)
    ng = mat != GAP
    for s in range(0, len(i), chunk):
        a, b = i[s:s + chunk], j[s:s + chunk]
        both = ng[a] & ng[b]
        same = (mat[a] == mat[b]) & both
        n = both.sum(1)
        nb[s:s + chunk] = n
        with np.errstate(invalid="ignore", divide="ignore"):
            ident[s:s + chunk] = np.where(n > 0, same.sum(1) / np.maximum(n, 1), np.nan)
    return ident, nb


def matched_site_pairs(M: dict, i: np.ndarray, j: np.ndarray) -> pd.DataFrame:
    """For pairs (i, j): every SOG where both are P -> (pair_k, sog, si_i, si_j)."""
    if len(i) == 0:
        return pd.DataFrame({"k": [], "sog": [], "si_a": [], "si_b": []})
    P, SI = M["P"], M["SI"]

    def coo_sorted(X):
        X = X.tocsr()
        X.eliminate_zeros()
        X.sort_indices()
        c = X.tocoo()
        return c.row, c.col, c.data

    k, s, _ = coo_sorted(P[i].multiply(P[j]))
    ka, sa_, va = coo_sorted(SI[i].multiply(P[j]))   # SI is non-zero exactly where P is
    kb, sb_, vb = coo_sorted(SI[j].multiply(P[i]))
    if not (np.array_equal(k, ka) and np.array_equal(s, sa_) and np.array_equal(k, kb) and np.array_equal(s, sb_)):
        raise AssertionError("SI / P sparsity mismatch")
    va = va - 1
    vb = vb - 1
    return pd.DataFrame({"k": k.astype(np.int64), "sog": s.astype(np.int64), "si_a": va.astype(np.int64),
                         "si_b": vb.astype(np.int64)})


# ----------------------------------------------------------------------------------------------- statistics
def bspline_basis(x: np.ndarray, knots_inner: np.ndarray, lo: float, hi: float, degree: int = 3) -> np.ndarray:
    """B-spline design (without the first column, to keep the intercept identifiable)."""
    t = np.r_[[lo] * (degree + 1), knots_inner, [hi] * (degree + 1)]
    nb = len(t) - degree - 1
    xx = np.clip(x, lo, hi - 1e-9)
    B = np.column_stack([BSpline(t, np.eye(nb)[k], degree, extrapolate=False)(xx) for k in range(nb)])
    B = np.nan_to_num(B)
    return B[:, 1:]


def fit_logistic(X: np.ndarray, succ: np.ndarray, trials: np.ndarray, w: np.ndarray | None = None,
                 n_iter: int = 50, ridge: float = 1e-6) -> np.ndarray:
    """Binomial GLM (logit) by IRLS on aggregated cells: succ / trials with optional cell weights w."""
    w = np.ones(len(succ)) if w is None else w
    keep = (trials > 0) & (w > 0)
    X, s, n, w = X[keep], succ[keep], trials[keep], w[keep]
    beta = np.zeros(X.shape[1])
    p0 = np.clip(s.sum() / n.sum(), 1e-6, 1 - 1e-6)
    beta[0] = np.log(p0 / (1 - p0))
    for _ in range(n_iter):
        eta = np.clip(X @ beta, -30, 30)
        mu = 1 / (1 + np.exp(-eta))
        W = w * n * mu * (1 - mu) + 1e-12
        z = eta + (s / n - mu) / (mu * (1 - mu) + 1e-12)
        XtW = X.T * W
        new = np.linalg.solve(XtW @ X + ridge * np.eye(X.shape[1]), XtW @ z)
        if np.max(np.abs(new - beta)) < 1e-8:
            beta = new
            break
        beta = new
    return beta


def expit(x):
    return 1 / (1 + np.exp(-np.clip(x, -30, 30)))


def ortholog_effect(cells: pd.DataFrame, knots: np.ndarray, w: np.ndarray | None = None, ref_mask=None) -> dict:
    """Logistic regression  logit P(match) = spline(identity) + beta * ortholog  on aggregated cells.

    cells: identity (cell midpoint), rel (1/2), succ, trials. Returns log-OR beta and the standardised marginal
    difference (ortholog minus paralog) averaged over the identity distribution of the paralog cells (trials-weighted).
    """
    x = cells["identity"].to_numpy(float)
    B = bspline_basis(x, knots, 0.0, 1.0)
    orth = (cells["rel"].to_numpy() == REL_ORTHO).astype(float)
    X = np.column_stack([np.ones(len(x)), B, orth])
    beta = fit_logistic(X, cells["succ"].to_numpy(float), cells["trials"].to_numpy(float), w)
    ww = cells["trials"].to_numpy(float) * (np.ones(len(x)) if w is None else w)
    ref = (orth == 0) if ref_mask is None else ref_mask
    X1, X0 = X.copy(), X.copy()
    X1[:, -1], X0[:, -1] = 1, 0
    wr = ww * ref
    d = float(np.sum(wr * (expit(X1 @ beta) - expit(X0 @ beta))) / max(wr.sum(), 1e-12))
    return {"log_or": float(beta[-1]), "or": float(np.exp(beta[-1])), "marginal_diff": d, "beta": beta}


def binned_effect(cells: pd.DataFrame, edges: np.ndarray, w: np.ndarray | None = None) -> dict:
    """Per identity bin: pooled match rate by relation; standardised difference with paralog bin weights."""
    w = np.ones(len(cells)) if w is None else w
    b = np.clip(np.digitize(cells["identity"].to_numpy(float), edges) - 1, 0, len(edges) - 2)
    rel = cells["rel"].to_numpy()
    s = cells["succ"].to_numpy(float) * w
    n = cells["trials"].to_numpy(float) * w
    nb = len(edges) - 1
    S = np.zeros((2, nb))
    N = np.zeros((2, nb))
    for r, ri in ((REL_ORTHO, 0), (REL_PARA, 1)):
        m = rel == r
        np.add.at(S[ri], b[m], s[m])
        np.add.at(N[ri], b[m], n[m])
    with np.errstate(invalid="ignore", divide="ignore"):
        R = S / N
    ok = (N[0] > 0) & (N[1] > 0)
    wt = N[1] * ok
    diff = float(np.nansum((R[0] - R[1]) * wt) / wt.sum()) if wt.sum() > 0 else float("nan")
    return {"rate_ortho": R[0], "rate_para": R[1], "n_ortho": N[0], "n_para": N[1], "std_diff": diff}


def bootstrap_weights(og_codes: np.ndarray, n_og: int, rep: int, seed: int) -> np.ndarray:
    """Per-row weights = multiplicity of its OG in a bootstrap resample of OGs (seed rule seed*1000+rep)."""
    rng = np.random.default_rng(seed * 1000 + rep)
    mult = np.bincount(rng.integers(0, n_og, n_og), minlength=n_og).astype(float)
    return mult[og_codes]


def pooled_ratio(num: np.ndarray, den: np.ndarray, w: np.ndarray | None = None) -> float:
    w = np.ones(len(num)) if w is None else w
    d = float(np.sum(den * w))
    return float(np.sum(num * w)) / d if d > 0 else float("nan")
