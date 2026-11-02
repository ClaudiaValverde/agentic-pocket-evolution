"""Part A estimands (DESIGN §4.2; round-2 amendments 1, 3, 4; reviews/round2/pi_statistics.md §2).

Inputs
  sig      significant site pairs: columns i, j, cls (final class name), gi, gj (root groups), and optionally
           any weighting columns; H0 and 'unknown' pairs may be present (they are excluded from every F).
  totals   exact class totals: H1, H2, H12, H3j, rest  (homology.exact_class_totals)
  split    fraction of H1/H2/H12 pairs that are same-site (positional homologs, 's'), per class; the rest are
           'x' (offsite). From SOGs when available, else the placeholder estimate (flagged).
  sample   the seeded uniform sample of 'rest' pairs with search-resolved class in {H3s, H4, N, unknown}
           and their groups gi, gj (bootstrap reweighting).
  q0       per-pair probability that a non-homologous pair is called significant (calibrated null).

Definitions
  K_c   = Σ_{significant pairs in c} weight            (per pair: weight 1)
  n_c   = exact (H1/H2/H12 × split, H3j) or p̂_c · rest  (H3s, H4, N, unknown)
  E0_c  = q0 · n_c                                       (expected false positives)
  r_c   = (K_c − E0_c) / n_c ; raw rate K_c / n_c also reported
  F_anc = Σ_{H1s,H2s,H12s,H3j,H3s} (K−E0)_+ / Σ_{all but H0, unknown} (K−E0)_+ ; F_offsite (H1x,H2x,H12x),
  F_H4, F_analog (N) likewise; the four sum to 1.
  Weightings: per_pair; per_site (pair weight = (1/deg_i + 1/deg_j)/2); per_family_pair (each unordered root-group
  cell with ≥1 significant pair has total weight 1, split equally over its pairs). For the two weighted versions the
  expected-FP subtraction uses weight 1 per expected FP pair (the maximum a pair can get; conservative for N).
"""
from __future__ import annotations

import numpy as np
import pandas as pd

ANC = ["H1s", "H2s", "H12s", "H3j", "H3s"]
OFFSITE = ["H1x", "H2x", "H12x"]
H4 = ["H4"]
ANALOG = ["N"]
ALL = ANC + OFFSITE + H4 + ANALOG
REPORTED = ALL + ["unknown"]
SAMPLE_CLASSES = ["H3s", "H4", "N", "unknown"]


def pair_group_weight(w: np.ndarray, gi: np.ndarray, gj: np.ndarray) -> np.ndarray:
    """Bootstrap weight of a pair: w_g w_h across groups, w_g within a group (see homology.exact_class_totals)."""
    return np.where(gi == gj, w[gi], w[gi] * w[gj])


def pair_weights(sig: pd.DataFrame, scheme: str) -> np.ndarray:
    keep = sig.cls.isin(ALL).to_numpy()
    w = np.zeros(len(sig))
    if scheme == "per_pair":
        w[keep] = 1.0
    elif scheme == "per_site":
        ii, jj = sig.i.to_numpy()[keep], sig.j.to_numpy()[keep]
        allv = np.concatenate([ii, jj])
        u, inv, cnt = np.unique(allv, return_inverse=True, return_counts=True)
        deg = cnt[inv]
        w[keep] = (1 / deg[:len(ii)] + 1 / deg[len(ii):]) / 2
    elif scheme == "per_family_pair":
        a = np.minimum(sig.gi.to_numpy(), sig.gj.to_numpy())[keep]
        b = np.maximum(sig.gi.to_numpy(), sig.gj.to_numpy())[keep]
        cell = pd.Series(a.astype(np.int64)) * (1 << 31) + b
        cnt = cell.map(cell.value_counts()).to_numpy()
        w[keep] = 1.0 / cnt
    else:
        raise ValueError(scheme)
    if "search_weight" in sig:      # Horvitz-Thompson weight of sampled searches (06_classify tiered search)
        w = w * sig.search_weight.to_numpy()
    if "pair_w" in sig:             # optional site-weight product w_i w_j (redundancy-reweighting sensitivity)
        w = w * sig.pair_w.to_numpy()
    return w


def class_n(totals: dict, split: dict, sample: pd.DataFrame, sample_w: np.ndarray | None = None) -> dict:
    n = {}
    for c in ("H1", "H2", "H12"):
        f = split.get(c, 0.0)
        n[c + "s"] = totals[c] * f
        n[c + "x"] = totals[c] * (1 - f)
    n["H3j"] = totals["H3j"]
    sw = np.ones(len(sample)) if sample_w is None else sample_w
    if "pair_w" in sample:
        sw = sw * sample.pair_w.to_numpy()
    tot = sw.sum()
    for c in SAMPLE_CLASSES:
        p = sw[(sample.cls == c).to_numpy()].sum() / tot if tot > 0 else 0.0
        n[c] = p * totals["rest"]
    return n


def estimands(sig: pd.DataFrame, totals: dict, split: dict, sample: pd.DataFrame, q0: float,
              group_w: np.ndarray | None = None, schemes=("per_pair", "per_site", "per_family_pair"),
              base_weights: dict | None = None) -> dict:
    """One evaluation (point estimate, or one bootstrap replicate if group_w is given)."""
    if group_w is not None:
        gw_sig = pair_group_weight(group_w, sig.gi.to_numpy(), sig.gj.to_numpy())
        gw_smp = pair_group_weight(group_w, sample.gi.to_numpy(), sample.gj.to_numpy())
    else:
        gw_sig = np.ones(len(sig))
        gw_smp = None
    n = class_n(totals, split, sample, gw_smp)
    out = {"n": n, "q0": q0}
    cls = sig.cls.to_numpy()
    for sch in schemes:
        bw = base_weights[sch] if base_weights and sch in base_weights else pair_weights(sig, sch)
        w = bw * gw_sig
        K = {c: float(w[cls == c].sum()) for c in REPORTED}
        E0 = {c: q0 * n.get(c, 0.0) for c in REPORTED}
        adj = {c: max(K[c] - E0[c], 0.0) for c in ALL}
        tot_adj = sum(adj.values())
        tot_raw = sum(K[c] for c in ALL)
        res = {"K": K, "E0": E0,
               "F_adj": _shares(adj, tot_adj), "F_raw": _shares({c: K[c] for c in ALL}, tot_raw)}
        if sch == "per_pair":
            res["r_raw"] = {c: (K[c] / n[c] if n.get(c, 0) > 0 else np.nan) for c in REPORTED}
            res["r_adj"] = {c: ((K[c] - E0[c]) / n[c] if n.get(c, 0) > 0 else np.nan) for c in REPORTED}
        out[sch] = res
    return out


def _shares(x: dict, tot: float) -> dict:
    g = {"F_anc": ANC, "F_offsite": OFFSITE, "F_H4": H4, "F_analog": ANALOG}
    out = {k: (sum(x[c] for c in v) / tot if tot > 0 else np.nan) for k, v in g.items()}
    out.update({f"share_{c}": (x[c] / tot if tot > 0 else np.nan) for c in ALL})
    return out


def flatten(d: dict, prefix="") -> dict:
    out = {}
    for k, v in d.items():
        key = f"{prefix}{k}"
        if isinstance(v, dict):
            out.update(flatten(v, key + "."))
        else:
            out[key] = v
    return out


def bootstrap(sig, totals_fn, split, sample, q0, n_groups: int, n_rep: int, seed: int, schemes, log=None):
    """Two-way (crossed, pigeonhole) bootstrap over root groups with Poisson(1) multiplicities.

    totals_fn(w) -> exact totals under group weights w (recomputed jointly with K_c and n_c, as the review asks).
    Seeds: seed*1000 + replicate (prereg meta.bootstrap_seed_rule).
    """
    base = {s: pair_weights(sig, s) for s in schemes}
    reps = []
    for r in range(n_rep):
        rng = np.random.default_rng(seed * 1000 + r)
        w = rng.poisson(1.0, n_groups).astype(np.float64)
        tot = totals_fn(w)
        est = estimands(sig, tot, split, sample, q0, group_w=w, schemes=schemes, base_weights=base)
        reps.append(flatten(est))
        if log and (r + 1) % 20 == 0:
            log(f"bootstrap {r + 1}/{n_rep}")
    return pd.DataFrame(reps)


def ci_table(point: dict, reps: pd.DataFrame, alpha=0.05) -> pd.DataFrame:
    flat = flatten(point)
    rows = []
    for k, v in flat.items():
        if k in reps:
            col = pd.to_numeric(reps[k], errors="coerce").dropna()
            lo, hi = (np.quantile(col, [alpha / 2, 1 - alpha / 2]) if len(col) else (np.nan, np.nan))
        else:
            lo = hi = np.nan
        rows.append({"key": k, "value": v, "ci_lo": lo, "ci_hi": hi})
    return pd.DataFrame(rows)


def wilson(k: float, n: float, z: float = 1.96) -> tuple[float, float]:
    if n <= 0:
        return (np.nan, np.nan)
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return (c - h, c + h)
