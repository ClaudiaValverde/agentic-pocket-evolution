"""Covariate-conditioned null for site-pair similarity (DESIGN §4.2 + round-2 amendment 2).

Step 1  covariates per site f (6): log site volume (sph8 site_volume_A3 of the rep pocket), burrow score of the rep,
        lining polarity, lining aromatic fraction (site residues -> amino acids of the UniProt sequence),
        mean pLDDT/100 of the rep, model source (1 = Swiss-Prot v6 tar [tier 2], 0 = proteome v6 tar [tier 1]).
Step 2  regression on 10^6 seeded random non-homologous pairs (seed 0):
            mu(f_i, f_j)        = c + a.(f_i + f_j) + b.(f_i^2 + f_j^2) + f_i^T W f_j      (W symmetric)
            log|s - mu|         = same design   ->  sigma = exp(pred + 0.6352)  (Gaussian E log|e| = log sd - 0.6352)
        The design is symmetric in (i, j) and low rank, so it is evaluated as matrix products on GPU.
Step 3  per query i: z_ij = (s_ij - mu_ij) / sigma_ij against a fixed seeded background B of random sites, with
        homologs of i masked (same protein, same root group, any shared protein-level label); threshold u_i = 95 %
        quantile; GPD fitted to the exceedances by probability-weighted moments (Hosking & Wallis 1987).
        p_i(z) = zeta * (1 + xi (z - u)/sc)^(-1/xi) for z > u (zeta = exceedance fraction), else empirical;
        E_ij = N_db * p_i(z_ij). Pair significant iff max(E_ij, E_ji) < threshold.
Tail variants (tail_variant): 'pooled' (primary, declared after the toy test showed per-query shape noise biases
          extreme tails, before any real-data result), 'perq' (prereg literal), 'xi0' (conservative);
          score 'cos' (primary) / 'csls' (sensitivity: s' = 2 s - r_i - r_j, r = mean cos of the 10 NN).
"""
from __future__ import annotations

import numpy as np

POLAR = set("DEHKNQRSTY")
AROMATIC = set("FWY")
COV_NAMES = ["log_vol", "burrow_score", "polar_frac", "arom_frac", "plddt", "src_sprot"]
GAUSS_LOGABS = 0.6352  # -E[log|Z|] for Z ~ N(0,1) = (gamma + log 2)/2


def lining_fractions(residues: list, seq: str | None) -> tuple[float, float, bool]:
    """(polar fraction, aromatic fraction, ok) of the site residues (1-based UniProt positions)."""
    if seq is None or len(residues) == 0:
        return np.nan, np.nan, False
    aa = [seq[r - 1] for r in residues if 1 <= r <= len(seq)]
    if len(aa) < len(residues) or not aa:
        return np.nan, np.nan, False
    return sum(a in POLAR for a in aa) / len(aa), sum(a in AROMATIC for a in aa) / len(aa), True


# ----------------------------------------------------------------------------------------------- design
def n_design(d: int) -> int:
    return 1 + d + d + d * (d + 1) // 2


def pair_design(fi: np.ndarray, fj: np.ndarray) -> np.ndarray:
    """Explicit symmetric design matrix rows for pairs (used for fitting only)."""
    d = fi.shape[1]
    cols = [np.ones(len(fi)), *(fi + fj).T, *(fi ** 2 + fj ** 2).T]
    for k in range(d):
        for l in range(k, d):
            cols.append(fi[:, k] * fj[:, l] + fi[:, l] * fj[:, k] if k != l else fi[:, k] * fj[:, k])
    return np.stack(cols, 1)


def unpack(beta: np.ndarray, d: int):
    """beta -> (c, a, b, W) with mu = c + a.(fi+fj) + b.(fi^2+fj^2) + fi^T W fj."""
    c = beta[0]
    a = beta[1:1 + d]
    b = beta[1 + d:1 + 2 * d]
    W = np.zeros((d, d))
    k = 1 + 2 * d
    for r in range(d):
        for s in range(r, d):
            if r == s:
                W[r, r] = beta[k]
            else:
                W[r, s] = W[s, r] = beta[k]
            k += 1
    return c, a, b, W


class CovModel:
    """Fitted mean / scale model; standardise() works on numpy or torch arrays."""

    def __init__(self, beta_mu, beta_ls, mean, sd, d):
        self.beta_mu, self.beta_ls = np.asarray(beta_mu), np.asarray(beta_ls)
        self.mean, self.sd, self.d = np.asarray(mean), np.asarray(sd), d

    def to_dict(self):
        return {"beta_mu": self.beta_mu.tolist(), "beta_logsd": self.beta_ls.tolist(), "cov_mean": self.mean.tolist(),
                "cov_sd": self.sd.tolist(), "d": self.d, "cov_names": COV_NAMES[:self.d]}

    @staticmethod
    def from_dict(r):
        return CovModel(r["beta_mu"], r["beta_logsd"], r["cov_mean"], r["cov_sd"], r["d"])

    def norm(self, F):
        return (F - self.mean) / self.sd

    @staticmethod
    def fit(Fi, Fj, s, mean=None, sd=None, ridge=1e-6):
        d = Fi.shape[1]
        if mean is None:
            allf = np.concatenate([Fi, Fj])
            mean, sd = allf.mean(0), allf.std(0) + 1e-9
        fi, fj = (Fi - mean) / sd, (Fj - mean) / sd
        X = pair_design(fi, fj)
        A = X.T @ X + ridge * np.eye(X.shape[1])
        bmu = np.linalg.solve(A, X.T @ s)
        r = s - X @ bmu
        bls = np.linalg.solve(A, X.T @ np.log(np.abs(r) + 1e-6))
        return CovModel(bmu, bls, mean, sd, d)

    def mu_logsd_pairs(self, fi, fj):
        """Elementwise for aligned pairs (normalised covariates). numpy or torch."""
        out = []
        for beta in (self.beta_mu, self.beta_ls):
            c, a, b, W = unpack(beta, self.d)
            lib = _lib(fi)
            a_, b_, W_ = (lib_as(fi, a), lib_as(fi, b), lib_as(fi, W))
            v = c + (fi + fj) @ a_ + (fi ** 2 + fj ** 2) @ b_ + ((fi @ W_) * fj).sum(-1)
            out.append(v)
        return out[0], out[1]

    def mu_logsd_block(self, fq, fb):
        """Dense (Q, B) matrices for query block × background (torch)."""
        out = []
        for beta in (self.beta_mu, self.beta_ls):
            c, a, b, W = unpack(beta, self.d)
            a_, b_, W_ = lib_as(fq, a), lib_as(fq, b), lib_as(fq, W)
            uq = fq @ a_ + (fq ** 2) @ b_
            ub = fb @ a_ + (fb ** 2) @ b_
            out.append(c + uq[:, None] + ub[None, :] + (fq @ W_) @ fb.T)
        return out[0], out[1]

    def z(self, s, mu, logsd):
        lib = _lib(s)
        return (s - mu) / lib.exp(logsd + GAUSS_LOGABS)


def _lib(x):
    import numpy
    if type(x).__module__.startswith("torch"):
        import torch
        return torch
    return numpy


def lib_as(ref, arr):
    if type(ref).__module__.startswith("torch"):
        import torch
        return torch.as_tensor(np.asarray(arr), dtype=ref.dtype, device=ref.device)
    return np.asarray(arr, dtype=ref.dtype)


# ----------------------------------------------------------------------------------------------- GPD
def gpd_pwm_rows(Z, valid, q=0.95):
    """Row-wise GPD fit (torch, GPU). Z (Q, B) z-scores; valid (Q, B) bool mask of usable background.

    Returns dict of (Q,) tensors: u (threshold), sc, xi, zeta (exceedance fraction), n (valid count), n_exc.
    PWM (Hosking & Wallis 1987): a0 = mean(y), a1 = mean((1 - p_k) y_(k)), p_k = (k - 0.35)/m (ascending);
    xi = 2 - a0/(a0 - 2 a1), sc = 2 a0 a1/(a0 - 2 a1).
    """
    import torch
    Zm = torch.where(valid, Z, torch.full_like(Z, -torch.inf))
    n = valid.sum(1)
    Zs, _ = torch.sort(Zm, dim=1, descending=True)          # valid first
    m = torch.clamp((n.float() * (1 - q)).floor().long(), min=20)
    mmax = int(m.max())
    top = Zs[:, :mmax + 1]                                   # top m+1 values
    u = top.gather(1, m[:, None]).squeeze(1)                 # (m+1)-th largest = threshold
    ar = torch.arange(mmax, device=Z.device)[None, :]
    inm = ar < m[:, None]
    y = (top[:, :mmax] - u[:, None]).clamp(min=0) * inm      # descending exceedances
    mf = m.float()
    a0 = y.sum(1) / mf
    # ascending rank k = m - idx  (idx 0 = largest) ; p_k = (k - 0.35)/m ; weight (1 - p_k)
    k_asc = (mf[:, None] - ar.float())
    wgt = (1 - (k_asc - 0.35) / mf[:, None]) * inm
    a1 = (wgt * y).sum(1) / mf
    den = a0 - 2 * a1
    xi = 2 - a0 / den
    sc = 2 * a0 * a1 / den
    zeta = mf / n.float()
    return {"u": u, "sc": sc, "xi": xi, "zeta": zeta, "n": n, "n_exc": m, "ybar": a0}


def tail_variant(fit: dict, variant: str, xi_pooled: float | None = None):
    """(sc, xi) for a tail variant.
    'perq'   per-query PWM shape and scale (prereg literal; noisy shape -> light-tail bias when extrapolating)
    'pooled' shape shared by all queries (median of per-query PWM shapes), per-query scale from the
             exceedance mean given that shape: E[y] = sc / (1 - xi)  ->  sc = ybar (1 - xi)
    'xi0'    shape clipped at >= 0; clipped queries use the exponential tail with sc = ybar (conservative for
             bounded tails)
    """
    import torch
    if variant == "perq":
        return fit["sc"], fit["xi"]
    if variant == "pooled":
        xi = torch.full_like(fit["xi"], float(xi_pooled))
        return fit["ybar"] * (1 - xi), xi
    if variant == "xi0":
        neg = fit["xi"] < 0
        return torch.where(neg, fit["ybar"], fit["sc"]), fit["xi"].clamp(min=0)
    raise ValueError(variant)


def gpd_sf(z, u, sc, xi, zeta, xi_floor=None):
    """Tail probability P(Z >= z) for z > u under the fitted GPD (torch or numpy, broadcastable).
    For z <= u returns zeta (caller substitutes the empirical value when needed)."""
    lib = _lib(z)
    if xi_floor is not None:
        xi = lib.clip(xi, xi_floor, None) if lib is np else xi.clamp(min=xi_floor)
    y = (z - u)
    y = y.clamp(min=0) if lib is not np else np.clip(y, 0, None)
    small = lib.abs(xi) < 1e-6
    t = 1 + xi * y / sc
    if lib is np:
        t = np.clip(t, 0, None)
        with np.errstate(divide="ignore", invalid="ignore"):
            pw = np.where(small, np.exp(-y / sc), np.where(t > 0, t ** (-1 / np.where(small, 1, xi)), 0.0))
    else:
        import torch
        t = t.clamp(min=0)
        xs = torch.where(small, torch.ones_like(xi), xi)
        pw = torch.where(small, torch.exp(-y / sc), torch.where(t > 0, t.clamp(min=1e-30) ** (-1 / xs),
                                                              torch.zeros_like(t)))
    return zeta * pw
