#!/usr/bin/env python
"""Part A stage 4: covariate-conditioned null (round-2 amendment 2) + calibration.

  1. 10^6 random non-homologous site pairs (seed 0) -> fit CovModel (mean and log-scale); size-only baseline.
  2. Background B: seeded random sites (prereg n_background_pockets_per_query x 1.2 drawn, homologs masked per
     query -> >= 10^4 usable for almost all queries; shortfalls counted).
  3. Per query: z against B, GPD (PWM) on the top 5 % -> per-site tail parameters.
  4. Calibration: held-out 10^6 random non-homologous pairs (seed 1, B excluded): QQ data, observed vs expected
     counts at p and E thresholds, per tail variant and direction; extended tail check on 10^9 random pairs (seed 2).
Score variants: --score cos (primary) | csls (sensitivity).
Outputs: <work>/null_<score>/{covmodel.json, tail.npz, calib_*.parquet}; results/partA/<snap>/null_<score>_*.csv
"""
import json
import time

import numpy as np
import pandas as pd
import torch

from _common import parser, resolve, logger, pairs
from pocketevo import nulls

ap = parser(__doc__)
ap.add_argument("--score", choices=["cos", "csls"], default="cos")
ap.add_argument("--qblock", type=int, default=8192)
ap.add_argument("--n-fit", type=int, default=1_000_000)
ap.add_argument("--n-heldout", type=int, default=1_000_000)
ap.add_argument("--n-extended", type=int, default=1_000_000_000)
a = resolve(ap.parse_args())
log = logger(a.work / "partA.log")
out = a.work / f"null_{a.score}"
out.mkdir(exist_ok=True)
dev = torch.device("cuda:0")
T0 = time.time()
pre = a.prereg["pairs"].get("null", a.prereg["pairs"].get(None))   # YAML parses the key `null:` as None
NBG = int(pre["n_background_pockets_per_query"]["value"])
QTH = float(pre["tail_fit"]["threshold_quantile"]["value"])

X = np.load(a.work / "fp.f16.npy", mmap_mode="r")
N = X.shape[0]
Xg = torch.empty((N, X.shape[1]), dtype=torch.float16, device=dev)
for s in range(0, N, 500_000):
    Xg[s:s + 500_000] = torch.from_numpy(np.ascontiguousarray(X[s:s + 500_000]))
ctx = np.load(a.work / "context.npz")
site_prot, prot_group, prot_lab, cov_ok = ctx["site_prot"], ctx["prot_group"], ctx["prot_lab"], ctx["cov_ok"]
cov = np.load(a.work / "covariates.npy")
cov = np.where(np.isfinite(cov), cov, np.nanmean(cov, 0)).astype(np.float64)   # imputed rows flagged by cov_ok
site_group = prot_group[site_prot]
site_plab = prot_lab[site_prot]
sp_g = torch.from_numpy(site_prot).to(dev)
sg_g = torch.from_numpy(site_group).to(dev)
spl_g = torch.from_numpy(site_plab.astype(np.int32)).to(dev)
N_DB = N
csls_r = None
if a.score == "csls":
    D = np.load(a.work / "knn" / "D.f16.npy", mmap_mode="r")
    I = np.load(a.work / "knn" / "I.i32.npy", mmap_mode="r")
    csls_r = np.empty(N, np.float32)
    for s in range(0, N, 200_000):
        d = D[s:s + 200_000, :11].astype(np.float32)
        selfcol = I[s:s + 200_000, 0] == np.arange(s, min(N, s + 200_000))
        r10 = np.where(selfcol, d[:, 1:11].mean(1), d[:, :10].mean(1))
        csls_r[s:s + 200_000] = r10
    csls_rg = torch.from_numpy(csls_r).to(dev)


def score_pairs(i, j):
    """similarity for aligned pair index tensors (GPU)."""
    s = (Xg[i].float() * Xg[j].float()).sum(1)
    if csls_r is not None:
        s = 2 * s - csls_rg[i] - csls_rg[j]
    return s


def nonhom_pairs(i, j):
    """numpy: True where (i, j) are on different proteins, different root groups and share no protein label."""
    ok = (site_prot[i] != site_prot[j]) & (site_group[i] != site_group[j])
    A, B = site_plab[i], site_plab[j]
    for k in range(A.shape[1]):
        col = A[:, k][:, None]
        ok &= ~(((col == B) & (col >= 0)).any(1))
    return ok


def draw_pairs(n, seed, exclude=None):
    rng = np.random.default_rng(seed)
    pool = np.flatnonzero(cov_ok if exclude is None else (cov_ok & ~exclude))
    got_i, got_j, drawn = [], [], 0
    while sum(len(x) for x in got_i) < n:
        m = int((n - sum(len(x) for x in got_i)) * 1.3) + 1000
        i, j = pool[rng.integers(0, len(pool), m)], pool[rng.integers(0, len(pool), m)]
        drawn += m
        ok = nonhom_pairs(i, j)
        got_i.append(i[ok]); got_j.append(j[ok])
    i, j = np.concatenate(got_i)[:n], np.concatenate(got_j)[:n]
    return i, j, drawn


def sims_np(i, j):
    out = np.empty(len(i), np.float64)
    for s in range(0, len(i), 2_000_000):
        ii = torch.from_numpy(i[s:s + 2_000_000]).to(dev)
        jj = torch.from_numpy(j[s:s + 2_000_000]).to(dev)
        out[s:s + 2_000_000] = score_pairs(ii, jj).double().cpu().numpy()
    return out


summary = {"score": a.score, "N_db": N_DB}
# ---------------------------------------------------------------- 1. regression fit
mpath = out / "covmodel.json"
if mpath.exists():
    model = nulls.CovModel.from_dict(json.loads(mpath.read_text()))
    log("covmodel loaded")
else:
    t = time.time()
    fi, fj, drawn = draw_pairs(a.n_fit, seed=0)
    s = sims_np(fi, fj)
    model = nulls.CovModel.fit(cov[fi], cov[fj], s)
    mu, ls = model.mu_logsd_pairs(model.norm(cov[fi]), model.norm(cov[fj]))
    r2 = 1 - np.var(s - mu) / np.var(s)
    # size-only baseline: log volume only
    sm = nulls.CovModel.fit(cov[fi][:, :1], cov[fj][:, :1], s)
    mu1, _ = sm.mu_logsd_pairs(sm.norm(cov[fi][:, :1]), sm.norm(cov[fj][:, :1]))
    r2_size = 1 - np.var(s - mu1) / np.var(s)
    z = model.z(s, mu, ls)
    summary["fit"] = {"n": a.n_fit, "drawn": drawn, "accept_rate": a.n_fit / drawn, "r2_cov": r2, "r2_size_only": r2_size,
                      "s_mean": s.mean(), "s_p99": np.quantile(s, 0.99), "z_sd": float(np.std(z)), "wall_s": time.time() - t}
    mpath.write_text(json.dumps({**model.to_dict(), "size_only": sm.to_dict(), "fit": summary["fit"]}, default=float))
    log("covmodel", json.dumps(summary["fit"], default=float))

# ---------------------------------------------------------------- 2. background
rng = np.random.default_rng(0)
B = np.sort(rng.choice(np.flatnonzero(cov_ok), int(NBG * 1.2), replace=False))
in_B = np.zeros(N, bool); in_B[B] = True
Fn = torch.from_numpy(model.norm(cov).astype(np.float32)).to(dev)
XB = Xg[B]
FB = Fn[B]
gB = sg_g[B]; pB = sp_g[B]
nlab = int(site_plab.max()) + 2
LB = torch.zeros((nlab, len(B)), dtype=torch.bool, device=dev)      # label incidence of background sites
plB = spl_g[B]
for k in range(plB.shape[1]):
    c = plB[:, k]
    m = c >= 0
    LB[c[m].long(), torch.nonzero(m).squeeze(1)] = True
LB[nlab - 1] = False   # code -1 -> row nlab-1 (empty)
rB = csls_rg[B] if csls_r is not None else None


def block_z(q):
    """(Q, |B|) z-scores and validity mask for query index tensor q."""
    S = Xg[q].float() @ XB.float().T
    if rB is not None:
        S = 2 * S - csls_rg[q][:, None] - rB[None, :]
    mu, ls = model.mu_logsd_block(Fn[q], FB)
    Z = model.z(S, mu, ls)
    valid = (sp_g[q][:, None] != pB[None, :]) & (sg_g[q][:, None] != gB[None, :])
    ql = spl_g[q].long()
    ql = torch.where(ql < 0, torch.full_like(ql, nlab - 1), ql)
    for k in range(ql.shape[1]):
        valid &= ~LB[ql[:, k]]
    return Z, valid


tail_path = out / "tail.npz"
if tail_path.exists():
    tl = dict(np.load(tail_path))
    log("tail loaded")
else:
    t = time.time()
    keys = ["u", "sc", "xi", "zeta", "n", "n_exc", "ybar"]
    tl = {k: np.empty(N, np.float32 if k not in ("n", "n_exc") else np.int32) for k in keys}
    for s in range(0, N, a.qblock):
        q = torch.arange(s, min(N, s + a.qblock), device=dev)
        Z, valid = block_z(q)
        fit = nulls.gpd_pwm_rows(Z, valid, q=QTH)
        for k in keys:
            tl[k][s:s + len(q)] = fit[k].cpu().numpy()
        if (s // a.qblock) % 50 == 0:
            log(f"tail {s}/{N} {time.time() - t:.0f}s")
    xi_pool = float(np.median(tl["xi"][cov_ok]))
    tl["xi_pooled"] = np.array(xi_pool)
    np.savez(tail_path, **tl)
    summary["tail"] = {"wall_s": time.time() - t, "xi_pooled": xi_pool,
                       "xi_q": np.quantile(tl["xi"], [0.05, 0.25, 0.5, 0.75, 0.95]).tolist(),
                       "n_bg_q": np.quantile(tl["n"], [0, 0.001, 0.01, 0.5]).tolist(),
                       "n_queries_bg_lt_10000": int((tl["n"] < NBG).sum())}
    log("tail", json.dumps(summary["tail"]))
xi_pool = float(tl["xi_pooled"])
TL = {k: torch.from_numpy(np.asarray(tl[k], np.float32)).to(dev) for k in ("u", "sc", "xi", "zeta", "ybar")}
VARIANTS = ("pooled", "perq", "xi0")


def tail_p(q, z, variant):
    fit = {k: TL[k][q] for k in TL}
    sc, xi = nulls.tail_variant(fit, variant, xi_pool)
    return nulls.gpd_sf(z, fit["u"], sc, xi, fit["zeta"])


# ---------------------------------------------------------------- 3. held-out calibration
cal_path = out / "calib_heldout.parquet"
if not cal_path.exists():
    t = time.time()
    hi, hj, drawn = draw_pairs(a.n_heldout, seed=1, exclude=in_B)
    rec = {"i": hi, "j": hj}
    s = sims_np(hi, hj)
    mu, ls = model.mu_logsd_pairs(model.norm(cov[hi]), model.norm(cov[hj]))
    z = model.z(s, mu, ls)
    rec["s"], rec["z"] = s, z
    # empirical p (both directions) against the query's own masked background, and GPD tail p
    for dname, (qa, za) in {"ij": (hi, z), "ji": (hj, z)}.items():
        pe = np.empty(len(qa)); pt = {v: np.empty(len(qa)) for v in VARIANTS}
        for s0 in range(0, len(qa), a.qblock):
            q = torch.from_numpy(qa[s0:s0 + a.qblock]).to(dev)
            zz = torch.from_numpy(za[s0:s0 + a.qblock]).to(dev).float()
            Z, valid = block_z(q)
            ge = ((Z >= zz[:, None]) & valid).sum(1).float()
            pe[s0:s0 + len(q)] = ((ge + 1) / (valid.sum(1).float() + 1)).cpu().numpy()
            for v in VARIANTS:
                pt[v][s0:s0 + len(q)] = tail_p(q, zz, v).cpu().numpy()
        rec[f"pemp_{dname}"] = pe
        for v in VARIANTS:
            # combined p: empirical below the threshold, GPD above (as used for significance)
            above = za > tl["u"][qa]
            rec[f"p_{v}_{dname}"] = np.where(above, pt[v], pe)
    cal = pd.DataFrame(rec)
    cal.to_parquet(cal_path, index=False)
    log(f"held-out calibration {time.time() - t:.0f}s")
cal = pd.read_parquet(cal_path)
rows = []
for v in VARIANTS:
    for t in (1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6):
        n = len(cal)
        o_ij = int((cal[f"p_{v}_ij"] < t).sum()); o_ji = int((cal[f"p_{v}_ji"] < t).sum())
        o_both = int(((cal[f"p_{v}_ij"] < t) & (cal[f"p_{v}_ji"] < t)).sum())
        rows.append({"variant": v, "p_threshold": t, "E_threshold": t * N_DB, "n_pairs": n, "expected_one_dir": n * t,
                     "obs_ij": o_ij, "obs_ji": o_ji, "obs_both": o_both,
                     "ratio_ij": o_ij / (n * t), "ratio_ji": o_ji / (n * t)})
calt = pd.DataFrame(rows)
calt.to_csv(a.results / f"null_{a.score}_calibration_heldout.csv", index=False)
log("held-out calibration table\n" + calt.to_string())

# ---------------------------------------------------------------- 4. extended tail check (10^9 pairs)
ext_path = a.results / f"null_{a.score}_calibration_extended.csv"
if not ext_path.exists() and a.n_extended > 0:
    t = time.time()
    rng2 = np.random.default_rng(2)
    pool = np.flatnonzero(cov_ok & ~in_B)
    pthr = np.array(sorted({1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1.0 / N_DB, 1e-2 / N_DB, 1e-4 / N_DB},
                           reverse=True))
    cnt = {v: {d: np.zeros(len(pthr), np.int64) for d in ("ij", "ji", "both")} for v in VARIANTS}
    n_done, n_drawn, chunk = 0, 0, 4_000_000
    while n_done < a.n_extended:
        i = pool[rng2.integers(0, len(pool), chunk)]; j = pool[rng2.integers(0, len(pool), chunk)]
        n_drawn += chunk
        ok = nonhom_pairs(i, j)
        i, j = i[ok], j[ok]
        i, j = i[:a.n_extended - n_done], j[:a.n_extended - n_done]
        ig, jg = torch.from_numpy(i).to(dev), torch.from_numpy(j).to(dev)
        s = score_pairs(ig, jg)
        mu, ls = model.mu_logsd_pairs(Fn[ig], Fn[jg])
        z = model.z(s, mu, ls)
        for v in VARIANTS:
            pij = torch.where(z > TL["u"][ig], tail_p(ig, z, v), torch.ones_like(z))
            pji = torch.where(z > TL["u"][jg], tail_p(jg, z, v), torch.ones_like(z))
            pm = torch.maximum(pij, pji)
            for k, th in enumerate(pthr):
                cnt[v]["ij"][k] += int((pij < th).sum()); cnt[v]["ji"][k] += int((pji < th).sum())
                cnt[v]["both"][k] += int((pm < th).sum())
        n_done += len(i)
        if (n_done // len(i)) % 50 == 0:
            log(f"extended {n_done:.3g} {time.time() - t:.0f}s")
    rows = []
    for v in VARIANTS:
        for k, th in enumerate(pthr):
            rows.append({"variant": v, "p_threshold": th, "E_threshold": th * N_DB, "n_pairs": n_done,
                         "expected_one_dir": n_done * th, "obs_ij": cnt[v]["ij"][k], "obs_ji": cnt[v]["ji"][k],
                         "obs_both": cnt[v]["both"][k], "ratio_ij": cnt[v]["ij"][k] / (n_done * th),
                         "ratio_both_vs_one_dir": cnt[v]["both"][k] / (n_done * th)})
    ext = pd.DataFrame(rows)
    ext.to_csv(ext_path, index=False)
    summary["extended"] = {"n_pairs": n_done, "drawn": n_drawn, "wall_s": time.time() - t}
    log("extended calibration\n" + ext.to_string())

summary["wall_s"] = time.time() - T0
prev = json.loads((out / "summary.json").read_text()) if (out / "summary.json").exists() else {}
prev.update(summary)
(out / "summary.json").write_text(json.dumps(prev, indent=1, default=float))
for k in ("snapshot_id", "prereg_sha256", "git_commit", "seed", "code_sha256_16", "written_utc"):
    prev.pop(k, None)              # provenance fields are rewritten by this run, never inherited
pairs.write_sidecar(out / "summary.json", a.snapshot, [0, 1, 2], {"stage": f"null_{a.score}", **prev})
(a.results / f"null_{a.score}_summary.json").write_text((out / "summary.json").read_text())
log("null done", a.score, f"{time.time() - T0:.0f}s")
