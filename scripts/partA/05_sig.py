#!/usr/bin/env python
"""Part A stage 5: significant pairs from the k = 1000 lists (round-2 amendments 1, 4).

For every kNN entry (i -> j, different proteins): z_ij from the covariate model, E_ij = N_db p_i(z_ij) with the
query's tail, E_ji = N_db p_j(z_ij) with the neighbour's tail (the model is symmetric, so z_ji = z_ij); a pair is
kept if max(E_ij, E_ji) < --keep-E under any tail variant (default 1, for the E-threshold sensitivity), then
symmetrised (unordered, deduplicated). Saturation: queries whose k-th neighbour still has E_ij < 0.01 (primary).
Also writes a seeded query sample with its top-200 neighbour lists for the descriptive top-k composition.
Outputs: <work>/sig_<score>/pairs.parquet (i, j, s, z, Emax_<variant>), saturation.npy, topk_sample.npz
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
ap.add_argument("--keep-E", type=float, default=1.0)
ap.add_argument("--block", type=int, default=16384)
ap.add_argument("--topk-sample", type=int, default=100_000)
a = resolve(ap.parse_args())
log = logger(a.work / "partA.log")
dev = torch.device("cuda:0")
T0 = time.time()
nd = a.work / f"null_{a.score}"
out = a.work / f"sig_{a.score}"
out.mkdir(exist_ok=True)
model = nulls.CovModel.from_dict(json.loads((nd / "covmodel.json").read_text()))
tl = dict(np.load(nd / "tail.npz"))
xi_pool = float(tl["xi_pooled"])
TL = {k: torch.from_numpy(np.asarray(tl[k], np.float32)).to(dev) for k in ("u", "sc", "xi", "zeta", "ybar")}
VARIANTS = ("pooled", "perq", "xi0")
I = np.load(a.work / "knn" / "I.i32.npy", mmap_mode="r")
D = np.load(a.work / "knn" / "D.f16.npy", mmap_mode="r")
N, K = I.shape
N_DB = N
ctx = np.load(a.work / "context.npz")
site_prot = torch.from_numpy(ctx["site_prot"]).to(dev)
cov = np.load(a.work / "covariates.npy")
cov = np.where(np.isfinite(cov), cov, np.nanmean(cov, 0))
Fn = torch.from_numpy(model.norm(cov).astype(np.float32)).to(dev)
csls = None
if a.score == "csls":
    Dn = D
    r = np.empty(N, np.float32)
    for s in range(0, N, 200_000):
        d = Dn[s:s + 200_000, :11].astype(np.float32)
        selfcol = I[s:s + 200_000, 0] == np.arange(s, min(N, s + 200_000))
        r[s:s + 200_000] = np.where(selfcol, d[:, 1:11].mean(1), d[:, :10].mean(1))
    csls = torch.from_numpy(r).to(dev)

rng = np.random.default_rng(0)
topk_q = np.sort(rng.choice(N, min(a.topk_sample, N), replace=False))
topk_I = np.full((len(topk_q), 200), -1, np.int32)
topk_E = np.full((len(topk_q), 200), np.inf, np.float32)
topk_pos = {int(q): k for k, q in enumerate(topk_q)}

sat = np.zeros(N, np.float32)       # E_ij (pooled) of the last neighbour
n_h0 = 0
chunks = []
for s in range(0, N, a.block):
    e = min(N, s + a.block)
    q = torch.arange(s, e, device=dev)
    Ib = torch.from_numpy(np.ascontiguousarray(I[s:e]).astype(np.int64)).to(dev)
    Sb = torch.from_numpy(np.ascontiguousarray(D[s:e])).to(dev).float()
    if csls is not None:
        Sb = 2 * Sb - csls[q][:, None] - csls[Ib]
    fq = Fn[q][:, None, :].expand(-1, K, -1)
    fj = Fn[Ib]
    mu, ls = model.mu_logsd_pairs(fq.reshape(-1, fq.shape[-1]), fj.reshape(-1, fj.shape[-1]))
    Z = model.z(Sb.reshape(-1), mu, ls).reshape(Sb.shape)
    diffp = site_prot[q][:, None] != site_prot[Ib]
    n_h0 += int((~diffp).sum()) - int((Ib == q[:, None]).sum())
    Emax = {}
    for v in VARIANTS:
        fit_q = {k: TL[k][q][:, None] for k in TL}
        fit_j = {k: TL[k][Ib] for k in TL}
        scq, xiq = nulls.tail_variant(fit_q, v, xi_pool)
        scj, xij = nulls.tail_variant(fit_j, v, xi_pool)
        pij = torch.where(Z > fit_q["u"], nulls.gpd_sf(Z, fit_q["u"], scq, xiq, fit_q["zeta"]), torch.ones_like(Z))
        pji = torch.where(Z > fit_j["u"], nulls.gpd_sf(Z, fit_j["u"], scj, xij, fit_j["zeta"]), torch.ones_like(Z))
        Eij = N_DB * pij
        Emax[v] = torch.maximum(Eij, N_DB * pji)
        if v == "pooled":
            sat[s:e] = Eij[:, -1].cpu().numpy()
            Eij_pooled = Eij
    keep = diffp & (torch.stack([Emax[v] for v in VARIANTS]).min(0).values < a.keep_E)
    r_, c_ = torch.nonzero(keep, as_tuple=True)
    ii, jj = q[r_], Ib[r_, c_]
    lo, hi = torch.minimum(ii, jj), torch.maximum(ii, jj)
    rec = {"i": lo.cpu().numpy().astype(np.int32), "j": hi.cpu().numpy().astype(np.int32),
           "s": Sb[r_, c_].cpu().numpy().astype(np.float32), "z": Z[r_, c_].cpu().numpy().astype(np.float32)}
    for v in VARIANTS:
        rec[f"Emax_{v}"] = Emax[v][r_, c_].cpu().numpy().astype(np.float32)
    chunks.append(pd.DataFrame(rec))
    # top-200 sample (excluding self and same-protein entries), pooled one-direction E
    sel = topk_q[(topk_q >= s) & (topk_q < e)]
    if len(sel):
        li = torch.from_numpy(sel - s).to(dev)
        for row, qq in zip(li.tolist(), sel.tolist()):
            m = diffp[row]
            idx = Ib[row][m][:200]
            topk_I[topk_pos[qq], :len(idx)] = idx.cpu().numpy()
            topk_E[topk_pos[qq], :len(idx)] = Eij_pooled[row][m][:200].cpu().numpy()
    if (s // a.block) % 20 == 0:
        log(f"sig block {s}/{N} kept {sum(len(c) for c in chunks)} {time.time() - T0:.0f}s")

df = pd.concat(chunks, ignore_index=True)
n_entries = len(df)
df = df.sort_values(["i", "j"]).drop_duplicates(["i", "j"], keep="first").reset_index(drop=True)
df.to_parquet(out / "pairs.parquet", index=False)
np.save(out / "saturation.npy", sat)
np.savez(out / "topk_sample.npz", q=topk_q, I=topk_I, E=topk_E)
summ = {"score": a.score, "N": N, "k": K, "keep_E": a.keep_E, "n_entries_kept": n_entries, "n_pairs_unique": len(df),
        "n_same_protein_entries_excl_self": n_h0,
        "n_sig_pairs_E0.01": {v: int((df[f"Emax_{v}"] < 0.01).sum()) for v in VARIANTS},
        "saturated_queries_E0.01": int((sat < 0.01).sum()), "saturated_frac": float((sat < 0.01).mean()),
        "wall_s": time.time() - T0}
(out / "summary.json").write_text(json.dumps(summ, indent=1))
(a.results / f"sig_{a.score}_summary.json").write_text(json.dumps(summ, indent=1))
pairs.write_sidecar(out / "pairs.parquet", a.snapshot, 0, {"stage": f"sig_{a.score}", **summ})
log("sig done", json.dumps(summ))
