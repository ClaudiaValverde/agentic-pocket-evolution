#!/usr/bin/env python
"""Part A stage 8 (descriptive): raw top-k class composition (k = 10/50/200) of the fingerprint neighbour lists,
and the size-only / covariate-only baselines: neighbours ranked by the similarity *predicted from covariates alone*
(size-only = log volume model; covariate-only = full CovModel mean) for the same seeded query sample, against all
sites. Classes are join/forest classes (H0 removed; 'rest' = not resolved by joins, not searched: H3s/H4/N).
Output: results/partA/<snap>/topk_composition.csv
"""
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from _common import parser, resolve, logger, pairs
from pocketevo import homology as hom
from pocketevo import nulls

ap = parser(__doc__)
ap.add_argument("--n-baseline-queries", type=int, default=20_000)
a = resolve(ap.parse_args())
log = logger(a.work / "partA.log")
T0 = time.time()
dev = torch.device("cuda:0")
c = np.load(a.work / "context.npz")
forest = hom.parse_orthoxml(str(Path(a.orthology) / "hogs.orthoxml.gz"))
ctx = hom.PairContext(c["site_prot"], c["site_lab"], c["prot_group"], c["prot_gene"], forest)
tk = np.load(a.work / "sig_cos" / "topk_sample.npz")
Q, I, E = tk["q"], tk["I"], tk["E"]
rows = []


def composition(q, nb, label, Emat=None):
    for k in (10, 50, 200):
        qq = np.repeat(q, k)
        jj = nb[:, :k].ravel()
        ok = jj >= 0
        cls = hom.classify_pairs(ctx, qq[ok].astype(np.int64), jj[ok].astype(np.int64))
        names = pd.Series(cls).map(hom.CLASS_NAMES)
        vc = names.value_counts(normalize=True)
        r = {"ranking": label, "k": k, "n_queries": len(q), "n_pairs": int(ok.sum())}
        r.update({f"frac_{n}": float(vc.get(n, 0.0)) for n in ("H1", "H2", "H12", "H3j", "rest")})
        if Emat is not None:
            r["frac_E_lt_0.01_onedir"] = float((Emat[:, :k].ravel()[ok] < 0.01).mean())
        rows.append(r)


composition(Q, I, "fingerprint_cos", E)
# baselines on a sub-sample of the same queries
model = json.loads((a.work / "null_cos" / "covmodel.json").read_text())
full = nulls.CovModel.from_dict(model)
size = nulls.CovModel.from_dict(model["size_only"])
cov = np.load(a.work / "covariates.npy")
cov = np.where(np.isfinite(cov), cov, np.nanmean(cov, 0))
qb = Q[:a.n_baseline_queries]
site_prot = torch.from_numpy(ctx.site_prot).to(dev)
for name, m, cols in (("covariate_only", full, slice(None)), ("size_only", size, slice(0, 1))):
    Fn = torch.from_numpy(m.norm(cov[:, cols]).astype(np.float32)).to(dev)
    nb = np.full((len(qb), 200), -1, np.int32)
    for s in range(0, len(qb), 256):
        q = torch.from_numpy(qb[s:s + 256]).to(dev)
        best_v = torch.full((len(q), 400), -torch.inf, device=dev)
        best_i = torch.full((len(q), 400), -1, dtype=torch.long, device=dev)
        for t0 in range(0, len(cov), 1_000_000):
            tt = torch.arange(t0, min(len(cov), t0 + 1_000_000), device=dev)
            mu, _ = m.mu_logsd_block(Fn[q], Fn[tt])
            mu = torch.where(site_prot[q][:, None] == site_prot[tt][None, :], torch.full_like(mu, -torch.inf), mu)
            v, ix = torch.topk(mu, 400, dim=1)
            cat_v = torch.cat([best_v, v], 1); cat_i = torch.cat([best_i, tt[ix]], 1)
            sv, si = torch.topk(cat_v, 400, dim=1)
            best_v, best_i = sv, cat_i.gather(1, si)
        nb[s:s + len(q)] = best_i[:, :200].cpu().numpy()
    composition(qb, nb, name)
    log(f"baseline {name} done {time.time() - T0:.0f}s")
# random-pair reference composition
rng = np.random.default_rng(0)
rnd = rng.integers(0, len(cov), (len(qb), 200)).astype(np.int32)
composition(qb, rnd, "random")
df = pd.DataFrame(rows)
df.to_csv(a.results / "topk_composition.csv", index=False)
pairs.write_sidecar(a.results / "topk_composition.csv", a.snapshot, 0, {"stage": "topk_baselines",
                    "wall_s": time.time() - T0})
log("topk\n" + df.to_string())
