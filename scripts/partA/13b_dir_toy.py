#!/usr/bin/env python
"""Part A v2.1 stage 13b: directional-test toy worlds over a fixed seed list (verifier item 7 / PI v2.1 item 3).

No-convergence BM worlds (src/pocketevo/directional.py::toy_world, 20 families x 32 species, the production toy
parameters) for seeds 0..49 (prereg seeds [0, 1, 2] are the first three): false-pass rate = passes / tested pairs,
with a Wilson CI over the pooled tests. Power for planted derived-clade convergence at frac 0.5 / 0.2 / 0.1 on seeds
0..9. The test and its cutoff are exactly those used on the real data (directional.decide).
Output: results/partA/<snap>/<run>/directional_toy_seeds.csv
"""
import time

import numpy as np
import pandas as pd

from _common import parser, resolve, logger, pairs
from pocketevo import directional as dr
from pocketevo.estimands import wilson

ap = parser(__doc__)
ap.add_argument("--n-seeds", type=int, default=50)
a = resolve(ap.parse_args())
log = logger(a.work / "partA.log")
T0 = time.time()
from multiprocessing import Pool


def job(arg):
    world, s, frac = arg
    r = dr.toy_world(s, planted=world == "planted", frac=frac if frac else 0.5)
    return {"world": world, "seed": s, "frac": frac, **{k: float(v) for k, v in r.items()}}


jobs = [("no_convergence", s, 0.0) for s in range(a.n_seeds)] + \
       [("planted", s, f) for f in (0.5, 0.2, 0.1) for s in range(10)]
with Pool(16) as pool:
    rows = pool.map(job, jobs)
df = pd.DataFrame(rows)
nc = df[df.world == "no_convergence"]
k = float((nc.false_pass * nc.n).sum()); n = float(nc.n.sum())
lo, hi = wilson(k, n)
summary = [{"world": "no_convergence_pooled", "seed": -1, "frac": 0.0, "n": n, "false_pass": k / n,
            "false_pass_wilson_lo": lo, "false_pass_wilson_hi": hi, "n_seeds_with_any_false_pass":
            int((nc.false_pass > 0).sum()), "n_seeds": len(nc),
            "prereg_seeds_0_1_2_false_pass": float((nc[nc.seed < 3].false_pass * nc[nc.seed < 3].n).sum()
                                                   / nc[nc.seed < 3].n.sum())}]
for frac in (0.5, 0.2, 0.1):
    p = df[(df.world == "planted") & (df.frac == frac)]
    summary.append({"world": "planted_pooled", "seed": -1, "frac": frac, "power": float(p.power.mean()),
                    "false_pass": float((p.false_pass * p.n).sum() / p.n.sum())})
out = pd.concat([pd.DataFrame(summary), df], ignore_index=True)
out.to_csv(a.results / "directional_toy_seeds.csv", index=False)
pairs.write_sidecar(a.results / "directional_toy_seeds.csv", a.snapshot, list(range(a.n_seeds)),
                    {"stage": "directional_toy_seeds", "wall_s": time.time() - T0})
log("toy seeds\n" + pd.DataFrame(summary).to_string())
