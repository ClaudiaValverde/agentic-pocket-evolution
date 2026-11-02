#!/usr/bin/env python
"""Part A figures: null calibration (required): (a) QQ of held-out one-direction p-values (10^6 pairs, seed 1) per
tail variant; (b) observed / expected counts vs p threshold on held-out (10^6) and extended (10^9, seed 2) random
non-homologous pairs. Each PNG gets a sidecar (generating script + input hashes)."""
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from _common import parser, resolve, pairs

ap = parser(__doc__)
ap.add_argument("--score", default="cos")
a = resolve(ap.parse_args())
COL = {"pooled": "#2a78d6", "perq": "#eb6834", "xi0": "#1baf7a"}
LAB = {"pooled": "pooled shape (primary)", "perq": "per-query shape (prereg literal)", "xi0": "shape clipped ≥ 0"}
INK, INK2 = "#0b0b0b", "#52514e"
cal_p = a.work / f"null_{a.score}" / "calib_heldout.parquet"
cal = pd.read_parquet(cal_p)
ho = pd.read_csv(a.results / f"null_{a.score}_calibration_heldout.csv")
ext_p = a.results / f"null_{a.score}_calibration_extended.csv"
ext = pd.read_csv(ext_p)
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False})
fig, ax = plt.subplots(1, 2, figsize=(10, 4.2), facecolor="#fcfcfb")
n = len(cal)
exp_ = -np.log10((np.arange(1, n + 1) - 0.5) / n)
for v in ("pooled", "perq", "xi0"):
    obs = -np.log10(np.sort(np.clip(cal[f"p_{v}_ij"].to_numpy(), 1e-30, 1)))
    ax[0].plot(exp_, np.minimum(obs, 12), color=COL[v], lw=2, label=LAB[v])
ax[0].plot([0, 6.2], [0, 6.2], color=INK2, lw=1, ls="--", label="y = x")
ax[0].set_xlabel("expected −log10 p (uniform)")
ax[0].set_ylabel("observed −log10 p (capped at 12)")
ax[0].set_title(f"a  held-out QQ, {n:,} random non-homologous pairs", loc="left", color=INK)
ax[0].legend(frameon=False, fontsize=8)
for v in ("pooled", "perq", "xi0"):
    d = ho[ho.variant == v]
    ax[1].plot(d.p_threshold, np.maximum(d.ratio_ij, 1e-3), color=COL[v], lw=2, marker="o", ms=5,
               label=f"{LAB[v]}, held-out 10^6")
    e = ext[ext.variant == v]
    ax[1].plot(e.p_threshold, np.maximum(e.ratio_ij, 1e-3), color=COL[v], lw=1.2, ls=":", marker="s", ms=4,
               label=f"{LAB[v]}, extended {int(e.n_pairs.iloc[0]):.0e}")
ax[1].axhline(1, color=INK2, lw=1, ls="--")
ax[1].axvline(0.01 / ext.E_threshold.iloc[0] * ext.p_threshold.iloc[0], color=INK2, lw=0.8)
ax[1].text(0.01 / ext.E_threshold.iloc[0] * ext.p_threshold.iloc[0] * 1.5, 2e-3, "E = 0.01", color=INK2, fontsize=8)
ax[1].set_xscale("log"); ax[1].set_yscale("log"); ax[1].invert_xaxis()
ax[1].set_xlabel("p threshold (one direction)")
ax[1].set_ylabel("observed / expected count")
ax[1].set_title("b  observed vs expected false positives", loc="left", color=INK)
ax[1].legend(frameon=False, fontsize=6.5)
fig.tight_layout()
out = a.results / f"fig_null_calibration_{a.score}.png"
fig.savefig(out, dpi=150, facecolor="#fcfcfb")
h = {str(p): hashlib.sha256(Path(p).read_bytes()).hexdigest()[:16] for p in (cal_p, ext_p)}
pairs.write_sidecar(out, a.snapshot, [1, 2], {"script": "scripts/partA/10_figures.py", "inputs_sha256_16": h})
print(out)
