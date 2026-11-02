#!/usr/bin/env python
"""Part A v2 stage 14: undetected-homology check on significant N pairs (prereg undetected_homology_correction:
200 pairs, 'HHsearch-level check or Foldseek TM-align mode'; HHsearch is not installed -> Foldseek TM-align mode).

200 seeded significant N pairs (pooled, E < 0.01) and, as reference, 200 seeded N pairs of the random 'rest' sample
are re-aligned with `foldseek structurealign --alignment-type 1` (TM-align) on the full chains. Hidden homolog
(lenient, conservative for N): max(qTM, tTM) >= 0.5; strict: min(qTM, tTM) >= 0.5. epsilon = hidden fraction (Wilson
CI); corrected F_analog = F_analog x (1 - epsilon) with the difference moved to H4.
"""
import json

import numpy as np
import pandas as pd

from _common import parser, resolve, logger, pairs
from pocketevo import homology as hom
from pocketevo.estimands import wilson

ap = parser(__doc__)
ap.add_argument("--n", type=int, default=int(200))
ap.add_argument("--foldseek", default=str(pairs.REPO / "envs/ancestra/bin/foldseek"))
ap.add_argument("--foldseek-db", default="/tmp/pe_big/partA/foldseek/afdb_v6")
a = resolve(ap.parse_args())
log = logger(a.work / "partA.log")
n_pre = int(a.prereg["pairs"]["homology_classes"]["undetected_homology_correction"]["n_pairs_sampled"]["value"])
n = a.n or n_pre
uni = pd.read_parquet(a.work / "universe.parquet", columns=["accession"])
ACC = uni.accession.to_numpy()
sig = pd.read_parquet(a.work / f"cls_cos{a.tag}" / "pairs_cls.parquet", columns=["i", "j", "cls", "Emax_pooled"])
sigN = sig[(sig.cls == "N") & (sig.Emax_pooled < 0.01)].sample(n, random_state=a.seed + 50)
smp = pd.read_parquet(a.work / f"cls_cos{a.tag}" / "sample_rest.parquet")
rndN = smp[smp.cls == "N"].sample(n, random_state=a.seed + 51)
lk = pd.read_csv(a.foldseek_db + ".lookup", sep="\t", header=None, names=["key", "name", "file"])
name = {x.replace("AF-", "").split("-F1")[0]: x for x in lk.name}
n2k = dict(zip(lk.name, lk.key))
rows = []
for tag, d in (("sig_N", sigN), ("random_N", rndN)):
    q = pd.DataFrame({"qname": pd.Series(ACC[d.i.values]).map(name).fillna("__m__").values,
                      "tname": pd.Series(ACC[d.j.values]).map(name).fillna("__m__").values})
    r = hom.run_pair_search(q, a.foldseek_db, n2k, a.work / "foldseek_aln" / f"tmalign_{tag}", "foldseek", a.foldseek,
                            threads=16, evalue=1e6, extra=("--alignment-type", "1", "--tmscore-threshold", "0", "--exact-tmscore", "1"),
                            log=log)
    tmax = np.fmax(r.qtmscore.values, r.ttmscore.values)
    tmin = np.fmin(r.qtmscore.values, r.ttmscore.values)
    k_l, k_s = int(np.nansum(tmax >= 0.5)), int(np.nansum(tmin >= 0.5))
    rows.append({"set": tag, "n": n, "aligned": int(np.isfinite(tmax).sum()), "hidden_lenient": k_l,
                 "eps_lenient": k_l / n, "eps_lenient_wilson": wilson(k_l, n), "hidden_strict": k_s,
                 "eps_strict": k_s / n, "eps_strict_wilson": wilson(k_s, n),
                 "tm_max_median": float(np.nanmedian(tmax))})
    d = d.assign(qtm=r.qtmscore.values, ttm=r.ttmscore.values)
    d.to_parquet(a.work / f"hidden_homology_{tag}.parquet", index=False)
df = pd.DataFrame(rows)
df.to_csv(a.results / "hidden_homology.csv", index=False)
pairs.write_sidecar(a.results / "hidden_homology.csv", a.snapshot, [a.seed + 50, a.seed + 51],
                    {"stage": "hidden_homology", "method": "foldseek structurealign --alignment-type 1 (TM-align)"})
log("hidden homology\n" + df.to_string())
