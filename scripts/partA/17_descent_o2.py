#!/usr/bin/env python
"""Part A v3: run the descent pair builder (scripts/partA_descent/build_descent_pairs.py, unmodified) on O2.

O2 SOGs contain tier-3 / tier-L sites that are not in the tier-1/2 fingerprint universe, so the builder's
M2 site matching (descent.matched_site_pairs) asserts. This wrapper substitutes a version that keeps only SOG
matches where both sites are in the universe (M2 then covers tier-1/2 sites only). M1 counts are not affected
(they do not use matched_site_pairs). All arguments are passed through to the builder.
"""
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))
from pocketevo import descent  # noqa: E402


def matched_site_pairs_universe_only(M, i, j):
    if len(i) == 0:
        return pd.DataFrame({"k": [], "sog": [], "si_a": [], "si_b": []})
    SI = M["SI"].tocsr()
    S1 = SI[i].tocsr(); S2 = SI[j].tocsr()
    mask = (S1 > 0).multiply(S2 > 0)

    def coo_sorted(X):
        X = X.tocsr(); X.eliminate_zeros(); X.sort_indices()
        c = X.tocoo()
        return c.row, c.col, c.data

    k, s, va = coo_sorted(S1.multiply(mask))
    kb, sb, vb = coo_sorted(S2.multiply(mask))
    assert np.array_equal(k, kb) and np.array_equal(s, sb)
    return pd.DataFrame({"k": k.astype(np.int64), "sog": s.astype(np.int64),
                         "si_a": (va - 1).astype(np.int64), "si_b": (vb - 1).astype(np.int64)})


descent.matched_site_pairs = matched_site_pairs_universe_only
spec = importlib.util.spec_from_file_location("bdp", REPO / "scripts/partA_descent/build_descent_pairs.py")
mod = importlib.util.module_from_spec(spec)
sys.modules["bdp"] = mod
sys.argv = [str(spec.origin)] + sys.argv[1:]
spec.loader.exec_module(mod)
mod.descent.matched_site_pairs = matched_site_pairs_universe_only
mod.main()
