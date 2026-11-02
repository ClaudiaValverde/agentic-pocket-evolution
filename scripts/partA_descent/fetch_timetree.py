#!/usr/bin/env python
"""Fetch TimeTree pairwise divergence times for all tier-1 species pairs (DESIGN §2; fetched once, DATA.md).

Writes data/timetree/<snapshot>_<orthology>_pairwise.csv with the raw TimeTree response per pair (empty = TimeTree
had no estimate), plus a sidecar. Endpoint: https://timetree.org/api/pairwise/<taxid_a>/<taxid_b> (CSV).
"""
import argparse
import csv
import io
import itertools
import json
import sys
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from pocketevo.pairs import write_sidecar  # noqa: E402

URL = "https://timetree.org/api/pairwise/{a}/{b}"


def fetch(a, b, retries=3):
    for k in range(retries):
        try:
            with urllib.request.urlopen(URL.format(a=a, b=b), timeout=30) as r:
                txt = r.read().decode()
            rows = list(csv.DictReader(io.StringIO(txt)))
            return rows[0] if rows else {}
        except Exception as e:  # network: retry, then record failure (counted, never silent)
            err = repr(e)
            time.sleep(1 + k)
    return {"error": err}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", default="results/sogs/20261030T085136Z/O1/web/species.parquet")
    ap.add_argument("--snapshot", default="20261030T085136Z")
    ap.add_argument("--orthology-version", default="O1")
    ap.add_argument("--out", default="data/timetree")
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args()
    sp = pd.read_parquet(a.species)
    tax = sorted(int(t) for t in sp.taxid)
    pairs = list(itertools.combinations(tax, 2))
    t0 = time.time()
    with ThreadPoolExecutor(a.threads) as ex:
        res = list(ex.map(lambda p: fetch(*p), pairs))
    rows = []
    for (x, y), r in zip(pairs, res):
        age = r.get("precomputed_age") or ""
        rows.append({"taxid_a": x, "taxid_b": y, "tt_taxon_a": r.get("taxon_a_id", ""), "tt_taxon_b": r.get("taxon_b_id", ""),
                     "n_studies": r.get("all_total", ""), "age_mya": age, "ci_low": r.get("precomputed_ci_low", ""),
                     "ci_high": r.get("precomputed_ci_high", ""), "adjusted_age": r.get("adjusted_age", ""),
                     "error": r.get("error", "")})
    out = Path(a.out) / f"{a.snapshot}_{a.orthology_version}_pairwise.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows)
    df.to_csv(out, index=False)
    n_ok = int((pd.to_numeric(df.age_mya, errors="coerce") > 0).sum())
    write_sidecar(out, a.snapshot, None, {"source": URL, "retrieved_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                                          "n_pairs": len(df), "n_with_age": n_ok, "runtime_s": time.time() - t0})
    print(out, len(df), "with age:", n_ok, f"{time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
