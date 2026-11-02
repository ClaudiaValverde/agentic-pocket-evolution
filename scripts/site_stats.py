#!/usr/bin/env python
"""Summary statistics of a site build (scripts/build_sites.py output) -> CSVs + SITES.md + checksums.

Reads only site metadata (never fingerprints) with DuckDB, plus the snapshot structure tables for
denominators (proteins with status ok). Also runs a threshold sensitivity on a deterministic unit sample
(re-merging with other Jaccard thresholds, in memory, no outputs) and records the within-model pairwise
residue-Jaccard distribution, which shows what the merge step can do at all.

    python scripts/site_stats.py --sites /tmp/pe_big/sites/20261030T085136Z \
        --manifest data/snapshot/20261030T085136Z/manifest.parquet --out results/sites/20261030T085136Z
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
import pyarrow.compute as pc
import pyarrow.parquet as pq

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from pocketevo import sites as S  # noqa: E402

JBINS = [0, 0.1, 0.2, 0.25, 0.3, 0.5, 1.0001]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while b := fh.read(1 << 20):
            h.update(b)
    return h.hexdigest()


def sens_unit(args):
    path, thresholds = args
    t = pq.read_table(path, columns=S.POCKET_COLUMNS)
    out = {"n_pockets": t.num_rows}
    for th in thresholds:
        out[f"sites_J{th}"] = S.build_sites(t, jaccard=th, with_fingerprints=False)[0].num_rows if t.num_rows else 0
    # within-model (same accession AND same fragment) pairwise Jaccard histogram
    t = t.sort_by([("accession", "ascending"), ("fragment", "ascending"), ("pocket_id", "ascending")])
    acc = pc.dictionary_encode(t["accession"].combine_chunks()).indices.to_numpy().astype(np.int64)
    frag = t["fragment"].to_numpy().astype(np.int64)
    row, res = S.parse_residues(t["residues"])
    key = (((acc[row] << 10) | frag[row]) << 20) | res      # model-local keys
    rk = np.unique(np.stack([row, key], 1), axis=0)
    i, j, jac = S.jaccard_edges(rk[:, 0], rk[:, 1], t.num_rows, 1e-9)
    out["pairs_hist"] = np.histogram(jac, JBINS)[0].tolist()
    out["max_within_model_J"] = float(jac.max()) if len(jac) else 0.0
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sites", required=True)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--sens-units", type=int, default=40, help="units per tier for the threshold sensitivity")
    ap.add_argument("--sens-thresholds", default="0.1,0.2,0.3,0.5")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--workers", type=int, default=16)
    a = ap.parse_args()
    O, out = Path(a.sites), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    mj = json.loads((Path(a.manifest).parent / "manifest.json").read_text())
    atlas = Path(mj["atlas_root"])
    man = pd.read_parquet(a.manifest)
    c = duckdb.connect()
    c.execute("SET threads TO 16")
    tiers = sorted(int(p.name[4:]) for p in (O / "sites").glob("tier*"))
    rows, hist_rows, spp_rows = [], [], []
    for T in tiers:
        files = sorted((O / "sites" / f"tier{T}").glob("*.parquet"))
        units = [int(f.stem) for f in files]
        mt = man[man.tier == T]
        sample = len(units) < len(mt)
        glob = str(O / "sites" / f"tier{T}" / "*.parquet")
        c.execute(f"CREATE OR REPLACE VIEW s AS SELECT * FROM read_parquet('{glob}')")
        g = c.sql("""SELECT count(*) n_sites, sum(n_pockets) n_pockets, count(DISTINCT accession) n_prot,
                     count(DISTINCT taxid) n_taxa,
                     avg((n_pockets>1)::INT) frac_multi_pocket, avg((n_fragments>1)::INT) frac_multi_fragment,
                     sum(CASE WHEN n_pockets>1 THEN n_pockets ELSE 0 END)/sum(n_pockets) frac_pockets_in_multi,
                     median(n_res) med_n_res, avg(n_res) mean_n_res, max(n_pockets) max_pockets_per_site
                     FROM s""").df().iloc[0].to_dict()
        # denominators: ok structures -> proteins (accessions) in these units
        spaths = [str(atlas / p) for p in mt[mt.unit_id.isin(units)].structures_path]
        den = c.sql(f"""SELECT count(DISTINCT accession) n_prot_ok,
                        count(DISTINCT CASE WHEN fragment>1 THEN accession END) n_prot_fragmented
                        FROM read_parquet({spaths!r}, union_by_name=true) WHERE status='ok'""").df().iloc[0]
        spp = c.sql("SELECT accession, count(*) k, sum(n_pockets) np, max(n_fragments) nf FROM s GROUP BY 1").df()
        k_all = np.r_[spp.k.to_numpy(), np.zeros(int(den.n_prot_ok) - len(spp), int)]
        frag = spp[spp.nf > 1]
        g.update(tier=T, n_units=len(units), units_in_snapshot=len(mt), is_sample=sample,
                 n_prot_ok=int(den.n_prot_ok), frac_prot_with_site=len(spp) / max(int(den.n_prot_ok), 1),
                 sites_per_prot_mean=float(spp.k.mean()), sites_per_prot_median=float(spp.k.median()),
                 sites_per_prot_p90=float(spp.k.quantile(.9)), sites_per_prot_p99=float(spp.k.quantile(.99)),
                 sites_per_prot_max=int(spp.k.max()), sites_per_prot_mean_incl0=float(k_all.mean()),
                 n_prot_fragmented_with_sites=len(frag),
                 fragmented_pockets=int(frag.np.sum()), fragmented_sites=int(frag.k.sum()))
        rows.append(g)
        h = c.sql("SELECT least(n_pockets, 6) b, count(*) n FROM s GROUP BY 1 ORDER BY 1").df()
        for b, n in zip(h.b, h.n):
            hist_rows.append(dict(tier=T, pockets_per_site=("6+" if b == 6 else str(int(b))), n_sites=int(n)))
        h2 = pd.Series(np.minimum(k_all, 21)).value_counts().sort_index()
        for b, n in h2.items():
            spp_rows.append(dict(tier=T, sites_per_protein=("21+" if b == 21 else str(int(b))), n_prot=int(n)))
    summ = pd.DataFrame(rows)
    summ.to_csv(out / "site_summary_by_tier.csv", index=False)
    pd.DataFrame(hist_rows).to_csv(out / "pockets_per_site_hist.csv", index=False)
    pd.DataFrame(spp_rows).to_csv(out / "sites_per_protein_hist.csv", index=False)

    # threshold sensitivity + within-model Jaccard on a deterministic sample of units per tier
    th = [float(x) for x in a.sens_thresholds.split(",")]
    sens = []
    rng = np.random.default_rng(a.seed)
    for T in tiers:
        mt = man[man.tier == T].sort_values("unit_id")
        pick = mt.iloc[np.sort(rng.choice(len(mt), min(a.sens_units, len(mt)), replace=False))]
        with ProcessPoolExecutor(a.workers) as ex:
            res = list(ex.map(sens_unit, [(str(atlas / p), th) for p in pick.pockets_path]))
        d = dict(tier=T, n_units=len(pick), n_pockets=sum(r["n_pockets"] for r in res))
        for x in th:
            d[f"sites_J{x}"] = sum(r[f"sites_J{x}"] for r in res)
            d[f"pockets_per_site_J{x}"] = d["n_pockets"] / max(d[f"sites_J{x}"], 1)
        hist = np.sum([r["pairs_hist"] for r in res], axis=0)
        for lo, hi, n in zip(JBINS[:-1], JBINS[1:], hist):
            d[f"pairs_J[{lo},{min(hi, 1)})"] = int(n)
        d["max_within_model_J"] = max(r["max_within_model_J"] for r in res)
        sens.append(d)
    sens = pd.DataFrame(sens)
    sens.to_csv(out / "threshold_sensitivity.csv", index=False)

    # checksums of all outputs
    with open(out / "checksums.sha256", "w") as fh:
        for f in sorted(O.rglob("*.parquet")):
            fh.write(f"{sha256(f)}  {f.relative_to(O)}\n")
    print(summ.T.to_string())
    print(sens.T.to_string())


if __name__ == "__main__":
    main()
