#!/usr/bin/env python
"""Seeded, stratified sample of tier-1 proteins for the re-detection experiment (prereg
states.measured_error_rates.redetection_sample_n_proteins = 2000) and a seeds subset (300).

Eligible: tier-1 structure rows of the snapshot manifest with status ok, single-fragment accession (F1 only),
>= 1 kept atlas pocket, and a Foldcomp afdb_uniprot_v4 entry for the accession (v4 model exists).
Strata: kingdom (Eukaryota / Bacteria / Archaea, from NCBI rankedlineage) x structure pLDDT band (<70, 70-90,
>=90): equal allocation n/9 per cell, cells with too few proteins take all and the remainder is spread
proportionally over the other cells. Seed = prereg meta.seeds[0]. The 300-protein seeds subset is a stratified
sub-sample (same rule) of the 2,000.

  python scripts/observation/sample_redetection.py --out /tmp/pe_big/obs/redetect_sample.parquet
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

P = Path(__file__).resolve().parents[2]
ATLAS = Path("/cluster/scratch/gmagnusson/pocketcompare_afdb")


def kingdoms(taxids) -> dict:
    merged = {}
    with open(ATLAS / "ncbi/merged.dmp") as fh:
        for ln in fh:
            f = [x.strip() for x in ln.split("|")]
            merged[int(f[0])] = int(f[1])
    tx = {int(t): merged.get(int(t), int(t)) for t in taxids}
    want = set(tx.values())
    out = {}
    with open(ATLAS / "ncbi/rankedlineage.dmp") as fh:
        for ln in fh:
            f = [x.strip() for x in ln.split("|")]
            t = int(f[0])
            if t in want:
                out[t] = f[9] or f[8]          # domain/superkingdom column (new_taxdump: col 10 = superkingdom)
    return {t: out.get(m, "unknown") for t, m in tx.items()}


def allocate(sizes: pd.Series, n: int) -> pd.Series:
    alloc = pd.Series(0, index=sizes.index)
    left, open_ = n, list(sizes.index)
    while left > 0 and open_:
        share = left / len(open_)
        changed = False
        for c in list(open_):
            room = sizes[c] - alloc[c]
            if room <= share:
                alloc[c] += room
                left -= room
                open_.remove(c)
                changed = True
        if not changed:
            base = int(left // len(open_))
            for c in open_:
                alloc[c] += base
            left -= base * len(open_)
            for c in open_[:left]:
                alloc[c] += 1
            left = 0
    return alloc.astype(int)


def strat_sample(df: pd.DataFrame, n: int, rng: np.random.Generator) -> pd.DataFrame:
    sizes = df.groupby("stratum").size()
    alloc = allocate(sizes, n)
    parts = []
    for c, k in alloc.items():
        sub = df[df.stratum == c]
        parts.append(sub.iloc[rng.permutation(len(sub))[:k]])
    return pd.concat(parts).sort_values("accession").reset_index(drop=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", default=str(P / "data/snapshot/20261030T085136Z/manifest.parquet"))
    ap.add_argument("--fc", default="/tmp/pe_big/obs/fc_t12_entries.parquet")
    ap.add_argument("--n", type=int, default=None)
    ap.add_argument("--n-seeds", type=int, default=300)
    ap.add_argument("--out", default="/tmp/pe_big/obs/redetect_sample.parquet")
    a = ap.parse_args()
    pre = yaml.safe_load(open(P / "config/prereg.yaml"))
    seed = int(pre["meta"]["seeds"][0])
    n = a.n or int(pre["states"]["measured_error_rates"]["redetection_sample_n_proteins"]["value"])
    man = pd.read_parquet(a.manifest)
    man = man[man.tier == 1]
    S = pd.concat([pd.read_parquet(ATLAS / p, columns=["unit_id", "structure_stem", "accession", "fragment", "taxid",
                                                       "organism", "source", "n_res", "plddt_mean", "status",
                                                       "n_pockets_all", "n_pockets_kept"])
                   for p in man.structures_path], ignore_index=True)
    nfrag = S.groupby("accession").fragment.max()
    S = S[(S.status == "ok") & (S.fragment == 1) & S.accession.map(nfrag).eq(1) & (S.n_pockets_kept >= 1)]
    fc = pd.read_parquet(a.fc)
    fc = fc[fc.name.str.endswith("-F1-model_v4")].drop_duplicates("acc")
    S = S.merge(fc[["acc", "name", "off", "len"]].rename(columns={"acc": "accession", "name": "fc_name"}),
                on="accession", how="inner")
    kd = kingdoms(S.taxid.unique())
    S["kingdom"] = S.taxid.map(kd).fillna("unknown")
    S["band"] = pd.cut(S.plddt_mean, [-1, 70, 90, 101], right=False, labels=["<70", "70-90", ">=90"]).astype(str)
    S["stratum"] = S.kingdom + "|" + S.band
    rng = np.random.default_rng(seed)
    samp = strat_sample(S, n, rng)
    rng2 = np.random.default_rng(seed + 1)
    seeds = strat_sample(samp, a.n_seeds, rng2)
    samp["in_seed_subset"] = samp.accession.isin(seeds.accession)
    samp["member"] = samp.structure_stem + ".cif.gz"
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    samp.to_parquet(a.out, index=False)
    info = dict(seed=seed, n=int(len(samp)), n_seed_subset=int(samp.in_seed_subset.sum()), n_eligible=int(len(S)),
                eligible_by_stratum=S.groupby("stratum").size().to_dict(),
                sample_by_stratum=samp.groupby("stratum").size().to_dict())
    Path(a.out).with_suffix(".json").write_text(json.dumps(info, indent=1))
    print(json.dumps(info, indent=1))


if __name__ == "__main__":
    main()
