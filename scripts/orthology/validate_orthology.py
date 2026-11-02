#!/usr/bin/env python
"""Validation of the O2 HOG tables (sanity checks with known answers + PANTHER agreement).

1. Known pairs (human HSP90AA1/HSP90AB1 paralogs, human-mouse/yeast orthologs, EGFR/ERBB2, globins ...).
2. Human-mouse: #ortholog pairs, fraction of human genes in HOGs with a mouse ortholog, 1:1 fraction.
3. PANTHER agreement (UniProt xref_panther, current release) on pairs sampled with a fixed seed from
   species with PANTHER annotations fetched per proteome: n_pairs split into ortholog / paralog / cross-root
   ('unrelated': proteins of the same species pair in different root HOGs). Reports for each relation the
   fraction of annotated pairs in the same PANTHER family (PTHRnnnnn) and same subfamily (PTHRnnnnn:SFn).
   Also the reverse view (recall): of human-mouse protein pairs sharing a PANTHER subfamily, the fraction
   that are O2 orthologs (same root HOG / ortholog relation).
4. Agreement of pocketevo.orthology.relation with Ancestra's own hogs_to_pairs on a sample (if importable).
Output: <out_dir>/validation/validation.json, panther_pairs.parquet
Usage: PYTHONPATH=src envs/ancestra/bin/python scripts/orthology/validate_orthology.py   (ancestra env for step 4;
       .venv works and skips step 4)
"""
from __future__ import annotations

import argparse
import io
import json
import os
import random
import sys
import time
import urllib.parse
import urllib.request

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import fetch, load_config, log  # noqa: E402
from pocketevo.orthology import HogIndex  # noqa: E402

KNOWN = [  # (a, b, expected relation, note)
    ("P07900", "P08238", "paralog", "HSP90AA1 vs HSP90AB1 (human)"),
    ("P07900", "P07901", "ortholog", "HSP90AA1 human vs mouse"),
    ("P08238", "P11499", "ortholog", "HSP90AB1 human vs mouse"),
    ("P07900", "P02829", "ortholog", "HSP90AA1 human vs yeast HSP82"),
    ("P00533", "P04626", "paralog", "EGFR vs ERBB2 (human)"),
    ("P00533", "Q01279", "ortholog", "EGFR human vs mouse"),
    ("P68871", "P02088", "ortholog", "HBB human vs mouse Hbb-b1"),
    ("P69905", "P68871", "paralog", "HBA vs HBB (human)"),
    ("P04637", "P02340", "ortholog", "TP53 human vs mouse"),
    ("P60709", "P63261", "paralog", "ACTB vs ACTG1 (human)"),
    ("P38398", "P48754", "ortholog", "BRCA1 human vs mouse"),
    ("P00533", "P07900", "unrelated", "EGFR vs HSP90AA1 (negative control)"),
    ("P68871", "P07900", "unrelated", "HBB vs HSP90AA1 (negative control)"),
]
PANTHER_SPECIES = ["UP000005640_9606", "UP000000589_10090", "UP000000437_7955", "UP000000803_7227",
                   "UP000001940_6239", "UP000002311_559292", "UP000006548_3702", "UP000000625_83333"]


def panther_for_proteome(up: str, dest: str) -> pd.DataFrame:
    url = "https://rest.uniprot.org/uniprotkb/stream?" + urllib.parse.urlencode(
        {"query": f"proteome:{up}", "fields": "accession,xref_panther", "format": "tsv"})
    fetch(url, dest)
    d = pd.read_csv(dest, sep="\t", dtype=str).rename(columns={"Entry": "protein", "PANTHER": "panther"})
    return d


def fam_sub(x):
    if not isinstance(x, str) or not x:
        return None, None
    ids = [i for i in x.strip(";").split(";") if i]
    fam = next((i for i in ids if ":" not in i), None) or ids[0].split(":")[0]
    sub = next((i for i in ids if ":" in i), None)
    return fam, sub


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config")
    ap.add_argument("--tables")
    ap.add_argument("--n-pairs", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    cfg = load_config(a.config)
    tdir = a.tables or os.path.join(cfg["out_dir"], "tables")
    vdir = os.path.join(os.path.dirname(tdir), "validation")
    os.makedirs(vdir, exist_ok=True)
    t0 = time.time()
    h = pd.read_parquet(os.path.join(tdir, "hogs.parquet"))
    p = pd.read_parquet(os.path.join(tdir, "protein_hogs.parquet"))
    idx = HogIndex(h, p)
    res = {"tables": tdir}

    # 1. known pairs
    kn = []
    for x, y, exp, note in KNOWN:
        r = idx.relation(x, y)
        kn.append(dict(a=x, b=y, note=note, expected=exp, got=r["relation"], level=r["level"], ok=r["relation"] == exp))
    res["known_pairs"] = kn
    res["known_pairs_ok"] = f"{sum(k['ok'] for k in kn)}/{len(kn)}"
    log(pd.DataFrame(kn).to_string())

    # 2. human-mouse orthologs
    hs, mm = "UP000005640_9606", "UP000000589_10090"
    sub = p[p.species.isin([hs, mm]) & (p.root >= 0)]
    pairs = []
    for r, g in sub.groupby("root"):
        H = g[g.species == hs].protein.tolist()
        M = g[g.species == mm].protein.tolist()
        for x in H:
            for y in M:
                if idx.relation(x, y)["relation"] == "ortholog":
                    pairs.append((x, y))
    hm = pd.DataFrame(pairs, columns=["h", "m"])
    nh = hm.groupby("h").m.nunique()
    nm = hm.groupby("m").h.nunique()
    one = hm[hm.h.map(nh).eq(1) & hm.m.map(nm).eq(1)]
    n_h = int((p.species == hs).sum())
    res["human_mouse"] = dict(n_ortholog_pairs=len(hm), n_1to1=len(one), n_human_proteins=n_h,
                              frac_human_with_mouse_ortholog=round(hm.h.nunique() / n_h, 4),
                              frac_human_1to1=round(len(one) / n_h, 4),
                              frac_pairs_1to1=round(len(one) / max(len(hm), 1), 4))
    log(res["human_mouse"])

    # 3. PANTHER
    pan = []
    for sp in PANTHER_SPECIES:
        up = sp.split("_")[0]
        d = panther_for_proteome(up, os.path.join(cfg["work_dir"], "panther", f"{up}.tsv"))
        d["species"] = sp
        pan.append(d)
    pan = pd.concat(pan)
    pan[["fam", "sub"]] = pan.panther.apply(lambda x: pd.Series(fam_sub(x)))
    fam = dict(zip(pan.protein, pan.fam))
    sfam = dict(zip(pan.protein, pan["sub"]))
    res["panther_coverage"] = {sp: round(float(pan[pan.species == sp].fam.notna().mean()), 3) for sp in PANTHER_SPECIES}
    rng = random.Random(a.seed)
    pp = p[p.species.isin(PANTHER_SPECIES) & (p.root >= 0) & p.protein.map(lambda x: fam.get(x) is not None)]
    by_root = {r: g.protein.tolist() for r, g in pp.groupby("root") if len(g) >= 2}
    roots = sorted(by_root)
    want = {"ortholog": int(a.n_pairs * 0.4), "paralog": int(a.n_pairs * 0.3)}
    want["unrelated"] = a.n_pairs - sum(want.values())
    got = {k: [] for k in want}
    tries = 0
    while any(len(got[k]) < want[k] for k in ("ortholog", "paralog")) and tries < 2_000_000:
        tries += 1
        r = rng.choice(roots)
        x, y = rng.sample(by_root[r], 2)
        rel = idx.relation(x, y)
        k = rel["relation"]
        if k in ("ortholog", "paralog") and len(got[k]) < want[k]:
            got[k].append((x, y, k, rel["level"]))
    allp = pp.protein.tolist()
    root_of = dict(zip(pp.protein, pp.root))
    while len(got["unrelated"]) < want["unrelated"]:
        x, y = rng.sample(allp, 2)
        if root_of[x] != root_of[y]:
            got["unrelated"].append((x, y, "unrelated", None))
    df = pd.DataFrame([t for k in got for t in got[k]], columns=["a", "b", "relation", "level"])
    df["fam_a"], df["fam_b"] = df.a.map(fam), df.b.map(fam)
    df["same_family"] = df.fam_a == df.fam_b
    df["sub_a"], df["sub_b"] = df.a.map(sfam), df.b.map(sfam)
    both_sub = df.sub_a.notna() & df.sub_b.notna()
    df["same_subfamily"] = (df.sub_a == df.sub_b).where(both_sub)
    df.to_parquet(os.path.join(vdir, "panther_pairs.parquet"), index=False)
    agg = df.groupby("relation").agg(n=("a", "size"), same_family=("same_family", "mean"),
                                     same_subfamily=("same_subfamily", "mean"),
                                     n_with_sub=("same_subfamily", "count")).round(4)
    res["panther_pairs"] = agg.reset_index().to_dict(orient="records")
    res["panther_pairs_tries"] = tries
    log(agg.to_string())
    # recall view: human-mouse pairs sharing a PANTHER subfamily
    hp = pan[(pan.species == hs) & pan["sub"].notna()][["protein", "sub"]]
    mp = pan[(pan.species == mm) & pan["sub"].notna()][["protein", "sub"]]
    hp = hp[hp.protein.isin(idx.path)]
    mp = mp[mp.protein.isin(idx.path)]
    j = hp.merge(mp, on="sub", suffixes=("_h", "_m"))
    j = j.sample(n=min(20000, len(j)), random_state=a.seed)
    rels = [idx.relation(x, y)["relation"] for x, y in zip(j.protein_h, j.protein_m)]
    vc = pd.Series(rels).value_counts(normalize=True).round(4).to_dict()
    res["human_mouse_same_panther_subfamily_pairs"] = dict(n=len(j), relation_fractions=vc)
    log(res["human_mouse_same_panther_subfamily_pairs"])

    # 4. cross-check with ancestra's flattening
    try:
        from ancestra.orthoxml import hogs_to_pairs, read_orthoxml
        ox_path = os.path.join(os.path.dirname(tdir), "accurate", "hogs.orthoxml.gz")
        if os.path.exists(ox_path):
            ox = read_orthoxml(ox_path)
            orth, para = hogs_to_pairs(ox, with_paralogs=True)
            rs = random.Random(a.seed)
            so = rs.sample(sorted(orth), 20000)
            sp_ = rs.sample(sorted(para), 20000)
            ok_o = sum(idx.relation(x, y)["relation"] == "ortholog" for x, y in so)
            ok_p = sum(idx.relation(x, y)["relation"] == "paralog" for x, y in sp_)
            res["ancestra_crosscheck"] = dict(n_orth_pairs_total=len(orth), n_para_pairs_total=len(para),
                                              orth_agree=f"{ok_o}/20000", para_agree=f"{ok_p}/20000")
            log(res["ancestra_crosscheck"])
    except ImportError:
        res["ancestra_crosscheck"] = "skipped (ancestra not importable)"
    res["seconds"] = round(time.time() - t0, 1)
    json.dump(res, open(os.path.join(vdir, "validation.json"), "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
