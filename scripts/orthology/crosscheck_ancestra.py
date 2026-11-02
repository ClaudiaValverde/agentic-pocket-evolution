#!/usr/bin/env python
"""Agreement of pocketevo.orthology.HogIndex.relation with Ancestra's own hogs_to_pairs flattening on random
samples of ortholog and paralog pairs (fixed seed). Run in the ancestra env (needs lxml):
  PYTHONPATH=src envs/ancestra/bin/python scripts/orthology/crosscheck_ancestra.py [--dir data/orthology/ancestra_O2]
"""
import argparse, json, os, random, sys, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "src"))
from ancestra.orthoxml import hogs_to_pairs, read_orthoxml
from pocketevo.orthology import HogIndex

ap = argparse.ArgumentParser()
ap.add_argument("--dir", default="data/orthology/ancestra_O2")
ap.add_argument("--n", type=int, default=50000)
ap.add_argument("--seed", type=int, default=0)
a = ap.parse_args()
t0 = time.time()
ox = read_orthoxml(os.path.join(a.dir, "accurate", "hogs.orthoxml.gz"))
orth, para = hogs_to_pairs(ox, with_paralogs=True)
idx = HogIndex.load(os.path.join(a.dir, "tables"))
rs = random.Random(a.seed)
so = rs.sample(sorted(orth), a.n)
sp = rs.sample(sorted(para), a.n)
res = dict(n_orth_pairs_total=len(orth), n_para_pairs_total=len(para),
           orth_agree=sum(idx.relation(x, y)["relation"] == "ortholog" for x, y in so),
           para_agree=sum(idx.relation(x, y)["relation"] == "paralog" for x, y in sp), n_sampled_each=a.n,
           seconds=round(time.time() - t0, 1))
os.makedirs(os.path.join(a.dir, "validation"), exist_ok=True)
json.dump(res, open(os.path.join(a.dir, "validation", "crosscheck_ancestra.json"), "w"), indent=1)
print(res)
