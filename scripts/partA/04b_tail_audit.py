#!/usr/bin/env python
"""Part A stage 4b: what is in the deep tail of the held-out null? (calibration diagnostic)

Held-out random 'non-homologous' pairs are non-homologous only by joins (root group, Pfam clan, TED CATH H,
PANTHER family). Pairs with p < 1e-4 in either direction, plus a seeded random 3,000 held-out pairs as reference,
are resolved by the H3/H4 search rules (MMseqs2 E < 1e-3; Foldseek E < 1e-3 or max TM >= 0.5). The calibration
table is then recomputed for pairs resolved to N only: obs_N(t) = obs(t) * frac_N(tail bin) vs
expected_N(t) = n * frac_N(reference) * t.
Output: results/partA/<snap>/null_<score>_tail_audit.csv
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

from _common import parser, resolve, logger, pairs, build_seqdb
from pocketevo import homology as hom

ap = parser(__doc__)
ap.add_argument("--score", default="cos")
ap.add_argument("--variant", default="pooled")
ap.add_argument("--threads", type=int, default=24)
ap.add_argument("--mmseqs", default=str(pairs.REPO / "envs/ancestra/bin/mmseqs"))
ap.add_argument("--foldseek", default=str(pairs.REPO / "envs/ancestra/bin/foldseek"))
ap.add_argument("--foldseek-db", default="/tmp/pe_big/partA/foldseek/afdb_v6")
ap.add_argument("--sequences", default="/tmp/pe_big/work/annotations/uniprot_sequences.parquet")
a = resolve(ap.parse_args())
log = logger(a.work / "partA.log")
cal = pd.read_parquet(a.work / f"null_{a.score}" / "calib_heldout.parquet")
uni = pd.read_parquet(a.work / "universe.parquet", columns=["accession"])
ACC = uni.accession.to_numpy()
v = a.variant
pmin = np.minimum(cal[f"p_{v}_ij"], cal[f"p_{v}_ji"])
tail = cal[pmin < 1e-4].copy()
ref = cal.sample(3000, random_state=0)
sel = pd.concat([tail.assign(set="tail"), ref.assign(set="reference")])
pa, pb = ACC[sel.i.values], ACC[sel.j.values]
sel["qacc"], sel["tacc"] = np.where(pa < pb, pa, pb), np.where(pa < pb, pb, pa)

seqdb = a.work / "mmseqs" / "seqdb"
build_seqdb(a, pd.read_parquet(a.work / "proteins.parquet", columns=["accession"]).accession.to_numpy(), log)


def lookup(db, strip):
    lk = pd.read_csv(str(db) + ".lookup", sep="\t", header=None, names=["key", "name", "file"])
    nm = lk.name.astype(str)
    if strip:
        nm = nm.str.replace(r"^AF-", "", regex=True).str.replace(r"-F1-model_v\d+.*$", "", regex=True)
    return dict(zip(nm, lk.key)), dict(zip(lk.key, lk.name.astype(str)))


upp = sel[["qacc", "tacc"]].drop_duplicates().reset_index(drop=True)
sl, _ = lookup(seqdb, False)
ms = hom.run_pair_search(pd.DataFrame({"qname": upp.qacc, "tname": upp.tacc}), str(seqdb), sl,
                         a.work / "mmseqs" / f"tail_{a.score}", "mmseqs", a.mmseqs, threads=a.threads, log=log)
fa, fk = lookup(a.foldseek_db, True)
fname = {x: fk[k] for x, k in fa.items()}
fs = hom.run_pair_search(pd.DataFrame({"qname": upp.qacc.map(fname).fillna("__m__"),
                                       "tname": upp.tacc.map(fname).fillna("__m__")}),
                         a.foldseek_db, {vv: k for k, vv in fk.items()}, a.work / "foldseek_aln" / f"tail_{a.score}",
                         "foldseek", a.foldseek, threads=a.threads, log=log)
h3 = (ms.status == "hit") & (ms.evalue < 1e-3)
h4 = (fs.status == "hit") & ((fs.evalue < 1e-3) | (np.fmax(fs.qtmscore, fs.ttmscore) >= 0.5))
unk = (ms.status == "missing") | (fs.status == "missing")
upp["cls"] = np.where(h3, "H3s", np.where(h4, "H4", np.where(unk, "unknown", "N")))
sel = sel.merge(upp, on=["qacc", "tacc"], how="left")
sel["pmin"] = np.minimum(sel[f"p_{v}_ij"], sel[f"p_{v}_ji"])
sel.to_parquet(a.work / f"null_{a.score}" / "tail_audit_pairs.parquet", index=False)
refN = float((sel[sel.set == "reference"].cls == "N").mean())
rows = []
n = len(cal)
for t in (1e-4, 1e-5, 1e-6):
    for d in ("ij", "ji"):
        hit = sel[(sel.set == "tail") & (sel[f"p_{v}_{d}"] < t)]
        k = len(hit)
        fN = float((hit.cls == "N").mean()) if k else np.nan
        comp = hit.cls.value_counts().to_dict()
        rows.append({"variant": v, "p_threshold": t, "direction": d, "obs_all": k, "expected_all": n * t,
                     "ratio_all": k / (n * t), "frac_N_in_tail": fN, "frac_N_reference": refN,
                     "obs_N": k * fN if k else 0.0, "expected_N": n * refN * t,
                     "ratio_N": (k * fN) / (n * refN * t) if k else 0.0, "composition": json.dumps(comp)})
df = pd.DataFrame(rows)
df.to_csv(a.results / f"null_{a.score}_tail_audit.csv", index=False)
pairs.write_sidecar(a.results / f"null_{a.score}_tail_audit.csv", a.snapshot, 0,
                    {"stage": "tail_audit", "reference_classes": sel[sel.set == "reference"].cls.value_counts().to_dict()})
log("tail audit\n" + df.to_string())
