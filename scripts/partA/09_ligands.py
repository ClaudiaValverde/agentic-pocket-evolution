#!/usr/bin/env python
"""Part A stage 9: functional validation — shared experimental ligand class (DESIGN §4.2, prereg
pairs.functional_validation).

Site ligand annotation: BioLiP2 site (data/annotations/biolip_sites.parquet, UniProt numbering) whose residues lie
>= 50 % inside the atlas site (prereg f2 biolip2_residue_overlap_min); regular ligands only (no metal / peptide /
nucleic acid). Shared class for a pair: same CCD id, or ECFP4 (Morgan r=2, 2048 bit, RDKit) Tanimoto >= 0.5 for
any ligand pair.
Test set: significant N pairs (primary null, E < 0.01 both directions) with both sites annotated; proteins in
pocketcompare's training/validation data (external/pocketcompare-agents/data/splits: BLiPS train/val sites,
PPI-L train/val pockets) excluded (also reported without the exclusion).
Control: random pairs of annotated sites, cross-group, no shared join label, resolved by the same MMseqs2/Foldseek
rules to N (seeded), reweighted to the test set's distribution of covariate-predicted similarity (deciles of mu_ij;
'covariate-matched'). Reference rows: significant H1s / H3j / H4 pairs.
Output: results/partA/<snap>/ligand_enrichment.csv (+ json)
"""
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from _common import parser, resolve, logger, pairs
from pocketevo import homology as hom
from pocketevo import nulls

ap = parser(__doc__)
ap.add_argument("--biolip", default=None)
ap.add_argument("--n-control", type=int, default=30_000)
ap.add_argument("--threads", type=int, default=32)
ap.add_argument("--mmseqs", default=str(pairs.REPO / "envs/ancestra/bin/mmseqs"))
ap.add_argument("--foldseek", default=str(pairs.REPO / "envs/ancestra/bin/foldseek"))
ap.add_argument("--foldseek-db", default="/tmp/pe_big/partA/foldseek/afdb_v6")
a = resolve(ap.parse_args())
log = logger(a.work / "partA.log")
T0 = time.time()
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator
RDLogger.DisableLog("rdApp.*")
TANI = float(a.prereg["pairs"]["functional_validation"]["ecfp4_tanimoto"]["value"])
OVL = float(a.prereg["part_b"]["evidence_filters"]["f2_experimental"]["biolip2_residue_overlap_min"]["value"])

uni = pd.read_parquet(a.work / "universe.parquet", columns=["site_id", "accession", "residues"])
bl = pd.read_parquet(a.biolip or Path(a.annotations) / "biolip_sites.parquet")
bl = bl[bl.ligand_class == "regular"]
bl = bl[bl.accession.isin(set(uni.accession))]
by_acc = {k: g for k, g in bl.groupby("accession")}
site_ligs = {}
for si, (acc, res) in enumerate(zip(uni.accession.values, uni.residues.values)):
    g = by_acc.get(acc)
    if g is None:
        continue
    rs = set(np.asarray(res).tolist())
    L = set()
    for lig, br in zip(g.ligand.values, g.residues_uniprot.values):
        br = np.asarray(br)
        if len(br) and sum(int(x) in rs for x in br) / len(br) >= OVL:
            L.add(str(lig))
    if L:
        site_ligs[si] = L
log(f"annotated sites: {len(site_ligs)} (BioLiP regular rows {len(bl)})")

# ECFP4 per CCD
import gzip
smiles, n_bad_rows = {}, 0
with gzip.open("/tmp/pe_big/raw/biolip2/ligand.tsv.gz", "rt") as fh:
    head = fh.readline().lstrip("#").rstrip("\n").split("\t")
    k_ccd, k_smi = head.index("CCD"), head.index("SMILES")
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) < len(head):
            n_bad_rows += 1
            continue
        smiles[f[k_ccd]] = f[k_smi]   # extra tab-separated fields only occur after SMILES (names)
gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
fps, n_fail = {}, 0
for ccd in sorted({x for v in site_ligs.values() for x in v}):
    fp = None
    for smi in str(smiles.get(ccd, "")).split("; "):
        m = Chem.MolFromSmiles(smi) if smi and smi != "nan" else None
        if m is not None:
            fp = gen.GetFingerprint(m)
            break
    if fp is None:
        n_fail += 1
    fps[ccd] = fp
tani_cache = {}


def shared(si, sj):
    A, B = site_ligs.get(si), site_ligs.get(sj)
    if not A or not B:
        return np.nan
    if A & B:
        return 1.0
    for x in A:
        for y in B:
            k = (x, y) if x < y else (y, x)
            if k not in tani_cache:
                fx, fy = fps.get(x), fps.get(y)
                tani_cache[k] = DataStructs.TanimotoSimilarity(fx, fy) if fx is not None and fy is not None else 0.0
            if tani_cache[k] >= TANI:
                return 1.0
    return 0.0


# pocketcompare training/validation proteins
sp = pairs.REPO / "external/pocketcompare-agents/data/splits"
train_acc = set()
bs = pd.read_parquet(sp / "blips_sites_v2.parquet")
part = "partition_v2" if "partition_v2" in bs.columns else "partition_v1"
train_acc |= set(bs[bs[part].astype(str).isin(["train", "val"])].uniprot.dropna().astype(str))
for f in ("ppi_l_pairs_train.parquet", "ppi_l_pairs_val.parquet"):
    p = sp / f
    if p.exists():
        d = pd.read_parquet(p)
        for col in d.columns:
            if "uniprot" in col.lower():
                train_acc |= set(d[col].dropna().astype(str))
log(f"pocketcompare train/val accessions: {len(train_acc)}")
acc_of = uni.accession.to_numpy()

sig = pd.read_parquet(a.work / "cls_cos" / "pairs_cls.parquet", columns=["i", "j", "cls", "Emax_pooled", "gi", "gj"])
sig = sig[sig.Emax_pooled < 0.01]
ann = np.zeros(len(uni), bool)
ann[list(site_ligs)] = True
sig = sig[ann[sig.i.values] & ann[sig.j.values]].copy()
sig["train"] = np.isin(acc_of[sig.i.values], list(train_acc)) | np.isin(acc_of[sig.j.values], list(train_acc))
sig["shared"] = [shared(x, y) for x, y in zip(sig.i.values, sig.j.values)]

# covariate-predicted similarity for matching
model = nulls.CovModel.from_dict(json.loads((a.work / "null_cos" / "covmodel.json").read_text()))
cov = np.load(a.work / "covariates.npy")
cov = np.where(np.isfinite(cov), cov, np.nanmean(cov, 0))
Fn = model.norm(cov)


def mu_of(i, j):
    return model.mu_logsd_pairs(Fn[i], Fn[j])[0]


# control: random annotated-site pairs resolved to N
c = np.load(a.work / "context.npz")
ctx = hom.PairContext(c["site_prot"], c["site_lab"], c["prot_group"], np.full(len(c["prot_group"]), -1), None)  # only REST needed
annot = np.flatnonzero(ann)
rng = np.random.default_rng(a.seed + 20)
ci, cj = [], []
while sum(map(len, ci)) < a.n_control:
    x, y = annot[rng.integers(0, len(annot), 200_000)], annot[rng.integers(0, len(annot), 200_000)]
    cl = hom.classify_pairs(ctx, x, y)
    ok = cl == hom.REST
    ci.append(x[ok]); cj.append(y[ok])
ctl = pd.DataFrame({"i": np.concatenate(ci)[:a.n_control], "j": np.concatenate(cj)[:a.n_control]})
pa, pb = acc_of[ctl.i.values], acc_of[ctl.j.values]
ctl["qacc"], ctl["tacc"] = np.where(pa < pb, pa, pb), np.where(pa < pb, pb, pa)
upp = ctl[["qacc", "tacc"]].drop_duplicates()


def lookup(db, strip):
    lk = pd.read_csv(str(db) + ".lookup", sep="\t", header=None, names=["key", "name", "file"])
    nm = lk.name.astype(str)
    if strip:
        nm = nm.str.replace(r"^AF-", "", regex=True).str.replace(r"-F1-model_v\d+.*$", "", regex=True)
    return dict(zip(nm, lk.key)), dict(zip(lk.key, lk.name.astype(str)))


seqdb = a.work / "mmseqs" / "seqdb"
sl, _ = lookup(seqdb, False)
ms = hom.run_pair_search(pd.DataFrame({"qname": upp.qacc.values, "tname": upp.tacc.values}), str(seqdb), sl,
                         a.work / "mmseqs" / "ligand_ctl", "mmseqs", a.mmseqs, threads=a.threads, log=log)
fa, fk = lookup(a.foldseek_db, True)
fname = {x: fk[k] for x, k in fa.items()}
fs = hom.run_pair_search(pd.DataFrame({"qname": upp.qacc.map(fname).fillna("__m__").values,
                                       "tname": upp.tacc.map(fname).fillna("__m__").values}),
                         a.foldseek_db, {v: k for k, v in fk.items()}, a.work / "foldseek_aln" / "ligand_ctl",
                         "foldseek", a.foldseek, threads=a.threads, log=log)
h3 = (ms.status == "hit") & (ms.evalue < 1e-3)
h4 = (fs.status == "hit") & ((fs.evalue < 1e-3) | (np.fmax(fs.qtmscore, fs.ttmscore) >= 0.5))
unk = (ms.status == "missing") | (fs.status == "missing")
upp["cls"] = np.where(h3.values, "H3s", np.where(h4.values, "H4", np.where(unk.values, "unknown", "N")))
ctl = ctl.merge(upp, on=["qacc", "tacc"])
ctl = ctl[ctl.cls == "N"].copy()
ctl["train"] = np.isin(ctl.qacc.values, list(train_acc)) | np.isin(ctl.tacc.values, list(train_acc))
ctl["shared"] = [shared(x, y) for x, y in zip(ctl.i.values, ctl.j.values)]
ctl["mu"] = mu_of(ctl.i.values, ctl.j.values)
sig["mu"] = mu_of(sig.i.values, sig.j.values)

rows = []
for excl in (True, False):
    T = sig[(sig.cls == "N") & (~sig.train if excl else True)]
    C = ctl[~ctl.train] if excl else ctl
    # decile matching on mu: control weights = test share / control share per bin
    edges = np.unique(np.quantile(T.mu, np.linspace(0, 1, 11))) if len(T) >= 10 else np.array([-np.inf, np.inf])
    tb = np.clip(np.searchsorted(edges, T.mu, side="right") - 1, 0, len(edges) - 2)
    cb = np.clip(np.searchsorted(edges, C.mu, side="right") - 1, 0, len(edges) - 2)
    tsh = np.bincount(tb, minlength=len(edges) - 1) / max(len(T), 1)
    csh = np.bincount(cb, minlength=len(edges) - 1) / max(len(C), 1)
    w = np.where(csh[cb] > 0, tsh[cb] / np.where(csh[cb] > 0, csh[cb], 1), 0.0)
    p_t = float(T.shared.mean()) if len(T) else np.nan
    p_c = float(np.sum(w * C.shared) / np.sum(w)) if np.sum(w) > 0 else np.nan
    p_c_raw = float(C.shared.mean()) if len(C) else np.nan
    k_t = float(T.shared.sum())
    lo, hi = (np.nan, np.nan)
    if len(T):
        from pocketevo.estimands import wilson
        lo, hi = wilson(k_t, len(T))
    # bootstrap CI of the enrichment ratio (pairs resampled within each set; seeded)
    rb = np.random.default_rng(7)
    ratios = []
    for _ in range(1000):
        if len(T) == 0 or len(C) == 0:
            break
        ti = rb.integers(0, len(T), len(T)); cii = rb.integers(0, len(C), len(C))
        pt = T.shared.values[ti].mean(); pc = np.sum(w[cii] * C.shared.values[cii]) / max(np.sum(w[cii]), 1e-12)
        ratios.append(pt / pc if pc > 0 else np.nan)
    ratios = np.array(ratios, float)
    rows.append({"exclude_pocketcompare_train": excl, "set": "sig_N_vs_matched_random_N", "n_test": len(T),
                 "n_control": len(C), "shared_test": p_t, "shared_test_wilson_lo": lo, "shared_test_wilson_hi": hi,
                 "shared_control_matched": p_c, "shared_control_unmatched": p_c_raw,
                 "enrichment": p_t / p_c if p_c and p_c > 0 else np.nan,
                 "enrichment_ci_lo": np.nanquantile(ratios, 0.025) if len(ratios) else np.nan,
                 "enrichment_ci_hi": np.nanquantile(ratios, 0.975) if len(ratios) else np.nan})
    for ref in ("H1s", "H2s", "H3j", "H3s", "H4", "H1x", "H2x"):
        R = sig[(sig.cls == ref) & (~sig.train if excl else True)]
        rows.append({"exclude_pocketcompare_train": excl, "set": f"sig_{ref}", "n_test": len(R),
                     "shared_test": float(R.shared.mean()) if len(R) else np.nan})
df = pd.DataFrame(rows)
df.to_csv(a.results / "ligand_enrichment.csv", index=False)
summ = {"annotated_sites": len(site_ligs), "ccd_fp_fail": n_fail, "n_ccd": len(fps), "train_accessions": len(train_acc),
        "control_resolved_classes": upp.cls.value_counts().to_dict(), "wall_s": time.time() - T0}
(a.results / "ligand_enrichment.json").write_text(json.dumps(summ, indent=1, default=str))
pairs.write_sidecar(a.results / "ligand_enrichment.csv", a.snapshot, a.seed + 20, {"stage": "ligands", **summ})
log("ligands\n" + df.to_string())
