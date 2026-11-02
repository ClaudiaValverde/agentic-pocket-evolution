#!/usr/bin/env python
"""Part A v2 stage 12: functional convergence estimate F_conv,func (PI request v2.2) and shared-ligand rates.

Experimental site labels (a label is attached to an atlas site when >= 50 % of the label's residues lie in the site;
prereg f2 biolip2_residue_overlap_min, used for all three sources):
  * BioLiP2 (PDB holo): CCD id; metal/peptide/nucleic acid classes as given (peptide, DNA, RNA dropped)
  * UniProt binding sites with experimental evidence only (ECO:0000269 / ECO:0007744): ChEBI ligand; residues of
    all experimental features of one (protein, ligand) are pooled; generic 'substrate' / DNA / RNA dropped;
    common ligands mapped to CCD ids for cross-source matching (config/partA_ligand_map.tsv)
  * M-CSA reference (literature-curated) catalytic residues: 4-digit EC of the entry
Categories: metal | cofactor | small_molecule | catalytic_EC.
Shared class of a pair: a common token (CCD / ChEBI / EC) or ECFP4 Tanimoto >= 0.5 between any two non-metal CCDs.

Definitions (all on the primary significant set: pooled tail, E < 0.01 both directions):
  s_sig   = share of annotated significant N pairs with a shared class
  s_ctl   = same share among random N pairs of annotated sites (MMseqs2/Foldseek-negative), reweighted to the
            test set's deciles of covariate-predicted similarity mu_ij (covariate-matched)
  delta   = s_sig - s_ctl, computed per pair and per family pair (each root-group cell of the test / control set
            weight 1; v2.1, verifier item 6)
  F_conv,func = F_analog x max(delta, 0) with matching weightings (headline: per-family-pair F_analog x
            per-family-pair delta; v2 used per-pair delta, kept as F_conv_func.per_family_pair_x_pair_delta_v2);
            assumes the annotated subset is representative of significant N pairs (stated)
  label variants (sensitivity): full | no_tanimoto (common tokens only) | no_metal (metal tokens dropped) |
            tokens_no_metal; category rows are EXPLORATORY (post hoc category scaling, ledger D18 correction)
  AF      = delta / s_sig      attributable fraction of the observed sharing
CIs: two-way bootstrap over root groups with the same Poisson multiplicities as 07_estimands (seed*1000 + r), so
F_analog and the shares are resampled jointly (pair weight w_g w_h, w_g within a group).
pocketcompare training/validation proteins are excluded (also reported without the exclusion).
"""
import gzip
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from _common import parser, resolve, logger, pairs
from pocketevo import homology as hom
from pocketevo import nulls
from pocketevo.estimands import pair_group_weight, wilson

ap = parser(__doc__)
ap.add_argument("--n-control", type=int, default=60_000)
ap.add_argument("--n-boot", type=int, default=200)
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
ANN = Path(a.annotations)
summ = {}

uni = pd.read_parquet(a.work / "universe.parquet", columns=["site_id", "accession", "residues", "plddt_mean_rep"])
acc_of = uni.accession.to_numpy()
site_rows = {}
for si, acc in enumerate(acc_of):
    site_rows.setdefault(acc, []).append(si)
RES = [set(np.asarray(r).tolist()) for r in uni.residues.values]

lm = pd.read_csv(pairs.REPO / "config/partA_ligand_map.tsv", sep="\t", comment="#", dtype=str).fillna("")
name2ccd = {r.uniprot_name: r.ccd for r in lm.itertuples() if r.uniprot_name}
cat_of_ccd = {r.ccd: r.category for r in lm.itertuples()}
site_tok: dict = {}


def attach(acc, residues, token):
    residues = [int(x) for x in residues]
    if not residues:
        return 0
    n = 0
    for si in site_rows.get(acc, []):
        if sum(r in RES[si] for r in residues) / len(residues) >= OVL:
            site_tok.setdefault(si, set()).add(token)
            n += 1
    return n


tok_cat = {}
# BioLiP2
bl = pd.read_parquet(ANN / "biolip_sites.parquet", columns=["accession", "ligand", "ligand_class", "residues_uniprot"])
bl = bl[bl.ligand_class.isin(["regular", "metal"]) & bl.accession.isin(site_rows.keys())]
n_bl = 0
for acc, lig, lc, res in zip(bl.accession.values, bl.ligand.values, bl.ligand_class.values, bl.residues_uniprot.values):
    t = f"CCD:{lig}"
    tok_cat[t] = "metal" if lc == "metal" else cat_of_ccd.get(str(lig), "small_molecule")
    n_bl += attach(acc, res, t)
# UniProt experimental binding sites
uf = pd.read_parquet(ANN / "uniprot_site_features.parquet",
                     columns=["uniprot_acc", "type", "start", "end", "ligand", "chebi", "experimental"])
uf = uf[(uf.type == "Binding site") & uf.experimental & uf.chebi.notna() & (uf.ligand != "substrate")
        & ~uf.ligand.isin(["DNA", "RNA"])]
upm = pd.read_parquet(ANN / "uniprot_proteins.parquet", columns=["accession", "uniprot_acc"])
uf = uf.merge(upm, on="uniprot_acc")
uf = uf[uf.accession.isin(site_rows.keys())]
n_up = 0
for (acc, lig, chebi), g in uf.groupby(["accession", "ligand", "chebi"]):
    res = sorted({p for s_, e_ in zip(g.start, g.end) if pd.notna(s_) and pd.notna(e_)
                  for p in range(int(s_), int(e_) + 1)})
    if lig in name2ccd:
        t = f"CCD:{name2ccd[lig]}"
        tok_cat[t] = cat_of_ccd.get(name2ccd[lig], "small_molecule")
    else:
        t = f"CHEBI:{chebi}"
        tok_cat[t] = "metal" if ("(2+)" in lig or "(3+)" in lig or "(+)" in lig and len(lig) <= 6
                                 or lig.endswith(" cation")) else "small_molecule"
    n_up += attach(acc, res, t)
# M-CSA reference catalytic residues
mr = pd.read_parquet(ANN / "mcsa_residues.parquet", columns=["mcsa_id", "accession", "resid", "is_reference"])
me = pd.read_parquet(ANN / "mcsa_entries.parquet", columns=["mcsa_id", "ec"])
mr = mr[mr.is_reference & mr.accession.isin(site_rows.keys())].merge(me, on="mcsa_id")
n_mc = 0
for (acc, ec), g in mr.groupby(["accession", "ec"]):
    if isinstance(ec, str) and ec.count(".") == 3:
        t = f"EC:{ec}"
        tok_cat[t] = "catalytic_EC"
        n_mc += attach(acc, g.resid.dropna().astype(int).tolist(), t)
summ["site_labels_attached"] = {"biolip2": n_bl, "uniprot_exp": n_up, "mcsa_ref": n_mc}
summ["annotated_sites"] = len(site_tok)
log("annotated sites", len(site_tok), json.dumps(summ["site_labels_attached"]))
pd.DataFrame({"si": list(site_tok), "tokens": [sorted(v) for v in site_tok.values()]}).to_parquet(
    a.work / "functional_site_labels.parquet", index=False)

# ECFP4 for CCD tokens (non-metal)
smiles = {}
with gzip.open("/tmp/pe_big/raw/biolip2/ligand.tsv.gz", "rt") as fh:
    head = fh.readline().lstrip("#").rstrip("\n").split("\t")
    k_ccd, k_smi = head.index("CCD"), head.index("SMILES")
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) >= len(head):
            smiles[f[k_ccd]] = f[k_smi]
gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
fps = {}
for t, c in tok_cat.items():
    if t.startswith("CCD:") and c != "metal":
        for smi in str(smiles.get(t[4:], "")).split("; "):
            m = Chem.MolFromSmiles(smi) if smi else None
            if m is not None:
                fps[t] = gen.GetFingerprint(m)
                break
tc = {}


def shared(si, sj, cat=None, tanimoto=True, metal=True):
    """1.0 shared / 0.0 not / nan if either side has no (eligible) label.
    tanimoto=False: common tokens only; metal=False: metal tokens dropped (generic-metal sensitivity)."""
    A, B = site_tok.get(si, set()), site_tok.get(sj, set())
    if cat is not None:
        A = {x for x in A if tok_cat[x] == cat}
        B = {x for x in B if tok_cat[x] == cat}
    if not metal:
        A = {x for x in A if tok_cat[x] != "metal"}
        B = {x for x in B if tok_cat[x] != "metal"}
    if not A or not B:
        return np.nan
    if A & B:
        return 1.0
    if not tanimoto:
        return 0.0
    for x in A:
        fx = fps.get(x)
        if fx is None:
            continue
        for y in B:
            fy = fps.get(y)
            if fy is None:
                continue
            k = (x, y) if x < y else (y, x)
            if k not in tc:
                tc[k] = DataStructs.TanimotoSimilarity(fx, fy)
            if tc[k] >= TANI:
                return 1.0
    return 0.0


# pocketcompare training/validation proteins
sp = pairs.REPO / "external/pocketcompare-agents/data/splits"
bs = pd.read_parquet(sp / "blips_sites_v2.parquet")
part = "partition_v2" if "partition_v2" in bs.columns else "partition_v1"
train_acc = set(bs[bs[part].astype(str).isin(["train", "val"])].uniprot.dropna().astype(str))
for f in ("ppi_l_pairs_train.parquet", "ppi_l_pairs_val.parquet"):
    if (sp / f).exists():
        d = pd.read_parquet(sp / f)
        for col in d.columns:
            if "uniprot" in col.lower():
                train_acc |= set(d[col].dropna().astype(str))
is_train = np.isin(acc_of, list(train_acc))

ann = np.zeros(len(uni), bool)
ann[list(site_tok)] = True
sig = pd.read_parquet(a.work / f"cls_cos{a.tag}" / "pairs_cls.parquet", columns=["i", "j", "cls", "Emax_pooled", "gi", "gj"])
sig = sig[(sig.Emax_pooled < 0.01) & (sig.cls != "H0")]
summ["sig_pairs_by_class"] = sig.cls.value_counts().to_dict()
sig = sig[ann[sig.i.values] & ann[sig.j.values]].copy()
sig["train"] = is_train[sig.i.values] | is_train[sig.j.values]
model = nulls.CovModel.from_dict(json.loads((a.work / "null_cos" / "covmodel.json").read_text()))
cov = np.load(a.work / "covariates.npy")
Fn = model.norm(np.where(np.isfinite(cov), cov, np.nanmean(cov, 0)))
CATS = [None, "metal", "cofactor", "small_molecule", "catalytic_EC"]


LABEL_VARIANTS = {"full": {}, "no_tanimoto": {"tanimoto": False}, "no_metal": {"metal": False},
                  "tokens_no_metal": {"tanimoto": False, "metal": False}}


def label_cols(df):
    for c in CATS:
        df["shared" if c is None else f"shared_{c}"] = [shared(x, y, c) for x, y in zip(df.i.values, df.j.values)]
    for v, kw in LABEL_VARIANTS.items():
        if v != "full":
            df[f"shared__{v}"] = [shared(x, y, None, **kw) for x, y in zip(df.i.values, df.j.values)]
    # mechanism of a 'shared' call: common non-metal token / common metal token only / ECFP4 Tanimoto only
    mech = []
    for x, y, sh in zip(df.i.values, df.j.values, df["shared"].values):
        if sh != 1.0:
            mech.append(None); continue
        com = site_tok[x] & site_tok[y]
        if any(tok_cat[t] != "metal" for t in com):
            mech.append("token_nonmetal")
        elif com:
            mech.append("token_metal_only")
        else:
            mech.append("tanimoto_only")
    df["shared_mechanism"] = mech
    df["mu"] = model.mu_logsd_pairs(Fn[df.i.values], Fn[df.j.values])[0]
    plddt = uni.plddt_mean_rep.to_numpy()
    df["plddt90"] = (plddt[df.i.values] >= 90) & (plddt[df.j.values] >= 90)
    return df


sig = label_cols(sig)

# control: random annotated-site pairs, resolved to N by the same searches
c = np.load(a.work / "context.npz")
ctx = hom.PairContext(c["site_prot"], c["site_lab"], c["prot_group"], np.full(len(c["prot_group"]), -1), None)
site_group = c["prot_group"][c["site_prot"]]
annot = np.flatnonzero(ann)
rng = np.random.default_rng(a.seed + 20)
ci, cj, drawn = [], [], 0
while sum(map(len, ci)) < a.n_control:
    x, y = annot[rng.integers(0, len(annot), 200_000)], annot[rng.integers(0, len(annot), 200_000)]
    drawn += 200_000
    ok = hom.classify_pairs(ctx, x, y) == hom.REST
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
                         a.work / "mmseqs" / "func_ctl", "mmseqs", a.mmseqs, threads=a.threads, log=log)
fa, fk = lookup(a.foldseek_db, True)
fname = {x: fk[k] for x, k in fa.items()}
fs = hom.run_pair_search(pd.DataFrame({"qname": upp.qacc.map(fname).fillna("__m__").values,
                                       "tname": upp.tacc.map(fname).fillna("__m__").values}),
                         a.foldseek_db, {v: k for k, v in fk.items()}, a.work / "foldseek_aln" / "func_ctl",
                         "foldseek", a.foldseek, threads=a.threads, log=log)
h3 = (ms.status == "hit") & (ms.evalue < 1e-3)
h4 = (fs.status == "hit") & ((fs.evalue < 1e-3) | (np.fmax(fs.qtmscore, fs.ttmscore) >= 0.5))
unk = (ms.status == "missing") | (fs.status == "missing")
upp["cls"] = np.where(h3.values, "H3s", np.where(h4.values, "H4", np.where(unk.values, "unknown", "N")))
summ["control_resolved_classes"] = upp.cls.value_counts().to_dict()
ctl = ctl.merge(upp, on=["qacc", "tacc"])
ctl = ctl[ctl.cls == "N"].copy()
ctl["train"] = is_train[ctl.i.values] | is_train[ctl.j.values]
ctl["gi"], ctl["gj"] = site_group[ctl.i.values], site_group[ctl.j.values]
ctl = label_cols(ctl)
summ["control_drawn"] = drawn


def match_w(T, C):
    if len(T) < 10 or len(C) == 0:
        return np.ones(len(C))
    edges = np.unique(np.quantile(T.mu, np.linspace(0, 1, 11)))
    tb = np.clip(np.searchsorted(edges, T.mu, side="right") - 1, 0, len(edges) - 2)
    cb = np.clip(np.searchsorted(edges, C.mu, side="right") - 1, 0, len(edges) - 2)
    tsh = np.bincount(tb, minlength=len(edges) - 1) / len(T)
    csh = np.bincount(cb, minlength=len(edges) - 1) / len(C)
    return np.where(csh[cb] > 0, tsh[cb] / np.where(csh[cb] > 0, csh[cb], 1), 0.0)


def share(df, col, w):
    v = df[col].to_numpy()
    ok = ~np.isnan(v)
    return float(np.sum(w[ok] * v[ok]) / np.sum(w[ok])) if np.sum(w[ok]) > 0 else np.nan


# bootstrap replicates of F_analog from 07 (same seeds); strata reps if present
G = int(c["prot_group"].max()) + 1
reps = pd.read_parquet(a.work / f"bootstrap_cos{a.tag}.parquet")
srep_p = a.work / f"bootstrap_strata_cos{a.tag}.parquet"
srep = pd.read_parquet(srep_p) if srep_p.exists() else None
est_dir = a.results
est = pd.read_csv(est_dir / f"estimands_cos{a.tag}.csv")
prim = est[(est.tag == "all") & (est.variant == "pooled") & (est.E == 0.01) & (est.q0kind == "nominal")]
strata = pd.read_csv(est_dir / f"strata_cos{a.tag}.csv")
F_point = {sch: float(prim[prim.key == f"{sch}.F_adj.F_analog"].value.iloc[0]) for sch in ("per_family_pair", "per_pair")}
F90 = {sch: float(strata[(strata.tag == "both_plddt_ge90") & (strata.key == f"{sch}.F_adj.F_analog")].value.iloc[0])
       for sch in ("per_family_pair", "per_pair")}
nb = min(a.n_boot, len(reps))
W = [np.random.default_rng(a.seed * 1000 + r).poisson(1.0, G).astype(float) for r in range(nb)]

site_half = np.load(a.work / "site_half.npy") if (a.work / "site_half.npy").exists() else None
Fhalf = None
if site_half is not None:
    Fhalf = {sch: float(strata[(strata.tag == "test_half") & (strata.key == f"{sch}.F_adj.F_analog")].value.iloc[0])
             for sch in ("per_family_pair", "per_pair")}


def cell_w(df):
    """per-family-pair weights: each unordered root-group cell of the set gets total weight 1."""
    if len(df) == 0:
        return np.zeros(0)
    cell = np.minimum(df.gi.values, df.gj.values).astype(np.int64) * (1 << 31) + np.maximum(df.gi.values, df.gj.values)
    cs = pd.Series(cell)
    return 1.0 / cs.map(cs.value_counts()).to_numpy()


def n_cells(df):
    if len(df) == 0:
        return 0
    return int(len(np.unique(np.minimum(df.gi.values, df.gj.values).astype(np.int64) * (1 << 31)
                             + np.maximum(df.gi.values, df.gj.values))))


def f_rep(subset, sch, k):
    if subset == "all":
        return reps[f"{sch}.F_adj.F_analog"].iloc[k]
    tag = {"plddt_ge90": "both_plddt_ge90", "test_half": "test_half"}[subset]
    if srep is None:
        return np.nan
    q = srep[(srep.tag == tag) & (srep.rep == k)]
    return q[f"{sch}.F_adj.F_analog"].iloc[0] if len(q) else np.nan


rows = []
SUBSETS = ["all", "plddt_ge90"] + (["test_half"] if site_half is not None else [])
for excl in (True, False):
    for subset in SUBSETS:
        T = sig[(sig.cls == "N") & (~sig.train if excl else True)]
        C = ctl[~ctl.train] if excl else ctl
        if subset == "plddt_ge90":
            T, C = T[T.plddt90], C[C.plddt90]
        if subset == "test_half":
            T = T[(site_half[T.i.values] == 1) & (site_half[T.j.values] == 1)]
            C = C[(site_half[C.i.values] == 1) & (site_half[C.j.values] == 1)]
        mw = match_w(T, C)
        n_any = int(T["shared"].notna().sum())
        Fsub = {"all": F_point, "plddt_ge90": F90, "test_half": Fhalf}[subset]
        specs = [(None, "full")] + [(None, v) for v in LABEL_VARIANTS if v != "full"] + \
                ([(c_, "full") for c_ in CATS[1:]] if subset == "all" else [])
        for cat, var in specs:
            col = ("shared" if var == "full" else f"shared__{var}") if cat is None else f"shared_{cat}"
            nT = int(T[col].notna().sum()); nC = int(C[col].notna().sum())
            frac_c = nT / n_any if (cat is not None and n_any) else 1.0
            r = {"exclude_train": excl, "subset": subset, "category": cat or "any", "label_variant": var,
                 "exploratory": cat is not None, "n_test": nT, "n_control": nC, "frac_of_annotated": frac_c,
                 "n_test_cells": n_cells(T[T[col].notna()])}
            wT_fam, wC_fam = cell_w(T), cell_w(C)
            for wname, wt0, wc0 in (("pair", np.ones(len(T)), mw), ("family_pair", wT_fam, wC_fam * mw)):
                sT, sC = share(T, col, wt0), share(C, col, wc0)
                d = sT - sC if np.isfinite(sT) and np.isfinite(sC) else np.nan
                r[f"s_sig.{wname}"], r[f"s_ctl_matched.{wname}"], r[f"delta.{wname}"] = sT, sC, d
                r[f"enrichment.{wname}"] = sT / sC if sC and sC > 0 else np.nan
                if wname == "pair":
                    lo, hi = wilson(np.nansum(T[col]), nT) if nT else (np.nan, np.nan)
                    r["s_sig.pair.wilson_lo"], r["s_sig.pair.wilson_hi"] = lo, hi
                    r["s_ctl_unmatched.pair"] = share(C, col, np.ones(len(C)))
                    r["AF.pair"] = d / sT if sT and sT > 0 else np.nan
            # F_conv,func: per-family-pair F_analog x per-family-pair delta (headline, v2.1); v2 form (x per-pair
            # delta) kept for comparison; per-pair F x per-pair delta
            d_fam, d_pair = r["delta.family_pair"], r["delta.pair"]
            r["F_conv_func.per_family_pair"] = Fsub["per_family_pair"] * frac_c * max(d_fam, 0.0) if np.isfinite(d_fam) else np.nan
            r["F_conv_func.per_family_pair_x_pair_delta_v2"] = Fsub["per_family_pair"] * frac_c * max(d_pair, 0.0) if np.isfinite(d_pair) else np.nan
            r["F_conv_func.per_pair"] = Fsub["per_pair"] * frac_c * max(d_pair, 0.0) if np.isfinite(d_pair) else np.nan
            if nT >= 5 and nC >= 5:
                bd = {"pair": [], "family_pair": []}
                bF = {"per_family_pair": [], "per_pair": []}
                for k, w in enumerate(W):
                    wt = pair_group_weight(w, T.gi.values, T.gj.values)
                    wc = pair_group_weight(w, C.gi.values, C.gj.values) * mw
                    dp = share(T, col, wt) - share(C, col, wc)
                    dfm = share(T, col, wt * wT_fam) - share(C, col, wc * wC_fam)
                    bd["pair"].append(dp); bd["family_pair"].append(dfm)
                    bF["per_family_pair"].append(f_rep(subset, "per_family_pair", k) * frac_c * max(dfm, 0.0))
                    bF["per_pair"].append(f_rep(subset, "per_pair", k) * frac_c * max(dp, 0.0))
                for kk, v in bd.items():
                    r[f"delta.{kk}.ci_lo"], r[f"delta.{kk}.ci_hi"] = np.nanquantile(v, [0.025, 0.975])
                for kk, v in bF.items():
                    v = np.array(v, float)
                    if np.isfinite(v).any():
                        r[f"F_conv_func.{kk}.ci_lo"], r[f"F_conv_func.{kk}.ci_hi"] = np.nanquantile(v, [0.025, 0.975])
            if cat is None and var == "full":
                sh = T[T["shared"] == 1.0]
                vc = sh.shared_mechanism.value_counts(normalize=True)
                for m in ("token_nonmetal", "token_metal_only", "tanimoto_only"):
                    r[f"mechanism_share.{m}"] = float(vc.get(m, 0.0))
            rows.append(r)
    for ref in ("H1s", "H1x", "H2s", "H2x", "H12s", "H3j", "H3s", "H4"):
        R_ = sig[(sig.cls == ref) & (~sig.train if excl else True)]
        rows.append({"exclude_train": excl, "subset": "all", "category": "any", "label_variant": "full",
                     "set": f"sig_{ref}", "n_test": int(R_.shared.notna().sum()),
                     "s_sig.pair": share(R_, "shared", np.ones(len(R_))),
                     "s_sig.pair.no_tanimoto": share(R_, "shared__no_tanimoto", np.ones(len(R_))),
                     "s_sig.pair.tokens_no_metal": share(R_, "shared__tokens_no_metal", np.ones(len(R_)))})
df = pd.DataFrame(rows)
df["set"] = df.get("set").fillna("sig_N_vs_matched_random_N") if "set" in df else "sig_N_vs_matched_random_N"
df.to_csv(a.results / f"functional_convergence{a.tag}.csv", index=False)
summ["wall_s"] = time.time() - T0
summ["n_boot"] = nb
summ["train_accessions"] = len(train_acc)
(a.results / f"functional_convergence{a.tag}.json").write_text(json.dumps(summ, indent=1, default=str))
pairs.write_sidecar(a.results / f"functional_convergence{a.tag}.csv", a.snapshot, [a.seed + 20, a.seed],
                    {"stage": "functional_convergence", **summ})
pd.set_option("display.width", 250)
log("functional\n" + df[["exclude_train", "subset", "category", "label_variant", "set", "n_test", "s_sig.pair",
                        "s_ctl_matched.pair", "delta.pair", "delta.family_pair", "F_conv_func.per_family_pair",
                        "F_conv_func.per_pair"]].to_string())
