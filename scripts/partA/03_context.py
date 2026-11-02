#!/usr/bin/env python
"""Part A stage 3: homology context and site covariates.

Outputs (<work>/):
  proteins.parquet   accession, pi, tier, gene (forest leaf or -1), group, group_kind, panther_fam, panther_sf, n_sites
  context.npz        site_prot, site_lab (join labels), prot_group, prot_gene, prot_lab (protein-level labels)
  covariates.npy     (N, 6) float32 site covariates (nulls.COV_NAMES); NaN rows counted (cov_ok False)
  labels_vocab.json  label code -> string
Inputs: universe, orthology run (--orthology/hogs.orthoxml.gz), optional placement table
(<orthology>/../placement*.parquet with columns accession, hog), UniProt proteins (PANTHER) + sequences,
Pfam regions + clans, TED domains.
"""
import glob
import json
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

from _common import parser, resolve, logger, pairs
from pocketevo import homology as hom
from pocketevo import nulls

ap = parser(__doc__)
ap.add_argument("--sequences", default="/tmp/pe_big/work/annotations/uniprot_sequences.parquet")
ap.add_argument("--placement", default=None, help="placement table (accession, hog); default: auto-detect")
ap.add_argument("--workers", type=int, default=32)
ap.add_argument("--idgap", default=None, help="O2 idgap_mapping.parquet (atlas_accession -> hog_member_accession)")
a = resolve(ap.parse_args())
log = logger(a.work / "partA.log")
t0 = time.time()
funnel = {}

uni = pd.read_parquet(a.work / "universe.parquet", columns=["si", "site_id", "accession", "tier", "residues",
                                                            "rep_site_volume_A3", "max_burrow_score",
                                                            "plddt_mean_rep"])
N = len(uni)
acc_u, site_prot = np.unique(uni.accession.to_numpy().astype(str), return_inverse=True)
P = len(acc_u)
prot = pd.DataFrame({"accession": acc_u, "pi": np.arange(P)})
prot["tier"] = uni.groupby(site_prot).tier.first().to_numpy()
prot["n_sites"] = np.bincount(site_prot, minlength=P)
log(f"context: {N} sites, {P} proteins")

# ---------------------------------------------------------------- annotations
ann = Path(a.annotations)
up_path = ann / "uniprot_proteins.parquet"
if up_path.exists():
    up = pd.read_parquet(up_path, columns=["accession", "panther"])
    pan = up.set_index("accession").panther.reindex(acc_u).fillna("")
    fam = pan.map(lambda s: sorted({x.split(":")[0] for x in s.split(";") if x}))
    sf = pan.map(lambda s: sorted({x for x in s.split(";") if ":SF" in x}))
    prot["panther_fam"] = fam.map(lambda v: v[0] if v else None).to_numpy()
    prot["panther_sf"] = sf.map(lambda v: v[0] if v else None).to_numpy()
    funnel["panther_fam_multi"] = int((fam.map(len) > 1).sum())
    funnel["panther_available"] = True
else:
    prot["panther_fam"] = None
    prot["panther_sf"] = None
    funnel["panther_available"] = False
    log("WARNING: uniprot_proteins.parquet missing -> no PANTHER labels/proxy (counted)")
funnel["prot_with_panther_fam"] = int(prot.panther_fam.notna().sum())

# ---------------------------------------------------------------- forest + groups
forest = hom.parse_orthoxml(str(Path(a.orthology) / "hogs.orthoxml.gz"))
g_of = pd.Series(np.arange(len(forest.gene_acc)), index=forest.gene_acc)
g_of = g_of[~g_of.index.duplicated()]
prot_gene = g_of.reindex(acc_u).fillna(-1).astype(np.int64).to_numpy().copy()
in_group = prot_gene >= 0
in_group[in_group] = forest.gene_parent[prot_gene[in_group]] >= 0
R = len(forest.root_names)
prot_group = np.full(P, -1, np.int64)
kind = np.array(["none"] * P, dtype=object)
prot_group[in_group] = forest.node_root_rank[forest.gene_parent[prot_gene[in_group]]]
kind[in_group] = "hog"
prot_gene[~in_group] = -1    # forest genes outside any group behave as singletons
funnel["prot_in_forest_groups"] = int(in_group.sum())
if a.idgap and Path(a.idgap).exists():
    # ID gap (O2): atlas accession -> the orthology member of the same gene (one member per model)
    ig = pd.read_parquet(a.idgap, columns=["atlas_accession", "hog_member_accession", "method"])
    ig = ig[ig.hog_member_accession.notna() & (ig.method != "direct")]
    mem = ig.set_index("atlas_accession").hog_member_accession.reindex(acc_u)
    gi_ = g_of.reindex(mem.values).fillna(-1).astype(np.int64).to_numpy()
    ok_ = (~in_group) & (gi_ >= 0)
    ok_[ok_] = forest.gene_parent[gi_[ok_]] >= 0
    # one atlas protein per member gene: keep the first, others -> not mapped (counted)
    dup = pd.Series(gi_).where(ok_).duplicated(keep="first").to_numpy() & ok_
    used = set(prot_gene[in_group].tolist())
    clash = ok_ & np.isin(gi_, list(used))
    ok_ &= ~dup & ~clash
    prot_gene[ok_] = gi_[ok_]
    prot_group[ok_] = forest.node_root_rank[forest.gene_parent[gi_[ok_]]]
    kind[ok_] = "idgap"
    in_group = in_group | ok_
    funnel["prot_idgap_mapped"] = int(ok_.sum())
    funnel["prot_idgap_dup_or_clash"] = int((dup | clash).sum())

placement = a.placement or next(iter(sorted(glob.glob(str(Path(a.orthology).parent / "placement*.parquet")))), None)
if placement:
    if Path(placement).is_dir():   # O2 layout: placement/tier2/part-*.parquet (accession, root_name, method, ...)
        pl = pd.concat([pd.read_parquet(f) for f in sorted(glob.glob(str(Path(placement) / "part-*.parquet")))])
        funnel["placement_methods"] = pl.method.value_counts().to_dict()
        pl = pl[pl.root_name.notna()].rename(columns={"root_name": "hog"}).drop_duplicates("accession")
    else:
        pl = pd.read_parquet(placement)
    rn = pd.Series(np.arange(R), index=forest.root_names)
    hog_root = pl.set_index("accession").hog.map(lambda h: rn.get(str(h).split(".")[0], -1))
    m = (prot_group < 0) & (hog_root.reindex(acc_u).fillna(-1).to_numpy() >= 0)
    prot_group[m] = hog_root.reindex(acc_u).to_numpy().astype(np.int64)[m]
    kind[m] = "placement"
    funnel["placement_table"] = placement
funnel["prot_placed"] = int((kind == "placement").sum())

# PANTHER-subfamily proxy: majority root HOG among forest proteins of the same subfamily
need = prot_group < 0
sfv = prot.panther_sf.to_numpy()
ref = pd.DataFrame({"sf": sfv[in_group], "g": prot_group[in_group]}).dropna()
maj = ref.groupby("sf").g.agg(lambda s: s.value_counts().index[0]) if len(ref) else pd.Series(dtype=np.int64)
m = need & pd.Series(sfv).map(maj).notna().to_numpy()
prot_group[m] = pd.Series(sfv[m]).map(maj).astype(np.int64).to_numpy()
kind[m] = "proxy_hog"
need = prot_group < 0
nxt = R
has_sf = need & pd.notna(sfv)
codes, uniq = pd.factorize(pd.Series(sfv[has_sf]))
prot_group[has_sf] = nxt + codes
kind[has_sf] = "proxy_sf"
nxt += len(uniq)
need = prot_group < 0
prot_group[need] = nxt + np.arange(need.sum())
kind[need] = "singleton"
prot["gene"], prot["group"], prot["group_kind"] = prot_gene, prot_group, kind
funnel["group_kind_proteins"] = prot.group_kind.value_counts().to_dict()
funnel["group_kind_sites"] = prot.groupby("group_kind").n_sites.sum().to_dict()
funnel["n_groups"] = int(prot_group.max() + 1)
log("groups", json.dumps(funnel["group_kind_proteins"]))

# ---------------------------------------------------------------- region labels (Pfam clan / TED CATH H)
min_ov = float(a.prereg["pairs"]["homology_classes"]["H3"]["region_overlap_min"]["value"])
pf = pd.read_parquet(ann / "pfam_regions.parquet", columns=["accession", "pfam", "clan", "seq_start", "seq_end"])
pf = pf[pf.accession.isin(acc_u)]
pf = pd.DataFrame({"accession": pf.accession.values, "rid": np.arange(len(pf)),
                   "label": np.where(pf.clan.notna(), pf.clan.fillna("").values, pf.pfam.values),
                   "start": pf.seq_start.values, "end": pf.seq_end.values})
ted = pd.read_parquet(ann / "ted_domains.parquet",
                      columns=["accession", "ted_id", "chopping", "cath_label", "cath_level", "af_fragment"])
ted = ted[ted.accession.isin(acc_u) & (ted.cath_level == "H")]
rows = []
for r in ted.itertuples(index=False):
    k = int(str(r.af_fragment).lstrip("F") or 1)
    off = 200 * (k - 1)
    for s, e in hom.parse_chopping(r.chopping):
        rows.append((r.accession, r.ted_id, r.cath_label, s + off, e + off))
ted_r = pd.DataFrame(rows, columns=["accession", "rid", "label", "start", "end"])
funnel["ted_H_domains"] = int(len(ted)); funnel["pfam_regions"] = int(len(pf))


def _chunk(args):
    accs, res, which = args
    regs = pf if which == "pfam" else ted_r
    regs = regs[regs.accession.isin(set(accs))]
    return hom.region_labels_for_sites(accs, res, regs, min_ov, "C:" if which == "pfam" else "T:")


def labels_parallel(which):
    accv = uni.accession.to_numpy().astype(str)
    resv = uni.residues.to_list()
    bounds = np.linspace(0, N, a.workers * 4 + 1).astype(int)
    # chunk boundaries on accession changes so no protein is split (sites of a protein are contiguous)
    jobs = [(accv[s:e], resv[s:e], which) for s, e in zip(bounds[:-1], bounds[1:])]
    with Pool(a.workers) as pool:
        outs = pool.map(_chunk, jobs)
    labs, hit = [], 0
    for l, h in outs:
        labs += l
        hit += h
    return labs, hit


lab_pf, hit_pf = labels_parallel("pfam")
lab_ted, hit_ted = labels_parallel("ted")
funnel["sites_with_pfam_region_label"] = hit_pf
funnel["sites_with_ted_H_region_label"] = hit_ted
log(f"region labels: pfam {hit_pf}, ted {hit_ted} of {N} sites ({time.time() - t0:.0f}s)")

vocab: dict = {}
fam_lab = [[f"F:{x}"] if x else [] for x in prot.panther_fam.to_numpy()]
site_lists = [lp + lt + fam_lab[p] for lp, lt, p in zip(lab_pf, lab_ted, site_prot)]
site_lab, over_s = hom.encode_label_slots(site_lists, vocab, hom.K_SITE + 1)
# protein-level labels: all Pfam clans/families, all TED H labels, PANTHER family
pl_pf = pf.groupby("accession").label.agg(lambda s: sorted(set("C:" + s)))
pl_te = ted.groupby("accession").cath_label.agg(lambda s: sorted(set("T:" + s)))
prot_lists = [list(pl_pf.get(x, [])) + list(pl_te.get(x, [])) + fam_lab[k] for k, x in enumerate(acc_u)]
prot_lab, over_p = hom.encode_label_slots(prot_lists, vocab, hom.K_PROT)
funnel["label_overflow_site"] = over_s; funnel["label_overflow_prot"] = over_p
funnel["n_labels"] = len(vocab)
funnel["ted_label_missing_sites"] = int(N - hit_ted)

# ---------------------------------------------------------------- covariates
def load_sequences(seq_path, up_path, accs, acc_table):
    """atlas accession -> UniProt sequence (via uniprot_acc); ok only if its length equals the model length."""
    sq = pd.read_parquet(seq_path)
    upm = pd.read_parquet(up_path, columns=["accession", "uniprot_acc"])
    m = upm.merge(sq, on="uniprot_acc", how="left").drop_duplicates("accession").set_index("accession")
    seqs = m.sequence.reindex(accs).to_numpy()
    nres = pd.read_parquet(acc_table, columns=["accession", "n_res"]).drop_duplicates("accession")
    nres = nres.set_index("accession").n_res.reindex(accs).to_numpy()
    lens = np.array([len(x) if isinstance(x, str) else -1 for x in seqs])
    # fragmented human models (> 2,700 aa) have a 1,400-residue F1; their numbering is UniProt (SITES.md)
    return seqs, (lens == nres) | ((nres == 1400) & (lens > 2700))

seq_ok = Path(a.sequences).exists() and up_path.exists()
if seq_ok:
    seqs, lens_ok = load_sequences(a.sequences, up_path, acc_u, ann / "accessions.parquet")
    funnel["prot_seq_missing"] = int(sum(1 for x in seqs if not isinstance(x, str)))
    funnel["prot_seq_len_mismatch_model"] = int((~lens_ok).sum() - funnel["prot_seq_missing"])
    seqs = np.array([x if ok_ else None for x, ok_ in zip(seqs, lens_ok)], dtype=object)
else:
    seqs = np.array([None] * P)
    log("WARNING: sequences missing -> polarity/aromatic NaN")
pol = np.full(N, np.nan); aro = np.full(N, np.nan)
res_list = uni.residues.to_list()
for k in range(N):
    s = seqs[site_prot[k]]
    if isinstance(s, str):
        pol[k], aro[k], ok = nulls.lining_fractions(res_list[k], s)
cov = np.stack([np.log(uni.rep_site_volume_A3.to_numpy()), uni.max_burrow_score.to_numpy(), pol, aro,
                uni.plddt_mean_rep.to_numpy() / 100.0, (uni.tier.to_numpy() == 2).astype(float)], 1).astype(np.float32)
cov_ok = np.isfinite(cov).all(1)
funnel["sites_cov_ok"] = int(cov_ok.sum())
funnel["sites_cov_missing_seq_or_range"] = int((~cov_ok).sum())
np.save(a.work / "covariates.npy", cov)

np.savez(a.work / "context.npz", site_prot=site_prot, site_lab=site_lab, prot_group=prot_group,
         prot_gene=prot_gene, prot_lab=prot_lab, cov_ok=cov_ok)
prot.to_parquet(a.work / "proteins.parquet", index=False)
(a.work / "labels_vocab.json").write_text(json.dumps({v: k for k, v in vocab.items()}))
funnel["wall_s"] = time.time() - t0
pairs.write_sidecar(a.work / "context.npz", a.snapshot, None, {"stage": "context", "orthology": a.orthology,
                    "funnel": funnel})
(a.results / "context_funnel.json").write_text(json.dumps(funnel, indent=1, default=str))
log("context done", json.dumps(funnel, default=str))
