#!/usr/bin/env python
"""Part A stage 6: homology classes for significant pairs and for the denominator samples.

  a. join/forest classes (H0/H1/H2/H12/H3j/rest) for every kept pair (max E < keep_E, any tail variant)
  b. same-site split of H1/H2/H12 (round-2 amendment 1): SOG membership if --sogs is given, else the declared
     PLACEHOLDER: MMseqs2 local alignment of the two proteins, site residues mapped through it, mapped residue
     Jaccard >= prereg tau -> same site ('s'), else offsite ('x'). Denominator split from a seeded sample of
     within-group site pairs (same rule).
  c. 'rest' pairs: MMseqs2 E < 1e-3 -> H3s; else Foldseek E < 1e-3 or max(qTM, tTM) >= 0.5 -> H4; else N; no
     sequence/structure -> unknown. The same searches label the seeded 10^5 'rest' sample (denominators).
Outputs: <work>/cls_<score>/{pairs_cls.parquet, sample_rest.parquet, sample_within.parquet, summary.json}
"""
import json
import os
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

from _common import parser, resolve, logger, pairs, build_seqdb
from pocketevo import homology as hom

ap = parser(__doc__)
ap.add_argument("--score", choices=["cos", "csls"], default="cos")
ap.add_argument("--sequences", default="/tmp/pe_big/work/annotations/uniprot_sequences.parquet")
ap.add_argument("--foldseek-db", default="/tmp/pe_big/partA/foldseek/afdb_v6")
ap.add_argument("--mmseqs", default=str(pairs.REPO / "envs/ancestra/bin/mmseqs"))
ap.add_argument("--foldseek", default=str(pairs.REPO / "envs/ancestra/bin/foldseek"))
ap.add_argument("--sogs", default=None, help="site_sog_all_tau.parquet (site_id, og_id, sog_t<tau>); placeholder "
                "mapped-Jaccard rule for pairs where either site has no SOG")
ap.add_argument("--tau-choice", default=None, help="tau_choice.json of the SOG run (sets tau)")
ap.add_argument("--tau", type=float, default=None, help="explicit tau (overrides --tau-choice; sensitivity runs)")
ap.add_argument("--n-rest-sample", type=int, default=100_000)
ap.add_argument("--n-within-sample", type=int, default=200_000)
ap.add_argument("--threads", type=int, default=32)
ap.add_argument("--search-E", type=float, default=None, help="only search pairs with max E below this")
ap.add_argument("--n-secondary", type=int, default=300_000,
                help="seeded sample size of kept 'rest' pairs outside the primary set (pooled E < 0.01)")
a = resolve(ap.parse_args())
log = logger(a.work / "partA.log")
T0 = time.time()
out = a.work / f"cls_{a.score}{a.tag}"
out.mkdir(exist_ok=True)
pre = a.prereg["pairs"]["homology_classes"]
EV_SEQ, EV_STR, TM_MIN = 1e-3, 1e-3, 0.5          # prereg H3 / H4 rules
TAU = float(a.prereg["sogs"]["tau_column_jaccard"]["value"])
if a.tau_choice and Path(a.tau_choice).exists():     # tau tuned on the tune half by the SOG builder (prereg: tunable)
    TAU = float(json.loads(Path(a.tau_choice).read_text())["tau"])
if a.tau is not None:
    TAU = a.tau
summ = {"score": a.score, "tau": TAU}

uni = pd.read_parquet(a.work / "universe.parquet", columns=["site_id", "accession", "residues"])
RES = uni.residues.to_numpy()
prot = pd.read_parquet(a.work / "proteins.parquet")
ACC = prot.accession.to_numpy()
c = np.load(a.work / "context.npz")
forest = hom.parse_orthoxml(str(Path(a.orthology) / "hogs.orthoxml.gz"))
ctx = hom.PairContext(c["site_prot"], c["site_lab"], c["prot_group"], c["prot_gene"], forest)
site_group = ctx.prot_group[ctx.site_prot]
N = len(ctx.site_prot)

# ---------------------------------------------------------------- a. join classes
sig = pd.read_parquet(a.work / f"sig_{a.score}" / "pairs.parquet")
if a.search_E is not None:
    sig = sig[sig[[c_ for c_ in sig.columns if c_.startswith("Emax_")]].min(1) < a.search_E].reset_index(drop=True)
i, j = sig.i.to_numpy().astype(np.int64), sig.j.to_numpy().astype(np.int64)
cls0 = np.empty(len(sig), np.int8)
for s in range(0, len(sig), 5_000_000):
    cls0[s:s + 5_000_000] = hom.classify_pairs(ctx, i[s:s + 5_000_000], j[s:s + 5_000_000])
sig["cls0"] = cls0
sig["gi"], sig["gj"] = site_group[i], site_group[j]
summ["n_pairs"] = len(sig)
summ["cls0_counts"] = {hom.CLASS_NAMES[k]: int(v) for k, v in zip(*np.unique(cls0, return_counts=True))}
log("join classes", json.dumps(summ["cls0_counts"]), f"{time.time() - T0:.0f}s")

# ---------------------------------------------------------------- sequences / DBs
seqdb = a.work / "mmseqs" / "seqdb"
summ.update(build_seqdb(a, ACC, log))   # shared builder (_common.build_seqdb)


def read_lookup(db, strip_afdb=False):
    lk = pd.read_csv(str(db) + ".lookup", sep="\t", header=None, names=["key", "name", "file"])
    names = lk.name.astype(str)
    if strip_afdb:   # AF-P12345-F1-model_v6 -> P12345
        names = names.str.replace(r"^AF-", "", regex=True).str.replace(r"-F1-model_v\d+.*$", "", regex=True)
    return dict(zip(names, lk.key)), dict(zip(lk.key, lk.name.astype(str)))


seq_lookup, _ = read_lookup(seqdb)
fs_lookup_acc, fs_key2name = read_lookup(a.foldseek_db, strip_afdb=True)
fs_lookup = {acc_: fs_key2name[k] for acc_, k in fs_lookup_acc.items()}   # acc -> entry name in foldseek db
fs_name2key = {v: k for k, v in fs_key2name.items()}
summ["seqdb_entries"] = len(seq_lookup)
summ["foldseek_entries"] = len(fs_lookup_acc)
summ["proteins_missing_seq"] = int(sum(1 for x in ACC if x not in seq_lookup))
summ["proteins_missing_structure"] = int(sum(1 for x in ACC if x not in fs_lookup_acc))


def search_mmseqs(pp: pd.DataFrame, tag: str) -> pd.DataFrame:
    """pp: unique (qacc, tacc). Returns evalue + alignment coordinates/cigar."""
    q = pd.DataFrame({"qname": pp.qacc.values, "tname": pp.tacc.values})
    return hom.run_pair_search(q, str(seqdb), seq_lookup, a.work / "mmseqs" / tag, "mmseqs", a.mmseqs,
                               threads=a.threads, evalue=10.0, log=log)


def search_foldseek(pp: pd.DataFrame, tag: str) -> pd.DataFrame:
    q = pd.DataFrame({"qname": pp.qacc.map(fs_lookup).values, "tname": pp.tacc.map(fs_lookup).values})
    q = q.fillna("__missing__")
    r = hom.run_pair_search(q, a.foldseek_db, fs_name2key, a.work / "foldseek_aln" / tag, "foldseek", a.foldseek,
                            threads=a.threads, evalue=10.0, log=log)
    r["qacc"], r["tacc"] = pp.qacc.values, pp.tacc.values
    return r


def resolve_rest(pp: pd.DataFrame, tag: str) -> pd.DataFrame:
    """pp unique protein pairs (qacc, tacc) -> cls in {H3s, H4, N, unknown} + search fields."""
    pp = pp.reset_index(drop=True)
    ms = search_mmseqs(pp, tag)
    h3 = (ms.status == "hit") & (ms.evalue < EV_SEQ)
    pp["mm_evalue"] = ms.evalue.values
    pp["mm_status"] = ms.status.values
    todo = pp[~h3.values].copy()
    fs = search_foldseek(todo[["qacc", "tacc"]], tag) if len(todo) else pd.DataFrame(columns=["status"])
    todo["fs_evalue"] = fs.evalue.values if len(todo) else []
    todo["fs_tm"] = np.fmax(fs.qtmscore.values, fs.ttmscore.values) if len(todo) else []
    todo["fs_status"] = fs.status.values if len(todo) else []
    pp = pp.merge(todo[["qacc", "tacc", "fs_evalue", "fs_tm", "fs_status"]], on=["qacc", "tacc"], how="left")
    h4 = (pp.fs_status == "hit") & ((pp.fs_evalue < EV_STR) | (pp.fs_tm >= TM_MIN))
    unknown = (pp.mm_status == "missing") | (pp.fs_status == "missing")
    cls = np.where(h3.values, "H3s", np.where(h4, "H4", np.where(unknown, "unknown", "N")))
    pp["cls"] = cls
    return pp


def prot_pairs(ii, jj):
    pa, pb = ACC[ctx.site_prot[ii]], ACC[ctx.site_prot[jj]]
    q = np.where(pa < pb, pa, pb)
    t = np.where(pa < pb, pb, pa)
    return q, t


# ---------------------------------------------------------------- c. rest: significant pairs + 10^5 sample
rest = (sig.cls0 == hom.REST).to_numpy()
primary = (sig.Emax_pooled < 0.01).to_numpy()
rng_s = np.random.default_rng(a.seed + 30)
sec = rest & ~primary
p_sec = min(1.0, a.n_secondary / max(int(sec.sum()), 1))
sampled_sec = sec & (rng_s.random(len(sig)) < p_sec)
searched = rest & (primary | sampled_sec)
sig["search_weight"] = 1.0
sig.loc[sec, "search_weight"] = 0.0
sig.loc[sampled_sec, "search_weight"] = 1.0 / p_sec
q, t = prot_pairs(i[searched], j[searched])
up = pd.DataFrame({"qacc": q, "tacc": t}).drop_duplicates()
summ["rest_sig_pairs_keptE1"] = int(rest.sum())
summ["rest_pairs_primary_exhaustive"] = int((rest & primary).sum())
summ["rest_pairs_secondary"] = int(sec.sum())
summ["rest_pairs_secondary_sampled"] = int(sampled_sec.sum())
summ["secondary_sampling_prob"] = p_sec
summ["rest_searched_protein_pairs"] = len(up)
log(f"rest: {int(rest.sum())} site pairs; searched {int(searched.sum())} ({len(up)} protein pairs)")

# seeded uniform sample of 'rest' site pairs (cross-group, no shared join label)
rng = np.random.default_rng(a.seed)
si_, sj_ = [], []
drawn = 0
while sum(map(len, si_)) < a.n_rest_sample:
    m = 400_000
    x, y = rng.integers(0, N, m), rng.integers(0, N, m)
    drawn += m
    cl = hom.classify_pairs(ctx, x, y)
    ok = cl == hom.REST
    si_.append(x[ok]); sj_.append(y[ok])
smp = pd.DataFrame({"i": np.concatenate(si_)[:a.n_rest_sample], "j": np.concatenate(sj_)[:a.n_rest_sample]})
summ["rest_sample_drawn"] = drawn
summ["rest_sample_accept"] = float(sum(map(len, si_)) / drawn)   # acceptance over all draws (v1 bug: n/drawn)
sq_, st_ = prot_pairs(smp.i.values, smp.j.values)
smp["qacc"], smp["tacc"] = sq_, st_
allpp = pd.concat([up, smp[["qacc", "tacc"]]]).drop_duplicates().reset_index(drop=True)
res_path = out / "rest_protein_pairs.parquet"
if res_path.exists():
    rp = pd.read_parquet(res_path)
else:
    rp = resolve_rest(allpp, "rest")
    rp.to_parquet(res_path, index=False)
key = dict(zip(zip(rp.qacc, rp.tacc), rp.cls))
sig["cls"] = pd.Series(sig.cls0.map(hom.CLASS_NAMES)).astype(object).values
cls_rest = [key[(x, y)] for x, y in zip(q, t)]
tmp = sig.cls.to_numpy()
tmp[rest] = "rest_unsearched"
tmp[searched] = cls_rest
sig["cls"] = tmp
smp["cls"] = [key[(x, y)] for x, y in zip(smp.qacc, smp.tacc)]
smp["gi"], smp["gj"] = site_group[smp.i.values], site_group[smp.j.values]
smp.to_parquet(out / "sample_rest.parquet", index=False)
summ["rest_sample_classes"] = smp.cls.value_counts().to_dict()
log("rest sample classes", json.dumps(summ["rest_sample_classes"]), f"{time.time() - T0:.0f}s")


# ---------------------------------------------------------------- b. same-site split (H1/H2/H12)
def within_sample(n, seed):
    """Uniform sample of within-group different-protein site pairs."""
    r = np.random.default_rng(seed)
    order = np.argsort(site_group, kind="stable")
    g_sorted = site_group[order]
    starts = np.searchsorted(g_sorted, np.arange(site_group.max() + 1))
    sizes = np.bincount(site_group, minlength=site_group.max() + 1).astype(np.float64)
    w = sizes * (sizes - 1)
    w = w / w.sum()
    xi, xj = [], []
    while sum(map(len, xi)) < n:
        g = r.choice(len(w), 200_000, p=w)
        a_ = order[starts[g] + (r.random(200_000) * sizes[g]).astype(np.int64)]
        b_ = order[starts[g] + (r.random(200_000) * sizes[g]).astype(np.int64)]
        ok = ctx.site_prot[a_] != ctx.site_prot[b_]
        xi.append(a_[ok]); xj.append(b_[ok])
    return np.concatenate(xi)[:n], np.concatenate(xj)[:n]


def _jac_chunk(rows):
    """rows: (qstart, tstart, cigar, [(site_q, site_t), ...]) -> list of (site_q, site_t, jaccard)."""
    out_ = []
    for qs_, ts_, cg, sp_ in rows:
        if not isinstance(cg, str):
            out_ += [(r, np.nan) for x, y, r in sp_]
            continue
        qpos, tpos = hom.cigar_map(cg, int(qs_), int(ts_))
        for x, y, r in sp_:
            out_.append((r, hom.mapped_jaccard(RES[x], RES[y], qpos, tpos)))
    return out_


def same_site(ii, jj, tag):
    """placeholder same-site call for site pairs (ii, jj) in the same root group -> (jaccard, same)."""
    pa, pb = ACC[ctx.site_prot[ii]], ACC[ctx.site_prot[jj]]
    swap = pa > pb
    qi, qj = np.where(swap, jj, ii), np.where(swap, ii, jj)
    qa, qb = np.where(swap, pb, pa), np.where(swap, pa, pb)
    sp = pd.DataFrame({"qacc": qa, "tacc": qb, "qi": qi, "qj": qj, "row": np.arange(len(ii))})
    upp = sp[["qacc", "tacc"]].drop_duplicates().reset_index(drop=True)
    ms = search_mmseqs(upp, tag)
    upp["qstart"], upp["tstart"], upp["cigar"] = ms.qstart.values, ms.tstart.values, ms.cigar.values
    sp = sp.merge(upp, on=["qacc", "tacc"], how="left")
    sp = sp.sort_values(["qacc", "tacc"], kind="stable").reset_index(drop=True)
    key_ = sp.qacc.values + "|" + sp.tacc.values
    starts = np.flatnonzero(np.r_[True, key_[1:] != key_[:-1]])
    ends = np.r_[starts[1:], len(sp)]
    QS, TS, CG = sp.qstart.values, sp.tstart.values, sp.cigar.values
    QI, QJ, RW = sp.qi.values, sp.qj.values, sp.row.values
    rows = [(QS[b], TS[b], CG[b], list(zip(QI[b:e], QJ[b:e], RW[b:e]))) for b, e in zip(starts, ends)]
    chunks = [rows[k:k + 20_000] for k in range(0, len(rows), 20_000)]
    with Pool(a.threads) as pool:
        parts = pool.map(_jac_chunk, chunks)
    jac = np.full(len(ii), np.nan)
    for p_ in parts:
        for r, v in p_:
            jac[r] = v
    return jac, np.nan_to_num(jac, nan=0.0) >= TAU


use_sog = bool(a.sogs and Path(a.sogs).exists())
split_mode = "sog+placeholder_fallback" if use_sog else "placeholder_mmseqs_mapped_jaccard"
summ["split_mode"] = split_mode
homcls = sig.cls.isin(["H1", "H2", "H12"]).to_numpy()
if use_sog:
    sog = pd.read_parquet(a.sogs, columns=["site_id", "og_id", f"sog_t{TAU:g}"])
    sid = (sog.og_id + ":g" + sog[f"sog_t{TAU:g}"].astype(str)).values
    sog_of = pd.Series(sid, index=sog.site_id.values).reindex(uni.site_id.values).to_numpy()
    has_sog = pd.notna(sog_of)
else:
    sog_of = np.full(N, None, dtype=object)
    has_sog = np.zeros(N, bool)


def split_rule(ii, jj, tag):
    """same-site call: SOG identity when both sites have a SOG, else the mapped-Jaccard placeholder."""
    both = has_sog[ii] & has_sog[jj]
    same = np.zeros(len(ii), bool)
    jac = np.full(len(ii), np.nan)
    same[both] = sog_of[ii[both]] == sog_of[jj[both]]
    if (~both).any():
        jj_, ss_ = same_site(ii[~both], jj[~both], tag)
        jac[~both] = jj_
        same[~both] = ss_
    src = np.where(both, "sog", "placeholder")
    return jac, same, src


jac, same, src = split_rule(i[homcls], j[homcls], "within_sig")
sig["site_jaccard"] = np.nan
sig.loc[homcls, "site_jaccard"] = jac
sig["same_site"] = False
sig.loc[homcls, "same_site"] = same
sig["split_source"] = None
sig.loc[homcls, "split_source"] = src
summ["split_source_sig"] = pd.Series(src).value_counts().to_dict()
newc = sig.cls.to_numpy().astype(object)
m = homcls
newc[m] = np.where(sig.same_site.to_numpy()[m].astype(bool), sig.cls.to_numpy()[m] + "s", sig.cls.to_numpy()[m] + "x")
sig["cls"] = newc

# denominator split sample
wi, wj = within_sample(a.n_within_sample, a.seed + 10)
wcls = hom.classify_pairs(ctx, wi, wj)
ws = pd.DataFrame({"i": wi, "j": wj, "cls0": [hom.CLASS_NAMES[x] for x in wcls]})
wjac, wsame, wsrc = split_rule(wi, wj, "within_sample")
ws["site_jaccard"], ws["same_site"], ws["split_source"] = wjac, wsame, wsrc
summ["split_source_sample"] = ws.groupby(["cls0", "split_source"]).size().astype(int).to_dict().__repr__()
ws["gi"] = site_group[wi]
ws.to_parquet(out / "sample_within.parquet", index=False)
split = ws.groupby("cls0").same_site.mean().to_dict()
summ["split_fraction"] = split
summ["split_sample_n"] = ws.cls0.value_counts().to_dict()
log("split", json.dumps(split))

sig.to_parquet(out / "pairs_cls.parquet", index=False)
summ["cls_counts_keptE1"] = sig.cls.value_counts().to_dict()
summ["cls_counts_E0.01_pooled"] = sig[sig.Emax_pooled < 0.01].cls.value_counts().to_dict()
summ["wall_s"] = time.time() - T0
(out / "summary.json").write_text(json.dumps(summ, indent=1, default=str))
(a.results / f"cls_{a.score}{a.tag}_summary.json").write_text(json.dumps(summ, indent=1, default=str))
pairs.write_sidecar(out / "pairs_cls.parquet", a.snapshot, a.seed, {"stage": f"cls_{a.score}{a.tag}", **summ})
log("classify done", json.dumps(summ, default=str))
