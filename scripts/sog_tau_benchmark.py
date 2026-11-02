#!/usr/bin/env python
"""DESIGN.md §8.1 site-matching benchmark: tune tau (prereg sogs.tau_column_jaccard) on the TUNE half only.

Annotated sites (experimental only):
  * UniProt `Binding site` features with ECO:0000269 / ECO:0007744 (column `experimental`), grouped per
    (accession, ChEBI, ligand_label) -> one annotated site = union of the feature positions; ligand key = ChEBI.
  * BioLiP2 rows (residues mapped to UniProt numbering by scripts/annotations/biolip.py), ligand key = CCD id;
    per (accession, ligand) the PDB entries are clustered (residue overlap coefficient >= 0.5, union-find) and a
    cluster's site = residues seen in >= half of its entries (union if that is empty).
  Ligand classes peptide/dna/rna are excluded (not pockets). Metal ions are excluded from the primary benchmark
  (prereg: metal-only sites are U; BURROW drops small metal sites) and reported as a sensitivity row.
Detection (prereg benchmarks.detectability_preflight): an annotated site is detected if an atlas site of the same
  protein covers >= overlap_min (0.5) of its residues; the best-covering site is used.
Pairs (unit = pair of distinct atlas sites in members of different species of one root OG):
  positive = the two atlas sites carry annotations of the same ligand from the same source;
  negative = annotations of different ligands whose annotated residues map to disjoint MSA column sets (and the
             site pair is not also a positive pair).
  recovered = both atlas sites in the same SOG (Leiden membership); false merge = negative pair in the same SOG.
  recall = TP / positives (conditional on detection; unconditional recall also reported),
  precision = TP / (TP + FP), F1 = 2PR / (P + R). 95 % CIs: bootstrap over root OGs (200 replicates,
  seed rule prereg meta.bootstrap_seed_rule).
tau choice: max F1 on the TUNE half (both sources pooled, metals excluded); ties -> larger precision -> smaller tau.

    .venv/bin/python scripts/sog_tau_benchmark.py --sogs /tmp/pe_big/sogs/20261030T085136Z/O1 \\
        --annotations data/annotations --prereg config/prereg.yaml --results results/sogs/20261030T085136Z/O1
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

P = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(P / "src"))
sys.path.insert(0, str(P / "scripts"))
from pocketevo import sog  # noqa: E402
from build_sogs import git_commit, code_sha256  # noqa: E402

METAL_RE = re.compile(r"^(?:[A-Z][a-z]?)(?:\(\d?[+-]\))$")   # UniProt ligand names like Zn(2+), Cl(-)


def uniprot_sites(ann: Path) -> pd.DataFrame:
    f, fp = ann / "uniprot_site_features.parquet", ann / "uniprot_proteins.parquet"
    if not f.exists() or not fp.exists():
        return pd.DataFrame()
    S = pd.read_parquet(f)
    prot = pd.read_parquet(fp, columns=["accession", "uniprot_acc"])
    S = S[(S.type == "Binding site") & S.experimental.astype(bool)]
    S = S[S.start.notna() & S.end.notna()]
    for c in ["start_modifier", "end_modifier"]:
        if c in S:
            S = S[S[c].isna() | (S[c] == "EXACT")]
    S = S.merge(prot, on="uniprot_acc")
    S = S[S.chebi.notna()]   # prereg: "same ligand ChEBI" (generic names such as "substrate" carry no ChEBI)
    S["ligand_label"] = S.ligand_label.fillna("")
    S["res"] = [list(range(int(a), int(b) + 1)) for a, b in zip(S.start, S.end)]
    g = S.groupby(["accession", "chebi", "ligand_label"]).agg(residues=("res", lambda x: sorted(set().union(*map(set, x)))),
                                                             ligand=("ligand", "first")).reset_index()
    g["source"] = "uniprot"
    g["ligand_key"] = g.chebi
    g["metal"] = g.ligand.fillna("").map(lambda x: bool(METAL_RE.match(x)) or x.endswith(" cation"))
    g["ann_id"] = "U|" + g.accession + "|" + g.chebi + "|" + g.ligand_label
    return g[["ann_id", "source", "accession", "ligand_key", "ligand", "metal", "residues"]]


def _uf_clusters(sets: list, min_overlap: float) -> list:
    n = len(sets)
    par = list(range(n))

    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x
    for i in range(n):
        for j in range(i + 1, n):
            a, b = sets[i], sets[j]
            if a and b and len(a & b) / min(len(a), len(b)) >= min_overlap:
                par[find(i)] = find(j)
    out = {}
    for i in range(n):
        out.setdefault(find(i), []).append(i)
    return list(out.values())


def biolip_sites(ann: Path, overlap: float) -> pd.DataFrame:
    f = ann / "biolip_sites.parquet"
    if not f.exists():
        return pd.DataFrame()
    B = pd.read_parquet(f, columns=["accession", "pdb_id", "chain", "ligand", "ligand_class", "ligand_name", "residues_uniprot"])
    B = B[B.ligand_class.isin(["regular", "metal"])]
    B = B[B.residues_uniprot.map(len) > 0]
    rows = []
    for (acc, lig), g in B.groupby(["accession", "ligand"], sort=True):
        sets = [frozenset(int(x) for x in r) for r in g.residues_uniprot]
        for k, cl in enumerate(sorted(_uf_clusters(sets, overlap), key=lambda c: min(c))):
            cnt = {}
            for i in cl:
                for r in sets[i]:
                    cnt[r] = cnt.get(r, 0) + 1
            cons = sorted(r for r, c in cnt.items() if c >= len(cl) / 2) or sorted(cnt)
            rows.append({"ann_id": f"B|{acc}|{lig}|{k}", "source": "biolip", "accession": acc, "ligand_key": lig,
                         "ligand": g.ligand_name.iloc[0], "metal": g.ligand_class.iloc[0] == "metal", "residues": cons,
                         "n_entries": len(cl)})
    return pd.DataFrame(rows)


def detect(annots: pd.DataFrame, sites: pd.DataFrame, overlap_min: float) -> pd.DataFrame:
    by_acc = {a: g for a, g in sites.groupby("accession")}
    best_site, best_ov = [], []
    for acc, res in zip(annots.accession, annots.residues):
        g = by_acc.get(acc)
        bs, bo = None, 0.0
        if g is not None and len(res):
            rs = set(res)
            for sid, r in zip(g.site_id, g.residues):
                ov = len(rs & set(r.tolist() if hasattr(r, "tolist") else r)) / len(rs)
                if ov > bo + 1e-12:
                    bs, bo = sid, ov
        best_site.append(bs)
        best_ov.append(bo)
    out = annots.copy()
    out["site_id"] = best_site
    out["coverage"] = best_ov
    out["detected"] = out.coverage >= overlap_min - 1e-12
    return out


def make_pairs(A: pd.DataFrame, msa: dict, taxid_of: dict) -> pd.DataFrame:
    """A: annotations of one OG with og_id, site_id, detected. msa: {acc: (aligned, pos2col)}."""
    A = A.reset_index(drop=True)
    cols = []
    for acc, res in zip(A.accession, A.residues):
        if acc in msa:
            p2c = msa[acc]
            r = np.asarray(res, dtype=np.int64)
            r = r[(r >= 1) & (r <= len(p2c))]
            cols.append(frozenset(p2c[r - 1].tolist()))
        else:
            cols.append(None)
    A["cols"] = cols
    A["taxid"] = A.accession.map(taxid_of)
    pos, neg = {}, {}
    recs = [r[:3] + ((r[3] if isinstance(r[3], str) else None),) + r[4:]
            for r in A[["taxid", "accession", "detected", "site_id", "ann_id", "source", "ligand_key", "metal", "cols"]]
            .itertuples(index=False, name=None)]
    n = len(recs)
    for i in range(n):
        ti, acci, deti, sidi, aidi, srci, ligi, mi, ci = recs[i]
        for j in range(i + 1, n):
            tj, accj, detj, sidj, aidj, srcj, ligj, mj, cj = recs[j]
            if ti == tj or acci == accj:
                continue
            det = bool(deti and detj)
            key = tuple(sorted([sidi or f"undetected:{aidi}", sidj or f"undetected:{aidj}"]))
            if srci == srcj and ligi == ligj:
                d = pos.setdefault(key, {"detected": False, "sources": set(), "metal": True, "ligands": set()})
                d["sources"].add(srci)
                d["metal"] &= bool(mi and mj)
                d["ligands"].add(ligi)
                d["detected"] |= det
            elif ligi != ligj and det and sidi != sidj and ci is not None and cj is not None and not (ci & cj):
                d = neg.setdefault(key, {"detected": True, "sources": set(), "metal": True, "ligands": set()})
                d["sources"].add(srci + "/" + srcj)
                d["metal"] &= bool(mi or mj)
                d["ligands"].add(ligi + "|" + ligj)
    rows = []
    for kind, dd in (("pos", pos), ("neg", neg)):
        for (a, b), d in dd.items():
            if kind == "neg" and (a, b) in pos:
                continue
            rows.append({"kind": kind, "site_a": a, "site_b": b, "detected": d["detected"],
                         "sources": ";".join(sorted(d["sources"])), "metal_only": d["metal"],
                         "ligands": ";".join(sorted(d["ligands"]))[:200]})
    return pd.DataFrame(rows)


def metrics(df: pd.DataFrame, tau: float) -> dict:
    c = f"same_t{tau:g}"
    pos = df[df.kind == "pos"]
    posd = pos[pos.detected]
    neg = df[df.kind == "neg"]
    tp = int(posd[c].sum())
    fp = int(neg[c].sum())
    rec = tp / len(posd) if len(posd) else float("nan")
    prec = tp / (tp + fp) if tp + fp else float("nan")
    f1 = 2 * prec * rec / (prec + rec) if (prec + rec) and not np.isnan(prec + rec) else float("nan")
    return {"tau": tau, "n_pos": int(len(pos)), "n_pos_detected": int(len(posd)), "n_neg": int(len(neg)),
            "tp": tp, "fp": fp, "recall": rec, "precision": prec, "f1": f1,
            "recall_unconditional": tp / len(pos) if len(pos) else float("nan"),
            "false_merge_rate": fp / len(neg) if len(neg) else float("nan"),
            "n_ogs": int(df.og_id.nunique())}


def boot_ci(df, tau, n_rep, seed):
    ogs = np.array(sorted(df.og_id.unique()))
    if len(ogs) < 2:
        return (float("nan"),) * 2
    g = {o: x for o, x in df.groupby("og_id")}
    vals = []
    for r in range(n_rep):
        rng = np.random.default_rng(seed * 1000 + r)
        pick = rng.choice(ogs, len(ogs), replace=True)
        vals.append(metrics(pd.concat([g[o] for o in pick]), tau)["f1"])
    v = np.array([x for x in vals if not np.isnan(x)])
    return (float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))) if len(v) else (float("nan"),) * 2


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sogs", required=True, help="/tmp/pe_big/sogs/<snapshot>/<orthology_version>")
    ap.add_argument("--annotations", default=str(P / "data/annotations"))
    ap.add_argument("--prereg", default=str(P / "config/prereg.yaml"))
    ap.add_argument("--results", required=True)
    ap.add_argument("--n-boot", type=int, default=200)
    ap.add_argument("--no-choice", action="store_true", help="report only (e.g. functional slice), no tau_choice.json")
    ap.add_argument("--label", default="")
    a = ap.parse_args()
    t0 = time.time()
    prereg = yaml.safe_load(open(a.prereg))
    params = sog.prereg_params(prereg)
    psha = sog.prereg_sha256(a.prereg)
    seed = int(params["seeds"][0])
    root = Path(a.sogs)
    ann = Path(a.annotations)
    taus = params["tau_grid"]
    members = pd.read_parquet(root / "inputs/members.parquet", columns=["og_id", "accession", "taxid", "shard"])
    ogs = pd.read_parquet(root / "inputs/ogs.parquet")
    ss = pd.concat([pd.read_parquet(f) for f in sorted(glob.glob(str(root / "sogs/sites_shard=*.parquet")))], ignore_index=True)
    ran = set(ss.og_id)
    stats = pd.concat([pd.read_parquet(f) for f in sorted(glob.glob(str(root / "sogs/ogstats_shard=*.parquet")))], ignore_index=True)
    U = uniprot_sites(ann)
    B = biolip_sites(ann, params["detect_overlap_min"])
    have = {"uniprot": len(U) > 0, "biolip": len(B) > 0}
    annots = pd.concat([x for x in (U, B) if len(x)], ignore_index=True) if any(have.values()) else pd.DataFrame()
    if annots.empty:
        raise SystemExit("no annotations available in " + str(ann))
    annots = annots.merge(members[["og_id", "accession", "shard"]], on="accession")
    # annotations are in current-UniProt numbering: keep only proteins whose model sequence equals it
    mcrc = pd.concat([pd.read_parquet(f, columns=["accession", "crc64"])
                      for f in sorted(glob.glob(str(root / "inputs/shard_models/*/*.parquet")))]).drop_duplicates("accession")
    up = pd.read_parquet(ann / "uniprot_proteins.parquet", columns=["accession", "seq_crc64"]) \
        if (ann / "uniprot_proteins.parquet").exists() else pd.DataFrame(columns=["accession", "seq_crc64"])
    chk = mcrc.merge(up, on="accession", how="left")
    seq_ok = set(chk[chk.seq_crc64.isna() | (chk.seq_crc64 == chk.crc64)].accession)
    n_before = len(annots)
    annots = annots[annots.accession.isin(seq_ok)]
    n_seq_dropped = n_before - len(annots)
    annots = annots[annots.og_id.isin(ran)]
    sites = pd.concat([pd.read_parquet(f, columns=["site_id", "accession", "residues"])
                       for f in sorted(glob.glob(str(root / "inputs/shard_sites/*/*.parquet")))], ignore_index=True)
    sites = sites[sites.accession.isin(set(annots.accession))]
    annots = detect(annots, sites, params["detect_overlap_min"])
    # per OG: pairs
    taxid_of = dict(zip(members.accession, members.taxid))
    pair_frames = []
    for shard, ga in annots.groupby("shard"):
        msa = pd.read_parquet(root / "msa" / f"shard={int(shard):04d}.parquet")
        msa = msa[msa.og_id.isin(set(ga.og_id))]
        for og_id, g in ga.groupby("og_id"):
            if g.accession.nunique() < 2:
                continue
            mm = msa[msa.og_id == og_id]
            p2c = {acc: np.flatnonzero(np.frombuffer(s.encode(), np.uint8) != ord("-")) for acc, s in zip(mm.accession, mm.aligned)}
            pr = make_pairs(g, p2c, taxid_of)
            if len(pr):
                pr.insert(0, "og_id", og_id)
                pair_frames.append(pr)
    pairs = pd.concat(pair_frames, ignore_index=True)
    lab = ss.set_index("site_id")
    for t in taus:
        c = f"sog_t{t:g}"
        la = pairs.site_a.map(lab[c])
        lb = pairs.site_b.map(lab[c])
        pairs[f"same_t{t:g}"] = (la.notna() & lb.notna() & (la == lb)).astype(bool)
    pairs = pairs.merge(ogs[["og_id", "half", "forced_test"]], on="og_id")
    idb = stats.set_index("og_id").identity
    pairs["identity"] = pairs.og_id.map(idb)
    bins = params["identity_bins"]
    pairs["identity_bin"] = pd.cut(pairs.identity, [b[0] for b in bins] + [bins[-1][1] + 1e-9], right=False,
                                   labels=[f"[{b[0]:g},{b[1]:g})" for b in bins]).astype(str)
    pairs["primary"] = ~pairs.metal_only
    # ---- metrics
    rows = []
    for half in ["tune", "test"]:
        for subset, m in [("primary", pairs.primary), ("with_metals", pd.Series(True, index=pairs.index)),
                          ("uniprot", pairs.primary & pairs.sources.str.contains("uniprot")),
                          ("biolip", pairs.primary & pairs.sources.str.contains("biolip"))]:
            d = pairs[m & (pairs.half == half)]
            for t in taus:
                r = metrics(d, t)
                r.update({"half": half, "subset": subset})
                if subset == "primary":
                    r["f1_ci_lo"], r["f1_ci_hi"] = boot_ci(d, t, a.n_boot, seed)
                rows.append(r)
        for b in pairs.identity_bin.unique():
            d = pairs[pairs.primary & (pairs.half == half) & (pairs.identity_bin == b)]
            for t in taus:
                r = metrics(d, t)
                r.update({"half": half, "subset": f"identity {b}"})
                rows.append(r)
    R = pd.DataFrame(rows)
    out = Path(a.results)
    out.mkdir(parents=True, exist_ok=True)
    tag = a.label
    R.to_csv(out / f"tau_benchmark{tag}.csv", index=False)
    pairs.drop(columns=[]).to_parquet(out / f"tau_benchmark_pairs{tag}.parquet", index=False)
    # annotation disjointness: annotations of TUNE-half OGs used for tuning (exclude from evidence filter 2)
    used = annots.merge(ogs[["og_id", "half"]], on="og_id")
    used = used[used.half == "tune"][["ann_id", "source", "accession", "og_id", "ligand_key", "site_id", "detected"]]
    used.to_parquet(out / f"tau_tuning_annotations{tag}.parquet", index=False)
    side = {"snapshot_id": json.load(open(P / "data/snapshot" / root.parent.name / "manifest.json"))["snapshot_id"]
            if (P / "data/snapshot" / root.parent.name / "manifest.json").exists() else root.parent.name,
            "prereg_sha256": psha, "git_commit": git_commit(), "seed": seed, "code_version": git_commit(),
            "annotations_available": have, "orthology_version": root.name}
    for f in [f"tau_benchmark{tag}.csv", f"tau_benchmark_pairs{tag}.parquet", f"tau_tuning_annotations{tag}.parquet"]:
        json.dump(side, open(out / (f + ".json"), "w"), indent=1)
    tune = R[(R.half == "tune") & (R.subset == "primary")].copy()
    choice = None
    if not a.no_choice:
        ok = tune[tune.f1.notna()]
        if len(ok) == 0:
            raise SystemExit("no tune-half pairs: cannot choose tau")
        ok = ok.sort_values(["f1", "precision", "tau"], ascending=[False, False, True])
        choice = float(ok.tau.iloc[0])
        tc = dict(side)
        tc.update({"tau": choice, "objective": prereg["sogs"]["tau_column_jaccard"]["objective"],
                   "grid": taus, "tie_break": "max F1, then max precision, then smaller tau",
                   "tune_f1": {f"{t:g}": float(v) for t, v in zip(tune.tau, tune.f1)},
                   "msa_method": "FAMSA", "leiden_resolution": params["leiden_resolution"],
                   "n_pos_tune": int(tune.n_pos.iloc[0]), "n_neg_tune": int(tune.n_neg.iloc[0]),
                   "created_utc": pd.Timestamp.now("UTC").isoformat()})
        json.dump(tc, open(out / "tau_choice.json", "w"), indent=1)
    # summary of annotation funnel
    funnel = {"n_annot_sites": int(len(annots)), "n_annot_dropped_model_seq_mismatch": int(n_seq_dropped), "by_source": annots.source.value_counts().to_dict(),
              "n_detected": int(annots.detected.sum()),
              "detected_frac_by_source": annots.groupby("source").detected.mean().round(4).to_dict(),
              "n_ogs_with_pairs": int(pairs.og_id.nunique()),
              "pairs": pairs.groupby(["half", "kind"]).size().to_dict(), "wall_s": time.time() - t0}
    funnel = {k: ({str(kk): vv for kk, vv in v.items()} if isinstance(v, dict) else v) for k, v in funnel.items()}
    json.dump(funnel, open(out / f"tau_benchmark_funnel{tag}.json", "w"), indent=1)
    print(R[R.subset == "primary"].to_string())
    print(json.dumps(funnel))
    print("tau choice:", choice)


if __name__ == "__main__":
    main()
