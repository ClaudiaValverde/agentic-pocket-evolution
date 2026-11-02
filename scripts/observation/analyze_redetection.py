#!/usr/bin/env python
"""Re-detection rates and technical noise floor (DESIGN §3.3 (i); prereg states.measured_error_rates).

Reference = atlas v6 sites (tier 1, snapshot manifest) of the sampled proteins. Comparisons (same protein):
  foldcomp_v4_kept   Foldcomp v4 model, pocketcompare describe --af-mode, atlas filters (score >= 0.03, lining
                     pLDDT >= 70)                                    -> "model source" (v6 file vs Foldcomp v4)
  foldcomp_v4_all    same, no filters (what a W call could use)
  seed_s0 / s1 / s2  v6 model, BURROW training seed swapped (s0 = control, must equal atlas), atlas filters
Sites = pocketevo.sites.build_sites (Jaccard 0.3 merge) on both sides. A reference site is re-detected if its best
comparison site of the same protein has residue Jaccard >= tau (UniProt positions). Matched pairs give
fingerprint cosines (rep pockets; fp16 -> float32, renormalised).

Only proteins whose v4 sequence equals the v6 sequence are used (mismatches counted). Every tau in --tau is
computed (cheap; tau final comes from the SOG benchmark).

  python scripts/observation/analyze_redetection.py --tau 0.2 0.3 0.4 0.5
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

P = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(P / "src"))
from pocketevo import observation as O  # noqa: E402
from pocketevo import sites as SI  # noqa: E402

ATLAS = Path("/cluster/scratch/gmagnusson/pocketcompare_afdb")
SNAP = "20261030T085136Z"
MIN_SCORE, MIN_PLDDT = 0.03, 70.0
SCORE_BINS = [0.0, 0.1, 0.3, 1.01]
SCORE_LABELS = ["<0.1", "0.1-0.3", ">=0.3"]


def atlas_pockets(units, accs) -> pa.Table:
    man = pd.read_parquet(P / f"data/snapshot/{SNAP}/manifest.parquet").set_index("unit_id")
    cols = SI.POCKET_COLUMNS + ["fingerprint"]
    parts = []
    aset = pa.array(sorted(accs))
    for u in sorted(set(units)):
        t = pq.read_table(ATLAS / man.at[u, "pockets_path"], columns=cols)
        parts.append(t.filter(pc.is_in(t["accession"], aset)))
    return pa.concat_tables(parts, promote_options="permissive")


def described_pockets(path: Path, acc_of_stem: dict, model_version: int, filt: bool) -> tuple[pa.Table, dict]:
    if path.is_dir() or not path.exists():
        raise FileNotFoundError(path)
    t = pq.read_table(path)
    rec = json.loads(t.schema.metadata[b"pocketcompare"])
    d = t.to_pandas()
    d = d[d["ok"].astype(bool)]
    d["stem"] = d.structure_stem.astype(str)
    d = d[d.stem.isin(acc_of_stem)]
    if filt:
        d = d[(d.burrow_score >= MIN_SCORE) & (np.nan_to_num(d.plddt_mean.to_numpy(float), nan=-1) >= MIN_PLDDT)]
    d = d.assign(unit_id=-1, tier=1, accession=d.stem.map(acc_of_stem), fragment=1, model_version=model_version,
                 taxid=0, organism="")
    V = np.stack(d.fingerprint.to_numpy()) if len(d) else np.zeros((0, 384), np.float32)
    d = d[[c for c in SI.POCKET_COLUMNS]]
    tb = pa.Table.from_pandas(d, preserve_index=False)
    fp = pa.FixedSizeListArray.from_arrays(pa.array(V.astype(np.float16).reshape(-1), pa.float16()), 384)
    tb = tb.append_column("fingerprint", fp)
    return tb.cast(pa.schema([SI.POCKET_SCHEMA.field(c) for c in tb.column_names if c in SI.POCKET_SCHEMA.names])), rec


def sites_df(pockets: pa.Table):
    s, f, _ = SI.build_sites(pockets, jaccard=SI.DEFAULT_JACCARD, with_fingerprints=True)
    s = s.to_pandas()
    V = np.asarray(f.column("fingerprint").combine_chunks().flatten(), np.float32).reshape(-1, 384)
    V /= np.maximum(np.linalg.norm(V, axis=1, keepdims=True), 1e-12)
    s["fp_row"] = np.arange(len(s))
    return s, V


def match(ref: pd.DataFrame, Vr, qry: pd.DataFrame, Vq, accs) -> pd.DataFrame:
    """Best-Jaccard comparison site per reference site (same accession)."""
    rows = []
    qg = {a: g for a, g in qry.groupby("accession")}
    for a, g in ref.groupby("accession"):
        if a not in accs:
            continue
        q = qg.get(a)
        qsets = [frozenset(x) for x in q.residues] if q is not None else []
        for r in g.itertuples():
            rs = frozenset(r.residues)
            bi, bj = -1, 0.0
            for j, qs in enumerate(qsets):
                if rs & qs:
                    v = O.jaccard(rs, qs)
                    if v > bj:
                        bi, bj = j, v
            cos = float(Vr[r.fp_row] @ Vq[q.fp_row.iloc[bi]]) if bi >= 0 else np.nan
            cov = O.coverage(rs, qsets[bi]) if bi >= 0 else 0.0
            rows.append(dict(site_id=r.site_id, accession=a, n_res=r.n_res, plddt=r.plddt_mean_rep,
                             score=r.max_burrow_score, best_jaccard=bj, best_cov=cov, cos_best=cos,
                             best_site=(q.site_id.iloc[bi] if bi >= 0 else None)))
    return pd.DataFrame(rows)


def reverse_rate(qry, ref, accs, tau):
    """Fraction of comparison sites with a reference counterpart (Jaccard >= tau)."""
    rg = {a: [frozenset(x) for x in g.residues] for a, g in ref.groupby("accession")}
    hit = tot = 0
    for r in qry.itertuples():
        if r.accession not in accs:
            continue
        tot += 1
        rs = frozenset(r.residues)
        hit += any(O.jaccard(rs, x) >= tau for x in rg.get(r.accession, []))
    return hit / tot if tot else np.nan, tot


def q(x, ps=(0.05, 0.10, 0.25, 0.50)):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return {f"p{int(p * 100):02d}": (float(np.quantile(x, p)) if len(x) else None) for p in ps} | {"n": int(len(x))}


def git_commit():
    try:
        return subprocess.run(["git", "-C", str(P), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    except OSError:
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", default="/tmp/pe_big/obs")
    ap.add_argument("--tau", type=float, nargs="+", default=[0.3])
    ap.add_argument("--primary-tau", type=float, default=0.3)
    ap.add_argument("--out", default=str(P / f"results/observation/{SNAP}"))
    a = ap.parse_args()
    t0 = time.time()
    W, out = Path(a.work), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    M = pd.read_parquet(W / "redetect_models.parquet")
    ok = M[M.seq_identical]
    accs = set(ok.accession)
    meta = M.set_index("accession")
    ref_p = atlas_pockets(M.unit_id, M.accession)
    atlas_hash = {O.parquet_method_hash(ATLAS / p) for p in
                  pd.read_parquet(P / f"data/snapshot/{SNAP}/manifest.parquet").set_index("unit_id")
                  .loc[sorted(set(M.unit_id)), "pockets_path"]}
    ref, Vr = sites_df(ref_p)
    comps = {}
    stem4 = {n: acc for n, acc in zip(ok.fc_name, ok.accession)}
    stem6 = {s: acc for s, acc in zip(ok.structure_stem, ok.accession)}
    hashes = {}
    t4, rec = described_pockets(W / "v4_describe.parquet", stem4, 4, True)
    hashes["foldcomp_v4"] = O.method_record_hash(rec)
    comps["foldcomp_v4_kept"] = (t4, accs)
    t4a, _ = described_pockets(W / "v4_describe.parquet", stem4, 4, False)
    comps["foldcomp_v4_all"] = (t4a, accs)
    seed_accs = set(ok[ok.in_seed_subset].accession)
    for s in ("s0", "s1", "s2"):
        parts = sorted(p for p in W.glob(f"seed_{s}_part*.parquet") if not p.name.endswith(".structures.parquet"))
        if not parts:
            continue
        tabs = []
        for p_ in parts:
            tb, rec = described_pockets(p_, stem6, 6, True)
            tabs.append(tb)
            hashes[f"seed_{s}"] = O.method_record_hash(rec)
        comps[f"seed_{s}"] = (pa.concat_tables(tabs), seed_accs)
    per_site, summary = [], []
    for name, (tb, acs) in comps.items():
        qs, Vq = sites_df(tb)
        m = match(ref, Vr, qs, Vq, acs)
        m["comparison"] = name
        per_site.append(m)
        for tau in a.tau:
            rr, nq = reverse_rate(qs, ref, acs, tau)
            summary.append(dict(comparison=name, tau=tau, n_proteins=len(acs & set(ref.accession)),
                                n_ref_sites=int(len(m)), redetect=float((m.best_jaccard >= tau).mean()),
                                n_cmp_sites=nq, reverse_redetect=rr))
    S = pd.concat(per_site, ignore_index=True)
    S["kingdom"] = S.accession.map(meta.kingdom)
    S["protein_band"] = S.accession.map(meta.band)
    S["band"] = [O.plddt_band(x) for x in S.plddt]
    S["score_bin"] = pd.cut(S.score, SCORE_BINS, right=False, labels=SCORE_LABELS).astype(str)
    S.to_parquet(out / "redetection_sites.parquet", index=False)
    # stratified table, every tau
    strata = []
    for tau in a.tau:
        S["hit"] = S.best_jaccard >= tau
        for by in (["comparison"], ["comparison", "band"], ["comparison", "kingdom"], ["comparison", "score_bin"],
                   ["comparison", "band", "kingdom"]):
            g = S.groupby(by, observed=True)
            for k, x in g:
                k = k if isinstance(k, tuple) else (k,)
                cos = x.loc[x.hit, "cos_best"]
                strata.append(dict(tau=tau, stratum="+".join(by[1:]) or "all", **dict(zip(by, k)),
                                   n_sites=int(len(x)), n_proteins=int(x.accession.nunique()),
                                   redetect=float(x.hit.mean()),
                                   **{f"cos_{kk}": v for kk, v in q(cos).items() if kk != "n"}, n_matched=int(x.hit.sum())))
    T = pd.DataFrame(strata)
    T.to_csv(out / "redetection_strata.csv", index=False)
    # headline numbers at the primary tau
    S["hit"] = S.best_jaccard >= a.primary_tau
    head = {}
    for name in comps:
        x = S[S.comparison == name]
        head[name] = dict(redetect=float(x.hit.mean()), n_sites=int(len(x)), n_proteins=int(x.accession.nunique()),
                          redetect_by_band={b: float(y.hit.mean()) for b, y in x.groupby("band")},
                          n_by_band={b: int(len(y)) for b, y in x.groupby("band")},
                          cos_matched=q(x.loc[x.hit, "cos_best"]),
                          cos_matched_by_band={b: q(y.loc[y.hit, "cos_best"]) for b, y in x.groupby("band")})
    res = dict(snapshot_id=SNAP, prereg_sha256=hashlib.sha256((P / "config/prereg.yaml").read_bytes()).hexdigest(),
               git_commit=git_commit(), seed=0, primary_tau=a.primary_tau, taus=a.tau,
               n_sample=int(len(M)), n_seq_identical=int(M.seq_identical.sum()),
               n_seq_mismatch=int((~M.seq_identical).sum()),
               method_hash_atlas_files=sorted(atlas_hash), method_hash_expected=O.ATLAS_METHOD_HASH,
               method_hash_comparisons=hashes,
               method_hash_ok=bool(atlas_hash == {O.ATLAS_METHOD_HASH} and all(h == O.ATLAS_METHOD_HASH for h in hashes.values())),
               summary=pd.DataFrame(summary).to_dict("records"), headline=head, runtime_s=round(time.time() - t0, 1))
    (out / "redetection.json").write_text(json.dumps(res, indent=1, default=float))
    print(json.dumps({k: v for k, v in res.items() if k not in ("headline",)}, indent=1, default=float))
    for name, h in head.items():
        print(name, round(h["redetect"], 3), h["n_sites"], {b: round(v, 3) for b, v in h["redetect_by_band"].items()},
              h["cos_matched"])


if __name__ == "__main__":
    main()
