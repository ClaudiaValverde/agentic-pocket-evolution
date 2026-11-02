#!/usr/bin/env python
"""Specificity-change events on real SOGs (DESIGN §5.2; prereg part_b.specificity_change, pi_freeze.specificity_null).

Residue-first: at every gene-tree node where the SOG is present throughout (>= 80 % of informative members P/W
on both sides), derived clade D vs its sister S:
  (1) >= 2 footprint columns with clade-consistent substitutions, BLOSUM62 < 0 (residues from the FAMSA MSA at
      the SOG footprint columns);
  (2) fingerprint mean shift 1 - cos(mean_D, mean_S), after regressing out site volume, site pLDDT (and model
      source; constant on tier 1), above the technical noise floor = q95 of (1 - same-site replicate cosine)
      across BURROW seeds s1/s2 (results/observation/<snap>/redetection_sites.parquet);
  (3) empirical null: shifts of the other sister pairs without residue signal, matched on depth and clade size;
      BH (alpha 0.05) across all residue-gated pairs. Clade permutation (secondary) for called pairs.
Null input: pseudo-clades (member identities shuffled among the P/W leaves of each SOG tree; residues and
fingerprints travel with the member), prereg seeds; funnel real vs null.

    python scripts/partB/run_specificity.py --release O1_r3 --tau 0.2
    python scripts/partB/run_specificity.py --orthology O2 --release O2_r1 --run-id O2_r1_spec1 --n-null-seeds 3
O2: states = rescue overlay (/tmp/pe_big/partB/<snap>/O2_r1/states_resc_t<tau>, as the O2 events stage), O2
orthoxml + guide tree, site fingerprints of tiers 1, 2, 4 (site_fp/) and 3 (fingerprints_O2/, atlas rows of the
rep pocket, scripts/partB/extract_site_fp.py); covariate 3 (model source) = 1 for AFDB v4 / Foldcomp models
(tiers 3, L), 0 for v6 (tiers 1, 2). Output: results/partB/<snap>/<run-id>/specificity/tau<tau> (run-id
defaults to the release, so the O1 path is unchanged).
"""
from __future__ import annotations

import argparse
import glob
import json
import re
import zlib
import sys
import time
from multiprocessing import get_context
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.dataset as ds

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from pocketevo import events, recon  # noqa: E402

SNAP = recon.cli_value("--snapshot") or recon.default_snapshot()   # --snapshot / $PE_SNAPSHOT / prereg
_G = {}


def _og(args):
    og, sog_ids, codes, seed = args
    t = _G["trees"][og]
    names = t.leaf_names()
    rows = []
    for j, s in enumerate(sog_ids):
        res, fp, cov = _G["res"].get(s, {}), _G["fp"].get(s, {}), _G["cov"].get(s, {})
        c = codes[:, j]
        if seed is not None:       # pseudo-clade null: permute member identities among P/W leaves
            rng = np.random.default_rng([seed, zlib.crc32(s.encode())])
            pw = np.flatnonzero(np.isin(c, (0, 1)))
            perm = pw[rng.permutation(len(pw))]
            mp = {names[a]: names[b] for a, b in zip(pw, perm)}
            res = {a: res[b] for a, b in mp.items() if b in res}
            fp = {a: fp[b] for a, b in mp.items() if b in fp}
            cov = {a: cov[b] for a, b in mp.items() if b in cov}
        rows += events.specificity_tests(t, c, res, fp, cov, _G["floor"], min_subst=_G["min_sub"],
                                         species_tree=_G["st"], sog_id=s)
    return rows


def funnel_row(df, n_tested, inp, seed):
    return dict(input=inp, seed=seed, n_pairs=len(df), residue_gate=int(df["residue_pass"].sum()),
                residue_and_shift=int((df["residue_pass"] & df["shift_pass"]).sum()), bh_tested=n_tested,
                called=int(df["called"].sum()),
                called_f1_ge3_species=int((df["called"] & (df["n_species_d"] >= 3)).sum()))


def d34_columns(df, rows, alpha):
    """D34 (post hoc sensitivity): tail-fitted p (normal = D17 rule; GPD), BH over the same tested set, same call
    rule (residue gate & shift above the noise floor & q <= alpha). Adds columns to df; returns funnel extras."""
    tp = events.specificity_tail_pvalues(rows)
    tested = df["tested"].to_numpy()
    ext = {}
    if tested.any():
        assert np.allclose(tp["p_emp"][tested], df["p_empirical"].to_numpy()[tested])
        ext["min_p_empirical"] = float(np.nanmin(df["p_empirical"]))
        ext["min_attainable_p_empirical"] = float(np.nanmin(1.0 / (1.0 + tp["n_pool"][tested])))
        ext["bh_threshold_single_call"] = alpha / int(tested.sum())
    for k in ("n_pool", "pool_max", "gpd_xi"):
        df[k] = tp[k]
    gate = (df["residue_pass"] & df["shift_pass"]).to_numpy()
    for v in ("normal", "gpd"):
        p = tp[f"p_{v}"]
        q = np.full(len(df), np.nan)
        if tested.any():
            q[tested] = events.bh(p[tested], alpha)[0]
        df[f"p_tail_{v}"] = p
        df[f"q_tail_{v}"] = q
        df[f"called_tail_{v}"] = gate & (np.nan_to_num(q, nan=2.0) <= alpha)
        ext[f"called_tail_{v}"] = int(df[f"called_tail_{v}"].sum())
        ext[f"min_p_tail_{v}"] = float(np.nanmin(p)) if tested.any() else np.nan
    return ext


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--release", "--sog-run", dest="release", default="O1_r3")
    recon.add_partb_args(ap, skip=("sog_run",))
    ap.add_argument("--tau", type=float, default=0.2)
    ap.add_argument("--min-p-members", type=int, default=3)
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--perm-n", type=int, default=None, help="default prereg permutation n (10000)")
    ap.add_argument("--n-null-seeds", type=int, default=20, help="pseudo-clade null seeds (D24: >= 20)")
    ap.add_argument("--orthoxml", default=None, help="default: O1 ancestra_tier1_v0 / O2 ancestra_O2")
    ap.add_argument("--guide", default=None)
    ap.add_argument("--states-glob", default=None, help="default: O1 bulk states_t<tau>; O2 rescue overlay")
    ap.add_argument("--fp-glob", action="append", default=None, help="site fingerprint parquet globs (repeatable)")
    ap.add_argument("--sites-glob", action="append", default=None, help="site metadata globs (covariates)")
    ap.add_argument("--rows-dir", default=None, help="D34: write all rows (real + each null seed) with tail-fit "
                    "p-values here (bulk, e.g. /tmp/pe_big/partB/...)")
    a = ap.parse_args()
    o2 = a.orthology == "O2"
    od = "ancestra_O2" if o2 else "ancestra_tier1_v0"
    a.orthoxml = a.orthoxml or str(ROOT / f"data/orthology/{od}/accurate/hogs.orthoxml.gz")
    a.guide = a.guide or str(ROOT / "data/orthology" / od / ("o2_guide_tree.nwk" if o2 else "tier1_guide_tree.nwk"))
    sb = f"/tmp/pe_big/sites/{SNAP}"
    a.fp_glob = a.fp_glob or ([f"{sb}/site_fp/tier{t}/*.parquet" for t in (1, 2, 4)]
                              + [f"{sb}/fingerprints_O2/tier3/*.parquet"] if o2 else [f"{sb}/site_fp/tier1/*.parquet"])
    a.sites_glob = a.sites_glob or ([f"{sb}/sites/tier{t}/*.parquet" for t in (1, 2, 3, 4)] if o2
                                    else [f"{sb}/sites/tier1/*.parquet"])
    run_id = a.run_id or a.release
    t0 = time.time()
    pr = recon.load_prereg()
    sc = pr["part_b"]["specificity_change"]
    seeds = pr["meta"]["seeds"]
    alpha = pr["meta"]["multiple_testing"]["alpha"]
    bulk = Path(f"/tmp/pe_big/sogs/{SNAP}/{a.release}")
    out = ROOT / f"results/partB/{SNAP}/{run_id}/specificity/tau{a.tau}"
    out.mkdir(parents=True, exist_ok=True)
    # ---- noise floor (same-site replicate cosines, BURROW seeds, matched at tau)
    rd = pd.read_parquet(ROOT / f"results/observation/{SNAP}/redetection_sites.parquet")
    rep = rd[rd["comparison"].isin(["seed_s1", "seed_s2"]) & (rd["best_jaccard"] >= a.tau) & rd["cos_best"].notna()]
    q = sc["noise_floor_quantile"]["value"]
    floor = events.noise_floor(rep["cos_best"].to_numpy(), q)
    # ---- states, SOG selection
    sg = a.states_glob or (f"/tmp/pe_big/partB/{SNAP}/{a.release}/states_resc_t{a.tau}/shard=*.parquet" if o2
                           else str(bulk / f"states_t{a.tau}/shard=*.parquet"))
    sds = ds.dataset(sorted(glob.glob(sg)), format="parquet")
    st = sds.to_table(columns=["og_id", "sog_id", "accession", "state", "site_id", "taxid"]
                      + (["model_tier"] if "model_tier" in sds.schema.names else [])).to_pandas()
    if "model_tier" not in st:      # O1 states: tier-1 models only
        st["model_tier"] = 1
    meta = pd.read_parquet(ROOT / f"results/sogs/{SNAP}/{a.release}/sogs_all_tau.parquet")
    meta = meta[np.isclose(meta["tau"], a.tau)]
    bad = set(meta.loc[meta["small_footprint"].astype(bool) | meta["empty_footprint"].astype(bool), "sog_id"])
    cnt = {"sogs_total": int(st["sog_id"].nunique()), "noise_floor": floor, "noise_floor_n": len(rep)}
    st = st[~st["sog_id"].isin(bad)]
    cnt["after_D8_flags"] = int(st["sog_id"].nunique())
    pcount = st[st["state"] == "P"].groupby("sog_id").agg(nP=("accession", "size"), nsp=("taxid", "nunique"))
    keep = set(pcount.index[(pcount["nP"] >= a.min_p_members) & (pcount["nsp"] >= 2)])
    st = st[st["sog_id"].isin(keep)]
    cnt["eligible_ge3P_ge2species"] = len(keep)
    print(f"[{time.time() - t0:.0f}s] {cnt}", flush=True)
    # ---- residues at footprint columns (P/W members)
    fpc = ds.dataset(sorted(glob.glob(str(bulk / "sogs/shard=*.parquet"))), format="parquet").to_table(
        columns=["og_id", "sog_id", "tau", "footprint_columns"]).to_pandas()
    fpc = fpc[np.isclose(fpc["tau"], a.tau) & fpc["sog_id"].isin(keep)]
    cols_of = {s: np.asarray(c, dtype=int) for s, c in zip(fpc["sog_id"], fpc["footprint_columns"])}
    ogs = set(fpc["og_id"])
    msa = ds.dataset(sorted(glob.glob(str(bulk / "msa/shard=*.parquet"))), format="parquet").to_table().to_pandas()
    msa = msa[msa["og_id"].isin(ogs)]
    aln = {(o, x): s for o, x, s in zip(msa["og_id"], msa["accession"], msa["aligned"])}
    del msa
    pw = st[st["state"].isin(["P", "W"])]
    res = {}
    for og, s, acc in zip(pw["og_id"], pw["sog_id"], pw["accession"]):
        al, c = aln.get((og, acc)), cols_of.get(s)
        if al is None or c is None:
            continue
        res.setdefault(s, {})[acc] = "".join(al[k] if k < len(al) else "-" for k in c)
    # ---- fingerprints + covariates of the matched sites
    # site_ids are unique only within a tier (tier 3 and tier L share accessions): key by (model tier, site_id);
    # the tier of a fingerprint / site file is the 'tier<T>' component of its path
    sid = pa.array(sorted(set(pw["site_id"].dropna())))
    def _load(globs, cols):
        out = []
        for g in globs:
            tier = int(re.search(r"tier(\d+)", g).group(1))
            files = sorted(glob.glob(g))
            if files:
                t_ = ds.dataset(files, format="parquet").to_table(filter=pc.is_in(pc.field("site_id"), value_set=sid),
                                                                  columns=cols)
                out.append((tier, t_))
        return out
    fps, nrow = {}, 0
    for tier, t_ in _load(a.fp_glob, ["site_id", "fingerprint"]):
        nrow += t_.num_rows
        fps.update(zip(((tier, x) for x in t_["site_id"].to_pylist()),
                       (np.asarray(x, dtype=np.float32) for x in t_["fingerprint"].to_numpy(zero_copy_only=False))))
    cnt["fingerprint_rows"] = nrow
    cnt["fingerprint_duplicate_keys"] = nrow - len(fps)
    covs = {}
    for tier, t_ in _load(a.sites_glob, ["site_id", "rep_site_volume_A3", "plddt_mean_rep", "model_version"]):
        d_ = t_.to_pandas()
        # model source covariate: constant 0 on tier 1 (O1 unchanged); 1 = AFDB v4 / Foldcomp (tiers 3, L) on O2
        covs.update({(tier, s_): np.array([v, p_, float(o2 and mv == 4)]) for s_, v, p_, mv in
                     zip(d_["site_id"], d_["rep_site_volume_A3"], d_["plddt_mean_rep"], d_["model_version"])})
    fpd, covd = {}, {}
    nofp = 0
    for s, acc, site, mt in zip(pw["sog_id"], pw["accession"], pw["site_id"], pw["model_tier"]):
        k = (int(mt), site)
        if k in fps:
            fpd.setdefault(s, {})[acc] = fps[k]
            covd.setdefault(s, {})[acc] = covs.get(k, np.array([np.nan, np.nan, 0.0]))
        else:
            nofp += 1
    cnt["PW_members_without_fingerprint"] = nofp
    del fps
    print(f"[{time.time() - t0:.0f}s] residues {len(res)} SOGs, fingerprints {len(fpd)} SOGs", flush=True)
    # ---- trees + codes
    trees = {t.tree_id: t for t in recon.iter_orthoxml_trees(a.orthoxml,
                                                             keep=ogs)}
    sttree = recon.SpeciesTree.from_file(a.guide)
    asm = events.assemble(trees, st.assign(footprint_plddt=None, reason=None), {})
    for s in covd:   # covariates with NaN -> drop covariate regression for that SOG (counted)
        if any(np.isnan(v).any() for v in covd[s].values()):
            covd[s] = {}
    _G.update(trees=trees, res=res, fp=fpd, cov=covd, floor=floor,
              min_sub=sc["min_clade_consistent_substitutions"]["value"], st=sttree)
    fun, outdf = [], None
    for seed in [None] + list(range(a.n_null_seeds)):
        with get_context("fork").Pool(a.workers) as pool:
            rows = [r for rr in pool.imap_unordered(_og, [(og, s, c, seed) for og, s, c, _ in asm], chunksize=16)
                    for r in rr]
        called, n_tested = events.call_specificity(rows, alpha=alpha)
        df = pd.DataFrame(called)
        fun.append(funnel_row(df, n_tested, "real" if seed is None else "null_pseudo_clade", seed))
        if a.rows_dir:      # D34 sensitivity columns + all rows to bulk
            fun[-1].update(d34_columns(df, rows, alpha))
            rd_ = Path(a.rows_dir); rd_.mkdir(parents=True, exist_ok=True)
            og_of_ = {s: og for og, ss, _, _ in asm for s in ss}
            df.assign(og_id=df["sog_id"].map(og_of_), input="real" if seed is None else "null_pseudo_clade",
                      seed=-1 if seed is None else seed).to_parquet(
                rd_ / f"rows_{'real' if seed is None else f'null{seed}'}.parquet", index=False)
        print(f"[{time.time() - t0:.0f}s] {fun[-1]}", flush=True)
        if seed is None:
            outdf = df
    # secondary null (clade permutation) for called pairs
    n_perm = a.perm_n or int(sc["permutation"]["n"]["value"])
    rng = np.random.default_rng(seeds[0])
    outdf["p_permutation"] = np.nan
    og_of = {s: og for og, ss, _, _ in asm for s in ss}
    for i in np.flatnonzero(outdf["called"].to_numpy()):
        r = outdf.iloc[i]
        t = trees[og_of[r["sog_id"]]]
        names = t.leaf_names()
        accs = [x for x in names[t.leaf_lo[r["node"]]:t.leaf_hi[r["node"]]] if x in fpd.get(r["sog_id"], {})]
        dset = set(names[t.leaf_lo[r["d_node"]]:t.leaf_hi[r["d_node"]]])
        F = np.array([fpd[r["sog_id"]][x] for x in accs], dtype=float)
        if covd.get(r["sog_id"]):
            F = events.residualize(F, np.array([covd[r["sog_id"]][x] for x in accs]))
        isd = np.array([x in dset for x in accs])
        if isd.sum() >= 2 and (~isd).sum() >= 1:
            outdf.iat[i, outdf.columns.get_loc("p_permutation")] = events.permutation_pvalue(F, isd, n_perm, rng)
    outdf["og_id"] = outdf["sog_id"].map(og_of)
    outdf["d_node_name"] = [trees[o].node_name[d] for o, d in zip(outdf["og_id"], outdf["d_node"])]
    outdf["d_species"] = [";".join(sorted({trees[o].species[trees[o].leaves[k]] for k in
                                           range(trees[o].leaf_lo[d], trees[o].leaf_hi[d])}))
                          for o, d in zip(outdf["og_id"], outdf["d_node"])]
    files = []
    f = out / "specificity_tested.parquet"
    outdf[outdf["tested"]].to_parquet(f, index=False); files.append(f)
    f = out / "specificity_called.parquet"
    outdf[outdf["called"]].to_parquet(f, index=False); files.append(f)
    ex = outdf[outdf["residue_pass"] & outdf["shift_pass"]].sort_values(["q_bh", "shift"], ascending=[True, False]).copy()
    ex["label"] = "exploratory, not significant"
    ex["exploratory_top5"] = False
    ex.iloc[:5, ex.columns.get_loc("exploratory_top5")] = True
    f = out / "specificity_exploratory.parquet"
    ex.to_parquet(f, index=False); files.append(f)
    f = out / "funnel_real_vs_null.csv"
    pd.DataFrame(fun).to_csv(f, index=False); files.append(f)
    cnt.update(n_pairs=len(outdf), runtime_s=time.time() - t0, alpha=alpha, permutation_n=n_perm,
               prereg_sha256=recon.prereg_sha256(), release=a.release, tau=a.tau, run_id=run_id,
               orthology=a.orthology, states_glob=sg, fp_glob=a.fp_glob, orthoxml=a.orthoxml, guide=a.guide,
               n_null_seeds=a.n_null_seeds)
    f = out / "manifest.json"
    f.write_text(json.dumps(cnt, indent=1, default=float)); files.append(f)
    for f in files:
        recon.write_sidecar(f, SNAP, seeds, {"script": "scripts/partB/run_specificity.py", "argv": sys.argv[1:]})
    print(pd.DataFrame(fun).to_string(index=False))


if __name__ == "__main__":
    main()
