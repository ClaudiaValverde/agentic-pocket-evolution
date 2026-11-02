#!/usr/bin/env python
"""Part A M1/M2 pair table (DESIGN §4.1). One pass over the SOG state shards; per root OG:

  * member x SOG indicator matrices (P, informative per variant) -> exact ordered counts for every member pair;
  * candidate pairs = pairs with >= 1 informative observation in any variant;
  * 'cap' sample: <= --cap pairs per OG (seeded), each pair weight 1 (big families cannot dominate);
  * 'unc' sample: <= --unc-max pairs per OG (seeded), Horvitz-Thompson weight n_candidates / n_sampled
    (estimates the all-pairs, uncapped rate);  control OGs (--control-ogs) additionally get every pair (in_all);
  * global identity from the OG MSA, LCA relation from the HOG tree, species of both members;
  * M2: every P/P observation of a 'cap' (or control in_all) pair -> fingerprint cosine of the two sites.

Outputs (bulk, regenerable) under --out/<snapshot>/<orthology>/: pairs.parquet, m2_matches.parquet (+ sidecars).
"""
import argparse
import json
import multiprocessing as mp
import sys
import time
import zlib
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))
from pocketevo import descent  # noqa: E402
from pocketevo.homology import lca_type, parse_orthoxml  # noqa: E402
from pocketevo.pairs import check_prereg, load_prereg, write_sidecar  # noqa: E402

COLS = ["og_id", "sog_id", "sog_index", "accession", "state", "state_s70", "site_id", "footprint_plddt",
        "aligned_frac", "frac_lt70", "len_ratio", "taxid", "n_filtered_pockets"]
G = {}


def init_globals(a):
    pre = load_prereg()
    ap = pre["states"]["A_primary"]
    G["crit"] = {"plddt": ap["footprint_mean_plddt"]["value"],
                 "plddt_s70": pre["states"]["A_sensitivity_stratum"]["footprint_mean_plddt"],
                 "frac_lt70": ap["footprint_frac_residues_below_70"]["value"],
                 "aligned": ap["footprint_aligned_non_gap_frac"]["value"],
                 "len_lo": ap["protein_length_vs_og_median"]["value"][0],
                 "len_hi": ap["protein_length_vs_og_median"]["value"][1]}
    G["seed"] = int(pre["meta"]["seeds"][0])
    G["ceiling"] = descent.load_ceiling(pd.read_parquet(a.error_rates), a.ceiling_source, a.tau)
    t0 = time.time()
    f = parse_orthoxml(str(Path(a.orthology) / "hogs.orthoxml.gz"))
    G["forest"] = f
    G["gene_of"] = {acc: k for k, acc in enumerate(f.gene_acc)}
    print(f"forest {len(f.gene_acc)} leaves, {len(f.root_names)} roots, {time.time()-t0:.0f}s", flush=True)
    u = pq.read_table(a.universe, columns=["si", "site_id", "n_res", "tier"]).to_pandas()
    u = u[u.tier == 1]
    G["si_of"] = pd.Series(u.si.to_numpy(np.int64), index=u.site_id.to_numpy())
    G["nres_of"] = pd.Series(u.n_res.to_numpy(np.int64), index=u.site_id.to_numpy())
    sg = pd.read_parquet(a.sogs_web, columns=["sog_id", "footprint_size"])
    G["big_fp"] = pd.Series(sg.footprint_size.to_numpy() > a.small_footprint_max, index=sg.sog_id.to_numpy())
    G["fp_size"] = pd.Series(sg.footprint_size.to_numpy(float), index=sg.sog_id.to_numpy())   # for qm90d11
    G["fp_path"] = a.fingerprints
    G["args"] = a
    G["controls"] = set(a.control_ogs.split(",")) if a.control_ogs else set()


def og_rng(og_id: str, salt: int) -> np.random.Generator:
    return np.random.default_rng([G["seed"], zlib.crc32(og_id.encode()), salt])


def process_og(og_id, rows, aln, fp):
    a = G["args"]
    rows = rows.reset_index(drop=True)
    acc = np.array(sorted(set(rows.accession)), dtype=object)
    midx = pd.Series(np.arange(len(acc)), index=acc)
    sogs = np.array(sorted(set(rows.sog_index)))
    sidx = pd.Series(np.arange(len(sogs)), index=sogs)
    rows["big_fp"] = G["big_fp"].reindex(rows.sog_id.to_numpy()).fillna(True).to_numpy(bool)
    rows["fp_res"] = np.round(G["fp_size"].reindex(rows.sog_id.to_numpy()).to_numpy()
                              * rows.aligned_frac.to_numpy(float))   # footprint residues (qm90d11, D11 symmetric)
    masks = descent.informative_masks(rows, G["crit"])
    nres = G["nres_of"].reindex(rows.site_id.to_numpy()).to_numpy(float)
    si = G["si_of"].reindex(rows.site_id.to_numpy()).fillna(-1).to_numpy(np.int64)
    small = np.nan_to_num(nres, nan=99) <= 10
    M = descent.og_matrices(midx[rows.accession.to_numpy()].to_numpy(), sidx[rows.sog_index.to_numpy()].to_numpy(),
                            len(acc), len(sogs), masks, small, si)
    C = descent.pair_counts(M, G["ceiling"])
    tot = C["lit"]["n_inf"] + C["s70"]["n_inf"]
    ii, jj = np.nonzero(np.triu(tot, 1))
    n_cand = len(ii)
    if n_cand == 0:
        return None, None
    ctrl = og_id in G["controls"]
    kc = descent.sample_pairs(ii, jj, a.cap, og_rng(og_id, 1))
    ku = descent.sample_pairs(ii, jj, a.unc_max, og_rng(og_id, 2))
    ka = np.arange(n_cand) if ctrl else np.zeros(0, np.int64)   # control OGs: every pair, flagged in_all
    keep = np.union1d(np.union1d(kc, ku), ka)
    pi, pj = ii[keep], jj[keep]
    # identity from the MSA
    ids = aln.accession.to_numpy()
    mrow = pd.Series(np.arange(len(ids)), index=ids)
    L = len(aln.aligned.iloc[0])
    mat = np.frombuffer("".join(aln.aligned).encode(), np.uint8).reshape(len(ids), L)
    ra = mrow.reindex(acc[pi]).to_numpy()
    rb = mrow.reindex(acc[pj]).to_numpy()
    okm = ~(pd.isna(ra) | pd.isna(rb))
    ident = np.full(len(pi), np.nan)
    ncol = np.zeros(len(pi), np.int64)
    if okm.any():
        ident[okm], ncol[okm] = descent.pair_identity(mat, ra[okm].astype(np.int64), rb[okm].astype(np.int64))
    # relation from the HOG tree
    ga = np.array([G["gene_of"].get(x, -1) for x in acc[pi]])
    gb = np.array([G["gene_of"].get(x, -1) for x in acc[pj]])
    rel = np.zeros(len(pi), np.int8)
    okg = (ga >= 0) & (gb >= 0)
    rel[okg] = lca_type(G["forest"], ga[okg], gb[okg])
    tax = rows.groupby("accession").taxid.first()
    df = pd.DataFrame({"og_id": og_id, "acc_a": acc[pi], "acc_b": acc[pj], "rel": rel,
                       "taxid_a": tax[acc[pi]].to_numpy(), "taxid_b": tax[acc[pj]].to_numpy(),
                       "identity": ident, "n_aligned_cols": ncol,
                       "in_cap": np.isin(keep, kc), "in_unc": np.isin(keep, ku),
                       "w_unc": np.where(np.isin(keep, ku), n_cand / len(ku), 0.0),
                       "in_all": np.isin(keep, ka), "n_cand": n_cand, "control": ctrl})
    for v in descent.VARIANTS:
        for q in ("n_inf", "n_match", "ceil_sum"):
            df[f"{q}_{v}"] = C[v][q][pi, pj].astype(np.float32)
        for k, st in enumerate(descent.STRATA):   # informative obs per (band_Y, size_X) stratum -> any ceiling later
            df[f"n_inf_{v}_s{k}"] = C[v]["strata"][st][pi, pj].astype(np.float32)
    # M2: P/P observations of cap pairs
    cap_rows = np.flatnonzero(df.in_cap.to_numpy() | df.in_all.to_numpy())
    mm = descent.matched_site_pairs(M, pi[cap_rows], pj[cap_rows])
    m2 = None
    if len(mm):
        mm = mm[(mm.si_a >= 0) & (mm.si_b >= 0)]
        if len(mm):
            fa = fp[mm.si_a.to_numpy()].astype(np.float32)
            fb = fp[mm.si_b.to_numpy()].astype(np.float32)
            cos = (fa * fb).sum(1) / (np.linalg.norm(fa, axis=1) * np.linalg.norm(fb, axis=1) + 1e-12)
            m2 = pd.DataFrame({"og_id": og_id, "pair_row": cap_rows[mm.k.to_numpy()],
                               "sog_index": sogs[mm.sog.to_numpy()], "si_a": mm.si_a.to_numpy(),
                               "si_b": mm.si_b.to_numpy(), "cosine": cos.astype(np.float32)})
    return df, m2


def process_shard(shard):
    a = G["args"]
    t0 = time.time()
    sp_ = Path(a.states) / shard
    st = pq.read_table(sp_, columns=COLS).to_pandas()
    st = st[(st.state.isin(["P", "A"])) | (st.state_s70 == "A")]
    msa = pq.read_table(Path(a.msa) / shard).to_pandas()
    msa_by = dict(tuple(msa.groupby("og_id")))
    fp = np.load(G["fp_path"], mmap_mode="r")
    outs, m2s, n_og, n_nomsa = [], [], 0, 0
    for og_id, rows in st.groupby("og_id", sort=True):
        if (rows.state == "P").sum() == 0:
            continue
        aln = msa_by.get(og_id)
        if aln is None:
            n_nomsa += 1
            continue
        df, m2 = process_og(og_id, rows, aln, fp)
        n_og += 1
        if df is not None:
            outs.append(df)
            if m2 is not None:
                m2s.append(m2)
    pairs = pd.concat(outs, ignore_index=True) if outs else None
    m2 = pd.concat(m2s, ignore_index=True) if m2s else None
    return shard, pairs, m2, {"n_og": n_og, "n_nomsa": n_nomsa, "s": time.time() - t0}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--snapshot", default="20261030T085136Z")
    ap.add_argument("--orthology-version", default="O1")
    ap.add_argument("--orthology", default="data/orthology/ancestra_tier1_v0/accurate")
    ap.add_argument("--sogs-bulk", default="/tmp/pe_big/sogs/20261030T085136Z/O1")
    ap.add_argument("--tau", type=float, default=None, help="default: results/sogs/.../tau_choice.json")
    ap.add_argument("--sogs-web", default="results/sogs/20261030T085136Z/O1/web/sogs.parquet")
    ap.add_argument("--universe", default="/tmp/pe_big/partA/20261030T085136Z/universe.parquet")
    ap.add_argument("--fingerprints", default="/tmp/pe_big/partA/20261030T085136Z/fp.f16.npy")
    ap.add_argument("--error-rates", default="results/observation/20261030T085136Z/error_rates.parquet")
    ap.add_argument("--ceiling-source", default="foldcomp_v4_kept")
    ap.add_argument("--cap", type=int, default=200)
    ap.add_argument("--unc-max", type=int, default=5000)
    ap.add_argument("--small-footprint-max", type=int, default=3)
    ap.add_argument("--control-ogs", default="HOG0027286,HOG0025274,HOG0035618")
    ap.add_argument("--workers", type=int, default=24)
    ap.add_argument("--limit-shards", type=int, default=0)
    ap.add_argument("--out", default="/tmp/pe_big/partA_descent")
    a = ap.parse_args()
    t0 = time.time()
    check_prereg()
    if a.tau is None:
        a.tau = float(json.load(open(REPO / f"results/sogs/{a.snapshot}/{a.orthology_version}/tau_choice.json"))["tau"])
    a.states = str(Path(a.sogs_bulk) / f"states_t{a.tau:g}")
    a.msa = str(Path(a.sogs_bulk) / "msa")
    init_globals(a)
    shards = sorted(p.name for p in Path(a.states).glob("shard=*.parquet"))
    if a.limit_shards:
        shards = shards[:a.limit_shards]
    out = Path(a.out) / a.snapshot / a.orthology_version
    out.mkdir(parents=True, exist_ok=True)
    pairs, m2s, stats = [], [], {"n_og": 0, "n_nomsa": 0}
    ctx = mp.get_context("fork")
    with ctx.Pool(a.workers) as pool:
        for k, (sh, p, m, s) in enumerate(pool.imap_unordered(process_shard, shards, chunksize=1)):
            if p is not None:
                pairs.append(p)
            if m is not None:
                m2s.append(m)
            stats["n_og"] += s["n_og"]
            stats["n_nomsa"] += s["n_nomsa"]
            if k % 64 == 0:
                print(f"{k+1}/{len(shards)} shards, {time.time()-t0:.0f}s", flush=True)
    P = pd.concat(pairs, ignore_index=True)
    M2 = pd.concat(m2s, ignore_index=True)
    # make pair_row global
    P["pair_id"] = np.arange(len(P), dtype=np.int64)
    key = P.groupby("og_id").pair_id.min()
    M2["pair_id"] = key.reindex(M2.og_id).to_numpy() + M2.pair_row.to_numpy()
    # pair rows are per-OG in the same order as written (concat preserves order within an OG)
    P.to_parquet(out / "pairs.parquet", index=False)
    M2.to_parquet(out / "m2_matches.parquet", index=False)
    rt = time.time() - t0
    extra = {"stage": "partA_descent_pairs", "orthology": a.orthology, "tau": a.tau, "cap": a.cap,
             "unc_max": a.unc_max, "ceiling": {f"{k[0]}|{k[1]}": v for k, v in G["ceiling"].table.items()},
             "ceiling_source": a.ceiling_source, "a_crit": G["crit"], "n_shards": len(shards),
             "n_ogs_processed": stats["n_og"], "n_ogs_no_msa": stats["n_nomsa"], "n_pairs_rows": len(P),
             "n_m2": len(M2), "runtime_s": rt, "workers": a.workers, "control_ogs": sorted(G["controls"])}
    write_sidecar(out / "pairs.parquet", a.snapshot, G["seed"], extra)
    write_sidecar(out / "m2_matches.parquet", a.snapshot, G["seed"], extra)
    print(json.dumps(extra, indent=1, default=str))


if __name__ == "__main__":
    main()
