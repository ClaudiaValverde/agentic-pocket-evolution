#!/usr/bin/env python
"""Part B on a SOG state table: Sankoff reconstruction, events, gain/loss candidates, specificity test, funnel,
yield preflight (DESIGN §5.1–5.3; prereg part_b, pi_freeze.yield_preflight).

Input: a SOG directory in the sogs-web-v1 layout (webapp/SCHEMA.md), e.g. results/sogs/<snapshot>/<orthology>/
  sog_states/bucket=*.parquet  (required)   species_tree.nwk (optional; else the tier-1 guide tree)
  footprint_residues.parquet   (optional: sog_id, accession, residues = aligned footprint residues, '-' gap)
  fingerprints.parquet         (optional: sog_id, accession, fingerprint list, site_volume_A3, footprint_plddt,
                                model_source)
  noise floor: --noise-floor <float> or --replicate-cosines <parquet with column cosine> (prereg q95)
plus the orthology OrthoXML whose root-group ids are the og_id values.
Calibrated thresholds: results/partB/sim/thresholds.json (written by run_sim.py; prereg sha256 checked).

Outputs (--out, default results/partB/events/<snapshot>/<orthology>[/<stratum>]):
  events.parquet       every MPR event with clade size >= frozen minimum, all cost settings, rule flags
  candidates.parquet   events with called == True
  node_states.parquet  ancestral state per SOG x gene-tree node (primary costs), unambiguity, conditional costs
  specificity.parquet  sister-pair tests (residue-first, shift, empirical-null p, BH q, call)
  funnel.csv           counts per stage and cost setting;  yield_preflight.json;  manifest.json (+ sidecars)

    python scripts/partB/run_events.py --sogs results/sogs/<snap>/O2 --orthoxml <O2 hogs.orthoxml.gz> \
        --snapshot <snap> --workers 16
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from pocketevo import events, recon, simulate  # noqa: E402

O1 = ROOT / "data/orthology/ancestra_tier1_v0/accurate/hogs.orthoxml.gz"
GUIDE = ROOT / "data/orthology/ancestra_tier1_v0/tier1_guide_tree.nwk"
THRESH = ROOT / "results/partB/sim/thresholds.json"


def load_inputs(sog_dir: Path, orthoxml: str, guide: str | None, states_glob: str | None = None,
                state_col: str = "state", sog_meta: str | None = None, tau: float | None = None,
                min_footprint: int = 4, prefilter: bool = False, counters: dict | None = None):
    """States (+ D8 small_footprint exclusion when sog_meta/tau are given; optional event-eligibility prefilter),
    gene trees of the OGs present, species tree. Counters record SOGs dropped by reason."""
    cnt = counters if counters is not None else {}
    states = events.load_sog_states(sog_dir, states_glob, state_col)
    cnt["sogs_total"] = int(states["sog_id"].nunique())
    fpsize = None
    if sog_meta:
        m = pd.read_parquet(sog_meta)
        if "tau" in m.columns and tau is not None:
            m = m[np.isclose(m["tau"], tau)]
        fpsize = dict(zip(m["sog_id"], m["footprint_size"]))
        small = {s for s, v in fpsize.items() if v < min_footprint}
        for flag in ("small_footprint", "empty_footprint"):      # SOG-agent flags (O1_r2+), D8
            if flag in m.columns:
                small |= set(m.loc[m[flag].astype(bool), "sog_id"])
        n0 = states["sog_id"].nunique()
        states = states[~states["sog_id"].isin(small)]
        cnt["dropped_small_footprint_D8"] = int(n0 - states["sog_id"].nunique())
    if prefilter:
        el = events.eligible_sogs(states)
        n0 = states["sog_id"].nunique()
        states = states[states["sog_id"].isin(el)]
        cnt["dropped_lt2P_or_lt2A"] = int(n0 - states["sog_id"].nunique())
    cnt["sogs_kept"] = int(states["sog_id"].nunique())
    keep = set(states["og_id"].unique())
    trees = {t.tree_id: t for t in recon.iter_orthoxml_trees(orthoxml, keep=keep)}
    nwk = Path(guide) if guide else (sog_dir / "species_tree.nwk" if (sog_dir / "species_tree.nwk").exists()
                                     else GUIDE)
    st = recon.SpeciesTree.from_file(nwk)
    return states, trees, st, nwk, fpsize


def split_by_class(assembled, fpsize, cut: int = 10):
    """Split assembled OG blocks into footprint classes (<= cut -> 'le10', > cut -> 'gt10')."""
    out = {"le10": [], "gt10": []}
    for og, sog_ids, codes, conf in assembled:
        for cls in out:
            idx = [j for j, s in enumerate(sog_ids)
                   if ((fpsize or {}).get(s, cut + 1) <= cut) == (cls == "le10")]
            if idx:
                out[cls].append((og, [sog_ids[j] for j in idx], codes[:, idx], conf[:, idx]))
    return out


def settings_by_class(args_thr: dict, grid, pr):
    """{class: [(costs, Rules)]}; a type with no FDR-passing cell (chosen None) is blocked (no calls)."""
    out = {}
    for cls, path in args_thr.items():
        th = load_thresholds(path)
        st = events.rules_by_cost(th, grid, pr)
        if th:
            fixed = []
            for (g, l), r in st:
                ch = th["chosen"].get(f"{int(g)}:{int(l)}", {})
                for typ in ("gain", "loss"):
                    if ch.get(typ) is None:
                        r = r.with_cal(**{f"{typ}_min_species_cal": 10**6})
                fixed.append(((g, l), r))
            st = fixed
        out[cls] = st
    return out


def load_thresholds(path) -> dict | None:
    p = Path(path)
    if not p.exists():
        return None
    th = json.loads(p.read_text())
    if th.get("prereg_sha256") != recon.prereg_sha256():
        raise SystemExit(f"thresholds {p} were calibrated under a different prereg sha256")
    return th


def load_spec_inputs(sog_dir: Path):
    """Footprint residues / fingerprints + covariates as {sog_id: {accession: value}} (None if absent)."""
    rpath, fpath = sog_dir / "footprint_residues.parquet", sog_dir / "fingerprints.parquet"
    if not rpath.exists():
        return None
    res = pd.read_parquet(rpath)
    resd = {s: dict(zip(g["accession"], g["residues"])) for s, g in res.groupby("sog_id")}
    fpd, covd = {}, {}
    if fpath.exists():
        fp = pd.read_parquet(fpath)
        for s, g in fp.groupby("sog_id"):
            fpd[s] = {a: np.asarray(v, dtype=float) for a, v in zip(g["accession"], g["fingerprint"])}
            covd[s] = {a: np.array([vv, pp, ms], dtype=float) for a, vv, pp, ms in
                       zip(g["accession"], g["site_volume_A3"], g["footprint_plddt"], g["model_source"])}
    return resd, fpd, covd


def specificity(assembled, trees, st, spec_inputs, floor: float, pr: dict, seed: int, permutation: bool = True):
    if spec_inputs is None:
        return None, 0
    resd, fpd, covd = spec_inputs
    sc = pr["part_b"]["specificity_change"]
    min_sub = sc["min_clade_consistent_substitutions"]["value"]
    rows = []
    for og, sog_ids, codes, _conf in assembled:
        t = trees[og]
        for j, s in enumerate(sog_ids):
            rows += events.specificity_tests(t, codes[:, j], resd.get(s), fpd.get(s), covd.get(s), floor,
                                             min_subst=min_sub, species_tree=st, sog_id=s)
    if not rows:
        return pd.DataFrame(), 0
    alpha = pr["meta"]["multiple_testing"]["alpha"]
    called, n_tested = events.call_specificity(rows, alpha=alpha)
    df = pd.DataFrame(called)
    df["p_permutation"] = np.nan
    if not permutation:
        return df, n_tested
    # secondary null: clade-level permutation for gate-passing rows (prereg n = 10000)
    n_perm = int(sc["permutation"]["n"]["value"])
    rng = np.random.default_rng(seed)
    og_of = {s: og for og, sog_ids, _, _ in assembled for s in sog_ids}
    for i in np.flatnonzero((df["residue_pass"] & df["shift_pass"]).to_numpy()):
        r = df.iloc[i]
        t = trees[og_of[r["sog_id"]]]
        names = t.leaf_names()
        lo, hi = t.leaf_lo[r["d_node"]], t.leaf_hi[r["d_node"]]
        vlo, vhi = t.leaf_lo[r["node"]], t.leaf_hi[r["node"]]
        fpm = fpd.get(r["sog_id"], {})
        accs = [x for x in names[vlo:vhi] if x in fpm]
        dset = set(names[lo:hi])
        F = np.array([fpm[x] for x in accs])
        if covd.get(r["sog_id"]):
            F = events.residualize(F, np.array([covd[r["sog_id"]][x] for x in accs]))
        isd = np.array([x in dset for x in accs])
        if isd.sum() >= 2 and (~isd).sum() >= 1:
            df.iat[i, df.columns.get_loc("p_permutation")] = events.permutation_pvalue(F, isd, n_perm, rng)
    return df, n_tested


def spec_funnel_row(spec, n_tested, n_chars, **key):
    n = len(spec) if spec is not None else 0
    return dict(key, type="specificity", n_sogs=n_chars, clade_size_ok=n,
                species=int(spec["residue_pass"].sum()) if n else 0,
                frac=int((spec["residue_pass"] & spec["shift_pass"]).sum()) if n else 0,
                candidate=int(n_tested), called=int(spec["called"].sum()) if n else 0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sogs", required=True)
    ap.add_argument("--orthoxml", default=str(O1))
    ap.add_argument("--guide", default=None)
    ap.add_argument("--thresholds", default=str(THRESH))
    ap.add_argument("--snapshot", required=True)
    ap.add_argument("--orthology", default=None)
    ap.add_argument("--stratum", default="primary", help="primary (pLDDT>=90 A) | secondary (pLDDT>=70)")
    ap.add_argument("--w-mode", default="present")
    ap.add_argument("--noise-floor", type=float, default=None)
    ap.add_argument("--replicate-cosines", default=None)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--workers", type=int, default=recon.n_workers())
    ap.add_argument("--out", default=None)
    ap.add_argument("--states-glob", default=None, help="bulk state shards instead of <sogs>/sog_states")
    ap.add_argument("--state-col", default="state", help="state (primary) | state_s70 (secondary stratum)")
    ap.add_argument("--sog-meta", default=None, help="sogs_all_tau.parquet (footprint_size -> D8 exclusion)")
    ap.add_argument("--tau", type=float, default=None)
    ap.add_argument("--thresholds-le10", default=None, help="class-specific thresholds (footprint <= 10)")
    ap.add_argument("--thresholds-gt10", default=None, help="class-specific thresholds (footprint > 10)")
    ap.add_argument("--prefilter", action="store_true", help="drop SOGs with < 2 P or < 2 A members (counted)")
    ap.add_argument("--all-node-states", action="store_true")
    recon.add_partb_args(ap, skip=("snapshot", "orthology"))
    a = ap.parse_args()
    t0 = time.time()
    sog_dir = Path(a.sogs)
    orth = a.orthology or sog_dir.name
    out = Path(a.out) if a.out else ROOT / "results/partB/events" / a.snapshot / orth / a.stratum
    out.mkdir(parents=True, exist_ok=True)
    pr = recon.load_prereg()
    th = load_thresholds(a.thresholds)
    cnt = {}
    states, trees, st, nwk, fpsize = load_inputs(sog_dir, a.orthoxml, a.guide, a.states_glob, a.state_col,
                                                 a.sog_meta, a.tau, prefilter=a.prefilter, counters=cnt)
    print(f"[{time.time() - t0:.0f}s] loaded {cnt}", flush=True)
    assembled = events.assemble(trees, states, cnt)
    del states
    print(f"[{time.time() - t0:.0f}s] assembled", flush=True)
    costs = recon.costs_from_prereg(pr)
    grid = [costs["primary"]] + [c for c in costs["grid"] if tuple(map(float, c)) != costs["primary"]]
    settings = events.rules_by_cost(th, grid, pr)
    if a.thresholds_le10 and a.thresholds_gt10:
        by = settings_by_class({"le10": a.thresholds_le10, "gt10": a.thresholds_gt10}, grid, pr)
        parts = split_by_class(assembled, fpsize)
        ev, node_states = [], []
        for cls, asm in parts.items():
            e_, n_ = events.run_assembled(asm, trees, by[cls], st, a.w_mode, a.workers, keep_states=True)
            for x in e_:
                x["footprint_class"] = cls
            ev += e_; node_states += n_
        settings = by["gt10"]
        th = {"chosen": {c: load_thresholds(p)["chosen"] for c, p in
                         (("le10", a.thresholds_le10), ("gt10", a.thresholds_gt10))}}
    else:
        ev, node_states = events.run_assembled(assembled, trees, settings, st, a.w_mode, a.workers,
                                               keep_states=True)
    print(f"[{time.time() - t0:.0f}s] reconstructed", flush=True)
    for e in ev:
        e["footprint_size"] = (fpsize or {}).get(e["char"])
        e["tau"] = a.tau
        e["stratum"] = a.stratum
    evdf = pd.DataFrame(ev)
    files = []
    f = out / "events.parquet"; evdf.to_parquet(f, index=False); files.append(f)
    cand = evdf[evdf["called"]] if len(evdf) else evdf
    f = out / "candidates.parquet"; cand.to_parquet(f, index=False); files.append(f)
    f = out / "node_states.parquet"
    nsd = pd.DataFrame(node_states, columns=events.NODE_STATE_COLUMNS)
    if not a.all_node_states and len(evdf):     # keep ancestral states of SOGs with a frozen-rule candidate
        nsd = nsd[nsd["sog_id"].isin(set(evdf.loc[evdf["candidate"], "char"]))]
    nsd.to_parquet(f, index=False); files.append(f)
    n_chars = sum(len(s) for _, s, _, _ in assembled)
    fun = []
    for (g, l), rules in settings:
        sub = [e for e in ev if (e["gain_cost"], e["loss_cost"]) == (g, l)]
        fr = events.funnel(sub, n_chars)
        for typ in ("gain", "loss"):
            fun.append(dict(input="real", cost=f"{int(g)}:{int(l)}", type=typ, n_sogs=n_chars, **fr[typ]))
    # specificity
    floor = a.noise_floor
    q = pr["part_b"]["specificity_change"]["noise_floor_quantile"]["value"]
    rc = a.replicate_cosines or (sog_dir / "replicate_cosines.parquet")
    if floor is None and Path(rc).exists():
        floor = events.noise_floor(pd.read_parquet(rc)["cosine"].to_numpy(), q)
    spec, n_tested = (None, 0)
    if floor is not None:
        spec, n_tested = specificity(assembled, trees, st, load_spec_inputs(sog_dir), floor, pr, a.seed)
    if spec is not None:
        f = out / "specificity.parquet"; spec.to_parquet(f, index=False); files.append(f)
        fun.append(spec_funnel_row(spec, n_tested, n_chars, input="real", cost="na"))
    fdf = pd.DataFrame(fun)
    f = out / "funnel.csv"; fdf.to_csv(f, index=False); files.append(f)
    prim = [e for e in ev if (e["gain_cost"], e["loss_cost"]) == settings[0][0]]
    yp = events.yield_preflight(prim)
    yp["stratum"] = a.stratum
    f = out / "yield_preflight.json"; f.write_text(json.dumps(yp, indent=1)); files.append(f)
    man = dict(stage="partB_events", snapshot_id=a.snapshot, orthology=orth, stratum=a.stratum,
               prereg_sha256=recon.prereg_sha256(), git_commit=recon.git_commit(), seeds=[a.seed],
               runtime_s=time.time() - t0, n_in=cnt["sogs_in"], n_out=n_chars,
               n_dropped_by_reason={k: v for k, v in cnt.items() if k.startswith("dropped") or
                                    k == "og_not_in_orthology"},
               tau=a.tau, state_col=a.state_col, states_glob=a.states_glob,
               member_counters=cnt, thresholds=str(a.thresholds) if th else None,
               thresholds_chosen=th["chosen"] if th else None, noise_floor=floor, w_mode=a.w_mode,
               species_tree=str(nwk), orthoxml=str(a.orthoxml), n_events=len(ev),
               n_called_primary=int(sum(e["called"] for e in prim)),
               code_sha256=recon.code_sha256(recon, events, simulate, __file__))
    f = out / "manifest.json"; f.write_text(json.dumps(man, indent=1, default=str)); files.append(f)
    for f in files:
        recon.write_sidecar(f, a.snapshot, a.seed, {"script": "scripts/partB/run_events.py", "argv": sys.argv[1:]})
    print(fdf.to_string(index=False))
    print(json.dumps(yp))


if __name__ == "__main__":
    main()
