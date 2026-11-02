#!/usr/bin/env python
"""Artifact cascade on every Part B candidate (DESIGN §5.4; prereg part_b.artifact_cascade_order, states.ted_domain_rule).

Members "used" by a call: gain = P members of the derived clade + A members of the outgroups; loss = A members of
the derived clade + P members of the outgroups; specificity = P/W members of D and S. One column per check;
nothing is dropped (a check without data is `na` and does not fail, but is listed in `artifact_na`).
  a_footprint_plddt   all A used have footprint mean pLDDT >= 90 and frac(<70) <= 0.2       (pass/fail/na)
  a_disorder          no A used has > 20 % of its footprint in a UniProt 'disordered' region
  a_pae_ted           prereg ted_domain_rule: an A whose footprint spans >= 2 TED domains needs PAE (not fetched)
                      -> that A is treated as U; fail if fewer than 2 A remain (gain: from >= 2 lineages not
                      re-derived; counted conservatively on members)
  a_fragment          no member used is a UniProt fragment; all A used have length ratio 0.7-1.3 vs OG median
  a_orthology_o2      >= 80 % of the call's core members (derived clade + supporting members of the 2 nearest
                      supporting outgroup lineages) fall in the modal O2 root HOG (O1 vs O2 agreement)
  a_orthology_panther same members: >= 80 % of PANTHER-annotated ones share the modal PANTHER family (>= 2 annotated)
  a_instability       gain/spec: median rescue hits of the derived-clade P members >= 1/6; loss: median rescue hits
                      of the outgroup P members >= 1/6 (P supported beyond the atlas run); A used are 0/6 after rescue
  a_fragment_present  (D24) no present-side protein (gain: derived-clade P; loss: outgroup P) has length ratio
                      outside 0.7-1.3 vs the OG median or is a UniProt fragment
  a_species_tree_lineages (D24) the supporting outgroup species (A for gain, P for loss) fall in >= 2 species-tree
                      lineages outside the clade's species-tree LCA (sister, then successive outgroups)
  flag_paralog_branch (D24) event directly below a duplication node; o2_paralog_split = O2 keeps the derived copy
                      in its own sub-HOG (no sister-copy member inside the derived clade's deepest common O2 HOG)
  flag_present_plddt_lt90 (D24) some present-side member has footprint pLDDT < 90
  deep                event at Eukaryota/Bacteria/root level (orthology-uncertain; reported, not a fail)
  orthology_split     (D28, revised) fail if < 50 % of the outgroup supporting members (A for gain, P for loss) are
                      reciprocal best hits of a derived-clade member (DIAMOND, rbh_check.py; whole proteomes of the
                      analysis); frac_arch_match = same ordered Pfam domains (where annotated). Separate columns:
                      o2_root_split (O1: > 50 % of supporters in another O2 root HOG) and panther_split (> 50 % of
                      supporters in another PANTHER family, only when all members are annotated; else unknown).
                      Neither separate column enters pass_artifacts.
pass_artifacts = no check fails; pass_artifacts_strict = also no flag.

    python scripts/partB/run_artifacts.py --release O1_r3
"""
from __future__ import annotations

import argparse
import glob
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.dataset as ds

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from pocketevo import recon  # noqa: E402
from run_rescue import member_footprints  # noqa: E402

SNAP = recon.cli_value("--snapshot") or recon.default_snapshot()   # --snapshot / $PE_SNAPSHOT / prereg
CHECKS = ["a_footprint_plddt", "a_disorder", "a_pae_ted", "a_fragment", "a_fragment_present", "a_orthology_o2",
          "a_orthology_panther", "a_instability", "a_species_tree_lineages", "orthology_split"]
FLAGS = ["flag_paralog_branch", "flag_present_plddt_lt90"]   # reported, not fails (D24: flag, do not drop)


def ranges(txt):
    out = []
    for tok in str(txt or "").replace(",", ";").split(";"):
        tok = tok.strip()
        if "-" in tok:
            a, b = tok.split("-")[:2]
            if a.isdigit() and b.isdigit():
                out.append((int(a), int(b)))
    return out


def ted_segments(chop):
    segs = []
    for part in str(chop).split("_"):
        segs += ranges(part)
    return segs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--release", "--sog-run", dest="release", default="O1_r3")
    recon.add_partb_args(ap, skip=("sog_run",))
    ap.add_argument("--tau", type=float, default=0.2)
    ap.add_argument("--results-root", default=str(ROOT / "results/partB"), help="Part B results root (tests: tmp)")
    a = ap.parse_args()
    import time as _time
    _t0 = _time.time()
    RB = Path(a.results_root)
    out = RB / f"{SNAP}/{a.release}"
    bulk = Path(f"/tmp/pe_big/sogs/{SNAP}/{a.release}")
    O2 = a.orthology == "O2"
    oxml = ROOT / ("data/orthology/ancestra_O2/accurate/hogs.orthoxml.gz" if O2 else
                   "data/orthology/ancestra_tier1_v0/accurate/hogs.orthoxml.gz")
    gtree = ROOT / ("data/orthology/ancestra_O2/o2_guide_tree.nwk" if O2 else
                    "data/orthology/ancestra_tier1_v0/tier1_guide_tree.nwk")
    cand = pd.read_parquet(out / "candidates.parquet")
    cand["kind"] = "gain_loss"
    spf = out / f"specificity/tau{a.tau}/specificity_called.parquet"
    spec = pd.read_parquet(spf) if spf.exists() else pd.DataFrame()
    exf = out / f"specificity/tau{a.tau}/specificity_exploratory.parquet"
    spec_called = len(spec) > 0
    if not spec_called and exf.exists():        # none BH-significant: audit the top-5 exploratory pairs instead
        ex = pd.read_parquet(exf)
        spec = ex[ex["exploratory_top5"]].copy()
    trees = {t.tree_id: t for t in recon.iter_orthoxml_trees(
        oxml,
        keep=set(cand["og_id"]) | set(spec.get("og_id", pd.Series(dtype=str))))}
    rows_all = []
    for tau in sorted(set(cand["tau"]) | ({a.tau} if len(spec) else set())):
        sogs = set(cand.loc[cand["tau"] == tau, "sog_id"]) | (set(spec["sog_id"]) if tau == a.tau and len(spec) else set())
        if O2:   # O2: states with the rescue overlays (prep_o2_states.py); A-only rescue -> P instability na
            post = ds.dataset(sorted(glob.glob(f"/tmp/pe_big/partB/{SNAP}/{a.release}/states_resc_t{tau}/shard=*.parquet")),
                              format="parquet").to_table(columns=["sog_id", "accession", "state"]).to_pandas()
            post = post[post["sog_id"].isin(sogs)]
            hits = {}
        else:
            evdir = RB / f"events/{SNAP}/{a.release}/tau{tau}/primary/rescue"
            post = pd.read_parquet(evdir / "sog_states/bucket=all.parquet")
            rh = pd.read_parquet(evdir / "rescue_hits.parquet")
            hits = dict(zip(zip(rh["sog_id"], rh["accession"]), rh["n_settings_hit"]))
        bs = ds.dataset(sorted(glob.glob(str(bulk / f"states_t{tau}/shard=*.parquet"))), format="parquet").to_table(
            columns=["sog_id", "accession", "state", "footprint_plddt", "frac_lt70", "len_ratio"]).to_pandas()
        bs = bs[bs["sog_id"].isin(sogs)]
        stt = {(s, x): st for s, x, st in zip(post["sog_id"], post["accession"], post["state"])}
        for s, x, st in zip(bs["sog_id"], bs["accession"], bs["state"]):
            stt.setdefault((s, x), st)
        q = {(s, x): (p, f, l) for s, x, p, f, l in zip(bs["sog_id"], bs["accession"], bs["footprint_plddt"],
                                                     bs["frac_lt70"], bs["len_ratio"])}
        accs_by_sog = bs.groupby("sog_id")["accession"].apply(set).to_dict()
        fpres = member_footprints(sogs, tau, bulk, accs_by_sog)
        rows_all.append((tau, stt, q, hits, fpres))
    allacc = set()
    for _, stt, *_ in rows_all:
        allacc |= {x for _, x in stt}
    u = pd.read_parquet(ROOT / "data/annotations/uniprot_proteins.parquet",
                        columns=["accession", "disordered", "panther", "fragment"])
    u = u[u["accession"].isin(allacc)].set_index("accession")
    dis = {x: ranges(v) for x, v in u["disordered"].items()}
    pan = {x: next((t for t in str(v).split(";") if t and ":" not in t), None) for x, v in u["panther"].items()}
    frag = u["fragment"].to_dict()
    ted = pd.read_parquet(ROOT / "data/annotations/ted_domains.parquet", columns=["accession", "chopping"])
    ted = ted[ted["accession"].isin(allacc)]
    tdom = ted.groupby("accession")["chopping"].apply(lambda s: [ted_segments(c) for c in s]).to_dict()
    o2t = pd.read_parquet(ROOT / "data/orthology/ancestra_O2/tables/protein_hogs.parquet",
                          columns=["protein", "root", "path"])
    o2t = o2t[o2t["protein"].isin(allacc)]
    o2 = dict(zip(o2t["protein"], o2t["root"]))
    o2path = {x: list(p) for x, p in zip(o2t["protein"], o2t["path"])}
    sptree = recon.SpeciesTree.from_file(gtree)
    rbf = out / "rbh_outgroups.parquet"
    rbh = {}
    if rbf.exists():
        rbt = pd.read_parquet(rbf)
        rbh = {(float(k[0]), k[1], int(k[2]), k[3]): g for k, g in rbt.groupby(["tau", "sog_id", "node", "type"])}

    def st_lineages(clade_sp, out_sp):
        """distinct species-tree lineages (outside the clade LCA) of the outgroup species."""
        tx = [recon.taxid_from_species_name(x) for x in clade_sp]
        c = sptree.lca_of_taxids(tx)
        T = sptree.tree
        anc = [c] + T.ancestors(c)
        lin = set()
        for x in out_sp:
            v = sptree.leaf_of_taxid.get(recon.taxid_from_species_name(x))
            if v is None or (T.leaf_lo[c] <= T.leaf_lo[v] < T.leaf_hi[c]):
                continue                 # inside the clade's species-tree LCA: not a species-tree outgroup
            prev = c
            for a_ in anc[1:]:
                if T.leaf_lo[a_] <= T.leaf_lo[v] < T.leaf_hi[a_]:
                    for ch in T.children(a_):
                        if ch != prev and T.leaf_lo[ch] <= T.leaf_lo[v] < T.leaf_hi[ch]:
                            lin.add(int(ch))
                    break
                prev = a_
        return len(lin)

    def evaluate(tau, sog, og, node, kind_type, d_node=None, s_node=None):
        R = next(r for r in rows_all if r[0] == tau)
        _, stt, qq, hits, fpres = R
        t = trees[og]
        names = t.leaf_names()
        if kind_type == "specificity":
            v = node
            lo, hi = t.leaf_lo[v], t.leaf_hi[v]
            dlo, dhi = t.leaf_lo[d_node], t.leaf_hi[d_node]
            used = [x for x in names[lo:hi] if stt.get((sog, x)) in ("P", "W")]
            a_used, p_core = [], [x for x in names[dlo:dhi] if stt.get((sog, x)) in ("P", "W")]
        else:
            lo, hi = t.leaf_lo[node], t.leaf_hi[node]
            inside, outside = names[lo:hi], names[:lo] + names[hi:]
            if kind_type == "gain":
                p_core = [x for x in inside if stt.get((sog, x)) == "P"]
                a_used = [x for x in outside if stt.get((sog, x)) == "A"]
                used = p_core + a_used
            else:
                a_used = [x for x in inside if stt.get((sog, x)) == "A"]
                p_core = [x for x in outside if stt.get((sog, x)) == "P"]
                used = a_used + p_core
        r = {}
        na = []
        # footprint pLDDT / frac<70 of A used
        if a_used:
            vals = [qq.get((sog, x), (np.nan, np.nan, np.nan)) for x in a_used]
            pl = np.array([v[0] for v in vals], float); fr = np.array([v[1] for v in vals], float)
            r["min_A_footprint_plddt"] = float(np.nanmin(pl)) if np.isfinite(pl).any() else np.nan
            r["a_footprint_plddt"] = "pass" if np.all((pl >= 90) & (fr <= 0.2)) else "fail"
        else:
            r["a_footprint_plddt"] = "na"; na.append("a_footprint_plddt")
        # disorder of A footprints
        nd = 0; nchk = 0
        for x in a_used:
            f = fpres.get((sog, x))
            if not f or x not in dis:
                continue
            nchk += 1
            inr = sum(any(lo_ <= p <= hi_ for lo_, hi_ in dis[x]) for p in f)
            nd += inr / len(f) > 0.2
        r["n_A_disordered"] = nd
        r["a_disorder"] = ("fail" if nd else "pass") if nchk else "na"
        if not nchk:
            na.append("a_disorder")
        # TED multi-domain footprints of A (no PAE fetched -> A treated as U)
        multi = 0; nt = 0
        for x in a_used:
            f = fpres.get((sog, x))
            if not f or x not in tdom:
                continue
            nt += 1
            doms = {i for i, segs in enumerate(tdom[x]) for p in f if any(lo_ <= p <= hi_ for lo_, hi_ in segs)}
            multi += len(doms) >= 2
        r["n_A_multi_ted"] = multi
        if kind_type == "specificity":
            r["a_pae_ted"] = "na"; na.append("a_pae_ted")
        else:
            r["a_pae_ted"] = "fail" if len(a_used) - multi < 2 else "pass"
            if not nt:
                na.append("a_pae_ted(no TED labels)")
        # fragments / length ratio
        fr_used = [x for x in used if frag.get(x) is True]
        lr = np.array([qq.get((sog, x), (np.nan,) * 3)[2] for x in a_used], float)
        r["n_uniprot_fragment"] = len(fr_used)
        r["a_fragment"] = "fail" if fr_used or (np.isfinite(lr).any() and np.any((lr < 0.7) | (lr > 1.3))) else "pass"
        # orthology: O2 root agreement + PANTHER, on the derived clade + the supporting members of the two
        # nearest outgroup lineages that support the call (the members the call rests on); all-used values kept
        near = list(used)
        if kind_type in ("gain", "loss"):
            from pocketevo.events import TreeIndex
            sup_state = "A" if kind_type == "gain" else "P"
            core = [x for x in names[lo:hi] if stt.get((sog, x)) == ("P" if kind_type == "gain" else "A")]
            sup, k_ = [], 0
            for (_a, _x, llo, lhi) in TreeIndex(t).lineages(node):
                m_ = [x for x in names[llo:lhi] if stt.get((sog, x)) == sup_state]
                if m_:
                    sup += m_; k_ += 1
                if k_ == 2:
                    break
            near = core + sup
        ra = [o2[x] for x in used if x in o2]
        r["o2_modal_root_frac_all_used"] = Counter(ra).most_common(1)[0][1] / len(ra) if ra else np.nan
        roots = [o2[x] for x in near if x in o2]
        if len(roots) >= 2:
            top = Counter(roots).most_common(1)[0][1] / len(roots)
            r["o2_modal_root_frac"] = top; r["n_used_in_O2"] = len(roots)
            r["a_orthology_o2"] = "pass" if top >= 0.8 else "fail"
        else:
            r["a_orthology_o2"] = "na"; na.append("a_orthology_o2")
        pf = [pan[x] for x in near if pan.get(x)]
        if len(pf) >= 2:
            top = Counter(pf).most_common(1)[0][1] / len(pf)
            r["panther_modal_frac"] = top; r["n_panther"] = len(pf)
            r["a_orthology_panther"] = "pass" if top >= 0.8 else "fail"
        else:
            r["a_orthology_panther"] = "na"; na.append("a_orthology_panther")
        # pocket-finder instability (rescue settings)
        hh = [hits[(sog, x)] for x in p_core if (sog, x) in hits]
        if hh:
            r["median_P_rescue_hits"] = float(np.median(hh))
            r["a_instability"] = "pass" if np.median(hh) >= 1 else "fail"
        else:
            r["a_instability"] = "na"; na.append("a_instability")
        # D24: present side (gain: derived P = p_core; loss: outgroup P = p_core)
        lrp = np.array([qq.get((sog, x), (np.nan,) * 3)[2] for x in p_core], float)
        plp = np.array([qq.get((sog, x), (np.nan,) * 3)[0] for x in p_core], float)
        nfp = int(np.sum((lrp < 0.7) | (lrp > 1.3)) + sum(frag.get(x) is True for x in p_core))
        r["n_present_len_outlier_or_fragment"] = nfp
        r["a_fragment_present"] = ("fail" if nfp else "pass") if len(p_core) else "na"
        r["present_min_footprint_plddt"] = float(np.nanmin(plp)) if np.isfinite(plp).any() else np.nan
        r["flag_present_plddt_lt90"] = bool(np.isfinite(plp).any() and np.nanmin(plp) < 90)
        if kind_type in ("gain", "loss"):
            sp_of = lambda xs: {t.species[t.leaves[names.index(x)]] for x in xs}  # noqa: E731
            clade_sp = {t.species[t.leaves[k]] for k in range(lo, hi)}
            out_support = a_used if kind_type == "gain" else p_core
            nl = st_lineages(clade_sp, sp_of(out_support))
            r["species_tree_out_lineages"] = nl
            r["a_species_tree_lineages"] = "pass" if nl >= 2 else "fail"
            par = int(t.parent[node])
            r["flag_paralog_branch"] = bool(t.kind[par] == recon.DUP)
            if r["flag_paralog_branch"]:
                inside_set = set(names[lo:hi])
                sis = [x for x in names[t.leaf_lo[par]:t.leaf_hi[par]] if x not in inside_set]
                ps = [o2path[x] for x in inside_set if x in o2path]
                if len(ps) >= 2 and any(x in o2path for x in sis):
                    common = set(ps[0]).intersection(*map(set, ps[1:]))
                    deepest = max(common, key=lambda h: ps[0].index(h)) if common else None
                    r["o2_paralog_split"] = bool(deepest is not None and
                                                 not any(deepest in o2path[x] for x in sis if x in o2path))
                else:
                    r["o2_paralog_split"] = None
        else:
            r["a_species_tree_lineages"] = "na"
            r["flag_paralog_branch"] = bool(t.kind[node] == recon.DUP)
        # D28 (revised): orthology_split = sequence RBH of the outgroup supporters (rbh_check.py); the O2-root and
        # PANTHER versions are kept as separate columns (o2_root_split, panther_split)
        r["orthology_split"] = "unknown"; r["o2_root_split"] = "unknown"; r["panther_split"] = "unknown"
        if kind_type in ("gain", "loss"):
            der = p_core if kind_type == "gain" else a_used
            outs = a_used if kind_type == "gain" else p_core
            dr = [o2[x] for x in der if x in o2]
            orr = [o2[x] for x in outs if x in o2]
            if dr and orr and not O2:
                mode = Counter(dr).most_common(1)[0][0]
                fr_ = sum(x != mode for x in orr) / len(orr)
                r["frac_outgroup_other_O2_root"] = fr_
                r["o2_root_split"] = "fail" if fr_ > 0.5 else "pass"
            if der and outs and all(pan.get(x) for x in der + outs):
                mode = Counter(pan[x] for x in der).most_common(1)[0][0]
                fr_ = sum(pan[x] != mode for x in outs) / len(outs)
                r["frac_outgroup_other_panther"] = fr_
                r["panther_split"] = "fail" if fr_ > 0.5 else "pass"
            rb = rbh.get((float(tau), sog, int(node), kind_type))
            if rb is not None and "testable" in rb:
                rb = rb[rb["testable"].astype(bool)]
            if rb is not None and len(rb):
                fr_ = float(rb["rbh"].astype(bool).mean())
                r["frac_outgroup_rbh"] = fr_
                r["n_rbh_tested"] = int(len(rb))
                am = rb["arch_match"].dropna()
                r["frac_arch_match"] = float(am.astype(bool).mean()) if len(am) else np.nan
                r["n_arch_compared"] = int(len(am))
                r["orthology_split"] = "fail" if fr_ < 0.5 else "pass"
        r["n_members_used"] = len(used); r["n_A_used"] = len(a_used)
        r["artifact_na"] = ";".join(na)
        r["pass_artifacts"] = not any(r[c] == "fail" for c in CHECKS)
        r["artifact_fails"] = ";".join(c for c in CHECKS if r[c] == "fail")
        r["pass_artifacts_strict"] = r["pass_artifacts"] and not any(r.get(f_, False) for f_ in FLAGS)
        return r

    res = pd.DataFrame([evaluate(r.tau, r.sog_id, r.og_id, int(r.node), r.type) for r in cand.itertuples()])
    # idempotent: a re-run replaces the cascade columns of an earlier run instead of duplicating them
    cand = pd.concat([cand.drop(columns=[c for c in res.columns if c in cand.columns]).reset_index(drop=True), res], axis=1)
    f = out / "candidates.parquet"
    cand.to_parquet(f, index=False)
    files = [f]
    if len(spec):
        sres = [evaluate(a.tau, r.sog_id, r.og_id, int(r.node), "specificity", d_node=int(r.d_node))
                for r in spec.itertuples()]
        sres = pd.DataFrame(sres)
        spec = pd.concat([spec.drop(columns=[c for c in sres.columns if c in spec.columns]).reset_index(drop=True), sres],
                         axis=1)
        spec["fdr_supported"] = spec_called
        if not spec_called:
            spec["label"] = "exploratory, not significant"
        f2 = out / "specificity_candidates.parquet"
        spec.to_parquet(f2, index=False)
        files.append(f2)
    summ = {}
    for name, df in (("gain_loss", cand), ("specificity", spec)):
        if not len(df):
            continue
        g = df[df.get("primary", pd.Series(True, index=df.index)).astype(bool)] if name == "gain_loss" else df
        by = g.groupby("type") if name == "gain_loss" else [("specificity", g)]
        for typ, x in by:
            summ[typ] = {"n": int(len(x)), "pass_artifacts": int(x["pass_artifacts"].sum()),
                         "pass_artifacts_strict": int(x["pass_artifacts_strict"].sum()),
                         **{f_: int(x[f_].sum()) for f_ in FLAGS},
                         "o2_root_split": x["o2_root_split"].value_counts().to_dict() if "o2_root_split" in x else {},
                         "panther_split": x["panther_split"].value_counts().to_dict() if "panther_split" in x else {},
                         "o2_paralog_split": x["o2_paralog_split"].value_counts(dropna=False).to_dict()
                         if "o2_paralog_split" in x else {},
                         **{c: x[c].value_counts().to_dict() for c in CHECKS}}
    f3 = out / "artifact_summary.json"
    f3.write_text(json.dumps(dict(summ, runtime_s=_time.time() - _t0), indent=1, default=int))
    files.append(f3)
    # append to CANDIDATES.md (written by write_candidates.py)
    mdp = out / "CANDIDATES.md"
    if mdp.exists():
        L = ["\n## Artifact cascade (run_artifacts.py; one column per check, nothing dropped)\n",
             "Checks (fail -> pass_artifacts false): " + ", ".join(CHECKS) + ". Flags (reported, not fails): " +
             ", ".join(FLAGS) + ". PAE was not fetched: the TED multi-domain rule turns such A calls into U. "
             "`na` = no data for that check (listed in `artifact_na`).\n",
             "| type | n | pass_artifacts | pass_strict | paralog_branch | O2 paralog split (T/F/na) | present pLDDT<90 | "
             "fail counts |", "|---|---|---|---|---|---|---|---|"]
        for typ, v in summ.items():
            if not isinstance(v, dict):
                continue
            fails = {c: v[c].get("fail", 0) for c in CHECKS if v[c].get("fail", 0)}
            o2s = v.get("o2_paralog_split", {})
            L.append(f"| {typ} | {v['n']} | {v['pass_artifacts']} | {v['pass_artifacts_strict']} | "
                     f"{v.get('flag_paralog_branch', 0)} | {o2s.get(True, 0)}/{o2s.get(False, 0)}/{o2s.get(None, 0)} | "
                     f"{v.get('flag_present_plddt_lt90', 0)} | {json.dumps(fails)} |")
        if len(spec) and not spec_called:
            sf = pd.read_csv(out / f"specificity/tau{a.tau}/funnel_real_vs_null.csv")
            L.append("\n## Specificity change (residue-first, prereg; empirical matched-depth null + BH)\n")
            L.append("| " + " | ".join(sf.columns) + " |\n|" + "---|" * len(sf.columns))
            for _, rr in sf.iterrows():
                L.append("| " + " | ".join(str(x) for x in rr.values) + " |")
            L.append("\n**No pair passes BH (negative result, 3 null seeds).** Pairs with >= 2 BLOSUM62<0 clade-consistent "
                     "substitutions AND a fingerprint shift above the noise floor are listed as **exploratory, not "
                     "significant** (`specificity/tau0.2/specificity_exploratory.parquet`); the top 5 (smallest q) are "
                     "in `specificity_candidates.parquet` with artifact columns for the case builder's audit. No claims.\n")
            cols_ = ["sog_id", "level", "n_d", "n_s", "n_species_d", "n_subst", "substitutions", "shift", "p_empirical",
                     "q_bh", "pass_artifacts", "artifact_fails"]
            L.append("| " + " | ".join(cols_) + " |\n|" + "---|" * len(cols_))
            for _, rr in spec.iterrows():
                L.append("| " + " | ".join(f"{rr[c]:.3g}" if isinstance(rr[c], float) else str(rr[c]) for c in cols_) + " |")
        L.append("\nPrimary calls only (τ 0.2, 2:1); specificity = BH-called pairs (or the exploratory top 5). Per-candidate columns are in "
                 "`candidates.parquet` / `specificity_candidates.parquet`.\n")
        mdp.write_text(mdp.read_text() + "\n".join(L) + "\n")
        files.append(mdp)
    for x in files:
        recon.write_sidecar(x, SNAP, 0, {"script": "scripts/partB/run_artifacts.py", "argv": sys.argv[1:]})
    print(json.dumps(summ, indent=1, default=int))


if __name__ == "__main__":
    main()
