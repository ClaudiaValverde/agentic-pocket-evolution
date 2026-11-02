#!/usr/bin/env python
"""M4: evidence filters + rule-based case selection -> frozen case list + cases/<case_id>/case.yaml (schema case-v1).

Re-reads the Part B candidate table every run (it grows: specificity rows, artifact-cascade columns), so the same
command works for O1 and, with other paths, for O2. Only derived counts and identifiers (PDB/UniProt/ChEMBL/M-CSA/
PubMed IDs) are written; no BioLiP rows (no licence) and no raw annotation records.

  .venv/bin/python scripts/cases/build_cases.py --orthology O1 --version v1            # O1 defaults
  .venv/bin/python scripts/cases/build_cases.py --orthology O2 --sog-release O2_r1 \
      --candidates results/partB/<snap>/O2/candidates.parquet --events-base <O2 pre-rescue events dir> \
      --hog-trees results/partB/hog_trees/O2/hog_trees.parquet --version v1

Rules: src/pocketevo/evidence.py (prereg values from config/prereg.yaml; M4 constants documented there and in
docs/EVIDENCE.md).
"""
from __future__ import annotations

import argparse
import glob
import json
import sys
import time
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.dataset as ds
import yaml

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from pocketevo import evidence as E  # noqa: E402
from pocketevo import recon  # noqa: E402

SNAP = "20261030T085136Z"
F3_INDEPENDENT = {"gain": ("conservation", "rescue_consistency", "source_consistency", "cofolding"),
                  "loss": ("conservation", "rescue_consistency", "source_consistency", "cofolding"),
                  # specificity is called on footprint substitutions, so conservation is not independent there
                  "specificity": ("rescue_consistency", "source_consistency", "cofolding"),
                  "convergence": ("source_consistency", "cofolding")}


def log(t0, msg):
    print(f"[{time.time() - t0:6.0f}s] {msg}", flush=True)


# --------------------------------------------------------------------------------------------- loaders
def load_states(sog_ids, events_base, sogs_web):
    """Post-rescue member states for candidate SOGs (rescue/sog_states), else the SOG release states."""
    parts = []
    rs = Path(events_base) / "rescue/sog_states/bucket=all.parquet" if events_base else None
    have = set()
    if rs and rs.exists():
        st = pd.read_parquet(rs)
        st = st[st.sog_id.isin(sog_ids)]
        st["state_source"] = "post-rescue"
        parts.append(st)
        have = set(st.sog_id)
    miss = set(sog_ids) - have
    if miss:
        f = sorted(glob.glob(str(Path(sogs_web) / "sog_states/bucket=*.parquet")))
        t = ds.dataset(f, format="parquet").to_table(
            filter=ds.field("sog_id").isin(sorted(miss)),
            columns=["og_id", "sog_id", "taxid", "accession", "state", "footprint_plddt", "reason"]).to_pandas()
        t["state_source"] = "sog-release (no rescue)"
        parts.append(t)
    st = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
    return st


def load_msas(og_ids, bulk):
    ogs = pd.read_parquet(Path(bulk) / "inputs/ogs.parquet", columns=["og_id", "shard"])
    ogs = ogs[ogs.og_id.isin(og_ids)]
    out = {}
    for sh, g in ogs.groupby("shard"):
        f = Path(bulk) / f"msa/shard={int(sh):04d}.parquet"
        m = pd.read_parquet(f)
        m = m[m.og_id.isin(set(g.og_id))]
        for og, mm in m.groupby("og_id"):
            out[og] = dict(zip(mm.accession, mm.aligned))
    return out


def load_annotations(accs, tune_excl):
    ann = ROOT / "data/annotations"
    accs = sorted({x for x in accs if isinstance(x, str)})
    bl = ds.dataset(str(ann / "biolip_sites.parquet")).to_table(
        filter=ds.field("accession").isin(accs),
        columns=["accession", "pdb_id", "chain", "ligand", "ligand_class", "residues_uniprot", "pubmed"]).to_pandas()
    up = ds.dataset(str(ann / "uniprot_site_features.parquet")).to_table(
        filter=ds.field("uniprot_acc").isin(accs)).to_pandas()
    up = up[up.experimental]
    mc = ds.dataset(str(ann / "mcsa_residues.parquet")).to_table(
        filter=ds.field("accession").isin(accs), columns=["accession", "mcsa_id", "resid", "source"]).to_pandas()
    ch = pd.read_parquet(ann / "chembl_by_accession.parquet")
    ch = ch[ch.accession.isin(accs)]
    prot = pd.read_parquet(ann / "uniprot_proteins.parquet",
                           columns=["accession", "seq_crc64", "length", "gene_primary", "protein_name", "organism_id",
                                    "reviewed", "panther", "pfam"])
    prot = prot[prot.accession.isin(accs)]
    # prereg tune_test.annotation_disjointness is applied per SOG at hit level (exclusion_keys / filter_hits)
    n_ex = {"rule": "per SOG (prereg literal); OG-wide reported as sensitivity"}
    return dict(biolip={a: g.to_dict("records") for a, g in bl.groupby("accession")},
                uniprot={a: g.to_dict("records") for a, g in up.groupby("uniprot_acc")},
                mcsa={a: g.to_dict("records") for a, g in mc.groupby("accession")},
                chembl=ch.set_index("accession").to_dict("index"),
                prot=prot.set_index("accession").to_dict("index"), n_excluded_tune=n_ex)


def load_manual(cdir: Path) -> list:
    """Reviewed post-hoc literature items from cases/<id>/manual.yaml (prereg process.case_provenance)."""
    f = Path(cdir) / "manual.yaml"
    if not f.exists():
        return []
    m = yaml.safe_load(f.read_text()) or {}
    return [x for x in m.get("literature_evidence", []) if x.get("about_change") and x.get("review_ok")]


def fetch_families(accs, cache: Path, batch=100, timeout=30) -> dict:
    """PANTHER/Pfam xrefs from UniProt REST for accessions outside the local annotation tables (O2 context species).
    Cached per accession (tsv rows); network failure leaves the accession out (-> side-family check pending)."""
    cache.mkdir(parents=True, exist_ok=True)
    out, todo = {}, []
    for x in sorted(set(accs)):
        f = cache / f"fam_{x}.json"
        if f.exists():
            out[x] = tuple(json.loads(f.read_text()))
        else:
            todo.append(x)
    for i in range(0, len(todo), batch):
        q = "+OR+".join(f"accession:{x}" for x in todo[i:i + batch])
        url = (f"https://rest.uniprot.org/uniprotkb/search?query={q}&fields=accession,xref_panther,xref_pfam"
               f"&format=tsv&size={batch}")
        try:
            with urllib.request.urlopen(url, timeout=timeout) as r:
                lines = r.read().decode().splitlines()[1:]
        except Exception:  # noqa: BLE001
            continue
        for ln in lines:
            acc, pan, pf = (ln.split("\t") + ["", ""])[:3]
            v = (pan.strip(";"), pf.strip(";"))
            out[acc] = v
            (cache / f"fam_{acc}.json").write_text(json.dumps(v))
    return out


def load_manual_artifacts(cdir: Path) -> list:
    f = Path(cdir) / "manual.yaml"
    if not f.exists():
        return []
    m = yaml.safe_load(f.read_text()) or {}
    return [dict(name=x["name"], status=x["status"], detail=x.get("detail", ""), source="manual.yaml")
            for x in m.get("artifact_checks", [])]


def rel_cases(ctx) -> Path:
    return ctx["cases_root"]


# --------------------------------------------------------------------------------------------- per-member evidence
def used_in_tuning(tune, pairs):
    """Restrict the tau-tuning annotation list (all annotations of tune-half OGs) to those that actually entered a
    tune-half benchmark pair (tau_benchmark_pairs.parquet; any kind/subset, detected or not: conservative).
    An annotation is used if its site key (site_id, or 'undetected:<ann_id>') is one side of a tune-half pair and its
    ligand key is among that pair's ligand tokens. (D25a: the tuning list alone over-excludes.)"""
    if tune is None or not len(tune):
        return tune
    if pairs is None:
        return tune
    pr = pairs[pairs.half == "tune"]
    tok = {}
    for a, b, lig in zip(pr.site_a, pr.site_b, pr.ligands.fillna("")):
        t = {x for part in str(lig).split(";") for x in part.split("|") if x}
        tok.setdefault(a, set()).update(t)
        tok.setdefault(b, set()).update(t)
    key = [sid if isinstance(sid, str) else f"undetected:{aid}" for sid, aid in zip(tune.site_id, tune.ann_id)]
    m = [lk in tok.get(k, ()) for k, lk in zip(key, tune.ligand_key)]
    return tune[np.array(m, dtype=bool)]


def exclusion_keys(tune, site_sog, tau_col="sog_t0.2"):
    """prereg tune_test.annotation_disjointness. `tune` = annotations actually used to tune tau (used_in_tuning).
    Returns (per_sog, og_wide) key sets: per_sog = {(sog_id, accession, source, ligand_key)} for used annotations whose
    matched site is a member of that SOG (prereg literal: 'excluded ... for the same SOG'); og_wide =
    {(og_id, accession, source, ligand_key)} for every used annotation of the OG (conservative sensitivity)."""
    if tune is None or not len(tune):
        return set(), set()
    # OG ids come from site_sog (the orthology being evaluated), so O1-tuned annotations map onto O2 SOGs
    t = tune.dropna(subset=["site_id"]).drop(columns=["og_id"]).merge(
        site_sog[["og_id", "site_id", tau_col]], on="site_id", how="left")
    acc_og = site_sog.drop_duplicates("site_id").assign(acc=lambda x: x.site_id.str.split(":").str[0])
    acc_og = acc_og.groupby("acc").og_id.apply(set).to_dict()
    og_wide = {(o, a, s_, k) for a, s_, k in zip(tune.accession, tune.source, tune.ligand_key)
               for o in acc_og.get(a, ())}
    t = t[t[tau_col].notna() & (t[tau_col] >= 0)]
    per_sog = {(f"{o}:g{int(g)}", a, s_, k) for o, g, a, s_, k in
               zip(t.og_id, t[tau_col], t.accession, t.source, t.ligand_key)}
    return per_sog, og_wide


def hit_key(h):
    if h["source"] == "biolip":
        return ("biolip", h["ligand"])
    if h["source"] == "uniprot":
        return ("uniprot", h.get("chebi") or "")
    return ("mcsa", str(h.get("mcsa_id")))


def filter_hits(hits, sog_id, og_id, acc, excl, mode):
    """Drop hits used for tau tuning. mode 'sog' (primary) or 'og' (sensitivity). Returns (kept, n_dropped)."""
    per_sog, og_wide = excl
    keep = []
    for h in hits:
        src, k = hit_key(h)
        key = (sog_id, acc, src, k) if mode == "sog" else (og_id, acc, src, k)
        if key in (per_sog if mode == "sog" else og_wide):
            continue
        keep.append(h)
    return keep, len(hits) - len(keep)


def member_evidence(acc, fp_res, A, rules, members_crc):
    """Site-level hits for one member; only if the model sequence equals the UniProt 2026_03 sequence."""
    p = A["prot"].get(acc)
    if p is None or members_crc.get(acc) != p["seq_crc64"]:
        return None
    fp = set(fp_res)
    return (E.biolip_hits(A["biolip"].get(acc, []), fp, rules) + E.uniprot_hits(A["uniprot"].get(acc, []), fp, rules)
            + E.mcsa_hits(A["mcsa"].get(acc, []), fp, rules))


def chembl_active(acc, A, rules):
    c = A["chembl"].get(acc)
    return bool(c and (c.get("n_cpd_pchembl6_single") or 0) > 0 and (c.get("max_pchembl") or 0) >= rules.pchembl_min)


# --------------------------------------------------------------------------------------------- artifacts
ART_PREFIX = ("art_", "artifact_", "a_", "flag_")   # Part B run_artifacts.py: a_* checks (pass/fail/na), flag_* (D24)


def artifact_status(row) -> tuple[str, list]:
    cols = [c for c in row.index if c.startswith(ART_PREFIX)]
    checks = []
    for c in cols:
        v = row[c]
        if c.startswith("flag_"):
            if v is not None and not pd.isna(v):
                checks.append(dict(name=c, status="na", detail=f"flag (reported, not a fail; D24): {bool(v)}"))
            continue
        if isinstance(v, str) and v in ("pass", "fail", "na", "pending"):
            checks.append(dict(name=c, status=v, detail=""))
        elif isinstance(v, (bool, np.bool_)):
            checks.append(dict(name=c, status="pass" if v else "fail", detail=""))
        elif v is not None and not (isinstance(v, float) and np.isnan(v)):
            checks.append(dict(name=c, status="na", detail=str(v)[:200]))
    if "pass_artifacts" in row.index and row["pass_artifacts"] is not None and not pd.isna(row["pass_artifacts"]):
        st = "pass" if bool(row["pass_artifacts"]) else "fail"
    else:
        st = "pending"
        checks.append(dict(name="artifact_cascade", status="pending",
                           detail="pass_artifacts not yet in the candidate table (Part B agent)"))
    return st, checks


# --------------------------------------------------------------------------------------------- one event candidate
def evaluate_event(row, ctx, rules):
    og, sog, typ = row.og_id, row.sog_id, row.type
    tr = ctx["trees"].get(og)
    st = ctx["states"][ctx["states"].sog_id == sog]
    fp_cols = ctx["fp"].get(sog, [])
    aln = ctx["msa"].get(og, {})
    rec = dict(case_id=f"{typ}_{sog.replace(':', '_')}_n{int(row.node)}_{ctx['orth']}", category=typ, og_id=og,
               sog_id=sog, node=int(row.node), level=row.level, fdr_supported=bool(row.fdr_supported),
               candidate_label=row.label, score_S_partial=float(row.score_S_partial), deep=bool(row.get("deep", False)))
    if tr is None or st.empty:
        rec.update(f1="fail", f2="fail", f3="fail", note="tree or states missing")
        return rec, None
    mask = E.subtree_mask(tr["parent"].values, int(row.node))
    leaves = set(tr.loc[mask & (tr["species"] != ""), "name"])
    st = st.assign(derived=st.accession.isin(leaves))
    d_st, o_st = st[st.derived], st[~st.derived]
    sp = lambda x: str(int(x))  # noqa: E731
    lca = row.get("clade_lca_species_tree")
    lca_taxa = ctx["sptree_taxa"].get(lca)
    lca_members = [(sp(t), s) for t, s in zip(st.taxid, st.state) if lca_taxa is None or int(t) in lca_taxa]
    if typ == "specificity" and "sister_node" in row.index and not pd.isna(row["sister_node"]):
        smask = E.subtree_mask(tr["parent"].values, int(row["sister_node"]))
        o_st = st[st.accession.isin(set(tr.loc[smask & (tr["species"] != ""), "name"]))]
    f1 = E.persistence(list(zip(d_st.taxid.map(sp), d_st.state)), list(zip(o_st.taxid.map(sp), o_st.state)), typ,
                       bool(row.f1_single_event), rules, lca_members)
    # sides
    dstate, ostate = E.change_state(typ)
    if typ == "gain":
        pres, absn = d_st[d_st.state == "P"], o_st[o_st.state != "G"]
    elif typ == "loss":
        pres, absn = o_st[o_st.state == "P"], d_st[d_st.state != "G"]
    else:
        pres, absn = d_st[d_st.state == "P"], o_st[o_st.state == "P"]
    hits_all = {}
    hits, hits_og, n_unmapped, fp_res_of, n_drop, n_drop_og = {}, {}, 0, {}, 0, 0
    for acc in set(pres.accession) | set(absn.accession):
        if acc not in aln:
            n_unmapped += 1
            continue
        fr = list(E.columns_to_residues(aln[acc], fp_cols).values())
        fp_res_of[acc] = fr
        h = member_evidence(acc, fr, ctx["ann"], rules, ctx["crc"])
        if h is None:
            n_unmapped += 1
            continue
        hits[acc], nd = filter_hits(h, sog, og, acc, ctx["excl"], "sog")
        n_drop += nd
        hits_og[acc], nd = filter_hits(h, sog, og, acc, ctx["excl"], "og")
        n_drop_og += nd
        hits_all[acc] = h
    ps = E.side_summary({a: hits[a] for a in pres.accession if a in hits})
    as_ = E.side_summary({a: hits[a] for a in absn.accession if a in hits})
    f2_og = E.experimental_filter(typ, E.side_summary({a: hits_og[a] for a in pres.accession if a in hits_og}),
                                  E.side_summary({a: hits_og[a] for a in absn.accession if a in hits_og}),
                                  len(pres), len(absn), 0)
    # diagnostic only (never used for labels): f2 if the tau-tuning annotations were not excluded
    f2_none = E.experimental_filter(typ, E.side_summary({a: hits_all[a] for a in pres.accession if a in hits_all}),
                                    E.side_summary({a: hits_all[a] for a in absn.accession if a in hits_all}),
                                    len(pres), len(absn), 0)
    n_ch = sum(chembl_active(a, ctx["ann"], rules) for a in pres.accession)
    n_ch_abs = sum(chembl_active(a, ctx["ann"], rules) for a in absn.accession)
    f2 = E.experimental_filter(typ, ps, as_, len(pres), len(absn), n_ch)
    f2["present"], f2["absent"] = ps, as_
    f2["chembl_present_proteins"], f2["chembl_absent_proteins"] = int(n_ch), int(n_ch_abs)
    f2["n_members_numbering_unverified"] = int(n_unmapped)
    f2["n_hits_excluded_tau_tuning"] = int(n_drop)
    f2["sensitivity_og_wide_exclusion"] = dict(status=f2_og["status"], n_hits_excluded=int(n_drop_og))
    # f3
    cons = E.conservation(aln, list(d_st[d_st.state.isin(E.INFORMATIVE)].accession),
                          list(o_st[o_st.state.isin(E.INFORMATIVE)].accession), fp_cols)
    rh = ctx["rescue_hits"]
    rh = rh[rh.sog_id == sog].set_index("accession")["n_settings_hit"] if len(rh) else pd.Series(dtype=float)
    abs_A = absn[absn.state == "A"].accession
    resc = E.rescue_consistency([rh.get(a, np.nan) for a in pres.accession], [rh.get(a, np.nan) for a in abs_A])
    rd = ctx["redetect"]
    src = E.source_consistency([rd.get(s, np.nan) for s in st[st.accession.isin(pres.accession)].get(
        "site_id", pd.Series(dtype=str)).dropna()] if "site_id" in st else [])
    lines = dict(conservation=cons, rescue_consistency=resc, source_consistency=src,
                 cofolding=dict(status="pending", detail="placeholder; filled by the co-folding agent"))
    f3 = E.computational_filter(lines, F3_INDEPENDENT[typ])
    art, art_checks = artifact_status(row)
    if typ in ("gain", "loss"):
        def fam(accs):
            out = []
            for x in accs:
                p = ctx["ann"]["prot"].get(x)
                out.append((p.get("panther"), p.get("pfam")) if p is not None else ctx["fam_extra"].get(x, ("", "")))
            return out
        sfc = E.side_family_check(fam(pres.accession), fam(absn[absn.state == "A"].accession))
        if ctx["side_family_mode"] == "report":
            sfc = dict(sfc, name="orthology_side_family_reported_only",
                       detail=sfc["detail"] + " [reported only; does not gate selection (PI D30)]")
        art_checks.append(sfc)
        if sfc["status"] == "fail" and ctx["side_family_mode"] == "gate":
            art = "fail"
    for m in load_manual_artifacts(rel_cases(ctx) / rec["case_id"]):
        art_checks.append(m)
        if m.get("status") == "fail":
            art = "fail"
    fl, lab = E.label(f1["status"], f2["status"], f3["status"], art)
    nev = E.evidence_count(f2, lines, F3_INDEPENDENT[typ])
    rec.update(f1=f1["status"], f2=f2["status"], f3=f3["status"], f2_strength=f2["strength"],
               f3_passing=";".join(f3["passing"]), artifacts=art, filter_label=fl, label=lab, evidence_count=nev,
               S=float(row.score_S_partial) * nev, n_species_change=f1["n_species_change"],
               clade_shape=f1["shape"], n_out_support=f1["n_out_support"],
               n_present_with_site_evidence=ps["n_members_with_evidence"],
               n_absent_with_site_evidence=as_["n_members_with_evidence"],
               present_accs=";".join(sorted(pres.accession)), absent_accs=";".join(sorted(absn[absn.state == "A"].accession)),
               present_pdb_ids=";".join(ps["pdb_ids"]), present_ligands=";".join(ps["ligands"]),
               f2_ogwide_exclusion=f2_og["status"], f2_diag_no_exclusion=f2_none["status"],
               f2_diag_no_exclusion_strength=f2_none["strength"], n_hits_excluded_tune=int(n_drop),
               artifact_fail_reasons=";".join(f"{c['name']}: {c.get('detail', '')[:120]}" for c in art_checks
                                              if c.get("status") == "fail"),
               cons_n_specific=cons.get("n_specific"), rescue_frac=resc.get("frac_present_majority"))
    detail = dict(f1=f1, f2=f2, f3=f3, lines=lines, art_checks=art_checks, members=st, pres=set(pres.accession),
                  absn=set(absn.accession), hits=hits, fp_res=fp_res_of, rh=rh, row=row)
    return rec, detail


# --------------------------------------------------------------------------------------------- convergence
def evaluate_convergence(row, ctx, rules):
    A = ctx["ann"]
    rec = dict(case_id=f"convergence_{row.site_a.replace(':', '_')}__{row.site_b.replace(':', '_')}_{ctx['orth']}",
               category="convergence", og_id=row.sog_b.split(":")[0] if isinstance(row.sog_b, str) else None,
               sog_id=row.sog_b, fdr_supported=bool(row.passed_primary),
               candidate_label="directional pass (primary)" if row.passed_primary else "random-null-only (sensitivity)")
    sites = ctx["sites"]
    tok, hits = {}, {}
    for side, sid, acc in (("a", row.site_a, row.acc_a), ("b", row.site_b, row.acc_b)):
        s = sites.get(sid)
        h = member_evidence(acc, list(s["residues"]) if s is not None else [], A, rules, ctx["crc_all"]) or []
        hits[side] = h
        tok[side] = E.side_summary({acc: h})["ligand_tokens"]
    shared = sorted(set(tok["a"]) & set(tok["b"]))
    f2 = dict(status="pass" if shared else "fail", strength="site-level shared ligand" if shared else "none",
              detail=f"site-level ligand tokens a {tok['a'] or 'none'}, b {tok['b'] or 'none'}; shared {shared or 'none'}")
    st_b = ctx["sog_P_species"].get(row.sog_b, 0)
    f1s = "pass" if (row.passed_primary and st_b >= rules.f1_min_species) else "fail"
    f1 = dict(status=f1s, detail=f"directional test passed (primary nulls): {bool(row.passed_primary)}; derived SOG "
                                 f"{row.sog_b} P in {st_b} species")
    rd = ctx["redetect"]
    src = E.source_consistency([rd.get(x, np.nan) for x in (row.site_a, row.site_b) if x in rd])
    lines = dict(source_consistency=src, cofolding=dict(status="pending", detail="placeholder"))
    f3 = E.computational_filter(lines, F3_INDEPENDENT["convergence"])
    checks = []
    for side, sid in (("a", row.site_a), ("b", row.site_b)):
        s = sites.get(sid)
        if s is None:
            checks.append(dict(name=f"site_{side}_quality", status="pending", detail="site not in site table"))
            continue
        ok = s["plddt_frac_lt70_max"] <= 0.2 and s["n_res"] > 10
        checks.append(dict(name=f"site_{side}_quality", status="pass" if ok else "fail",
                           detail=f"{sid}: {s['n_res']} residues, pLDDT {s['plddt_mean_rep']:.1f}, frac<70 "
                                  f"{s['plddt_frac_lt70_max']:.2f}, volume {s['rep_site_volume_A3']:.0f} A3, burrow "
                                  f"{s['max_burrow_score']:.2f} (rule: frac<70 <= 0.2 and > 10 residues, as A/D11)"))
    art = "fail" if any(c["status"] == "fail" for c in checks) else (
        "pending" if any(c["status"] == "pending" for c in checks) else "pass")
    fl, lab = E.label(f1s, f2["status"], f3["status"], art)
    nev = E.evidence_count(dict(status=f2["status"], present={"n_biolip": int(bool(shared))}), lines,
                           F3_INDEPENDENT["convergence"])
    rec.update(f1=f1s, f2=f2["status"], f3=f3["status"], f2_strength=f2["strength"], f3_passing=";".join(f3["passing"]),
               artifacts=art, filter_label=fl, label=lab, evidence_count=nev, S=float(st_b) * nev,
               score_S_partial=float(st_b), n_species_change=st_b, present_ligands=";".join(shared), delta=row.delta,
               p=row.p)
    return rec, dict(f1=f1, f2=f2, f3=f3, lines=lines, art_checks=checks, row=row, hits=hits)


# --------------------------------------------------------------------------------------------- identifier checks
def verify_ids(kind, ids, cache_dir, timeout=8):
    """Resolve PDB / UniProt IDs against their APIs; cached. Returns {id: {resolved, title}}."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    out = {}
    for i in ids:
        cf = cache_dir / f"{kind}_{i}.json"
        if cf.exists():
            out[i] = json.loads(cf.read_text()); continue
        url = (f"https://data.rcsb.org/rest/v1/core/entry/{i}" if kind == "pdb" else
               f"https://rest.uniprot.org/uniprotkb/{i}.json?fields=accession,protein_name,organism_name")
        try:
            with urllib.request.urlopen(url, timeout=timeout) as r:
                j = json.loads(r.read())
            if kind == "pdb":
                rec = dict(resolved=True, title=j.get("struct", {}).get("title", ""),
                           pubmed=j.get("rcsb_primary_citation", {}).get("pdbx_database_id_pub_med"))
            else:
                pn = j.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "")
                rec = dict(resolved=j.get("primaryAccession") == i, title=pn,
                           organism=j.get("organism", {}).get("scientificName", ""))
        except urllib.error.HTTPError as e:
            rec = dict(resolved=False, error=f"HTTP {e.code}") if e.code == 404 else None
        except Exception as e:  # noqa: BLE001  network failure != unresolvable
            rec = None
            out[i] = dict(resolved=None, error=str(e)[:100])
            continue
        if rec is None:
            out[i] = dict(resolved=None, error="http error"); continue
        rec["retrieved_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        cf.write_text(json.dumps(rec))
        out[i] = rec
    return out


# --------------------------------------------------------------------------------------------- case.yaml
def clean(o):
    if isinstance(o, dict):
        return {str(k): clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple, set)):
        return [clean(v) for v in (sorted(o) if isinstance(o, set) else o)]
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return None if np.isnan(o) else float(o)
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, np.ndarray):
        return clean(o.tolist())
    if isinstance(o, float) and np.isnan(o):
        return None
    return o


def write_case(cdir: Path, rec, det, ctx, args, selection):
    cdir.mkdir(parents=True, exist_ok=True)
    (cdir / "evidence").mkdir(exist_ok=True)
    cat = rec["category"]
    prot = ctx["ann"]["prot"]
    proteins, evid, lit, pdb_all = [], [], [], set()
    if cat in ("gain", "loss", "specificity", "control", "duplication_gain_exploratory") and det is not None and "members" in det:
        m = det["members"]
        for r in m.itertuples():
            if r.state == "G":  # gene-absent rows are counted in the notes, not listed
                continue
            role = "derived" if r.derived else ("sister" if cat == "specificity" else "outgroup")
            p = prot.get(r.accession, {})
            hs = det["hits"].get(r.accession, [])
            side = "present" if r.accession in det["pres"] else ("absent" if r.accession in det["absn"] else "other")
            pre = ctx.get("pre_states", {}).get((rec["sog_id"], r.accession))
            proteins.append(dict(accession=r.accession, taxid=int(r.taxid), species=str(r.taxid), role=role,
                                 side=side, state=r.state,
                                 **({"state_pre_filter": pre} if pre is not None else {}), footprint_plddt=None if pd.isna(r.footprint_plddt) else float(r.footprint_plddt),
                                 reason=r.reason or "", gene=p.get("gene_primary") or "",
                                 protein_name=(p.get("protein_name") or "")[:80],
                                 rescue_settings_hit=None if r.accession not in det["rh"].index else int(det["rh"][r.accession]),
                                 residues=" ".join(f"A:{x}" for x in det["fp_res"].get(r.accession, [])) or None,
                                 n_site_evidence=len([h for h in hs if h.get("counted", True)])))
            for h in hs:
                side = "present" if r.accession in det["pres"] else "absent"
                about = side == "present" and cat != "control"
                if h["source"] == "biolip":
                    pdb_all.add(h["pdb_id"])
                    evid.append(dict(type="pdb", id=h["pdb_id"], about_change=about,
                                     detail=f"{r.accession} ({side} side): ligand {h['ligand']} ({h['ligand_class']}) "
                                            f"chain {h['chain']}, {h['overlap']:.2f} of BioLiP2 site residues in footprint"))
                elif h["source"] == "uniprot":
                    evid.append(dict(type="uniprot", id=r.accession, about_change=about,
                                     detail=f"{h['type']} {h['ligand']} {h['chebi']} (experimental ECO), "
                                            f"{h['overlap']:.2f} in footprint; refs {';'.join(h['refs'][:6])}"))
                    for ref in h["refs"]:
                        if ref.startswith("PubMed:"):
                            lit.append(dict(pmid=ref.split(":")[1], note=f"UniProt {r.accession} {h['type']}"))
                elif h["source"] == "mcsa":
                    evid.append(dict(type="mcsa", id=str(h["mcsa_id"]), about_change=about and h["counted"],
                                     detail=f"{r.accession}: M-CSA {h['origin']} catalytic residues, {h['overlap']:.2f} "
                                            f"in footprint" + ("" if h["counted"] else " (homology-transferred, not counted)")))
            if chembl_active(r.accession, ctx["ann"], ctx["rules"]):
                c = ctx["ann"]["chembl"][r.accession]
                evid.append(dict(type="chembl", id=str(c.get("single_protein_targets", "")).split(";")[0], about_change=False,
                                 detail=f"{r.accession}: {c['n_cpd_pchembl6_single']} compounds pChEMBL>=6 "
                                        f"(max {c['max_pchembl']}), protein-level only"))
    elif cat == "convergence":
        row = det["row"]
        for side, acc, sid in (("a", row.acc_a, row.site_a), ("b", row.acc_b, row.site_b)):
            p = prot.get(acc, {})
            proteins.append(dict(accession=acc, role="derived" if side == "b" else "outgroup", state="P", site_id=sid,
                                 gene=p.get("gene_primary") or "", protein_name=(p.get("protein_name") or "")[:80]))
            for h in det["hits"][side]:
                if h["source"] == "biolip":
                    pdb_all.add(h["pdb_id"])
                evid.append(dict(type=h["source"] if h["source"] != "biolip" else "pdb",
                                 id=h.get("pdb_id") or acc, about_change=False, detail=json.dumps(clean(h))[:200]))
    # identifier verification (capped; network failure -> unverified, a 404 -> fails f2 per prereg)
    vid = {}
    if args.verify_ids:
        vid["pdb"] = verify_ids("pdb", sorted(pdb_all)[:args.max_verify], ctx["id_cache"])
        vid["uniprot"] = verify_ids("uniprot", sorted({p["accession"] for p in proteins})[:args.max_verify],
                                    ctx["id_cache"])
    for e in evid:
        e["filter"] = "experimental"
    for x in load_manual(ctx["cases_root"] / rec["case_id"]):
        evid.append(dict(type="pmid", id=str(x["pmid"]), about_change=True, filter="experimental", source="literature",
                         post_hoc=True, side=x.get("side"), detail=f"{x.get('finding', '')} | quote: \"{x.get('quote', '')}\"",
                         author_agent=x.get("author_agent"), reviewed_by=x.get("reviewed_by"),
                         retrieved_at=x.get("retrieved_at")))
        lit.append(dict(pmid=str(x["pmid"]), note=f"post hoc literature f2 line (D25b): {x.get('finding', '')[:120]}"))
    unresolved = [i for k in vid for i, v in vid[k].items() if v.get("resolved") is False]
    for e in evid:
        v = vid.get("pdb" if e["type"] == "pdb" else "uniprot", {}).get(e["id"])
        e["verified"] = None if v is None else v.get("resolved")
    for k in ("pdb",):
        for i, v in vid.get(k, {}).items():
            if v.get("pubmed"):
                lit.append(dict(pmid=str(v["pubmed"]), note=f"primary citation of PDB {i}: {v.get('title', '')[:80]}"))
    f2 = det["f2"] if det else {"status": rec.get("f2"), "detail": ""}
    f2_status = f2["status"]
    if unresolved and f2_status == "pass":
        f2_status = "fail"
    filters = [dict(name="persistence", status=det["f1"]["status"], detail=det["f1"]["detail"]),
               dict(name="experimental", status=f2_status,
                    detail=f2["detail"] + (f"; unresolvable identifiers {unresolved}" if unresolved else "")),
               dict(name="computational", status=det["f3"]["status"], detail=det["f3"]["detail"])]
    lines = [dict(name=f"f3:{k}", status=v["status"], detail=v.get("detail", "")) for k, v in det["lines"].items()]
    uniq_lit, seen = [], set()
    for x in lit:
        if x["pmid"] not in seen:
            seen.add(x["pmid"]); uniq_lit.append(x)
    title = {"gain": "Pocket gain", "loss": "Pocket loss (audit-level, not FDR-calibrated)",
             "specificity": "Pocket specificity change", "convergence": "Convergent pocket pair (artifact-prone)",
             "control": "False-loss negative control",
             "duplication_gain_exploratory": "Duplication-associated gain (exploratory, not preregistered)"}[cat]
    lvl = rec.get("level") or ""
    summary = {"gain": f"SOG {rec['sog_id']}: pocket present in the {lvl} clade (node {rec.get('node')}) and absent-observed "
                       f"in outgroup lineages.",
               "loss": f"SOG {rec['sog_id']}: pocket present in outgroup lineages and absent-observed in the {lvl} clade. "
                       "Losses are NOT FDR-calibrated (D21): audit-level candidate only.",
               "specificity": f"SOG {rec['sog_id']}: pocket conserved in position, footprint substitutions on the {lvl} branch.",
               "convergence": f"Sites {rec['case_id']}: directional convergence candidate from Part A v2. The Part A verifier "
                              "(reviews/verify/verify_partA_v2.md §7) judged all 3 primary passes artifact-like.",
               "control": (f"SOG {rec['sog_id']}: loss called on pre-quality-filter states whose derived-clade A members are "
                           "all poor models (bad_model, D24); it disappears once A_primary + D11 apply. Negative control "
                           "for false losses: filters are expected to fail." if rec.get("control_kind") == "bad_model" else
                           f"SOG {rec['sog_id']}: loss called before rescue, overturned after rescue (A->W). Negative "
                           "control for false losses: filters are expected to fail."),
               "duplication_gain_exploratory": f"SOG {rec['sog_id']}: pocket gained on one copy after a duplication "
                          f"({lvl}). Not preregistered; the null is selection-conditioned. Not counted in headline numbers."}[cat]
    case = dict(schema="case-v1", case_id=rec["case_id"], title=f"{title}: {rec['sog_id']} ({lvl})",
                category=cat, track="discovery", label=rec["label"], toy=False, summary=summary,
                og_id=rec.get("og_id"), sog_id=rec.get("sog_id"), score=rec.get("S"),
                proteins=proteins, filters=filters + lines, artifact_checks=det["art_checks"],
                evidence=evid[:200],
                cofold=[dict(model_id="pending", status="pending", role="derived",
                             note="placeholder; co-folding agent fills models per DESIGN §6 and the D23 concordance "
                                  "rule (ligand pLDDT AND centroid distance or site-residue recall in >= 2/3 seeds)")],
                figures=[],
                literature=uniq_lit[:50],
                headline="convincing_only",
                filter_label_with_literature_post_hoc=rec.get("filter_label_with_literature"),
                notes=(f"filter_label={rec.get('filter_label')} (prereg: convincing = all three filters; supported = f1 + "
                       f"one of f2/f3); final label also needs the artifact cascade ({rec.get('artifacts')}). "
                       f"candidate label: {rec.get('candidate_label')}; fdr_supported={rec.get('fdr_supported')}. "
                       f"selection: {selection}. G (gene-absent) rows not listed: "
                       f"{int(det['members'].state.eq('G').sum()) if det and 'members' in det else 0}; every other member "
                       f"is listed with side (present/absent/other), state and reason"
                       f"{'; state_pre_filter = state in the pre-quality-filter run that called the false loss' if cat == 'control' else ''}. Evidence counts are derived from BioLiP2/UniProt 2026_03/M-CSA/ChEMBL 37 "
                       f"(DATA.md); BioLiP rows are not redistributed, only PDB IDs and counts."),
                provenance=dict(builder="scripts/cases/build_cases.py", orthology=ctx["orth"], snapshot=SNAP,
                                case_list=str(ctx["case_list_rel"]), prereg_sha256=ctx["prereg_sha"],
                                git_commit=ctx["git"], seed=args.seed, built_utc=ctx["now"],
                                n_excluded_tune_annotations=ctx["ann"]["n_excluded_tune"]))
    (cdir / "case.yaml").write_text(yaml.safe_dump(clean(case), sort_keys=False, allow_unicode=True, width=120))
    ev_json = dict(f1=det["f1"], f2={k: v for k, v in f2.items()}, f3=det["f3"], lines=det["lines"],
                   identifier_checks=vid)
    (cdir / "evidence/evidence_summary.json").write_text(json.dumps(clean(ev_json), indent=1, default=str))


# --------------------------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--orthology", default="O1")
    ap.add_argument("--sog-release", default="O1_r3")
    ap.add_argument("--sog-bulk", default=None, help="SOG bulk dir with msa/ inputs/ (default /tmp/pe_big/sogs/<snap>/<release>)")
    ap.add_argument("--candidates", default=f"results/partB/{SNAP}/O1_r3/candidates.parquet")
    ap.add_argument("--false-loss", default=f"results/partB/{SNAP}/O1_r3/false_loss_controls.parquet")
    ap.add_argument("--false-loss-kind", default=None, help="e.g. bad_model (D24); default all rows")
    ap.add_argument("--n-controls", type=int, default=0, help="seeded pick of n controls (0 = producer's examples)")
    ap.add_argument("--dup-gains", default=None, help="duplication_gains_exploratory.parquet (optional, exploratory)")
    ap.add_argument("--n-dup-gains", type=int, default=2)
    ap.add_argument("--pre-filter-states", default=f"results/sogs/{SNAP}/O1_r2/web/sog_states/bucket=*.parquet",
                    help="states of the pre-quality-filter run (false-loss controls: state that drove the call)")
    ap.add_argument("--side-family-mode", choices=["gate", "report"], default="report",
                    help="M4 PANTHER side-family check: gate the artifact status, or report only")
    ap.add_argument("--specificity", default=None, help="specificity candidates parquet (optional)")
    ap.add_argument("--events-base", default=f"results/partB/events/{SNAP}/O1_r3/tau0.2/primary",
                    help="pre-rescue events dir of the primary setting (contains rescue/)")
    ap.add_argument("--orthoxml", default=None, help="build HOG trees from this OrthoXML instead of --hog-trees "
                                                     "(same node indexing as scripts/partB/run_events.py)")
    ap.add_argument("--states-glob", default=None, help="member state shards (e.g. O2 rescued states)")
    ap.add_argument("--rescue-hits", default=None, help="per (sog_id, accession) n_settings_hit table")
    ap.add_argument("--hog-trees", default="results/partB/hog_trees/O1/hog_trees.parquet")
    ap.add_argument("--convergence", default=f"results/partA/{SNAP}/v2/convergence_candidates.parquet")
    ap.add_argument("--tau-tuning", default=f"results/sogs/{SNAP}/O1/tau_tuning_annotations.parquet")
    ap.add_argument("--tau-pairs", default=f"results/sogs/{SNAP}/O1/tau_benchmark_pairs.parquet",
                    help="benchmark pairs actually used to tune tau (restricts --tau-tuning, D25a)")
    ap.add_argument("--redetection", default=f"results/observation/{SNAP}/redetection_sites.parquet")
    ap.add_argument("--sites", default=f"/tmp/pe_big/sites/{SNAP}/sites")
    ap.add_argument("--out", default=f"results/cases/{SNAP}")
    ap.add_argument("--cases-dir", default="cases")
    ap.add_argument("--manual-dir", default="cases", help="where cases/<id>/manual.yaml inputs live")
    ap.add_argument("--version", default="v1")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--pool", choices=["convincing", "supported", "eligible"], default="convincing",
                    help="selection bar: prereg convincing (v1) or the D25c supported co-fold pool (v2)")
    ap.add_argument("--verify-ids", action="store_true")
    ap.add_argument("--max-verify", type=int, default=25)
    ap.add_argument("--no-case-files", action="store_true")
    ap.add_argument("--snapshot", default=SNAP, help="atlas snapshot id; default paths are rewritten to it")
    a = ap.parse_args()
    if a.snapshot != SNAP:
        # rewrite every default path that embeds the default snapshot (explicit arguments are left unchanged)
        for k, v in vars(a).items():
            if isinstance(v, str) and SNAP in v and v == ap.get_default(k):
                setattr(a, k, v.replace(SNAP, a.snapshot))
        globals()["SNAP"] = a.snapshot
    t0 = time.time()
    rel = lambda p: Path(p) if Path(p).is_absolute() else ROOT / p  # noqa: E731
    pr = recon.load_prereg()
    rules = E.rules_from_prereg(pr)
    bulk = Path(a.sog_bulk or f"/tmp/pe_big/sogs/{SNAP}/{a.sog_release}")
    if not (bulk / "msa").exists():
        raise SystemExit(f"SOG bulk dir {bulk} has no msa/ (regenerate the SOG release or pass --sog-bulk)")
    sogs_web = rel(f"results/sogs/{SNAP}/{a.sog_release}/web")

    cand = pd.read_parquet(rel(a.candidates))
    n_all = len(cand)
    cand = cand[cand["primary"]] if "primary" in cand else cand
    cand = cand.drop_duplicates(["sog_id", "type", "node"]).reset_index(drop=True)
    log(t0, f"candidates {n_all} rows, {len(cand)} primary unique: {cand.type.value_counts().to_dict()}; "
            f"columns: {len(cand.columns)} (pass_artifacts: {'pass_artifacts' in cand})")
    conv = pd.read_parquet(rel(a.convergence)) if a.convergence and rel(a.convergence).exists() else pd.DataFrame()
    conv_p = conv[conv.passed_primary] if len(conv) else conv
    fl = pd.read_parquet(rel(a.false_loss)) if a.false_loss and rel(a.false_loss).exists() else pd.DataFrame()
    n_fl_all = len(fl)
    if len(fl) and a.false_loss_kind and "kind" in fl:
        fl = fl[fl.kind == a.false_loss_kind]
    if len(fl) and a.n_controls:
        # seeded control pick (PI 16:00): uniform sample of n rows of the chosen kind, keys sorted for determinism
        fl = fl.sort_values(["sog_id", "node"]).reset_index(drop=True)
        rng = np.random.default_rng(a.seed)
        pick = rng.choice(len(fl), size=min(a.n_controls, len(fl)), replace=False)
        fl = fl.assign(producer_example=fl["negative_control_example"] if "negative_control_example" in fl else False)
        fl["negative_control_example"] = False
        fl.loc[pick, "negative_control_example"] = True
        fl = fl[fl.negative_control_example]
    log(t0, f"false-loss controls: {len(fl)} picked of {n_fl_all} (kind {a.false_loss_kind})")
    spec = pd.DataFrame()
    if a.specificity and rel(a.specificity).exists():
        spec = pd.read_parquet(rel(a.specificity))
        sig = None
        for c in ("bh_significant", "called", "significant"):
            if c in spec:
                sig = spec[c].fillna(False).astype(bool); break
        spec = spec.assign(type="specificity", exploratory=~sig if sig is not None else True)
        if "d_node" in spec:      # Part B specificity: node = parent, d_node = derived clade root
            spec["parent_node"], spec["node"] = spec["node"], spec["d_node"]
        if "n_species_d" in spec and "score_S_partial" not in spec:
            spec["score_S_partial"] = spec["n_species_d"].astype(float)
        for c, v in dict(fdr_supported=False, label="exploratory, not significant", score_S_partial=0.0,
                         f1_single_event=True, level="", clade_lca_species_tree=None).items():
            if c not in spec:
                spec[c] = v
        spec.loc[~spec.exploratory, "label"] = "BH-significant specificity change"
        spec["fdr_supported"] = ~spec.exploratory
        if "n_species_P" in spec and (spec.score_S_partial == 0).all():
            spec["score_S_partial"] = spec["n_species_P"].astype(float)
        log(t0, f"specificity rows {len(spec)} ({int(spec.exploratory.sum())} exploratory)")

    if len(spec):
        cand = pd.concat([cand, spec], ignore_index=True)
    dup = pd.DataFrame()
    if a.dup_gains and rel(a.dup_gains).exists():
        # PI 17:10: at most 2 duplication-associated gains passing the cascade except orthology_split, top by
        # pattern-species count (ties: sog_id); separate category, never headline
        dup = pd.read_parquet(rel(a.dup_gains))
        fails = dup.artifact_fails.fillna("").map(lambda x: {f for f in x.split(";") if f})
        dup = dup[[f <= {"orthology_split"} for f in fails]]
        dup = dup.sort_values(["n_species_pattern", "sog_id"], ascending=[False, True]).head(a.n_dup_gains).copy()
        dup = dup.assign(type="gain", primary=True, fdr_supported=False,
                         label="duplication-associated gain (exploratory, not preregistered; null is selection-conditioned)",
                         score_S_partial=dup.n_species_pattern.astype(float), f1_single_event=True,
                         clade_lca_species_tree=dup.level, dup_exploratory=True,
                         pass_artifacts=True, art_orthology_split_ignored="orthology_split ignored for this category (PI)")
        dup = dup.drop(columns=["artifact_fails"])
        log(t0, f"duplication-associated gains taken: {dup.sog_id.tolist()}")
        cand = pd.concat([cand, dup], ignore_index=True)
    sog_ids = set(cand.sog_id) | set(fl.get("sog_id", [])) | set(conv_p.get("sog_b", pd.Series(dtype=str)).dropna())
    og_ids = {s.split(":")[0] for s in sog_ids}
    if a.states_glob:
        f = sorted(glob.glob(a.states_glob))
        states = ds.dataset(f, format="parquet").to_table(
            filter=ds.field("sog_id").isin(sorted(sog_ids)),
            columns=["og_id", "sog_id", "taxid", "accession", "state", "footprint_plddt", "reason"]).to_pandas()
        states["state_source"] = f"states-glob {a.states_glob}"
    else:
        states = load_states(sog_ids, rel(a.events_base), sogs_web)
    # site ids of P members (for source consistency) from the release states
    f = sorted(glob.glob(str(sogs_web / "sog_states/bucket=*.parquet")))
    sid = ds.dataset(f, format="parquet").to_table(filter=ds.field("sog_id").isin(sorted(sog_ids)),
                                                   columns=["sog_id", "accession", "site_id"]).to_pandas()
    states = states.merge(sid.drop_duplicates(["sog_id", "accession"]), on=["sog_id", "accession"], how="left")
    log(t0, f"states {len(states)} rows for {states.sog_id.nunique()} SOGs")
    sg = ds.dataset(str(sogs_web / "sogs.parquet")).to_table(filter=ds.field("sog_id").isin(sorted(sog_ids)),
                                                             columns=["sog_id", "footprint_columns"]).to_pandas()
    fp = dict(zip(sg.sog_id, [list(map(int, x)) for x in sg.footprint_columns]))
    if a.orthoxml:
        rows = []
        for t in recon.iter_orthoxml_trees(rel(a.orthoxml), keep=og_ids):
            rows.append(pd.DataFrame(dict(og_id=t.tree_id, node=np.arange(t.n_nodes), parent=t.parent,
                                          name=t.node_name, species=t.species, taxid=t.taxid)))
        ht = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=["og_id", "node", "parent", "name",
                                                                                   "species"])
    else:
        ht = pd.read_parquet(rel(a.hog_trees))
    ht = ht[ht.og_id.isin(og_ids)]
    trees = {og: g.sort_values("node").reset_index(drop=True) for og, g in ht.groupby("og_id")}
    msa = load_msas(og_ids, bulk)
    log(t0, f"trees {len(trees)}, msas {len(msa)}, footprints {len(fp)}")
    mem = pd.read_parquet(bulk / "inputs/members.parquet", columns=["og_id", "accession", "seq_crc64"])
    mem = mem[mem.og_id.isin(og_ids)]
    crc = dict(zip(mem.accession, mem.seq_crc64))
    accs = set(states.accession) | set(conv_p.get("acc_a", [])) | set(conv_p.get("acc_b", []))
    tune = pd.read_parquet(rel(a.tau_tuning)) if a.tau_tuning and rel(a.tau_tuning).exists() else None
    if tune is not None:
        # by accession (not og_id): tau was tuned on O1, the OG ids of another orthology differ
        tune = tune[tune.accession.isin(set(states.accession.dropna()))]
    ann = load_annotations(accs, tune)
    ss = pd.read_parquet(rel(f"results/sogs/{SNAP}/{a.sog_release}/site_sog_all_tau.parquet"),
                         columns=["og_id", "site_id", "sog_t0.2"])
    tp = rel(a.tau_pairs) if a.tau_pairs else None
    pairs = pd.read_parquet(tp, columns=["og_id", "site_a", "site_b", "ligands", "half"]) if tp and tp.exists() else None
    n_listed = 0 if tune is None else len(tune)
    tune_used = used_in_tuning(tune, pairs)
    excl = exclusion_keys(tune_used, ss[ss.og_id.isin(og_ids)])  # ss = site_sog of the evaluated release
    ann["n_excluded_tune"]["n_tuning_list_rows"] = int(n_listed)
    ann["n_excluded_tune"]["n_used_in_tune_pairs"] = int(0 if tune_used is None else len(tune_used))
    ann["n_excluded_tune"]["pairs_file"] = str(a.tau_pairs) if pairs is not None else "missing (whole list used)"
    ann["n_excluded_tune"]["n_keys_per_sog"], ann["n_excluded_tune"]["n_keys_og_wide"] = len(excl[0]), len(excl[1])
    log(t0, f"annotations for {len(accs)} accessions; excluded (tau tuning): {ann['n_excluded_tune']}")
    # convergence sites: sequence identity check for numbering uses UniProt crc of the site protein's model
    sites = {}
    need = set(conv_p.get("site_a", [])) | set(conv_p.get("site_b", []))
    if need and Path(a.sites).exists():
        files = sorted(glob.glob(str(Path(a.sites) / "tier[12]/*.parquet")))
        t = ds.dataset(files, format="parquet").to_table(filter=ds.field("site_id").isin(sorted(need))).to_pandas()
        sites = {r["site_id"]: r for r in t.to_dict("records")}
    crc_all = dict(crc)
    for x in set(conv_p.get("acc_a", [])) | set(conv_p.get("acc_b", [])):
        p = ann["prot"].get(x)
        if p is not None and x not in crc_all:
            crc_all[x] = None   # model sequence not checked for non-SOG proteins -> no site-level evidence
    rd = pd.read_parquet(rel(a.redetection)) if a.redetection and rel(a.redetection).exists() else pd.DataFrame()
    rd = rd[rd.comparison == "foldcomp_v4_kept"] if len(rd) else rd
    redetect = dict(zip(rd.site_id, rd.best_jaccard)) if len(rd) else {}
    rh_f = rel(a.rescue_hits) if a.rescue_hits else rel(a.events_base) / "rescue/rescue_hits.parquet"
    rescue_hits = pd.read_parquet(rh_f) if rh_f.exists() else pd.DataFrame()
    spt = recon.SpeciesTree.from_file(sogs_web / "species_tree.nwk")
    T = spt.tree
    sptree_taxa = {nm: {int(T.taxid[v]) for v in T.leaves[T.leaf_lo[i]:T.leaf_hi[i]]}
                   for nm, i in spt.node_of_level.items()}
    sP = states[states.state == "P"].groupby("sog_id").taxid.nunique().to_dict()
    out_dir = rel(a.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    list_path = out_dir / f"case_list_{a.version}{'' if a.orthology == 'O1' else '_' + a.orthology}.parquet"
    missing = {x for x in states.loc[states.state.isin(["P", "A"]), "accession"].dropna() if x not in ann["prot"]}
    fam_extra = fetch_families(missing, Path("/tmp/pe_big/cases/id_cache")) if missing else {}
    log(t0, f"PANTHER/Pfam from UniProt REST for {len(missing)} non-atlas accessions: {len(fam_extra)} resolved")
    pre_states = {}
    if len(fl) and a.pre_filter_states:
        f = sorted(glob.glob(str(rel(a.pre_filter_states))))
        ps = ds.dataset(f, format="parquet").to_table(filter=ds.field("sog_id").isin(sorted(set(fl.sog_id))),
                                                      columns=["sog_id", "accession", "state", "state_s70"]).to_pandas()
        strat = dict(zip(fl.sog_id, fl.get("pre_stratum", pd.Series("primary", index=fl.index))))
        for r in ps.itertuples():
            st_ = r.state_s70 if "70" in str(strat.get(r.sog_id, "")) else r.state
            pre_states[(r.sog_id, r.accession)] = f"{st_} ({strat.get(r.sog_id, 'primary')}, pre-D11)"
    ctx = dict(pre_states=pre_states, side_family_mode=a.side_family_mode, fam_extra=fam_extra, cases_root=rel(a.manual_dir), excl=excl, trees=trees, states=states, fp=fp, msa=msa, ann=ann, crc=crc, crc_all=crc_all, redetect=redetect,
               rescue_hits=rescue_hits, sptree_taxa=sptree_taxa, orth=a.orthology, sites=sites, sog_P_species=sP,
               rules=rules, id_cache=Path(f"/tmp/pe_big/cases/id_cache"), prereg_sha=recon.prereg_sha256(),
               git=recon.git_commit(), now=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
               case_list_rel=list_path.relative_to(ROOT) if list_path.is_relative_to(ROOT) else list_path)

    recs, dets = [], {}
    for row in cand.itertuples(index=False):
        r = pd.Series(row._asdict())
        if r.get("dup_exploratory") is not True:
            r = r.drop(labels=[c for c in ("dup_exploratory", "art_orthology_split_ignored") if c in r.index])
        rec, det = evaluate_event(r, ctx, rules)
        if r.get("dup_exploratory") is True:
            rec["category"] = "duplication_gain_exploratory"
            rec["case_id"] = rec["case_id"].replace("gain_", "dupgain_", 1)
        recs.append(rec); dets[rec["case_id"]] = det
    # convergence: tau tuning annotations are SOG-internal ortholog pairs; not applied (sites are in different OGs)
    for row in conv_p.itertuples(index=False):
        rec, det = evaluate_convergence(pd.Series(row._asdict()), ctx, rules)
        recs.append(rec); dets[rec["case_id"]] = det
    # false-loss negative controls: evaluated as losses on the post-rescue states (node from the control table)
    pe_f = rel(a.events_base) / "rescue/events/events.parquet"
    post_f1 = {}
    if pe_f.exists():
        pe = pd.read_parquet(pe_f, columns=["char", "node", "type", "f1_single_event"])
        pe = pe[pe.type == "loss"]
        post_f1 = dict(zip(zip(pe["char"], pe["node"].astype(int)), pe["f1_single_event"]))
    for row in fl.itertuples(index=False):
        r = pd.Series(row._asdict())
        r["type"], r["fdr_supported"], r["label"], r["score_S_partial"] = "loss", False, "false-loss control", 0.0
        # single tree-consistent event at this node after rescue (False if the loss is no longer an event there)
        k = (r.sog_id, int(r.node))
        r["f1_single_event"], r["clade_lca_species_tree"] = bool(post_f1.get(k, False)), r.get("level")
        rec, det = evaluate_event(r, ctx, rules)
        rec["case_id"] = f"control_falseloss_{r.sog_id.replace(':', '_')}_n{int(r.node)}_{a.orthology}"
        rec["category"] = "control"
        rec["negative_control_example"] = bool(r.negative_control_example)
        rec["control_kind"] = r.get("kind", "rescue_instability")
        rec["control_expectation"] = ("fail (A side is a bad model; loss gone after quality filters)"
                                      if rec["control_kind"] == "bad_model" else "fail (loss overturned by rescue)")
        rec["producer_example"] = bool(r.get("producer_example", False))
        rec["control_ok"] = rec.get("filter_label") == "rejected"
        recs.append(rec); dets[rec["case_id"]] = det
    df = pd.DataFrame(recs)
    ev = df[df.category != "control"].copy()
    # literature f2 line (D25b, post hoc): cases/<id>/manual.yaml, reviewed items only; never used for the headline
    lit = {cid: load_manual(rel(a.manual_dir) / cid) for cid in ev.case_id}
    ev["n_lit_items"] = ev.case_id.map(lambda c: len(lit[c]))
    ev["f2_with_literature"] = [("pass" if (f == "pass" or len(lit[c])) else f) for c, f in zip(ev.case_id, ev.f2)]
    ev["filter_label_with_literature"] = [E.label(x, y, z)[0] for x, y, z in zip(ev.f1, ev.f2_with_literature, ev.f3)]
    if a.pool == "convincing":
        bar_desc = "convincing (all three filters) and artifact cascade pass; pending artifacts -> label pending"
        sel = E.select_cases(ev, rules, a.seed, bar="convincing", label_col="label") if len(ev) else ev
    else:
        # D25c: co-folding pool = prereg case_selection applied to the 'supported' pool (filter_label convincing or
        # supported); a case whose artifact cascade fails is excluded (logged); pending artifacts stay eligible
        bar_desc = ("D25c co-fold pool: filter_label in {convincing, supported}, artifacts != fail; "
                    "headline = convincing_only")
        has_pa = "pass_artifacts" in cand
        ev["artifact_excluded"] = (ev.artifacts != "pass") if has_pa else (ev.artifacts == "fail")
        if has_pa:
            bar_desc = ("D25c co-fold pool: filter_label in {convincing, supported} AND artifacts pass "
                        "(pass_artifacts and M4 orthology_side_family); no backfill; headline = convincing_only")
        expl = ev.get("exploratory", pd.Series(False, index=ev.index)).fillna(False).astype(bool)
        ev["exploratory"] = expl
        sel = E.select_cases(ev, rules, a.seed, bar=("convincing", "supported"), label_col="filter_label",
                             eligible=~ev["artifact_excluded"] & ~expl
                             & ~ev.category.eq("duplication_gain_exploratory")) if len(ev) else ev
        if a.pool == "eligible":
            # staging for scripts/cases/combine_pool.py: every pool-eligible row, selection done there
            sel["selected"] = (~ev["artifact_excluded"] & ~expl & ev.filter_label.isin(["convincing", "supported"])
                               & ~ev.category.eq("duplication_gain_exploratory"))
            sel["selection"] = np.where(sel["selected"], "eligible", "")
        # specificity exploratory rows (not BH-significant): audit rows only, never top-k
        # (only if the artifact cascade passes: the co-fold pool never takes artifact-failing rows)
        sel.loc[sel.category.eq("duplication_gain_exploratory") & ~ev["artifact_excluded"], ["selected", "selection"]] = \
            [True, "exploratory_duplication"]
        m_ex = expl & sel.category.eq("specificity") & ~ev["artifact_excluded"]
        sel.loc[m_ex, ["selected", "selection"]] = [True, "audit_exploratory"]
    keep = [c for c in ("n_lit_items", "f2_with_literature", "filter_label_with_literature", "artifact_excluded",
                        "exploratory")
            if c in ev]
    df = df.merge(ev[["case_id"] + keep], on="case_id", how="left")
    df = df.merge(sel[["case_id", "selected", "selection"]], on="case_id", how="left")
    df["selected"] = df["selected"].fillna(False).astype(bool)
    df["selection"] = df["selection"].fillna("")
    ctl = (df.category == "control") & df.get("negative_control_example", pd.Series(False, index=df.index)).fillna(False).astype(bool)
    df.loc[ctl, ["selected", "selection"]] = [True, "negative_control"]
    if "artifact_excluded" in df:
        df["artifact_excluded"] = df["artifact_excluded"].astype("boolean").fillna(False).astype(bool)
    df["orthology"] = a.orthology
    df["loss_note"] = np.where(df.category == "loss", "not FDR-calibrated (D21)", "")
    for c in df.columns:
        if df[c].dtype == object:
            df[c] = df[c].map(lambda v: v if v is None or isinstance(v, str) else str(v))
    df.to_parquet(list_path, index=False)
    counts = {cat: dict(n=int(len(g)), f1=int((g.f1 == "pass").sum()), f2=int((g.f2 == "pass").sum()),
                        f3=int((g.f3 == "pass").sum()),
                        filter_label=g.filter_label.value_counts().to_dict(), label=g.label.value_counts().to_dict(),
                        selected=int(g.selected.sum()),
                        convincing_filters_only=int((g.filter_label == "convincing").sum()),
                        convincing_headline=int(((g.filter_label == "convincing") & (g.artifacts != "fail")).sum()),
                        convincing_with_literature_post_hoc=int(((g.get("filter_label_with_literature",
                                                                        pd.Series(dtype=str)) == "convincing")
                                                                 & (g.artifacts != "fail")).sum()),
                        artifact_fail=int((g.artifacts == "fail").sum()))
              for cat, g in df.groupby("category")}
    recon.write_sidecar(list_path, SNAP, a.seed, {
        "script": "scripts/cases/build_cases.py", "argv": sys.argv[1:], "orthology": a.orthology,
        "inputs": {k: (str(v), recon.file_sha256(rel(v)) if v and rel(v).is_file() else None) for k, v in
                   dict(candidates=a.candidates, false_loss=a.false_loss, convergence=a.convergence,
                        hog_trees=a.hog_trees, tau_tuning=a.tau_tuning, redetection=a.redetection).items()},
        "sog_release": a.sog_release, "sog_bulk": str(bulk), "events_base": a.events_base,
        "rules": {**rules.__dict__, "m4_constants": {k: getattr(E, k) for k in dir(E) if k.isupper()}},
        "bar": bar_desc, "pool": a.pool, "headline": "convincing_only",
        "excluded_artifacts": ({c: r for c, r in zip(df.loc[df.artifact_excluded, "case_id"],
                                                     df.loc[df.artifact_excluded, "artifact_fail_reasons"])}
                               if "artifact_excluded" in df else {}),
        "headline_rule": "convincing = all three filters AND artifact cascade not failed; literature line post hoc",
        "counts": counts, "n_excluded_tune_annotations": ann["n_excluded_tune"],
        "code_sha256": recon.code_sha256(E.__file__, __file__), "runtime_s": round(time.time() - t0, 1)})
    log(t0, f"wrote {list_path} ({len(df)} rows)")
    print(json.dumps(counts, indent=1, default=str))
    if not a.no_case_files:
        n = 0
        for r in df[df.selected].to_dict("records"):
            write_case(rel(a.cases_dir) / r["case_id"], r, dets[r["case_id"]], ctx, a, r["selection"])
            n += 1
        log(t0, f"wrote {n} case dirs under {a.cases_dir}/")
    log(t0, "done")


if __name__ == "__main__":
    main()
