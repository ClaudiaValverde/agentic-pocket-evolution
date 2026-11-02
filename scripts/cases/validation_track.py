#!/usr/bin/env python
"""Annotation-first validation track (DESIGN §5.3): run the existing pipeline outputs on the control-list families and
report what they do. No tuning: SOGs, states, events, specificity tests and the artifact cascade are read as they are.

For each item and orthology (O1, O2):
- OG of the derived / ancestral proteins (same OG or not);
- site SOG: the SOG of the reference protein's OG whose footprint (mapped to the reference protein) overlaps most with
  the reference's experimental site residues (UniProt ECO:0000269/0007744 binding/active sites + BioLiP2 regular
  ligands); if the reference has no experimental site, the SOG in which it is P with the most P species;
- state pattern (counts and the state of every named protein);
- events at that SOG (gain/loss: Sankoff events, pre/post rescue, called, candidate, artifact cascade; specificity:
  residue-first tests on O1), and whether a called event of the expected type sits on a clade that contains a derived
  protein and no ancestral one (= recovered);
- negative controls: human state for uricase/GULO (must be G, never A); M. leprae vs M. tuberculosis A-call rates
  before/after filters and loss calls involving M. leprae; detection-limited convergence: catalytic-residue coverage by
  BURROW sites and whether any convergence pair links two of the proteins.

Writes results/cases/<snap>/validation_table.{parquet,md} (+ sidecar) and cases/val_<id>/case.yaml (category
validation_known; f2 = the published known answer, DOIs) for the known-answer items.
"""
from __future__ import annotations

import argparse
import glob
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.dataset as ds
import yaml

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts/cases"))
from pocketevo import evidence as E  # noqa: E402
from pocketevo import recon  # noqa: E402
import build_cases as B  # noqa: E402

SNAP = "20261030T085136Z"


def orth_cfg(snap):
    return {
        "O1": dict(bulk=f"/tmp/pe_big/sogs/{snap}/O1_r3", web=f"results/sogs/{snap}/O1_r3/web",
                   pre_states=f"results/sogs/{snap}/O1_r2/web/sog_states/bucket=*.parquet",
                   resc_states=f"results/partB/events/{snap}/O1_r3/tau0.2/primary/rescue/sog_states/bucket=all.parquet",
                   ev_pre=f"results/partB/events/{snap}/O1_r3/tau0.2/primary/events.parquet",
                   ev_post=f"results/partB/events/{snap}/O1_r3/tau0.2/primary/rescue/events/events.parquet",
                   cand=f"results/partB/{snap}/O1_r3/candidates.parquet",
                   spec=f"results/partB/{snap}/O1_r3/specificity/tau0.2/specificity_tested.parquet",
                   trees="results/partB/hog_trees/O1/hog_trees.parquet", orthoxml=None),
        "O2": dict(bulk=f"/tmp/pe_big/sogs/{snap}/O2_r1", web=f"results/sogs/{snap}/O2_r1/web",
                   pre_states=None, resc_states=None,
                   states_glob=f"/tmp/pe_big/partB/{snap}/O2_r1/states_resc_t0.2/shard=*.parquet",
                   ev_pre=None, ev_post=f"results/partB/events/{snap}/O2_r1/tau0.2/primary/events.parquet",
                   cand=f"results/partB/{snap}/O2_r1/candidates.parquet", spec=None, trees=None,
                   orthoxml="data/orthology/ancestra_O2/accurate/hogs.orthoxml.gz"),
    }


def rel(p):
    return None if p is None else (Path(p) if Path(p).is_absolute() else ROOT / p)


def read_states(cfg, sog_ids=None, og_ids=None):
    cols = ["og_id", "sog_id", "taxid", "accession", "state", "footprint_plddt", "reason"]
    if cfg.get("states_glob"):
        files = sorted(glob.glob(cfg["states_glob"]))
    else:
        files = sorted(glob.glob(str(rel(cfg["web"]) / "sog_states/bucket=*.parquet")))
    flt = ds.field("og_id").isin(sorted(og_ids)) if og_ids else ds.field("sog_id").isin(sorted(sog_ids))
    st = ds.dataset(files, format="parquet").to_table(filter=flt, columns=cols).to_pandas()
    if cfg.get("resc_states") and rel(cfg["resc_states"]).exists():
        rs = pd.read_parquet(rel(cfg["resc_states"]))
        rs = rs[rs.og_id.isin(og_ids)] if og_ids else rs[rs.sog_id.isin(sog_ids)]
        keep = st[~st.sog_id.isin(set(rs.sog_id))]
        st = pd.concat([keep.assign(src="release"), rs[cols].assign(src="post-rescue")], ignore_index=True)
    else:
        st["src"] = "release"
    return st


def load_trees(cfg, og_ids):
    if cfg["orthoxml"]:
        out = {}
        for t in recon.iter_orthoxml_trees(rel(cfg["orthoxml"]), keep=set(og_ids)):
            out[t.tree_id] = pd.DataFrame(dict(node=np.arange(t.n_nodes), parent=t.parent, name=t.node_name,
                                               species=t.species, level=t.level))
        return out
    ht = pd.read_parquet(rel(cfg["trees"]))
    ht = ht[ht.og_id.isin(og_ids)]
    return {og: g.sort_values("node").reset_index(drop=True) for og, g in ht.groupby("og_id")}


def clade_accs(tree, node):
    m = E.subtree_mask(tree["parent"].values, int(node))
    return set(tree.loc[m & (tree["species"] != ""), "name"])


def known_site(acc, ann_dir):
    up = ds.dataset(str(ann_dir / "uniprot_site_features.parquet")).to_table(
        filter=ds.field("uniprot_acc") == acc).to_pandas()
    up = up[up.experimental & up.type.isin(["Binding site", "Active site"])]
    res = set()
    for s, e in zip(up.start, up.end):
        if pd.notna(s) and pd.notna(e):
            res.update(range(int(s), int(e) + 1))
    bl = ds.dataset(str(ann_dir / "biolip_sites.parquet")).to_table(
        filter=ds.field("accession") == acc, columns=["ligand_class", "residues_uniprot"]).to_pandas()
    for r in bl[bl.ligand_class == "regular"].residues_uniprot:
        res.update(int(x) for x in r)
    return res


def eval_known(it, orth, cfg, ctx):
    der, anc = it.get("derived", []), it.get("ancestral", [])
    mem = ctx["members"][orth]
    og_of = mem.set_index("accession").og_id.to_dict()
    ogs = {a: og_of.get(a) for a in der + anc}
    ref = it.get("reference") or (der + anc)[0]
    og = og_of.get(ref)
    row = dict(item=it["id"], orthology=orth, role=it["role"], expected=";".join(it.get("expected_type", [])),
               reference=ref, og_id=og, ogs=json.dumps(ogs),
               same_og=bool(og) and all(ogs.get(a) == og for a in der + anc if ogs.get(a)),
               n_named_in_og=sum(ogs.get(a) == og for a in der + anc if og))
    if og is None:
        row.update(sog_id=None, sog_choice="reference not in orthology", recovered=False,
                   outcome="reference protein not in this orthology")
        return row, None
    sg = ctx["sogs"][orth]
    sg = sg[sg.og_id == og]
    msa = B.load_msas({og}, rel(cfg["bulk"])).get(og, {})
    st = read_states(cfg, og_ids={og})
    site = known_site(ref, ROOT / "data/annotations")
    best, how = None, ""
    if site and ref in msa:
        ov = []
        for s, cols in zip(sg.sog_id, sg.footprint_columns):
            r = set(E.columns_to_residues(msa[ref], list(cols)).values())
            ov.append((len(r & site), s))
        ov = [x for x in ov if x[0] > 0]
        if ov:
            n, best = max(ov)
            how = f"max overlap with {len(site)} experimental site residues of {ref}: {n}"
    if best is None:
        p = st[(st.accession == ref) & (st.state == "P")]
        if len(p):
            nsp = st[st.state == "P"].groupby("sog_id").taxid.nunique()
            best = max(set(p.sog_id), key=lambda s: (nsp.get(s, 0), s))
            how = f"{ref} is P; most P species ({nsp.get(best, 0)})" + ("" if site else "; no experimental site")
    if best is None:
        row.update(sog_id=None, sog_choice="no SOG with the reference site / reference never P", recovered=False,
                   outcome="no site SOG found")
        return row, None
    s = st[st.sog_id == best]
    cnt = s.state.value_counts().to_dict()
    named = {a: (s[s.accession == a].state.iloc[0] + (f"({s[s.accession == a].reason.iloc[0]})"
                                                       if s[s.accession == a].reason.iloc[0] else ""))
             if (s.accession == a).any() else "no row" for a in der + anc}
    small = bool(sg.loc[sg.sog_id == best, "small_footprint"].iloc[0]) if "small_footprint" in sg else False
    row.update(sog_id=best, sog_choice=how, small_footprint=small, n_species_P=int(s[s.state == "P"].taxid.nunique()),
               states=json.dumps({k: int(v) for k, v in cnt.items()}), named_states=json.dumps(named),
               state_source=";".join(sorted(set(s.src))))
    # events
    tree = ctx["trees"][orth].get(og)
    evs, rec_ev = [], None
    for phase in ("ev_pre", "ev_post"):
        f = rel(cfg.get(phase))
        if f is None or not f.exists():
            continue
        e = ctx["events"].setdefault((orth, phase), pd.read_parquet(f))
        e = e[e["char"] == best]
        for r in e.itertuples():
            ca = clade_accs(tree, r.node) if tree is not None else set()
            d_in, a_in = bool(ca & set(der)), bool(ca & set(anc))
            evs.append(dict(phase=phase, type=r.type, node=int(r.node), level=r.level, called=bool(r.called),
                            candidate=bool(r.candidate), derived_in_clade=d_in, ancestral_in_clade=a_in,
                            gain_cost=float(r.gain_cost)))
    ev = pd.DataFrame(evs)
    exp = set(it.get("expected_type", []))
    final_phase = "ev_post"
    if len(ev):
        ok = ev[(ev.phase == final_phase) & ev.called & ev.type.isin(exp) & ev.derived_in_clade & ~ev.ancestral_in_clade
                & (ev.gain_cost == 2.0)]
        if len(ok):
            rec_ev = ok.iloc[0].to_dict()
    # specificity tests (O1)
    sp_desc, sp_rec = "not tested (no specificity run for this orthology)", None
    if cfg.get("spec") and rel(cfg["spec"]).exists() and "specificity" in exp:
        t = ctx["spec"].setdefault(orth, pd.read_parquet(rel(cfg["spec"])))
        t = t[t.sog_id == best]
        rows = []
        for r in t.itertuples():
            ca = clade_accs(tree, r.d_node) if tree is not None else set()
            rows.append(dict(d_node=int(r.d_node), level=r.level, residue_pass=bool(r.residue_pass),
                             shift_pass=bool(r.shift_pass), called=bool(r.called), q=float(r.q_bh),
                             derived=bool(ca & set(der)), ancestral=bool(ca & set(anc))))
        tt = pd.DataFrame(rows)
        if len(tt):
            on = tt[tt.derived & ~tt.ancestral]
            sp_desc = (f"{len(tt)} branch tests in SOG; on a derived-only clade: {len(on)} "
                       f"(residue_pass {int(on.residue_pass.sum())}, shift_pass {int(on.shift_pass.sum())}, "
                       f"BH-called {int(on.called.sum())}; min q {on.q.min() if len(on) else float('nan'):.2f})")
            if len(on) and on.called.any():
                sp_rec = on[on.called].iloc[0].to_dict()
        else:
            sp_desc = "SOG not tested for specificity (not P throughout at any node with >= 2 members per side)"
    # artifact cascade of the matching candidate, if any
    cand = ctx["cand"].setdefault(orth, pd.read_parquet(rel(cfg["cand"])))
    cc = cand[(cand.sog_id == best) & cand.get("primary", True)]
    art = ";".join(f"{r.type}@{r.level}: pass_artifacts={r.pass_artifacts} ({r.artifact_fails or ''})"
                   for r in cc.itertuples()) if len(cc) and "pass_artifacts" in cc else "no candidate at this SOG"
    recovered = rec_ev is not None or sp_rec is not None
    called = ev[(ev.phase == final_phase) & ev.called & (ev.gain_cost == 2.0)] if len(ev) else ev
    ev_desc = ("; ".join(f"{r.type}@{r.level} n{r.node} (derived in clade {r.derived_in_clade}, ancestral "
                         f"{r.ancestral_in_clade})" for r in called.itertuples()) if len(called) else
               ("events but none called" if len(ev) else "no Sankoff event (SOG not event-eligible: < 2 P or < 2 A "
                                                          "members, or small footprint)"))
    elig = "event-eligible" if len(ev) else "not event-eligible"
    outcome = ("RECOVERED: " + (f"{rec_ev['type']} at {rec_ev['level']}" if rec_ev else
                                f"specificity at {sp_rec['level']}")) if recovered else "not recovered"
    why = ""
    if not recovered:
        if not row["same_og"]:
            why = "derived and ancestral proteins fall in different OGs (orthology split); "
        if "specificity" in exp:
            why += ("specificity: " + sp_desc + "; ") if cfg.get("spec") else "specificity not run on this orthology; "
        if exp & {"gain", "loss"}:
            why += ("no A members, SOG not event-eligible" if not len(ev) else
                    "no called event of the expected type on a derived-only clade")
    row.update(reason_not_recovered=why.strip("; "), events_called=ev_desc, sog_event_eligibility=elig, specificity=sp_desc, artifact_cascade=art,
               recovered=recovered, outcome=outcome)
    return row, dict(states=s, site=site, events=ev)


def eval_gene_loss(it, orth, cfg, ctx):
    """Gene loss must give G, never A. Human members of the OG (if any) are paralogs (the lost gene has no protein);
    the check is on the reference's site SOG (reference P, most P species): no human A there and no called loss on a
    clade containing a human member."""
    mem = ctx["members"][orth]
    og_of = mem.set_index("accession").og_id.to_dict()
    ref, tax = it["reference"], it["lost_in_taxid"]
    og = og_of.get(ref)
    row = dict(item=it["id"], orthology=orth, role=it["role"], expected=f"G (no gene) for taxid {tax}, never A",
               reference=ref, og_id=og)
    if og is None:
        row.update(recovered=False, outcome="reference not in orthology")
        return row, None
    st = read_states(cfg, og_ids={og})
    hm = mem[(mem.og_id == og) & (mem.taxid == tax)].accession.tolist()
    names = ctx["names"]
    p = st[(st.accession == ref) & (st.state == "P")]
    if not len(p):
        row.update(recovered=False, outcome=f"reference never P; human members in OG: {hm}")
        return row, None
    nsp = st[st.state == "P"].groupby("sog_id").taxid.nunique()
    sog = max(set(p.sog_id), key=lambda x: (nsp.get(x, 0), x))
    h = st[(st.sog_id == sog) & (st.taxid == tax)]
    hc = h.state.value_counts().to_dict()
    ev_bad = 0
    f = rel(cfg.get("ev_post"))
    if f is not None and f.exists() and hm:
        e = ctx["events"].setdefault((orth, "ev_post"), pd.read_parquet(f))
        e = e[(e["char"] == sog) & (e.type == "loss") & e.called]
        tree = ctx["trees"][orth].get(og)
        ev_bad = sum(bool(clade_accs(tree, n) & set(hm)) for n in e.node) if tree is not None else 0
    ok = hc.get("A", 0) == 0 and ev_bad == 0
    hm_desc = ", ".join(f"{a} {names.get(a, '')}" for a in hm) or "none"
    row.update(sog_id=sog, states=json.dumps({k: int(v) for k, v in hc.items()}), n_human_members=len(hm),
               recovered=ok,
               outcome=(("PASS: " if ok else "FAIL: ") + f"site SOG {sog}: human rows {hc or 'none (G)'}; "
                        f"called losses on a human clade {ev_bad}; human members of the OG (paralogs, the lost gene "
                        f"has no protein): {hm_desc}"))
    return row, None


def eval_truncation(it, ctx, snap):
    out = []
    t, r = it["taxid_test"], it["taxid_ref"]
    for name, glob_, col in (("O1_r2 pLDDT>=70 stratum (before quality filters)",
                              f"results/sogs/{snap}/O1_r2/web/sog_states/bucket=*.parquet", "state_s70"),
                             ("O1_r2 A_primary", f"results/sogs/{snap}/O1_r2/web/sog_states/bucket=*.parquet", "state"),
                             ("O1_r3 A_primary + D11", f"results/sogs/{snap}/O1_r3/web/sog_states/bucket=*.parquet",
                              "state")):
        files = sorted(glob.glob(str(ROOT / glob_)))
        x = ds.dataset(files, format="parquet").to_table(filter=ds.field("taxid").isin([t, r]),
                                                         columns=["taxid", col]).to_pandas()
        rate = x.groupby("taxid")[col].apply(lambda v: float((v == "A").mean())).to_dict()
        out.append(f"{name}: A-rate M. leprae {rate.get(t, float('nan')):.4f} vs M. tuberculosis "
                   f"{rate.get(r, float('nan')):.4f}")
    ev = pd.read_parquet(ROOT / f"results/partB/events/{snap}/O1_r3/tau0.2/primary/rescue/events/events.parquet")
    lo = ev[(ev.type == "loss") & ev.clade_species.str.contains(f"_{t}(?:;|$)", regex=True)]
    only = lo[~lo.clade_species.str.contains(f"_{r}(?:;|$)", regex=True)]
    out.append(f"post-rescue loss events whose clade contains M. leprae: {len(lo)} (called {int(lo.called.sum())}); "
               f"containing M. leprae but not M. tuberculosis: {len(only)} (called {int(only.called.sum())})")
    ok = int(only.called.sum()) == 0
    return dict(item=it["id"], orthology="O1", role=it["role"], expected="no M. leprae-only loss calls",
                outcome=("PASS: " if ok else "FAIL: ") + " | ".join(out), recovered=ok)


def eval_detection(it, orth, ctx, snap):
    accs = it["accessions"]
    mc = pd.read_parquet(ROOT / "data/annotations/mcsa_residues.parquet", columns=["accession", "resid", "source"])
    up = pd.read_parquet(ROOT / "data/annotations/uniprot_site_features.parquet",
                         columns=["uniprot_acc", "type", "start", "end", "experimental"])
    files = sorted(glob.glob(f"/tmp/pe_big/sites/{snap}/sites/tier[12]/*.parquet"))
    sites = ds.dataset(files, format="parquet").to_table(filter=ds.field("accession").isin(accs),
                                                        columns=["accession", "residues"]).to_pandas()
    desc = []
    for a in accs:
        cat = set(mc[(mc.accession == a) & (mc.source == "reference")].resid.astype(int))
        if not cat:
            u = up[(up.uniprot_acc == a) & (up.type == "Active site") & up.experimental]
            cat = {int(s) for s in u.start.dropna()}
        cov = set()
        for r in sites[sites.accession == a].residues:
            cov |= set(int(x) for x in r) & cat
        desc.append(f"{a}: {len(cov)}/{len(cat)} catalytic residues in a BURROW site" if cat else f"{a}: no catalytic annotation")
    conv = ctx["conv"]
    pairs = conv[conv.acc_a.isin(accs) & conv.acc_b.isin(accs)]
    desc.append(f"convergence candidates linking two of them: {len(pairs)} (primary passes {int(pairs.passed_primary.sum())})")
    return dict(item=it["id"], orthology="O1 sites / Part A v2", role=it["role"], expected="detection-limited",
                outcome="DETECTION-LIMITED: " + "; ".join(desc), recovered=False)


def write_case(it, rows, cases_dir: Path, ctx):
    cid = f"val_{it['id']}"
    d = cases_dir / cid
    d.mkdir(parents=True, exist_ok=True)
    rec = any(r.get("recovered") for r in rows)
    best = next((r for r in rows if r.get("recovered")), rows[0])
    filters = [
        dict(name="persistence", status="pass" if rec else "fail",
             detail="; ".join(f"{r['orthology']}: {r.get('outcome')} | SOG {r.get('sog_id')} ({r.get('sog_choice', '')}) | "
                              f"called events: {r.get('events_called', '')}" for r in rows)),
        dict(name="experimental", status="pass",
             detail=f"known answer from the literature: {it.get('known_answer', '')} (DOIs below)"),
        dict(name="computational", status="na",
             detail="; ".join(f"{r['orthology']}: specificity {r.get('specificity', '')}; artifact cascade "
                              f"{r.get('artifact_cascade', '')}" for r in rows))]
    evid = [dict(type="doi", id=x, about_change=True, filter="experimental", source="literature",
                 detail=f"known-answer reference for {it['id']}") for x in it.get("dois", [])]
    proteins = [dict(accession=a, role="derived") for a in it.get("derived", [])] + \
               [dict(accession=a, role="ancestral") for a in it.get("ancestral", [])]
    for p in proteins:
        for r in rows:
            ns = json.loads(r.get("named_states") or "{}")
            if p["accession"] in ns:
                p[f"state_{r['orthology']}"] = ns[p["accession"]]
    case = dict(schema="case-v1", case_id=cid, title=f"Validation (known answer): {it['title']}",
                category="validation_known", track="annotation-first",
                label="supported" if rec else "rejected", toy=False,
                summary=f"Known answer: {it.get('known_answer', '')}. Pipeline outcome: "
                        + "; ".join(f"{r['orthology']}: {r.get('outcome')}" for r in rows),
                og_id=best.get("og_id"), sog_id=best.get("sog_id"), proteins=proteins, filters=filters,
                artifact_checks=[dict(name=f"artifact_cascade_{r['orthology']}", status="na",
                                      detail=r.get("artifact_cascade", "")) for r in rows],
                evidence=evid, literature=[dict(doi=x) for x in it.get("dois", [])], cofold=[], figures=[],
                notes="label = whether the unchanged pipeline recovers the known event (supported) or not (rejected); "
                      "f2 is the published known answer, legitimately (annotation-first track). Not a discovery.",
                provenance=dict(builder="scripts/cases/validation_track.py", snapshot=ctx["snap"],
                                prereg_sha256=recon.prereg_sha256(), git_commit=recon.git_commit()))
    (d / "case.yaml").write_text(yaml.safe_dump(B.clean(case), sort_keys=False, allow_unicode=True, width=120))
    return cid


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--snapshot", default=SNAP)
    ap.add_argument("--controls", default="scripts/cases/validation_controls.yaml")
    ap.add_argument("--orthologies", default="O1,O2")
    ap.add_argument("--out", default=None)
    ap.add_argument("--cases-dir", default="cases")
    a = ap.parse_args()
    t0 = time.time()
    snap = a.snapshot
    cfgs = orth_cfg(snap)
    orths = a.orthologies.split(",")
    items = yaml.safe_load((ROOT / a.controls).read_text())["items"]
    ctx = dict(snap=snap, members={}, sogs={}, trees={}, events={}, spec={}, cand={})
    all_accs = {x for it in items for k in ("derived", "ancestral", "accessions") for x in it.get(k, [])}
    for o in orths:
        m = pd.read_parquet(rel(cfgs[o]["bulk"]) / "inputs/members.parquet", columns=["og_id", "accession", "taxid"])
        ctx["members"][o] = m
        og_ids = set(m[m.accession.isin(all_accs)].og_id)
        ctx["members"][o] = m[m.og_id.isin(og_ids)]
        sg = pd.read_parquet(rel(cfgs[o]["web"]) / "sogs.parquet")
        ctx["sogs"][o] = sg[sg.og_id.isin(og_ids)]
        ctx["trees"][o] = load_trees(cfgs[o], og_ids)
        print(f"[{time.time() - t0:.0f}s] {o}: {len(og_ids)} OGs", flush=True)
    u = pd.read_parquet(ROOT / "data/annotations/uniprot_proteins.parquet", columns=["accession", "gene_primary"])
    ctx["names"] = dict(zip(u.accession, u.gene_primary.fillna("")))
    ctx["conv"] = pd.read_parquet(ROOT / f"results/partA/{snap}/v2/convergence_candidates.parquet")
    rows, case_ids = [], []
    for it in items:
        if it["role"] == "known":
            rr = []
            for o in orths:
                r, _ = eval_known(it, o, cfgs[o], ctx)
                rr.append(r)
            rows += rr
            case_ids.append(write_case(it, rr, ROOT / a.cases_dir, ctx))
        elif it["role"] == "negative_gene_loss":
            rows += [eval_gene_loss(it, o, cfgs[o], ctx)[0] for o in orths]
        elif it["role"] == "negative_truncation":
            rows.append(eval_truncation(it, ctx, snap))
        elif it["role"] == "detection_limited":
            rows.append(eval_detection(it, "O1", ctx, snap))
        for it_r in rows[-len(orths):]:
            it_r.setdefault("title", it["title"])
            it_r.setdefault("dois", ";".join(it.get("dois", [])))
        print(f"[{time.time() - t0:.0f}s] {it['id']}: " + " | ".join(str(r.get('outcome'))[:90] for r in rows
                                                                       if r["item"] == it["id"]), flush=True)
    df = pd.DataFrame(rows)
    for c in df.columns:
        if df[c].dtype == object:
            df[c] = df[c].map(lambda v: v if v is None or isinstance(v, str) else str(v))
    out = rel(a.out or f"results/cases/{snap}/validation_table.parquet")
    df.to_parquet(out, index=False)
    known = df[df.role == "known"]
    n_known = known.item.nunique()
    k_any = known.groupby("item").recovered.apply(lambda v: v.astype(str).eq("True").any()).sum()
    k_by = {o: int(known[known.orthology == o].recovered.astype(str).eq("True").sum()) for o in orths}
    neg = df[df.role.str.startswith("negative")]
    md = [f"# Annotation-first validation track (DESIGN §5.3), snapshot {snap}", "",
          f"Pipeline outputs read unchanged (no re-tuning): SOGs O1_r3/O2_r1, states (post-rescue where run), Sankoff events "
          f"(gain cost 2:1, called), O1 residue-first specificity tests, Part B artifact cascade.",
          "", f"**Recall: of {n_known} known events, the pipeline recovers {int(k_any)}** (O1 {k_by.get('O1')}, O2 "
              f"{k_by.get('O2')}). Recovered = a called event of the expected type on a clade that contains a derived "
              f"protein and no ancestral one (specificity: BH-called branch test).",
          f"Negative controls: {int(neg.recovered.astype(str).eq('True').sum())}/{len(neg)} pass. "
          f"Detection-limited convergence controls are reported as such (not counted).", "",
          "| item | orth | role | expected | OG | same OG | site SOG | states | named states | events called | specificity | artifact cascade | outcome | why not recovered | DOIs |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in df.to_dict("records"):
        def cell(k):
            v = r.get(k)
            v = "" if v is None or str(v) in ("nan", "None") else str(v)
            return v.replace("|", "/")[: (2000 if k in ("outcome", "reason_not_recovered") else 300)]
        md.append("| " + " | ".join(cell(k) for k in
                                    ("item", "orthology", "role", "expected", "og_id", "same_og", "sog_id", "states",
                                     "named_states", "events_called", "specificity", "artifact_cascade", "outcome",
                                     "reason_not_recovered", "dois")) + " |")
    md += ["", "Columns `same_og`, `sog_choice`, `ogs` (OG of each named protein) are in validation_table.parquet.",
           f"Case files (category validation_known): {', '.join(case_ids)}."]
    out.with_suffix(".md").write_text("\n".join(md) + "\n")
    recon.write_sidecar(out, snap, 0, {"script": "scripts/cases/validation_track.py", "argv": sys.argv[1:],
                                       "recall": dict(n_known=int(n_known), recovered_any=int(k_any), by_orthology=k_by),
                                       "runtime_s": round(time.time() - t0, 1)})
    print("\n".join(md[:8]))


if __name__ == "__main__":
    main()
