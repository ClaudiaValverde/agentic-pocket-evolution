"""Evidence filters (DESIGN §5.5) and rule-based case selection (prereg `part_b.case_selection`).

Pure functions on small in-memory tables so that every rule is testable on toy data (tests/test_evidence.py).
I/O (reading candidates, SOG states, MSAs and annotation tables) lives in scripts/cases/build_cases.py.

Filters
- f1 persistence: number of derived species carrying the change, clade purity (informative derived members with the
  change state, and the same over all SOG members under the species-tree LCA of the derived clade: clean clade vs
  scattered), single tree-consistent event, outgroup lineages carrying the ancestral state.
- f2 experimental: site-level evidence mapped to each member's footprint residues (UniProt numbering, only where the
  model sequence equals the UniProt sequence): BioLiP2 holo ligands, UniProt binding/active sites with
  experimental ECO codes only, M-CSA reference catalytic residues; ChEMBL/BindingDB activity is protein-level only.
  Evidence must discriminate the two sides of the change; absence of annotation on the absent side is recorded as
  weak, a site-level ligand on the absent side contradicts the call.
- f3 computational: clade-specific footprint conservation vs a background of non-footprint columns, BURROW
  rescue-setting consistency of the present side, v6/v4 model-source re-detection, co-folding (placeholder).

M4 rule choices that the prereg leaves open are module constants below, written before any case result was seen
(recorded in the case-list sidecar and in docs/EVIDENCE.md). The prereg values (f1 species minimum, BioLiP overlap,
ECO codes, pChEMBL, top-k, audit count) are read from config/prereg.yaml, never hardcoded here.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

# --------------------------------------------------------------------------------------------- M4 rule constants
SITE_OVERLAP_MIN_DEFAULT = 0.5      # overwritten by prereg f2_experimental.biolip2_residue_overlap_min
F1_CLADE_PURITY_MIN = 0.8           # = prereg candidate rule "(>= 80 % of informative members)"
F1_MIN_OUT_LINEAGES = 2             # = prereg candidate rule ">= 2 sister/outgroup lineages"
CONS_CLADE_MIN = 0.8                # derived majority residue carried by >= 80 % of derived members
CONS_OUT_MAX = 0.5                  # ... and by < 50 % of outgroup members => clade-specific column
CONS_MIN_COLUMNS = 2                # = prereg specificity "min_clade_consistent_substitutions"
RESCUE_MAJORITY = 3                 # present-side member re-detected in >= 3 of 6 BURROW settings
RESCUE_FRAC_MIN = 0.8               # ... for >= 80 % of present-side members
REDETECT_JACCARD_MIN = 0.3          # v6 vs Foldcomp v4 same-protein re-detection hit (= site merge Jaccard)
SPEC_LIGAND_CLASSES = ("regular",)  # BioLiP ligand classes that count as small-molecule pocket evidence
INFORMATIVE = ("P", "A")            # W/U/G are not informative for the change


@dataclass
class Rules:
    f1_min_species: int = 3
    site_overlap_min: float = SITE_OVERLAP_MIN_DEFAULT
    eco_experimental: tuple = ("ECO:0000269", "ECO:0007744")
    pchembl_min: float = 6.0
    top_k: int = 5
    n_audit: int = 1
    extra: dict = field(default_factory=dict)


def rules_from_prereg(pr: dict) -> Rules:
    s = pr["part_b"] if "part_b" in pr else pr
    ef = s["evidence_filters"]
    f2 = ef["f2_experimental"]
    f1 = ef["f1_persistence"]
    import re
    m = re.search(r">=\s*(\d+)\s*species", f1)
    cs = s["case_selection"]
    return Rules(f1_min_species=int(m.group(1)) if m else 3,
                 site_overlap_min=float(f2["biolip2_residue_overlap_min"]["value"]),
                 eco_experimental=tuple(f2["uniprot_eco"]),
                 pchembl_min=float(f2["pchembl_min"]["value"]),
                 top_k=int(cs["top_k_per_category"]["value"]),
                 n_audit=int(cs["random_audit_cases_per_category"]["value"]))


# --------------------------------------------------------------------------------------------- tree helpers
def subtree_mask(parent: np.ndarray, node: int) -> np.ndarray:
    """Boolean mask of `node` and all its descendants; nodes are in preorder (parent index < child index)."""
    parent = np.asarray(parent)
    m = np.zeros(len(parent), dtype=bool)
    m[node] = True
    for v in range(node + 1, len(parent)):
        p = parent[v]
        if p >= 0 and m[p]:
            m[v] = True
    return m


def columns_to_residues(aligned: str, cols) -> dict:
    """MSA column (0-based) -> residue number (1-based, ungapped position) for non-gap columns of one row."""
    pos = np.cumsum([c not in "-." for c in aligned])
    out = {}
    for c in cols:
        c = int(c)
        if 0 <= c < len(aligned) and aligned[c] not in "-.":
            out[c] = int(pos[c])
    return out


def overlap_frac(site_residues, footprint) -> float:
    s = set(int(x) for x in site_residues)
    if not s:
        return 0.0
    return len(s & set(footprint)) / len(s)


# --------------------------------------------------------------------------------------------- f1 persistence
def change_state(event_type: str) -> tuple[str, str]:
    """(state on the derived side, state on the ancestral/outgroup side)."""
    if event_type == "gain":
        return "P", "A"
    if event_type == "loss":
        return "A", "P"
    return "P", "P"          # specificity: present on both sides


def persistence(derived: list, outgroup: list, event_type: str, single_event: bool, rules: Rules,
                lca_members: list | None = None) -> dict:
    """derived/outgroup/lca_members: lists of (species, state). Returns metrics, status and detail."""
    dstate, ostate = change_state(event_type)
    inf_d = [(s, x) for s, x in derived if x in INFORMATIVE]
    sp_change = sorted({s for s, x in derived if x == dstate})
    purity = (sum(x == dstate for _, x in inf_d) / len(inf_d)) if inf_d else 0.0
    out_lin = sorted({s for s, x in outgroup if x == ostate})
    lca_inf = [(s, x) for s, x in (lca_members or derived) if x in INFORMATIVE]
    lca_purity = (sum(x == dstate for _, x in lca_inf) / len(lca_inf)) if lca_inf else 0.0
    ok = (len(sp_change) >= rules.f1_min_species and purity >= F1_CLADE_PURITY_MIN and bool(single_event)
          and len(out_lin) >= F1_MIN_OUT_LINEAGES)
    shape = "clean clade" if (lca_purity >= F1_CLADE_PURITY_MIN and single_event) else "scattered"
    detail = (f"{dstate} in {len(sp_change)} derived species (purity {purity:.2f} of informative derived members; "
              f"{lca_purity:.2f} over all members under the species-tree LCA: {shape}); single tree-consistent event: "
              f"{bool(single_event)}; outgroup support: {ostate} in {len(out_lin)} outgroup species")
    return dict(status="pass" if ok else "fail", n_species_change=len(sp_change), clade_purity=round(purity, 3),
                lca_purity=round(lca_purity, 3), single_event=bool(single_event), n_out_support=len(out_lin),
                shape=shape, detail=detail)


# --------------------------------------------------------------------------------------------- f2 experimental
def biolip_hits(rows, footprint: set, rules: Rules) -> list:
    """rows: iterable of dicts with pdb_id, chain, ligand, ligand_class, residues_uniprot. Site-level hits."""
    out = []
    for r in rows:
        f = overlap_frac(r["residues_uniprot"], footprint)
        if f >= rules.site_overlap_min:
            out.append(dict(source="biolip", pdb_id=str(r["pdb_id"]).upper(), chain=r["chain"], ligand=r["ligand"],
                            ligand_class=r["ligand_class"], overlap=round(f, 3)))
    return out


def uniprot_hits(rows, footprint: set, rules: Rules) -> list:
    """rows: UniProt features (type, start, end, ligand, chebi, ligand_label, eco_codes, evidence_refs).
    Only experimental ECO codes count. Binding sites are grouped per ligand instance; a group hits if >= overlap
    of its residues lie in the footprint. Active sites hit if the residue is in the footprint."""
    groups = {}
    for r in rows:
        eco = set(str(r.get("eco_codes") or "").split(";"))
        if not eco & set(rules.eco_experimental):
            continue
        if r["type"] not in ("Binding site", "Active site"):
            continue
        try:
            s, e = int(r["start"]), int(r["end"])
        except (TypeError, ValueError):
            continue
        key = (r["type"], r.get("chebi") or r.get("ligand") or "", r.get("ligand_label") or "") \
            if r["type"] == "Binding site" else (r["type"], f"{s}", "")
        g = groups.setdefault(key, dict(res=set(), refs=set(), ligand=r.get("ligand") or "", chebi=r.get("chebi")))
        g["res"].update(range(s, e + 1))
        g["refs"].update(x for x in str(r.get("evidence_refs") or "").split(";") if x)
    out = []
    for (typ, lig, lab), g in groups.items():
        f = overlap_frac(g["res"], footprint)
        if f >= rules.site_overlap_min:
            out.append(dict(source="uniprot", type=typ, ligand=g["ligand"], chebi=g["chebi"] or "",
                            overlap=round(f, 3), refs=sorted(g["refs"])))
    return out


def mcsa_hits(rows, footprint: set, rules: Rules) -> list:
    """rows: M-CSA residues (mcsa_id, resid, source). Reference (curated) residues count; homologue rows are
    sequence-transferred and returned with counted=False."""
    by = {}
    for r in rows:
        by.setdefault((int(r["mcsa_id"]), r["source"]), set()).add(int(r["resid"]))
    out = []
    for (mid, src), res in by.items():
        f = overlap_frac(res, footprint)
        if f >= rules.site_overlap_min:
            out.append(dict(source="mcsa", mcsa_id=mid, origin=src, overlap=round(f, 3), counted=src == "reference"))
    return out


def side_summary(member_hits: dict) -> dict:
    """member_hits: accession -> list of hit dicts (biolip/uniprot/mcsa). Counted site-level hits only."""
    n_mem, pdbs, ligs, lig_regular, n_bl, n_up, n_mc = 0, set(), set(), set(), 0, 0, 0
    for acc, hits in member_hits.items():
        counted = [h for h in hits if h.get("counted", True)]
        if counted:
            n_mem += 1
        for h in counted:
            if h["source"] == "biolip":
                n_bl += 1
                pdbs.add(h["pdb_id"])
                ligs.add(h["ligand"])
                if h["ligand_class"] in SPEC_LIGAND_CLASSES:
                    lig_regular.add(h["ligand"])
            elif h["source"] == "uniprot":
                n_up += 1
                if h.get("chebi"):
                    ligs.add(h["chebi"])
                    if h["type"] == "Binding site":
                        lig_regular.add(h["chebi"])
            elif h["source"] == "mcsa":
                n_mc += 1
    return dict(n_members_with_evidence=n_mem, n_biolip=n_bl, n_uniprot_exp=n_up, n_mcsa_ref=n_mc,
                pdb_ids=sorted(pdbs), ligands=sorted(ligs), ligand_tokens=sorted(lig_regular))


def experimental_filter(event_type: str, present: dict, absent: dict, n_present: int, n_absent: int,
                        chembl_present: int = 0) -> dict:
    """present/absent: side_summary dicts. Gain: present = derived P, absent = outgroup members.
    Loss: present = outgroup P (ancestral state), absent = derived members. Specificity: present = derived,
    absent = sister (both have the pocket; evidence must show different ligand classes)."""
    if event_type == "specificity":
        a, b = set(present["ligand_tokens"]), set(absent["ligand_tokens"])
        if not a or not b:
            return dict(status="fail", strength="none", about_change=False,
                        detail=f"ligand classes at the footprint: derived {sorted(a) or 'none'}, sister "
                               f"{sorted(b) or 'none'}; both sides needed")
        jac = len(a & b) / len(a | b)
        ok = not (a & b)
        return dict(status="pass" if ok else "fail", strength="site-level" if ok else "shared-ligand",
                    about_change=ok, jaccard=round(jac, 3),
                    detail=f"derived ligands {sorted(a)} vs sister {sorted(b)} (Jaccard {jac:.2f}; pass needs disjoint)")
    pres = present["n_members_with_evidence"]
    contra = absent["n_members_with_evidence"]
    side_p = "derived (pocket gained)" if event_type == "gain" else "outgroup/ancestral (pocket present)"
    side_a = "outgroup (pocket absent)" if event_type == "gain" else "derived (pocket lost)"
    if contra:
        return dict(status="fail", strength="contradicted", about_change=False,
                    detail=f"site-level ligand/annotation in {contra} {side_a} member(s) "
                           f"(PDB {absent['pdb_ids'][:5]}): absent-side evidence contradicts the call")
    if pres:
        return dict(status="pass", strength="site-level present side; absent side: absence of annotation (weak)",
                    about_change=True,
                    detail=f"{side_p}: {pres}/{n_present} members with site-level evidence (BioLiP {present['n_biolip']}, "
                           f"UniProt exp {present['n_uniprot_exp']}, M-CSA {present['n_mcsa_ref']}; PDB "
                           f"{present['pdb_ids'][:8]}); {side_a}: 0/{n_absent} annotated (absence of annotation is weak)")
    if chembl_present:
        return dict(status="fail", strength="protein-level only", about_change=False,
                    detail=f"ChEMBL pChEMBL>=6 activity on {chembl_present} {side_p} protein(s): protein-level only, "
                           f"not mapped to the footprint")
    return dict(status="fail", strength="none", about_change=False, detail="no experimental evidence at the footprint")


# --------------------------------------------------------------------------------------------- f3 computational
def conservation(aligned: dict, derived: list, outgroup: list, fp_cols, min_members: int = 2) -> dict:
    """Clade-specific footprint columns: derived majority residue in >= CONS_CLADE_MIN of derived members and in
    < CONS_OUT_MAX of outgroup members. Background = same rate over non-footprint columns where both sides have
    >= min_members residues. Pass: >= CONS_MIN_COLUMNS clade-specific footprint columns and rate_fp > rate_bg."""
    d = [aligned[a] for a in derived if a in aligned]
    o = [aligned[a] for a in outgroup if a in aligned]
    if len(d) < min_members or len(o) < min_members:
        return dict(status="na", detail=f"too few aligned members (derived {len(d)}, outgroup {len(o)})")
    L = min(len(s) for s in d + o)
    D = np.array([list(s[:L]) for s in d])
    O = np.array([list(s[:L]) for s in o])
    fp = set(int(c) for c in fp_cols if int(c) < L)

    def col_stats(c):
        dc = D[:, c][~np.isin(D[:, c], ["-", "."])]
        oc = O[:, c][~np.isin(O[:, c], ["-", "."])]
        if len(dc) < min_members or len(oc) < min_members:
            return None
        vals, cnt = np.unique(dc, return_counts=True)
        maj = vals[np.argmax(cnt)]
        cd = cnt.max() / len(dc)
        co = float(np.mean(oc == maj))
        return cd, co

    spec_fp = n_fp = spec_bg = n_bg = 0
    cd_fp, co_fp = [], []
    for c in range(L):
        st = col_stats(c)
        if st is None:
            continue
        isspec = st[0] >= CONS_CLADE_MIN and st[1] < CONS_OUT_MAX
        if c in fp:
            n_fp += 1; spec_fp += isspec; cd_fp.append(st[0]); co_fp.append(st[1])
        else:
            n_bg += 1; spec_bg += isspec
    if n_fp == 0:
        return dict(status="na", detail="no footprint column with residues on both sides")
    r_fp, r_bg = spec_fp / n_fp, (spec_bg / n_bg if n_bg else 0.0)
    ok = spec_fp >= CONS_MIN_COLUMNS and r_fp > r_bg
    return dict(status="pass" if ok else "fail", n_specific=int(spec_fp), n_fp_cols=int(n_fp),
                rate_fp=round(r_fp, 3), rate_bg=round(r_bg, 3), mean_cons_derived=round(float(np.mean(cd_fp)), 3),
                mean_identity_out=round(float(np.mean(co_fp)), 3),
                detail=f"{spec_fp}/{n_fp} footprint columns clade-specific (rate {r_fp:.2f} vs background {r_bg:.2f}); "
                       f"derived conservation {np.mean(cd_fp):.2f}, outgroup identity to derived majority "
                       f"{np.mean(co_fp):.2f}")


def rescue_consistency(n_hits_present: list, n_hits_absent: list, n_settings: int = 6) -> dict:
    """Present-side members re-detected in >= RESCUE_MAJORITY settings for >= RESCUE_FRAC_MIN of them, and the
    absent-side A members in 0 settings (they are A after rescue by construction; reported for transparency)."""
    hp = [int(x) for x in n_hits_present if x is not None and not np.isnan(x)]
    ha = [int(x) for x in n_hits_absent if x is not None and not np.isnan(x)]
    if not hp:
        return dict(status="na", detail="no rescue run for present-side members")
    frac = float(np.mean([h >= RESCUE_MAJORITY for h in hp]))
    ok = frac >= RESCUE_FRAC_MIN and all(h == 0 for h in ha)
    return dict(status="pass" if ok else "fail", frac_present_majority=round(frac, 3),
                median_hits_present=float(np.median(hp)), max_hits_absent=max(ha) if ha else None,
                detail=f"present side re-detected in >= {RESCUE_MAJORITY}/{n_settings} BURROW settings for "
                       f"{frac:.2f} of {len(hp)} members (median {np.median(hp):.0f}); absent side max "
                       f"{max(ha) if ha else 'n/a'}/{n_settings} over {len(ha)} members")


def source_consistency(jaccards: list) -> dict:
    """v6 atlas vs Foldcomp v4 re-description of present-side sites (only proteins in the re-detection sample)."""
    j = [float(x) for x in jaccards if x is not None and not np.isnan(x)]
    if not j:
        return dict(status="na", detail="no present-side member in the v6/v4 re-detection sample")
    n_hit = sum(x >= REDETECT_JACCARD_MIN for x in j)
    ok = n_hit / len(j) >= RESCUE_FRAC_MIN
    return dict(status="pass" if ok else "fail", n=len(j), n_redetected=n_hit,
                detail=f"site re-detected in the v4 model for {n_hit}/{len(j)} sampled present-side members "
                       f"(Jaccard >= {REDETECT_JACCARD_MIN})")


def computational_filter(lines: dict, independent: tuple) -> dict:
    passing = [k for k in independent if lines.get(k, {}).get("status") == "pass"]
    pend = [k for k, v in lines.items() if v.get("status") == "pending"]
    st = "pass" if passing else ("pending" if pend and not any(lines.get(k, {}).get("status") == "fail"
                                                              for k in independent) else "fail")
    return dict(status=st, passing=passing,
                detail="; ".join(f"{k}: {v.get('status')} ({v.get('detail', '')})" for k, v in lines.items()))


# --------------------------------------------------------------------------------------------- artifact: orthology
def _fams(ids: str, strip_sub: bool = True) -> set:
    out = set()
    for x in str(ids or "").split(";"):
        x = x.strip()
        if x:
            out.add(x.split(":")[0] if strip_sub else x)
    return out


def side_family_check(present: list, absent: list) -> dict:
    """Artifact cascade 'orthology_misassignment' (DESIGN §5.4, PANTHER disagreement) for a gain/loss: present and
    absent members are (panther, pfam) string pairs. Fail if the two sides share no PANTHER family (subfamilies
    collapsed); if either side lacks PANTHER, Pfam is used; if neither is available on a side -> pending."""
    for k, name in ((0, "PANTHER"), (1, "Pfam")):
        a = set().union(*[_fams(m[k]) for m in present]) if present else set()
        b = set().union(*[_fams(m[k]) for m in absent]) if absent else set()
        if a and b:
            shared = sorted(a & b)
            ok = bool(shared)
            return dict(name="orthology_side_family", status="pass" if ok else "fail",
                        detail=f"{name}: present side {sorted(a)[:4]}, absent side {sorted(b)[:4]}; shared "
                               f"{shared[:4] or 'none'}" + ("" if ok else " -> different families on the two sides "
                                                             "(misassigned orthology / merged HOG)"))
    return dict(name="orthology_side_family", status="pending", detail="no PANTHER/Pfam family on one side")


# --------------------------------------------------------------------------------------------- labels, score, selection
def label(f1: str, f2: str, f3: str, artifacts: str = "pass") -> tuple[str, str]:
    """(filter_label, final label). filter_label per prereg: convincing = all three pass, supported = f1 + one of
    f2/f3, else rejected. Final label: an artifact failure rejects; artifacts still pending -> 'pending'."""
    if f1 == "pass" and f2 == "pass" and f3 == "pass":
        fl = "convincing"
    elif f1 == "pass" and (f2 == "pass" or f3 == "pass"):
        fl = "supported"
    else:
        fl = "rejected"
    if artifacts == "fail":
        return fl, "rejected"
    if artifacts == "pending" and fl != "rejected":
        return fl, "pending"
    return fl, fl


def evidence_count(f2_hits: dict, f3_lines: dict, independent: tuple) -> int:
    """Number of distinct passing evidence lines: site-level present-side kinds (BioLiP, UniProt exp, M-CSA) when f2
    passes, plus passing independent computational lines."""
    n = 0
    if f2_hits.get("status") == "pass":
        p = f2_hits.get("present", {})
        n += sum(int(p.get(k, 0) > 0) for k in ("n_biolip", "n_uniprot_exp", "n_mcsa_ref")) or 1
    n += sum(f3_lines.get(k, {}).get("status") == "pass" for k in independent)
    return int(n)


def select_cases(df, rules: Rules, seed: int, bar="convincing", score_col: str = "S",
                 category_col: str = "category", label_col: str = "filter_label", eligible=None):
    """Top-k by S per category among rows whose label meets the bar (a label or a tuple of labels) and, if given,
    the boolean `eligible` column/Series, plus n_audit seeded random passing rows not in the top-k. Ties broken by
    case_id. Returns df with `selected` (bool) and `selection` ('top_k' | 'audit' | '')."""
    out = df.copy()
    out["selected"] = False
    out["selection"] = ""
    rng = np.random.default_rng(seed)
    for cat in sorted(out[category_col].unique()):
        bars = (bar,) if isinstance(bar, str) else tuple(bar)
        m = (out[category_col] == cat) & out[label_col].isin(bars)
        if eligible is not None:
            m &= (out[eligible] if isinstance(eligible, str) else eligible).astype(bool)
        sub = out[m].sort_values([score_col, "case_id"], ascending=[False, True])
        top = sub.index[:rules.top_k]
        out.loc[top, ["selected", "selection"]] = [True, "top_k"]
        rest = sub.index[rules.top_k:]
        if len(rest) and rules.n_audit:
            pick = rng.choice(np.array(sorted(out.loc[rest, "case_id"])), size=min(rules.n_audit, len(rest)),
                              replace=False)
            out.loc[out["case_id"].isin(pick), ["selected", "selection"]] = [True, "audit"]
    return out
