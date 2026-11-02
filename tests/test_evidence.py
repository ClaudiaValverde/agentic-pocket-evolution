"""Known-answer toy tests for src/pocketevo/evidence.py (M4 evidence filters and case selection)."""
import numpy as np
import pandas as pd
import pytest

from pocketevo import evidence as E
from pocketevo import recon

R = E.Rules()


def test_rules_from_prereg_reads_frozen_values():
    r = E.rules_from_prereg(recon.load_prereg())
    assert r.f1_min_species == 3 and r.site_overlap_min == 0.5 and r.top_k == 5 and r.n_audit == 1
    assert set(r.eco_experimental) == {"ECO:0000269", "ECO:0007744"} and r.pchembl_min == 6


def test_subtree_mask_preorder():
    #        0
    #      1   4
    #     2 3   5
    parent = np.array([-1, 0, 1, 1, 0, 4])
    assert E.subtree_mask(parent, 1).tolist() == [False, True, True, True, False, False]
    assert E.subtree_mask(parent, 4).tolist() == [False, False, False, False, True, True]


def test_columns_to_residues_skips_gaps():
    assert E.columns_to_residues("--AB-C", [0, 2, 3, 4, 5]) == {2: 1, 3: 2, 5: 3}


def test_persistence_clean_vs_scattered_and_minimums():
    d = [("s1", "P"), ("s2", "P"), ("s3", "P"), ("s3", "U")]
    o = [("o1", "A"), ("o2", "A")]
    f = E.persistence(d, o, "gain", True, R)
    assert f["status"] == "pass" and f["n_species_change"] == 3 and f["shape"] == "clean clade"
    # only 2 species with the change
    assert E.persistence(d[:2], o, "gain", True, R)["status"] == "fail"
    # one outgroup lineage
    assert E.persistence(d, o[:1], "gain", True, R)["status"] == "fail"
    # not a single event
    assert E.persistence(d, o, "gain", False, R)["status"] == "fail"
    # scattered: under the LCA there are many A members outside the derived subtree
    lca = d + [("x1", "A"), ("x2", "A"), ("x3", "A")]
    assert E.persistence(d, o, "gain", True, R, lca)["shape"] == "scattered"
    # loss: A in derived, P in outgroup
    fl = E.persistence([("s1", "A"), ("s2", "A"), ("s3", "A")], [("o1", "P"), ("o2", "P")], "loss", True, R)
    assert fl["status"] == "pass"


def test_biolip_overlap_threshold():
    rows = [dict(pdb_id="1abc", chain="A", ligand="ATP", ligand_class="regular", residues_uniprot=[10, 11, 12, 13]),
            dict(pdb_id="2xyz", chain="B", ligand="HEM", ligand_class="regular", residues_uniprot=[10, 50, 51, 52])]
    h = E.biolip_hits(rows, {10, 11, 12}, R)
    assert [x["pdb_id"] for x in h] == ["1ABC"] and h[0]["overlap"] == 0.75


def test_uniprot_only_experimental_eco_counts():
    fp = {5, 6, 7}
    rows = [dict(type="Binding site", start=5, end=5, ligand="ATP", chebi="CHEBI:30616", ligand_label=None,
                 eco_codes="ECO:0000269", evidence_refs="PubMed:1"),
            dict(type="Binding site", start=6, end=6, ligand="ATP", chebi="CHEBI:30616", ligand_label=None,
                 eco_codes="ECO:0000269", evidence_refs="PDB:1ABC"),
            dict(type="Binding site", start=6, end=7, ligand="Zn", chebi="CHEBI:29105", ligand_label=None,
                 eco_codes="ECO:0000250", evidence_refs=""),          # by similarity: not experimental
            dict(type="Active site", start=7, end=7, ligand=None, chebi=None, ligand_label=None,
                 eco_codes="ECO:0007744", evidence_refs="")]
    h = E.uniprot_hits(rows, fp, R)
    assert sorted(x["type"] for x in h) == ["Active site", "Binding site"]
    b = [x for x in h if x["type"] == "Binding site"][0]
    assert b["chebi"] == "CHEBI:30616" and b["refs"] == ["PDB:1ABC", "PubMed:1"]


def test_mcsa_reference_counted_homologue_not():
    rows = [dict(mcsa_id=1, resid=5, source="reference"), dict(mcsa_id=1, resid=6, source="reference"),
            dict(mcsa_id=2, resid=5, source="homologue")]
    h = E.mcsa_hits(rows, {5, 6}, R)
    assert {(x["mcsa_id"], x["counted"]) for x in h} == {(1, True), (2, False)}
    s = E.side_summary({"P1": h})
    assert s["n_mcsa_ref"] == 1


def _side(n, pdbs=(), ligs=()):
    return dict(n_members_with_evidence=n, n_biolip=len(pdbs), n_uniprot_exp=0, n_mcsa_ref=0, pdb_ids=list(pdbs),
                ligands=list(ligs), ligand_tokens=list(ligs))


def test_experimental_discrimination_gain_loss():
    # present-side holo, absent side unannotated -> pass, absence recorded as weak
    f = E.experimental_filter("gain", _side(2, ["1ABC"], ["ATP"]), _side(0), 3, 2)
    assert f["status"] == "pass" and "weak" in f["strength"] and f["about_change"]
    # absent side has a holo ligand in the footprint -> contradicted
    f = E.experimental_filter("gain", _side(2, ["1ABC"], ["ATP"]), _side(1, ["9XYZ"], ["ATP"]), 3, 2)
    assert f["status"] == "fail" and f["strength"] == "contradicted"
    # protein-level ChEMBL only -> fail, recorded as protein-level only
    f = E.experimental_filter("loss", _side(0), _side(0), 3, 2, chembl_present=1)
    assert f["status"] == "fail" and f["strength"] == "protein-level only"
    assert E.experimental_filter("loss", _side(0), _side(0), 3, 2)["strength"] == "none"


def test_experimental_specificity_needs_different_ligands():
    assert E.experimental_filter("specificity", _side(1, ligs=["TES"]), _side(1, ligs=["EST"]), 2, 2)["status"] == "pass"
    assert E.experimental_filter("specificity", _side(1, ligs=["TES"]), _side(1, ligs=["TES", "EST"]), 2, 2)["status"] == "fail"
    assert E.experimental_filter("specificity", _side(1, ligs=["TES"]), _side(0), 2, 2)["status"] == "fail"


def test_conservation_clade_specific_columns():
    # columns 0-1 footprint: derived all K/W, outgroup all A/G -> clade-specific; background columns identical
    aln = {"d1": "KWLLLLLL", "d2": "KWLLLLLL", "d3": "KWLLLLLL", "o1": "AGLLLLLL", "o2": "AGLLLLLL"}
    c = E.conservation(aln, ["d1", "d2", "d3"], ["o1", "o2"], [0, 1])
    assert c["status"] == "pass" and c["n_specific"] == 2 and c["rate_bg"] == 0.0
    # same residues on both sides -> fail
    aln2 = {k: "KWLLLLLL" for k in aln}
    assert E.conservation(aln2, ["d1", "d2", "d3"], ["o1", "o2"], [0, 1])["status"] == "fail"
    # background equally divergent -> not enriched -> fail
    aln3 = {"d1": "KWKW", "d2": "KWKW", "o1": "AGAG", "o2": "AGAG"}
    assert E.conservation(aln3, ["d1", "d2"], ["o1", "o2"], [0, 1])["status"] == "fail"
    assert E.conservation(aln, ["d1"], ["o1", "o2"], [0, 1])["status"] == "na"


def test_rescue_and_source_consistency():
    assert E.rescue_consistency([6, 5, 4, 6, 6], [0, 0])["status"] == "pass"
    assert E.rescue_consistency([1, 1, 6], [0])["status"] == "fail"
    assert E.rescue_consistency([6, 6], [2])["status"] == "fail"
    assert E.rescue_consistency([], [0])["status"] == "na"
    assert E.source_consistency([0.8, 0.5, 0.9, 0.7, 0.6])["status"] == "pass"
    assert E.source_consistency([0.1, 0.9])["status"] == "fail"
    assert E.source_consistency([])["status"] == "na"


def test_computational_filter_independence_and_pending():
    lines = dict(conservation=dict(status="pass"), rescue_consistency=dict(status="fail"),
                 cofolding=dict(status="pending"))
    assert E.computational_filter(lines, ("conservation", "rescue_consistency", "cofolding"))["status"] == "pass"
    # for specificity conservation is not independent
    assert E.computational_filter(lines, ("rescue_consistency", "cofolding"))["status"] == "fail"
    assert E.computational_filter(dict(cofolding=dict(status="pending")), ("cofolding",))["status"] == "pending"


@pytest.mark.parametrize("f1,f2,f3,art,exp", [
    ("pass", "pass", "pass", "pass", ("convincing", "convincing")),
    ("pass", "fail", "pass", "pass", ("supported", "supported")),
    ("pass", "pass", "fail", "pass", ("supported", "supported")),
    ("fail", "pass", "pass", "pass", ("rejected", "rejected")),
    ("pass", "pass", "pass", "fail", ("convincing", "rejected")),
    ("pass", "pass", "pass", "pending", ("convincing", "pending")),
    ("pass", "fail", "fail", "pending", ("rejected", "rejected")),
])
def test_labels(f1, f2, f3, art, exp):
    assert E.label(f1, f2, f3, art) == exp


def test_select_cases_top_k_audit_seeded_and_bar_not_relaxed():
    rows = [dict(case_id=f"g{i:02d}", category="gain", label="convincing", S=float(i)) for i in range(10)]
    rows += [dict(case_id=f"l{i}", category="loss", label="supported", S=100.0) for i in range(3)]
    df = pd.DataFrame(rows)
    r = E.Rules(top_k=5, n_audit=1)
    a = E.select_cases(df, r, seed=0, label_col="label")
    b = E.select_cases(df, r, seed=0, label_col="label")
    top = set(a[a.selection == "top_k"].case_id)
    assert top == {"g09", "g08", "g07", "g06", "g05"}
    aud = a[a.selection == "audit"].case_id.tolist()
    assert len(aud) == 1 and aud[0] in {"g00", "g01", "g02", "g03", "g04"}
    assert aud == b[b.selection == "audit"].case_id.tolist()            # seeded
    assert not a[a.category == "loss"].selected.any()                   # supported never fills the convincing bar


def test_evidence_count():
    f2 = dict(status="pass", present=dict(n_biolip=2, n_uniprot_exp=1, n_mcsa_ref=0))
    lines = dict(conservation=dict(status="pass"), rescue_consistency=dict(status="fail"))
    assert E.evidence_count(f2, lines, ("conservation", "rescue_consistency")) == 3
    assert E.evidence_count(dict(status="fail"), lines, ("rescue_consistency",)) == 0


def _builder():
    import importlib.util
    from pathlib import Path
    p = Path(__file__).resolve().parents[1] / "scripts/cases/build_cases.py"
    spec = importlib.util.spec_from_file_location("build_cases", p)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_disjointness_excludes_only_annotations_used_in_tune_pairs():
    """D25a regression: the tuning list holds every tune-half annotation; only those in a tune-half benchmark pair
    were used to tune tau and may be excluded."""
    B = _builder()
    tune = pd.DataFrame(dict(
        ann_id=["B|X1|ATP|0", "B|X2|ATP|0", "B|Y1|HEM|0", "U|Z1|CHEBI:1|", "B|W1|ADP|0"],
        source=["biolip", "biolip", "biolip", "uniprot", "biolip"],
        accession=["X1", "X2", "Y1", "Z1", "W1"], og_id=["OG1", "OG1", "OG1", "OG2", "OG3"],
        ligand_key=["ATP", "ATP", "HEM", "CHEBI:1", "ADP"],
        site_id=["X1:s0", "X2:s0", "Y1:s0", None, "W1:s0"], detected=[True, True, True, False, True]))
    pairs = pd.DataFrame(dict(og_id=["OG1", "OG1", "OG2", "OG3"], site_a=["X1:s0", "X1:s0", "undetected:U|Z1|CHEBI:1|", "W1:s0"],
                              site_b=["X2:s0", "Y9:s0", "Q:s0", "Q2:s0"], ligands=["ATP", "ATP|GTP", "CHEBI:1", "ADP"],
                              half=["tune", "tune", "tune", "test"]))
    used = B.used_in_tuning(tune, pairs)
    # X1, X2 (pos pair) and the undetected Z1 are used; Y1 never paired; W1 only in a test-half pair
    assert sorted(used.accession) == ["X1", "X2", "Z1"]
    site_sog = pd.DataFrame(dict(og_id=["OG1"] * 3, site_id=["X1:s0", "X2:s0", "Y1:s0"], **{"sog_t0.2": [0, 0, 0]}))
    per_sog, og_wide = B.exclusion_keys(used, site_sog)
    assert per_sog == {("OG1:g0", "X1", "biolip", "ATP"), ("OG1:g0", "X2", "biolip", "ATP")}
    assert ("OG1", "Y1", "biolip", "HEM") not in og_wide
    h = [dict(source="biolip", ligand="ATP"), dict(source="biolip", ligand="HEM")]
    kept, n = B.filter_hits(h, "OG1:g0", "OG1", "X1", (per_sog, og_wide), "sog")
    assert n == 1 and kept[0]["ligand"] == "HEM"


def test_side_family_check():
    ok = E.side_family_check([("PTHR1;PTHR1:SF2", "PF1")], [("PTHR1:SF9", "PF2")])
    assert ok["status"] == "pass"
    bad = E.side_family_check([("PTHR1", "PF1")], [("PTHR2", "PF1")])
    assert bad["status"] == "fail"
    # PANTHER missing on one side -> Pfam decides
    assert E.side_family_check([("", "PF1")], [("PTHR2", "PF1;PF3")])["status"] == "pass"
    assert E.side_family_check([("", "")], [("PTHR2", "PF1")])["status"] == "pending"


def test_combine_pool_dedup_by_protein_set():
    import importlib.util
    from pathlib import Path
    p = Path(__file__).resolve().parents[1] / "scripts/cases/combine_pool.py"
    spec = importlib.util.spec_from_file_location("combine_pool", p)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    pool = pd.DataFrame(dict(case_id=["o1_a", "o2_a", "o2_b", "o1_l"], category=["gain", "gain", "gain", "loss"],
                             S=[5.0, 5.0, 3.0, 5.0], orth_rank=[0, 1, 1, 0],
                             present_accs=["P1;P2;P3", "P1;P2;P3;P4", "Q1;Q2", "P1;P2;P3"]))
    d = m.dedup(pool).set_index("case_id")
    assert d.loc["o1_a", "duplicate_of"] == "o2_a"        # same S, later list (O2) wins the tie
    assert pd.isna(d.loc["o2_a", "duplicate_of"]) and pd.isna(d.loc["o2_b", "duplicate_of"])
    assert pd.isna(d.loc["o1_l", "duplicate_of"])          # different category is never a duplicate
