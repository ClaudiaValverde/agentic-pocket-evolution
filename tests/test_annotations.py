"""Known-answer and consistency checks for data/annotations/*.parquet (scripts/annotations/).

Run:  .venv/bin/python -m pytest tests/test_annotations.py -v
Tables that have not been built are skipped (not failed). Unit tests of parsers use toy inputs.
"""
from __future__ import annotations

import hashlib
import json
import os
import sys

import pandas as pd
import pytest

P = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
A = os.path.join(P, "data", "annotations")
sys.path.insert(0, os.path.join(P, "scripts", "annotations"))


def load(name, **kw):
    p = os.path.join(A, name)
    if not os.path.exists(p):
        pytest.skip(f"{name} not built")
    return pd.read_parquet(p, **kw)


# ------------------------------------------------------------------ toy / unit tests
def test_uniprot_parse_entry_toy():
    import uniprot
    e = {"primaryAccession": "X1", "entryType": "UniProtKB reviewed (Swiss-Prot)",
         "proteinDescription": {"recommendedName": {"fullName": {"value": "toy"}}, "flag": "Fragment"},
         "sequence": {"value": "MAH", "length": 3, "crc64": "C", "md5": "M"},
         "features": [
             {"type": "Binding site", "location": {"start": {"value": 2}, "end": {"value": 2}},
              "ligand": {"name": "Zn(2+)", "id": "ChEBI:CHEBI:29105"},
              "evidences": [{"evidenceCode": "ECO:0000250"}]},
             {"type": "Binding site", "location": {"start": {"value": 3}, "end": {"value": 3}},
              "ligand": {"name": "Zn(2+)", "id": "ChEBI:CHEBI:29105"},
              "evidences": [{"evidenceCode": "ECO:0007744", "source": "PDB", "id": "1ABC"}]},
             {"type": "Region", "location": {"start": {"value": 1}, "end": {"value": 2}}, "description": "Disordered"}]}
    prot, sites = uniprot.parse_entry(e)
    assert prot["fragment"] and prot["reviewed"] and prot["disordered"] == "1-2"
    assert [s["experimental"] for s in sites] == [False, True]   # ECO:0000250 (by similarity) is not experimental
    assert sites[0]["chebi"] == "CHEBI:29105"


def test_biolip_mapping_toy():
    import biolip
    m, ident, cov = biolip.map_positions(("ACDEF", "MMACDEFGG"))
    assert m[1] == 3 and m[5] == 7 and ident == 1.0
    # receptor with an unobserved residue (gap) and one mutation
    m, ident, cov = biolip.map_positions(("MKTAYLWVNPQRS", "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILWVNPQRS"))
    assert m[13] == 44 and m[1] == 1


def test_ted_domain_md5_toy():
    import ted
    assert ted.domain_md5("ABCDEFGH", "2-3_6-7") == hashlib.md5(b"BCFG").hexdigest()


# ------------------------------------------------------------------ UniProt known answers
def test_uniprot_coverage():
    acc = load("accessions.parquet")
    pr = load("uniprot_proteins.parquet")
    assert pr.accession.is_unique
    assert set(pr.accession) == set(acc.accession)
    assert pr.uniprot_acc.notna().mean() > 0.99
    assert pr.uniprot_release.nunique() == 1


def test_ca2_zinc_binding():
    s = load("uniprot_site_features.parquet")
    z = s[(s.uniprot_acc == "P00918") & (s.type == "Binding site") & (s.chebi == "CHEBI:29105")]
    # UniProt numbering of human CA2 (initiator Met counted): His94, His96, His119 (= CA numbering His94/96/119)
    assert {94, 96, 119} <= set(z.start)
    assert z[z.start.isin([94, 96, 119])].experimental.all()


def test_cdk2_atp_binding():
    s = load("uniprot_site_features.parquet")
    atp = s[(s.uniprot_acc == "P24941") & (s.type == "Binding site") & (s.ligand == "ATP")]
    assert len(atp) >= 3
    assert (atp.chebi == "CHEBI:30616").all()
    assert 33 in set(atp.start)            # catalytic Lys33
    assert ((atp.start <= 10) & (atp.end >= 10)).any() or 10 in set(atp.start)   # glycine-rich loop G-x-G-x-x-G
    pr = load("uniprot_proteins.parquet")
    r = pr.set_index("accession").loc["P24941"]
    assert "PF00069" in r.pfam and "2.7.11.22" in r.ec and r.gene_primary == "CDK2"


def test_experimental_flag_definition():
    s = load("uniprot_site_features.parquet", columns=["eco_codes", "experimental"])
    s = s.sample(min(len(s), 200000), random_state=20261030)
    exp = s.eco_codes.fillna("").str.contains("ECO:0000269|ECO:0007744", regex=True)
    assert (exp == s.experimental).all()


# ------------------------------------------------------------------ M-CSA
def test_mcsa_chymotrypsin_trypsin_triad():
    r = load("mcsa_residues.parquet")
    chy = r[(r.uniprot_acc == "P00766") & (r.mcsa_id == 387) & r.is_reference]
    got = dict(zip(chy.resid, chy.code))
    # bovine chymotrypsinogen A: UniProt numbering == chymotrypsin numbering
    assert got.get(57) == "His" and got.get(102) == "Asp" and got.get(195) == "Ser"
    # bovine cationic trypsin P00760 (homologue of entry 387): UniProt 63/107/200 = chymotrypsin 57/102/195
    tr = r[(r.uniprot_acc == "P00760") & (r.mcsa_id == 387)]
    m = dict(zip(tr.resid, tr.ref_auth_resid))
    assert m.get(63) == 57 and m.get(107) == 102 and m.get(200) == 195


def test_riziotis_cases():
    c = load("mcsa_convergent_cases.parquet")
    assert c.case_id.nunique() == 34
    assert set(c.category) == {"structural+mechanistic", "mechanistic-only", "reaction-only"}
    ca = c[c.case_id == c[c.mcsa_id == 216].case_id.iloc[0]]
    assert {216, 516} <= set(ca.mcsa_id)      # alpha- vs gamma-carbonic anhydrase
    assert (c.groupby("case_id").size() >= 2).all()


# ------------------------------------------------------------------ Pfam / TED
def test_pfam_clans():
    c = load("pfam_clans.parquet").set_index("pfam")
    assert c.loc["PF00069", "clan"] == "CL0016"     # Pkinase in Protein kinase superfamily clan
    assert c.loc["PF00089", "clan"] == "CL0124"     # Trypsin in Peptidase PA clan


def test_pfam_regions_cdk2():
    r = load("pfam_regions.parquet")
    k = r[(r.accession == "P24941") & (r.pfam == "PF00069")]
    assert len(k) == 1 and k.seq_start.iloc[0] < 33 < k.seq_end.iloc[0]


def test_ted_cdk2():
    t = load("ted_domains.parquet")
    k = t[t.accession == "P24941"]
    assert len(k) >= 1
    assert k.cath_label.fillna("").str.match(r"^(3\.30\.200\.20|1\.10\.510\.10)").any()


# ------------------------------------------------------------------ BioLiP / ChEMBL
def test_biolip_ca2_zinc():
    b = load("biolip_sites.parquet")
    z = b[(b.accession == "P00918") & (b.ligand == "ZN")]
    assert len(z) > 50
    hits = z.residues_uniprot.map(lambda r: {94, 96, 119} <= set(r))
    assert hits.mean() > 0.8
    assert (b.n_mapped <= b.n_site_residues).all()


def test_chembl_cdk2():
    c = load("chembl_by_accession.parquet").set_index("accession")
    assert c.loc["P24941", "n_cpd_pchembl6_single"] > 300
    assert c.loc["P00918", "n_cpd_pchembl6_single"] > 300     # CA2: heavily profiled sulfonamides
    assert (c.n_cpd_pchembl6_single <= c.n_cpd_pchembl6_any).all()


def test_meta_files_record_versions():
    for name, key in [("uniprot_meta.json", "uniprot_release"), ("pfam_meta.json", "pfam_release"),
                      ("chembl_meta.json", "chembl_version"), ("ted_meta.json", "sha256"),
                      ("biolip_meta.json", "sha256"), ("mcsa_meta.json", "sha256")]:
        p = os.path.join(A, name)
        if os.path.exists(p):
            assert json.load(open(p)).get(key)
