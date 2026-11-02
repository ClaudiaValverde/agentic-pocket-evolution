# ANNOTATIONS — per-accession annotation tables (snapshot 20261030T085136Z, tiers 1+2)

Built by `scripts/annotations/run_all.sh <manifest.parquet>` (re-runnable; each script `--help`). Raw downloads in
`/tmp/pe_big/raw/<dataset>/` (not in git, re-downloaded if missing); outputs in `data/annotations/` (161 MB total).
Every table carries its source version column; `*_meta.json` hold sha256 of raw files + coverage. Sources,
versions, checksums and licences: DATA.md. Tests: `tests/test_annotations.py` (15 known-answer/toy tests).
Join key everywhere: `accession` = atlas accession (`structures/tier*/…parquet`). In UniProt 2026_03 all
989,696 atlas accessions are current primary accessions (0 secondary, 0 missing).

## Coverage (accessions; tier 1 = 552,810, tier 2 = 436,886, no overlap)
| annotation | tier 1 | tier 2 | total |
|---|---|---|---|
| UniProt entry found | 552,810 | 436,886 | 989,696 |
| Pfam xref (UniProt) | 415,860 | 426,497 | 842,357 |
| Pfam region coordinates (Pfam-A.regions) | 422,606 | 200,106 | 622,712 |
| TED domain (any) / with CATH H-level label | 464,280 / 291,925 | 421,805 / 361,138 | 886,085 / 653,063 |
| UniProt binding/active site (any evidence) | 84,102 | 213,708 | 297,810 |
| … with experimental ECO (0000269/0007744) | 4,156 | 3,359 | 7,515 |
| M-CSA catalytic residues (reference or homologue) | 15,816 | 105,039 | 120,855 |
| BioLiP2 site (any ligand) / regular small molecule | 13,646 / 7,874 | 9,148 / 6,015 | 22,794 / 13,889 |
| ChEMBL ≥1 compound pChEMBL ≥ 6 (single protein, conf 9) | 4,831 | 480 | 5,311 |
| PDB xref | 24,087 | 15,171 | 39,258 |

## Tables
**accessions.parquet** — accession, tier, taxid, organism, proteome (UP id, tier 1), n_res, fragment (AFDB fragment no).

**uniprot_proteins.parquet** (UniProt 2026_03; 1 row/accession) — accession, tiers, uniprot_acc, secondary_accs,
uniprot_id, reviewed, gene_primary, protein_name, organism_id, length, seq_crc64, seq_md5, sequence_version,
fragment (UniProt "Fragment" flag; 21,254 true), protein_existence, annotation_score, ec, cofactors (`name|CHEBI`),
pfam, panther, interpro, supfam, gene3d (`;`-joined ids), n_pdb, pdb_ids, signal, transmem, disordered
(`start-end;…`; disordered = UniProt Region "Disordered", MobiDB-lite), accession_status.
Full sequences are written outside git to `/tmp/pe_big/work/annotations/uniprot_sequences.parquet`.

**uniprot_site_features.parquet** (1.79 M rows) — uniprot_acc, type (Binding site 1.42 M / Active site 214 k /
Site 62 k / Mutagenesis 103 k), start, end, start/end_modifier, description, ligand, chebi (`CHEBI:n`),
ligand_label, ligand_note, ligand_part, alt_seq (mutagenesis), eco_codes, evidence_refs (`PubMed:id`, `PDB:id`),
experimental (= any ECO:0000269 or ECO:0007744; ECO:0000250/0000255/0000305/RuleBase are *not* experimental).
Experimental rows: Binding 50,391, Active 2,116, Site 3,024, Mutagenesis 101,925.

**mcsa_entries.parquet** (1,003) — mcsa_id, enzyme_name, reference_uniprot, uniprot_ids, ec, all_ecs,
reference_pdb, cath, n_residues, n_mechanisms, mechanism_rating, is_detailed_mechanism, cofactors, reference_in_atlas.
**mcsa_residues.parquet** (789 k) — mcsa_id, uniprot_acc, resid (UniProt numbering), code, is_reference,
source (reference/homologue), roles_summary, function_location (side chain / main-N …), ptm, ref_pdb, ref_chain,
ref_auth_resid (PDB numbering of the reference residue, e.g. chymotrypsin numbering), ref_cath, accession, in_atlas.
Homologue rows come from M-CSA `homologues_residues.json` (sequence-transferred, not individually curated).
**mcsa_convergent_cases.parquet** — Riziotis et al. 2025 (FEBS J, doi:10.1111/febs.17332) Tables 1–3: 34 cases,
72 member enzymes; case_id, category (structural+mechanistic 14 / mechanistic-only 11 / reaction-only 9), enzyme,
ec, cath, mcsa_id, pdb, residue_types, superposed_residues, reference_uniprot, reference_in_atlas,
n_atlas_accessions_with_mcsa_residues. No supplementary data file exists; parsed from the open-access article XML.

**pfam_clans.parquet** (Pfam 38.2; 30,134 families, 13,933 in a clan) — pfam, clan, clan_id, pfam_id, description.
**pfam_regions.parquet** (1.14 M) — accession, pfam, clan, seq_start/seq_end (envelope), ali_start/ali_end, crc64,
seq_matches_uniprot (96.95 % of regions are on the current UniProt sequence). Pfam-A.regions covers only
pfamseq (reference-proteome based, UniProt 2025_03), hence fewer tier-2 accessions than the UniProt Pfam xref.

**ted_domains.parquet** (1.67 M domains) — accession, ted_id, domain_no, consensus_level, chopping
(`s-e[_s-e]`), nres_domain, num_segments, plddt, cath_label, cath_level (H / T), cath_method, packing_density,
norm_rg, globular, af_fragment, seq_matches_uniprot (md5 of chopped current UniProt sequence == TED md5_domain;
99.40 % match). TED is on AFDB **v4** models; use only `seq_matches_uniprot` rows for residue-level overlap with v6.
Long proteins (AFDB multi-fragment) have fragment-local numbering (`af_fragment` ≠ F1).

**biolip_sites.parquet** (677 k sites, 358 k PDB chains) — accession, uniprot_acc, pdb_id, chain, resolution,
site_code, ligand (CCD / dna / rna / peptide), ligand_class (regular 322 k, metal 149 k, rna 141 k, dna 37 k,
peptide 29 k), ligand_name, ligand_inchikey, ligand_chembl, ligand_chain, ligand_serial, residues_pdb,
residues_uniprot (list), n_site_residues, n_mapped, aln_identity, aln_coverage, catalytic_uniprot, ec,
affinity_manual/moad/pdbbind/bindingdb, pubmed. UniProt numbering is ours: BioLiP receptor sequence aligned
(semi-global) to the UniProt sequence; residues kept only if amino acid identical (99.2 % of residues mapped,
96.3 % of sites fully mapped). Downloaded file is dated 2025-11-12 (989 k rows), although the site advertises
weekly updates; BioLiP states no licence (free download, cite) — keep the raw file out of git.

**chembl_targets.parquet** (ChEMBL 37) — accession × target: target_chembl_id, target_name, target_type,
target_tax_id, n_cpd_pchembl6 (distinct molregno, pChEMBL ≥ 6, validity comment null/“Manually validated”,
potential_duplicate = 0), n_cpd_pchembl6_conf9, n_act_pchembl6, max_pchembl.
**chembl_by_accession.parquet** (10,966) — n_chembl_targets, single_protein_targets, max_pchembl,
n_cpd_pchembl6_single (SINGLE PROTEIN target, assay confidence 9), n_cpd_pchembl6_any (any target incl. complexes).
Computed from the full ChEMBL SQLite (in /tmp, 29 GB) — the REST API cannot give distinct counts efficiently.

## Known answers checked (tests)
CA2 P00918 Zn (CHEBI:29105) binding at His94/His96/His119 — UniProt numbering equals CA numbering (experimental);
CDK2 P24941 ATP binding 10–18, K33, 81–83, 86, 129–132, D145; Pfam PF00069 (CL0016) spans K33; TED kinase CATH;
M-CSA chymotrypsinogen A P00766 His57/Asp102/Ser195 (UniProt = chymotrypsin numbering); bovine trypsin P00760
His63/Asp107/Ser200 map to reference 57/102/195; BioLiP CA2 ZN sites contain 94/96/119; ChEMBL CDK2 847 and
CA2 6,031 single-protein compounds with pChEMBL ≥ 6; Riziotis case α- vs γ-carbonic anhydrase (M-CSA 216/516).

## Caveats
- Experimental binding/active-site features exist for only 7.5 k accessions (mostly human, yeast, E. coli,
  model organisms); most UniProt sites are ECO:0000250/0000255 (by similarity / rule) and must not be used as
  experimental evidence (DESIGN §5.4) — they are, however, valid for homology-class checks.
- AFDB v6 model length == UniProt 2026_03 length for 99.55 % of F1 models (4,407 differ, mostly tier 1):
  treat UniProt numbering as model numbering only where lengths (ideally sequences) agree.
- M-CSA homologue residues, UniProt by-similarity sites and TED/Pfam labels are not independent of homology:
  do not use them as evidence for convergence.
