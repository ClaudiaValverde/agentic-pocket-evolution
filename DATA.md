# DATA — datasets used

Each entry: source, version, checksum, license. Restricted data are never committed. Large data stay at their source path (read-only mounts) and are referenced, not copied.

| dataset | source / path | version | checksum | license | notes |
|---|---|---|---|---|---|
| AFDB pocket atlas | /cluster/scratch/gmagnusson/pocketcompare_afdb/ | snapshot manifest (see data/snapshot/) | manifest sha256 in data/snapshot/*/manifest.json | derived from AFDB (CC BY 4.0) | live-growing; tier 3 partial |
| AlphaFold DB v6 models | /cluster/work/beltrao/shared/alphafold/v6/ | v6 (49 proteome tars) | per-tar size + mtime in data/snapshot/<id>/manifest.parquet; full sha256 not computed (read-only, ~TB) | CC BY 4.0 | |
| PDB mmCIF mirror | /cluster/scratch/jjaenes/pdb/ | mirror subdirs dated 2024-02-23 … 2026-10-09 (mtime) | not computed (read-only mirror); entries used are listed by PDB id in case/calibration files | CC0 | |
| AF3 cached MSAs (human, yeast) | /cluster/project/beltrao/shared/25.0{5,6}_alphafold3_msas_* | 2025.05 / 2025.06 | not computed; each used MSA's query sequence was checked against the folded sequence (docs/COFOLD.md) | internal | not committed |
| TimeTree 5 pairwise divergence times (tier-1 species, 1,035 pairs) | https://timetree.org/api/pairwise/<taxid_a>/<taxid_b>, fetched by scripts/partA_descent/fetch_timetree.py -> data/timetree/20261030T085136Z_O1_pairwise.csv | TimeTree API version 20260820, retrieved 2026-10-30 | sha256 b23c19fd6ddeb92e7652582e7852e67df3e0a8ba8f7dedfdb9d36b1feedf0f57 | TimeTree terms (free academic use; cite Kumar et al. 2022, Mol Biol Evol 39:msac174) | 552 pairs with an age; 483 HTTP 500 = all pairs whose MRCA is the root of cellular life (TimeTree: "No TimeTree study info available for this MRCA"): rank proxy 4,200 Mya (LUCA, Moody et al. 2024), flagged `rank_proxy`; 3 pairs node-median, 1 missing (Herpotrichiellaceae) |
| AF3 weights | /cluster/project/beltrao/shared/alphafold3/models/af3.bin.zst | file mtime 2026-06-04, 1,020,545,840 bytes | sha256 74d0258616917cd122f5eab6d076afe4a8930e96823851e65e4f777dfb1f33ff | AF3 model terms (restricted, non-commercial) | never committed |

## Code dependencies (pinned commits)
| ancestra-mvp-ai | https://github.com/evocellnet/ancestra-mvp-ai.git | ffb0bbef28bc31d05755264bbe0094e4c0d4a5fa |
| burrow | https://github.com/zoidpuff/burrow.git | b2da7abfbbc56d4be4850dfd44e02e55c7b5a648 |
| pocketcompare | https://github.com/zoidpuff/pocketcompare.git | d82357e977a08006506ebec59420eb58fb469dd0 |
| pocketcompare-agents | https://github.com/zoidpuff/pocketcompare-agents.git | 173daa9ed67166249489ae07dbd304c041416474 |
| 3Dmol.js (vendored, webapp/static/vendor/3Dmol-min.js) | https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.4.2/3Dmol-min.js | 2.4.2 | sha256 06b6d2fc7d418e8bef62a32cca44373ff21b2b787bd377613dd4039394a4e9ff | BSD-3-Clause (3Dmol.js, Rego & Koes 2015) |

## ColabFold MMseqs2 MSAs (server provenance)

Unpaired MSAs for proteins without a cached AF3 (jackhmmer) MSA. Fetched with `scripts/cofold/colabfold_msa.py`. The ColabFold server does not expose a version endpoint; the returned tar member names are the server's own database version strings and are recorded verbatim. License: ColabFold/MMseqs2 server is free for non-commercial and academic use (Mirdita et al. 2022, Nat Methods); underlying DBs UniRef30 (CC BY 4.0) and ColabFoldDB envdb.

| acc | len | seq sha256[:12] | server | mode | ticket | fetched (UTC) | db files returned (n seqs) | n seqs | depth/sqrt(L) |
|---|---|---|---|---|---|---|---|---|---|
| P24941 | 298 | fe22d91e18cd | https://api.colabfold.com | env | KcWlpZpQZF772Cnuy_-S6HtpFDN31dmbOwUdLA | 2026-10-30T09:45:23Z | `bfd.mgnify30.metaeuk30.smag30.a3m` (11192); `uniref.a3m` (8820) | 20013 | 1159.32 |

## Annotation datasets (scripts/annotations/, processed tables in data/annotations/, see docs/ANNOTATIONS.md)
Raw downloads live in /tmp/pe_big/raw/<dataset>/ (not committed); sha256 below is of the raw file as downloaded 2026-10-30.

| dataset | source / path | version | checksum (sha256) | license | notes |
|---|---|---|---|---|---|
| UniProtKB entries (tier-1 proteomes + Swiss-Prot) | https://rest.uniprot.org/uniprotkb/search (JSON, cursor pages) | release 2026_03 (X-UniProt-Release header; 2026-09-02) | per query file in data/annotations/uniprot_meta.json (`raw_files`) | CC BY 4.0 | proteins + binding/active-site/site/mutagenesis features with ECO codes |
| Pfam-A clans | https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz | Pfam 38.2 (2026-01; UniProtKB 2025_03) | 86062b7ef1a0e0caee0c28cef479ac0d294c80789c51cbf75e513b368ce3a6f6 | CC0 | md5 verified against FTP md5_checksums |
| Pfam-A regions | https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.regions.tsv.gz | Pfam 38.2 | 6be816409bc6d3f58005fbecaacfe1542046b87c64ec6e69636b864fc1a9dce5 | CC0 | 5.1 GB raw; filtered to atlas accessions |
| M-CSA entries/residues/homologue residues | https://www.ebi.ac.uk/thornton-srv/m-csa/api/ (entries, residues, homologues_residues.json) | live API, 1,003 entries, fetched 2026-10-30 | entries_all.json c99c0f05…589d; residues.json 92dc99ac…ca06; homologues_residues.json bdf77a8c…48b9 (full hashes in mcsa_meta.json) | CC BY 4.0 | catalytic residues in UniProt numbering |
| Riziotis et al. 2025 convergent enzymes (Tables 1-3) | Europe PMC full text XML PMC11796326 (doi:10.1111/febs.17332) | FEBS J 292:537 (2025) | eff77d51358dde36cdca3f6be543c0277099e11393d07de526ecf3bd1d7187b2 | CC BY 4.0 (open-access article) | no supplementary data file exists; 34 cases parsed from the article tables |
| BioLiP2 annotations | https://zhanggroup.org/BioLiP/download/BioLiP.txt.gz (+ data/ligand.tsv.gz) | file Last-Modified 2025-11-12 (989,058 rows; site header says weekly, 2026-08-30) | BioLiP.txt.gz c92229bbc8c55c3bd84a9813c3e278ba62f4cfa44e6315cc98d9bf63ed64b6ec; ligand.tsv.gz 7d107cfaf0cc873d393e9dcee7bbe16fc7dc024a9748495a92adf6e7658f61ec | no explicit licence ("all data freely available"; cite Yang 2013 NAR / Zhang 2024 NAR); PDBbind affinities removed by BioLiP | derived table only; do not redistribute raw |
| ChEMBL SQLite | https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_37_sqlite.tar.gz | ChEMBL 37 (2026-05-01) | 33c203740555f96067710cdfc1c3c55d890660e5908ec5cbf5817492c290d281 (matches FTP checksums.txt) | CC BY-SA 3.0 (attribution + ChEMBL IDs/release kept) | 5.8 GB raw, 29 GB extracted in /tmp/pe_big/work/chembl |
| TED domain summary | https://zenodo.org/records/13908086 ted_365m.domain_summary.cath.globularity.taxid.tsv.gz | Zenodo v5 (2024-10-31), AFDB v4 models | 0b0193460e8bbda890ffc1f6dab580db433260f71838ea94e30e505e9d38abde (md5 51468b8e… verified vs Zenodo) | CC BY 4.0 | 19.9 GB raw; filtered to atlas accessions |

## Orthology O2 inputs (scripts/orthology/, outputs in data/orthology/ancestra_O2/, see docs/ORTHOLOGY.md)
Raw downloads in /tmp/pe_big/ancestra/ (not committed; regenerable with the scripts). Checksums are md5 of the file as downloaded 2026-10-30.

| dataset | source / path | version | checksum (md5) | license | notes |
|---|---|---|---|---|---|
| UniProt reference proteome list | https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/README | Release 2026_03 (02-Sep-2026) | e442d51330199a402ca7f646da913bf5 | CC BY 4.0 | candidate list for context selection (37,190 proteomes) |
| O2 proteomes (46 tier-1 + 111 context) | `.../reference_proteomes/<Kingdom>/<UP>/<UP>_<taxid>.fasta.gz` (canonical, one protein per gene) | Release 2026_03 | per file in data/orthology/ancestra_O2/o2_proteomes.tsv (md5 of that list 6397650ef9bcfa446d4332c581e2ea5c) | CC BY 4.0 | proteome ID list = o2_proteomes.tsv; selection rule/coverage = context_species.tsv |
| NCBI taxonomy dump | https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz → /tmp/pe_big/ancestra/taxdump | downloaded 2026-10-30 08:57 UTC | taxdump.tar.gz a7931e8257ae5d55950eb5a43ac004d1 | public domain (NCBI) | guide tree, lineages, placement compatibility |
| Swiss-Prot FASTA | https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz | Release 2026_03 (02-Sep-2026) | bc9d398533e6df582b563c6c03093bd0 | CC BY 4.0 | tier-2 placement query sequences |
| UniProtKB / UniParc REST (ID gap, context placement) | rest.uniprot.org uniprotkb/accessions, uniprotkb/<acc>.json, uniprotkb/stream?query=organism_id:<taxid>, uniparc/stream?query=upi:/dbid: | Release 2026_03, retrieved 2026-10-30 | idgap caches: uniprot_missing.tsv b446970988076076b019e4249e62ffb4, uniprot_inactive.tsv 17f19d23c78fb1c5b06e5e10c43b9fb3 | CC BY 4.0 | sequences of non-member / deleted atlas accessions |
| PANTHER family xrefs (validation) | rest.uniprot.org/uniprotkb/stream?query=proteome:<UP>&fields=accession,xref_panther (8 proteomes) | UniProt 2026_03 xrefs (PANTHER version as cross-referenced by UniProt) | cached in /tmp/pe_big/ancestra/o2/panther/ | PANTHER: GPL/free for research; UniProt CC BY 4.0 | validation only |
| Ancestra O2 output | data/orthology/ancestra_O2/accurate/hogs.orthoxml.gz (+ fast/) | Ancestra ffb0bbe, accurate preset | accurate 8f61b5e80e05cbf59d3f24b1dba7e175; fast 505e82f6561023176ce44d02d9f279f2; tables/hogs.parquet e4dd12b81436de5762cbe984e9a085f1; tables/protein_hogs.parquet 2017a3f3ed09d0fdecadc0028e8f6f78; idgap_mapping.parquet a3fe68cb8f2af40928a4edc079fa6c4b | derived (CC BY 4.0 inputs) | generated here |
| CDK2_site_gly | 298 | 6cdf399ff732 | https://api.colabfold.com | env | MglOWdBBPj1w3OgNfgojzFAZadYkGBiav8YDcw | 2026-10-30T09:58:59Z | `bfd.mgnify30.metaeuk30.smag30.a3m` (8418); `uniref.a3m` (8428) | 16847 | 975.92 |
| A0A3Q0KPI7 | 1262 | 0111d2909038 | https://api.colabfold.com | env | okPymrAPzrgKWlUcTUeHSmYpX17xNhiiG0c-bg | 2026-10-30T14:59:25Z | `bfd.mgnify30.metaeuk30.smag30.a3m` (3202); `uniref.a3m` (2789) | 5992 | 168.67 |
| Q08623 | 228 | 012183f9b101 | https://api.colabfold.com | env | DxMcV71pdpFvTeXfyWJxfOoKfESiX1WWyrn-Qg | 2026-10-30T15:06:27Z | `bfd.mgnify30.metaeuk30.smag30.a3m` (5412); `uniref.a3m` (5865) | 11278 | 746.9 |
| A0A1D8PCA4 | 1300 | 0abe0d353bba | https://api.colabfold.com | env | dAIBvBHcyjwLINybAjghV8fzpgNbF568_jFf0A | 2026-10-30T14:43:18Z | `bfd.mgnify30.metaeuk30.smag30.a3m` (3642); `uniref.a3m` (3113) | 6756 | 187.38 |
| P46708 | 699 | 1ee3150b7272 | https://api.colabfold.com | env | gR7U5AVW_Ue8l3gD433HTNgBBgkkDiAD0AM3XA | 2026-10-30T14:47:23Z | `bfd.mgnify30.metaeuk30.smag30.a3m` (15325); `uniref.a3m` (8750) | 24076 | 910.64 |
| P07648 | 1122 | 2a14320720a5 | https://api.colabfold.com | env | prxHQ5GCic7dj58ILeVQCaXm4Eyo_L887HcMDQ | 2026-10-30T15:02:25Z | `bfd.mgnify30.metaeuk30.smag30.a3m` (6504); `uniref.a3m` (4973) | 11478 | 342.67 |
| Q4DVM5 | 242 | 2b9cfe2554f5 | https://api.colabfold.com | env | 5QJS44zbSxYiejD0VzUee1Yb7NWGZTjP5PDcPA | 2026-10-30T14:57:36Z | `bfd.mgnify30.metaeuk30.smag30.a3m` (3698); `uniref.a3m` (2459) | 6158 | 395.85 |
| A8Q9D7 | 923 | 425ee9f457c8 | https://api.colabfold.com | env | pv9ZztF0VrW9KaBFHqj0aZmi6YWAeJKkgPcOgw | 2026-10-30T14:47:48Z | `bfd.mgnify30.metaeuk30.smag30.a3m` (1441); `uniref.a3m` (2499) | 3941 | 129.72 |
| Q5F4Z2 | 515 | 45123f703b71 | https://api.colabfold.com | env | I9VEpVewIsW_ag-aWYKQhFmWBy0HEfTZvnTnOw | 2026-10-30T14:59:56Z | `bfd.mgnify30.metaeuk30.smag30.a3m` (15254); `uniref.a3m` (6848) | 22103 | 973.97 |
| G5EEY6 | 251 | 4b0046ee9116 | https://api.colabfold.com | env | 4oJzcMDc4LPUwuEly_BXA0YZAjasrtz4kABziw | 2026-10-30T14:57:35Z | `bfd.mgnify30.metaeuk30.smag30.a3m` (1144); `uniref.a3m` (2824) | 3969 | 250.52 |
| U7PTL3 | 266 | 4e12bb246f8d | https://api.colabfold.com | env | BL2q2mvRv9NhkP1_0_a2lywRawebKh0Wv-0agQ | 2026-10-30T15:11:23Z | `bfd.mgnify30.metaeuk30.smag30.a3m` (13917); `uniref.a3m` (9280) | 23198 | 1422.36 |
| A0A8M3B7A6 | 1034 | 570a8216c49a | https://api.colabfold.com | env | HqZaKnIGK4D23rEy4jGndS4b2es4iEW72kLzzg | 2026-10-30T14:37:43Z | `bfd.mgnify30.metaeuk30.smag30.a3m` (8441); `uniref.a3m` (6612) | 15054 | 468.16 |
| A0A0N4UML4 | 310 | 58871a831696 | https://api.colabfold.com | env | dfam3y7pOOE6k0JlOBd9P1Ptaho84YPzrS6AbQ | 2026-10-30T14:39:46Z | `bfd.mgnify30.metaeuk30.smag30.a3m` (2227); `uniref.a3m` (3839) | 6067 | 344.58 |
| A0A175VQU5 | 399 | 5d937e5ebfcc | https://api.colabfold.com | env | WVbFVt99CA9Z8giVSy9c6KYpbA2QGR54-5BoBQ | 2026-10-30T14:42:40Z | `bfd.mgnify30.metaeuk30.smag30.a3m` (11499); `uniref.a3m` (7418) | 18918 | 947.08 |
| B0V1B1 | 296 | 611b0510140a | https://api.colabfold.com | env | 2yb9ZPIAkMZWPMaljbTvnhv15aeahWg7EAVfqg | 2026-10-30T14:50:10Z | `bfd.mgnify30.metaeuk30.smag30.a3m` (1280); `uniref.a3m` (4482) | 5763 | 334.97 |
| A0A5K4FBA2 | 1154 | 619dc652fd34 | https://api.colabfold.com | env | _iVXI5yFj7DnrsBhQxwNcKD5HezlZzQIdi3kRg | 2026-10-30T14:33:39Z | `bfd.mgnify30.metaeuk30.smag30.a3m` (2896); `uniref.a3m` (3303) | 6200 | 182.51 |
| A0A0H3GQA5 | 427 | 636e9530469f | https://api.colabfold.com | env | FcMBFdamnj_nbSgnwVRmCheW8NrIkeGxenRcwA | 2026-10-30T14:35:12Z | `bfd.mgnify30.metaeuk30.smag30.a3m` (12199); `uniref.a3m` (6718) | 18918 | 915.51 |
| Q0JA69 | 326 | 664de37f2262 | https://api.colabfold.com | env | -8hkfSkwn0fvVXE_aPcPCtBvOlPbUCPw1Chhjw | 2026-10-30T14:50:43Z | `bfd.mgnify30.metaeuk30.smag30.a3m` (1595); `uniref.a3m` (1479) | 3075 | 170.31 |
| E7F1E4 | 800 | 672ce9080053 | https://api.colabfold.com | env | _WEjsjJvlXJT4DlHNdLuPm4YJTwDPTGygSFbkA | 2026-10-30T14:53:03Z | `bfd.mgnify30.metaeuk30.smag30.a3m` (414); `uniref.a3m` (1482) | 1897 | 67.07 |
| Q4DVM4 | 242 | 687c129b59ae | https://api.colabfold.com | env | APlFIxk2w2pgeKG4KA2LjCkq5n39BUV_F-jJPQ | 2026-10-30T14:56:21Z | `bfd.mgnify30.metaeuk30.smag30.a3m` (3685); `uniref.a3m` (2574) | 6260 | 402.41 |
| A0A0D2DCK2 | 945 | 6ccaafc6e036 | https://api.colabfold.com | env | FparbkOegFlzqQz9DBU2fo04Aw1cFBf-QyoQiQ | 2026-10-30T14:34:31Z | `bfd.mgnify30.metaeuk30.smag30.a3m` (2553); `uniref.a3m` (2794) | 5348 | 173.97 |
| C1HC76 | 402 | 731e231ad486 | https://api.colabfold.com | env | KtIMzsSJ6-1UEoywquoNd55rnl3LzXjVOBW_Cw | 2026-10-30T14:51:11Z | `bfd.mgnify30.metaeuk30.smag30.a3m` (742); `uniref.a3m` (1114) | 1857 | 92.62 |
| O75899 | 941 | 83e80395cf51 | https://api.colabfold.com | env | zgCvyPmF9F44ppgSiquKAV_xdDoUbVQsJprh9g | 2026-10-30T15:01:37Z | `bfd.mgnify30.metaeuk30.smag30.a3m` (1396); `uniref.a3m` (2452) | 3849 | 125.47 |
| C0NQ87 | 975 | 87f9eb14d93d | https://api.colabfold.com | env | uGBxhs4DKaAiITlo2uHz7EAbYrZnBIo0FBTAAA | 2026-10-30T14:51:53Z | `bfd.mgnify30.metaeuk30.smag30.a3m` (2136); `uniref.a3m` (3026) | 5163 | 165.35 |
| Q01915 | 508 | 890ce3d93c22 | https://api.colabfold.com | env | x1o2_dL6g1U5hl63VvdbmkBSflWOkMq0Uzn2Uw | 2026-10-30T14:49:25Z | `bfd.mgnify30.metaeuk30.smag30.a3m` (13976); `uniref.a3m` (6737) | 20714 | 919.04 |
| A0A0H3GPJ3 | 349 | 89fac1d210fc | https://api.colabfold.com | env | Jg_huYWuOAxxe7j10r_8LDz8MrgWlWtnfxisKA | 2026-10-30T14:37:14Z | `bfd.mgnify30.metaeuk30.smag30.a3m` (7767); `uniref.a3m` (8220) | 15988 | 855.82 |
| O68883 | 427 | 97c992c19d18 | https://api.colabfold.com | env | jaG_TQAqon5fVCRPZgeeOEPgTgoBhiqDzYf47Q | 2026-10-30T14:57:56Z | `bfd.mgnify30.metaeuk30.smag30.a3m` (11904); `uniref.a3m` (6739) | 18644 | 902.25 |
| Q9VLI0 | 982 | 9bedfef15325 | https://api.colabfold.com | env | vh3lal_uKExWFs2gUl8xv1J-c8F0LigZqCfjSg | 2026-10-30T15:06:31Z | `bfd.mgnify30.metaeuk30.smag30.a3m` (3140); `uniref.a3m` (4202) | 7343 | 234.32 |
| A0A0H3GUM2 | 249 | 9c54a22d618f | https://api.colabfold.com | env | rpaceIjT6M0PxcvKDevJJnw-e4wdGXtOLJJYmg | 2026-10-30T14:38:54Z | `bfd.mgnify30.metaeuk30.smag30.a3m` (10277); `uniref.a3m` (7789) | 18067 | 1144.95 |
| C0NAV3 | 394 | a36e1052dbc1 | https://api.colabfold.com | env | tx76mJits8cCYWJtrRwDOkYSew3qonxI5QoIhA | 2026-10-30T14:48:49Z | `bfd.mgnify30.metaeuk30.smag30.a3m` (710); `uniref.a3m` (989) | 1700 | 85.64 |
| A8E7K1 | 929 | a64126ae2cb2 | https://api.colabfold.com | env | e608YzYE2qv_TSe9nKbXFc4DkMmlEozpSGKfqQ | 2026-10-30T14:49:23Z | `bfd.mgnify30.metaeuk30.smag30.a3m` (1222); `uniref.a3m` (7904) | 9127 | 299.45 |
| A0A1D8PRG1 | 331 | a781c785dfdd | https://api.colabfold.com | env | bWzYRQE9SgJbQWhNlJ1eWOrNlie5crLFawHtBw | 2026-10-30T14:29:36Z | `bfd.mgnify30.metaeuk30.smag30.a3m` (1586); `uniref.a3m` (1847) | 3434 | 188.75 |
| A0A0H3GSY7 | 400 | ad53bb8e462c | https://api.colabfold.com | env | plf_9x895Vsg_yciPmx7c4ZC2rbfr5rAqeW_5w | 2026-10-30T14:38:33Z | `bfd.mgnify30.metaeuk30.smag30.a3m` (13811); `uniref.a3m` (8871) | 22683 | 1134.15 |
| G5EBZ1 | 935 | af78eb125bd9 | https://api.colabfold.com | env | FPoEZ55c2gJmh_Gd2jYiol1_W13Gh7w4bDdNCw | 2026-10-30T14:53:22Z | `bfd.mgnify30.metaeuk30.smag30.a3m` (1177); `uniref.a3m` (1539) | 2717 | 88.86 |
| A0A8R1TKT0 | 701 | b16552c143d4 | https://api.colabfold.com | env | s5O2PT6P0n-NC0gC3qY4xuzXf7t-MdtLHplbbw | 2026-10-30T14:45:59Z | `bfd.mgnify30.metaeuk30.smag30.a3m` (11832); `uniref.a3m` (6198) | 18031 | 681.02 |
| C0NJ85 | 231 | c3f7c22c8a1f | https://api.colabfold.com | env | dGb2cG5Scputp5JTjuURURi1zzVEwQMWHpvznA | 2026-10-30T14:50:48Z | `bfd.mgnify30.metaeuk30.smag30.a3m` (5592); `uniref.a3m` (3949) | 9542 | 627.82 |
| Q19673 | 693 | c41fbbe6e0ea | https://api.colabfold.com | env | YWH2raD3mVOdaRJya91obD6I8axQLIn2mLNRRw | 2026-10-30T15:09:31Z | `bfd.mgnify30.metaeuk30.smag30.a3m` (11834); `uniref.a3m` (6546) | 18381 | 698.24 |
| A0A0D2H481 | 410 | c6b34b21bba6 | https://api.colabfold.com | env | jxOEixCwkxuQ8oxGuZ066uVmAlSZxylrEA9xmg | 2026-10-30T14:27:54Z | `bfd.mgnify30.metaeuk30.smag30.a3m` (1606); `uniref.a3m` (1870) | 3477 | 171.72 |
| Q24451 | 1108 | c9adc660f0a8 | https://api.colabfold.com | env | lMiRgS7K9s27oLJ2aRCxM2gDvVk8ANBWvc-Nng | 2026-10-30T14:55:21Z | `bfd.mgnify30.metaeuk30.smag30.a3m` (2604); `uniref.a3m` (3157) | 5762 | 173.1 |
| A0A1D8PDL6 | 241 | c9f33c806deb | https://api.colabfold.com | env | DfZdYwU4_uYVcZPeAaQadxf4nDRPM1fWFq207g | 2026-10-30T14:28:35Z | `bfd.mgnify30.metaeuk30.smag30.a3m` (3106); `uniref.a3m` (2710) | 5817 | 374.71 |
| Q7CPE1 | 513 | cf0f126805ea | https://api.colabfold.com | env | FMjcjzyRMGqiC4xymOcg0bEaqYvxPHt3y-BI1A | 2026-10-30T15:01:56Z | `bfd.mgnify30.metaeuk30.smag30.a3m` (14183); `uniref.a3m` (5487) | 19671 | 868.5 |
| A0A8R1Y3G8 | 641 | d4543e5e0a4b | https://api.colabfold.com | env | BdwBRNqGnXvfmEvQw-pc0-PlEiG4FByAKOXp8Q | 2026-10-30T14:45:48Z | `bfd.mgnify30.metaeuk30.smag30.a3m` (8948); `uniref.a3m` (5931) | 14880 | 587.72 |
| A0A3Q0KJG9 | 371 | d6ba551627b4 | https://api.colabfold.com | env | RvicP7afHCeo_bJs3OSVdKhvdR2vR5tyx3YwsA | 2026-10-30T14:42:22Z | `bfd.mgnify30.metaeuk30.smag30.a3m` (959); `uniref.a3m` (843) | 1803 | 93.61 |
| P44945 | 1121 | d98e35d7eaa3 | https://api.colabfold.com | env | Ui0LfICFrUXDjOH07oP60Wsi1k2YGXwaBlaKNQ | 2026-10-30T15:05:40Z | `bfd.mgnify30.metaeuk30.smag30.a3m` (5551); `uniref.a3m` (3989) | 9541 | 284.96 |
| P55995 | 835 | da3d40051f55 | https://api.colabfold.com | env | AF60Z37spjjR5lOL_YWfVdyuCgDLwQ4A7OsQKw | 2026-10-30T15:08:33Z | `bfd.mgnify30.metaeuk30.smag30.a3m` (13801); `uniref.a3m` (8479) | 22281 | 771.07 |
| P0AFB5 | 349 | e4d22c03b5d3 | https://api.colabfold.com | env | 5RUEz4tDjxtPer5hSYBgVsg_YnoNAKFNiuq7dQ | 2026-10-30T15:03:57Z | `bfd.mgnify30.metaeuk30.smag30.a3m` (7702); `uniref.a3m` (7142) | 14845 | 794.63 |
| O00411 | 1230 | e79aa9b48007 | https://api.colabfold.com | env | ixyFTxREWw-SyaOMBfiAf8AmixBPMzV_6dme-A | 2026-10-30T15:00:17Z | `bfd.mgnify30.metaeuk30.smag30.a3m` (3155); `uniref.a3m` (3535) | 6691 | 190.78 |
| I1KY39 | 513 | e7ff353ac924 | https://api.colabfold.com | env | z9_eivlC9khw8MgVI3i0WI2fJbzoCHoS9SduPQ | 2026-10-30T14:55:54Z | `bfd.mgnify30.metaeuk30.smag30.a3m` (13167); `uniref.a3m` (7379) | 20547 | 907.17 |
| P33570 | 667 | ed002aca245e | https://api.colabfold.com | env | 4sayPppWlFBHhogrMV486dyZBUZ7jhUGzGgtiQ | 2026-10-30T14:42:17Z | `bfd.mgnify30.metaeuk30.smag30.a3m` (15815); `uniref.a3m` (8276) | 24092 | 932.85 |
| A0A0H3GR39 | 585 | ed980437f2b6 | https://api.colabfold.com | env | UmywatzpQr6DwqnrmEriWjaqnwk1xdZBgp-nEQ | 2026-10-30T14:34:52Z | `bfd.mgnify30.metaeuk30.smag30.a3m` (3912); `uniref.a3m` (4794) | 8707 | 359.99 |
| G5ECB1 | 253 | f12e17832e35 | https://api.colabfold.com | env | UmrRvGQ38Eg0faUhiaDiadTeHFnX-0rFy1JVvw | 2026-10-30T14:55:54Z | `bfd.mgnify30.metaeuk30.smag30.a3m` (873); `uniref.a3m` (2197) | 3071 | 193.07 |
| Q9JKY0 | 299 | f79d1fe6b57d | https://api.colabfold.com | env | L3q9F2gLBfOTmyG04ydK4bN8BLhGdp342Gj_rg | 2026-10-30T15:02:57Z | `bfd.mgnify30.metaeuk30.smag30.a3m` (1607); `uniref.a3m` (1487) | 3095 | 178.99 |
| A0A1C1CNY1 | 239 | f86e4f810d11 | https://api.colabfold.com | env | gaJ5hztl57dgHeYqREm7rF2_XwhDFWoIKZuVEw | 2026-10-30T14:40:54Z | `bfd.mgnify30.metaeuk30.smag30.a3m` (5277); `uniref.a3m` (4051) | 9329 | 603.44 |
| A0A158Q514 | 702 | fff70c95d58b | https://api.colabfold.com | env | mDM9K7TkC8OMRYVPc5LOEyOGlvTIQgOjtVk9Yg | 2026-10-30T14:45:21Z | `bfd.mgnify30.metaeuk30.smag30.a3m` (11768); `uniref.a3m` (6286) | 18055 | 681.44 |

## Cross-OG site matching (exploratory D31; CROSSOG.md)
| dataset | source / path | version | checksum | license | notes |
|---|---|---|---|---|---|
| Foldseek DB of AFDB v6 models (reused from Part A) | /tmp/pe_big/partA/foldseek/afdb_v6 (built with Foldseek 10.941cd33 `createdb`) | AFDB v6, tier 1/2 | afdb_v6.lookup sha256 8004ffe7c9321359… | CC BY 4.0 (AFDB) | regenerable, not committed; tier 3/L models absent (O2 reps restricted to tier 1/2) |
