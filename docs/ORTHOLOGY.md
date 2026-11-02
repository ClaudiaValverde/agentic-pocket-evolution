# Orthology O2: Ancestra on tier 1 + 111 context proteomes, ID-gap mapping, placement (2026-10-30)

P = <repo>. Config: `config/orthology.yaml`. Snapshot: `data/snapshot/20261030T085136Z`.
Outputs: `data/orthology/ancestra_O2/` (116 MB). Bulk / caches: `/tmp/pe_big/ancestra/o2/` (regenerable, not persistent).

## Pipeline (all re-runnable, parameterized by the snapshot manifest + config)
```bash
.venv/bin/python scripts/orthology/select_context_species.py        # 3 min  -> context_species.tsv, o2_proteomes.tsv
.venv/bin/python scripts/orthology/map_idgap.py                     # 20 min (UniProt REST) -> idgap_mapping.parquet
scripts/orthology/run_ancestra_o2.sh 40 accurate                    # 98.8 min -> accurate/hogs.orthoxml.gz (+verify)
PYTHONPATH=src envs/ancestra/bin/python scripts/orthology/export_hog_tables.py --orthoxml data/orthology/ancestra_O2/accurate/hogs.orthoxml.gz \
   --tree data/orthology/ancestra_O2/o2_guide_tree.nwk --taxdump /tmp/pe_big/ancestra/taxdump --out data/orthology/ancestra_O2   # 25 s
.venv/bin/python scripts/orthology/place_proteins.py --set tier2    # 3 min ; --set context 8 min ; --set tier3 --fasta <trembl.fa.gz>
PYTHONPATH=src:external/ancestra-mvp-ai .venv/bin/python scripts/orthology/validate_orthology.py   # 2 min
PYTHONPATH=src envs/ancestra/bin/python scripts/orthology/crosscheck_ancestra.py
.venv/bin/python -m pytest tests/test_orthology.py                  # toy OrthoXML with known answers (14 pass, 1 needs ancestra env)
```
Run-time helpers: `scripts/orthology/_common.py`, `src/pocketevo/taxonomy.py`, query API `src/pocketevo/orthology.py`.

## 1. Context species (`select_context_species.py`; rule in the docstring)
Rule revised by the PI before the run (atlas-complete eukaryotes are only protists/early fungi and cannot densify
clades near human/fly/plants). Two parts, total protein budget 1.55 M (tier 1 = 557,396):
- **(a) 40 eukaryotic densification proteomes, chosen regardless of atlas coverage** (`pockets_source=local` for
  39: pockets to be computed from the Foldcomp v4 DB). 27 clade slots in config (Primates, Laurasiatheria, Metatheria,
  Aves, Lepidosauria, Amphibia, Actinopterygii, Chondrichthyes, Cyclostomata, Insecta×4, Arachnida, Nematoda×3,
  Platyhelminthes, Mollusca, Echinodermata, Cnidaria, Choanoflagellata, Embryophyta×5, Chlorophyta, Saccharomycotina×2,
  Pezizomycotina×3, Basidiomycota×2, Mucoromycota, Kinetoplastea, Apicomplexa, Stramenopiles×2, Amoebozoa). In a slot:
  ≤1 per order and genus, ≤32k proteins, priority UniProt mnemonic > #isoforms > taxid; **accepted only if ≥0.8 of the
  proteome's accessions are in `fc_entries_afdb_uniprot_v4.parquet` or the atlas** (measured: min 0.83, median 0.94).
  Many recent proteomes fail this (new accessions after AFDB v4: chicken 0.29, Xenopus 0.42, salmon/carp/lamprey/
  hagfish ~0) and are skipped. Cyclostomata: 0/1 (no lamprey/hagfish proteome passes). Picks include capuchin, fruit bat,
  koala, zebra finch, wall lizard, caecilian, tilapia, ghost shark, Nasonia, D. pseudoobscura, Thrips, cowpea weevil,
  Varroa, Steinernema, Trichinella, Soboliphyme, Opisthorchis, Lottia, sea urchin, Nematostella, Salpingoeca,
  Physcomitrium, Aquilegia, peach, cacao, walnut, Chlamydomonas, Cryptococcus, Puccinia, Fusarium (G. moniliformis),
  Monilinia, Trichophyton, two Saccharomycotina, Phycomyces, L. major, P. knowlesi, Saprolegnia, Cafeteria, Acanthamoeba.
- **(b) 71 atlas-complete taxa** (`pockets_source=atlas`): ≥0.9 of the proteome's accessions present in the snapshot
  (measured overlap min 0.905, median 0.989); groups protists 8, bacteria 33 (of quota 60; budget-limited), archaea 30;
  round-robin over phyla, priority = #standard ranks shared with the nearest tier-1 species, ≤6 per class, ≤1 per genus.
- Candidates evaluated 171 (`context_considered.tsv`, with reason); 157 proteomes, 1,549,554 proteins in O2.

## 2. Tier-1 ID gap (`map_idgap.py` → `idgap_mapping.parquet`, every tier-1 atlas accession)
37,535 / 552,810 tier-1 atlas accessions are not proteome members (29,086 active TrEMBL, 8,449 deleted). Methods:
direct 515,275 · identical sequence in same proteome 25,135 · same unique gene name + DIAMOND ≥95 % id, ≥0.8 qcov 3,380 ·
DIAMOND only 3,997 · gene name only (low confidence, not used by `relation()`) 1,040 · unmapped 3,983 (rat 2,026,
zebrafish 1,737). High-confidence link: **99.1 %** of tier-1 atlas accessions (99.1 % of those with pockets). Atlas
proteins with pockets that sit in a multi-gene O2 HOG: 85.7 % direct → **91.7 %** with the mapping. 497 members receive
>1 atlas accession (max 14; redundant gene models — duplicates of one gene, not paralogs; dedupe per gene downstream).

## 3. Ancestra O2 run (accurate preset; 40 threads on a shared node, load 40–125)
| stage | wall | peak RSS |
|---|---|---|
| S0 families | 26.5 min | 11.3 GB |
| S2 scores | 38.2 min | 46.9 GB |
| S3 merge (accurate, FAMSA+FastTree local trees) | 32.8 min | 14.0 GB |
| output | 1.1 min | 9.2 GB |
| **total** | **98.8 min, 24.7 CPU-h (15× CPU)** | **48 GB** |
`fast` preset reusing S0/S2 checkpoints: merge 5.5 min (16 threads), saved in `fast/` for comparison only.
`verify_orthoxml.py`: PASS (XSD, geneRef integrity, each gene once, TaxRange = LCA, pairs count, pyham load).
Result: 57,332 root (top-level) groups (52,119 orthologGroups + 5,213 lineage-specific paralogGroups; pyham counts
56,120 top-level HOGs), 638,734 orthologGroups + 176,999 paralogGroups, 1,234,439 / 1,549,554 proteins (79.7 %) in
groups, 37.5 M ortholog pairs, largest root 5,324 genes. Per-species fraction in groups median 0.82; lowest:
Reticulomyxa filosa 0.41, Giardia muris 0.43, Steinernema 0.53, Tritrichomonas 0.53, Trichinella 0.53, Puccinia 0.53.
Fast preset: 58,065 roots, 673,703 S / 174,651 D groups (more, smaller groups).

## 4. Outputs (`data/orthology/ancestra_O2/`)
- `accurate/hogs.orthoxml.gz` (24 MB), `manifest.json`, `time.txt`, `verify.txt`; `fast/` same for the fast preset.
- `tables/hogs.parquet`: hog, hog_name, parent, root, kind (S/D), level, **level_taxid**, depth, n_genes (815,733 rows).
- `tables/protein_hogs.parquet`: protein → species, root HOG (−1 = singleton), **path** = nested HOG ids root→deepest
  (join hogs.parquet for levels/taxids). = deliverables (a) and (b).
- (c) `pocketevo.orthology.relation(a, b)` / `HogIndex.relation` → ortholog | paralog | unrelated | same | unknown, MRCA
  HOG id, level and NCBI taxid (deepest group containing both; QfO flattening identical to Ancestra's hogs_to_pairs).
  Atlas accessions are resolved through the ID-gap mapping automatically.
- (d) `o2_guide_tree.nwk` (names as used by Ancestra), `tables/guide_tree_taxids.nwk` (all nodes = NCBI taxids),
  `tables/levels.tsv` (level name → taxid, rank, parent).
- `context_species.tsv`, `context_considered.tsv`, `context_selection.json`, `o2_proteomes.tsv` (proteome IDs + md5).
- `idgap_mapping.parquet`, `idgap_summary.tsv`, `idgap_by_taxon.tsv`.
- `placement/{tier2,context}/part-*.parquet` + `placement.json` (provenance, counts, runtime).
- `validation/validation.json`, `panther_pairs.parquet`, `crosscheck_ancestra.json`.

## 5. Placement (`place_proteins.py`)
DIAMOND 2.2.8 `--fast`, e ≤ 1e-10, query cov ≥ 50 %, k = 25, vs all 1.55 M O2 proteins; best bitscore hit → root HOG;
deepest compatible HOG = deepest group on the hit's path whose level taxid is in the query's NCBI lineage;
`second_root_ratio` flags ambiguous family assignment. DIAMOND hits are cached by query content, so new HOG tables
only re-assign (1 min for tier 2). Tier 3: `--set tier3 --fasta <FASTA>` streams in chunks of 500k (resumable
parquet parts; tested on a 1.9k chunk). Sequences: Swiss-Prot 2026_03 FASTA; UniProt REST per context taxon; UniParc
(`dbid:`) for deleted/merged accessions (16 tier-2, 85.7k context accessions).
| set | atlas acc. | members | DIAMOND placed | no hit | in a HOG | compatible (of placed in HOG) | 2nd-root ratio > 0.9 | DIAMOND time |
|---|---|---|---|---|---|---|---|---|
| tier 2 (Swiss-Prot) | 436,886 | 9,479 | 386,639 (90.5 %) | 40,768 | 89.7 % | 98.9 % | 7.5 % | 99 s |
| context taxa | 466,282 | 363,763 | 89,821 | 12,698 | 66.4 % | 99.9 % | 7.4 % | 44 s |

## 6. Validation (`validation/validation.json`)
- Toy OrthoXML unit tests (`tests/test_orthology.py`): 15/15 pass in the ancestra env (relations, levels, roundtrip,
  aliases, agreement with Ancestra's hogs_to_pairs).
- Cross-check of `relation()` with Ancestra's flattening on 50k ortholog + 50k paralog pairs (of 37.5 M / 146.9 M): 50,000/50,000 and 50,000/50,000 agree.
- Known pairs 13/13: HSP90AA1–HSP90AB1 paralog (dup. at Gnathostomata), HSP90AA1 human–mouse/–yeast HSP82 ortholog,
  EGFR–ERBB2 paralog, EGFR/HBB/TP53/BRCA1/HSP90AB1 human–mouse ortholog, HBA–HBB paralog, ACTB–ACTG1 paralog,
  2 negative controls unrelated. Caveat: ACTB–ACTG1 duplication placed at Homo sapiens (should be vertebrate-old) and
  HBA–HBB at Deuterostomia: near-identical / fast-evolving families get wrong duplication levels.
- Human–mouse: 26,969 ortholog pairs, 16,585 1:1 (80.3 % of human proteins 1:1; 88.8 % have a mouse ortholog).
  O1 (tier 1 only): 23,233 / 16,647 — more context species add many-to-many pairs but leave 1:1s unchanged.
- PANTHER (UniProt xref_panther; 2,000 pairs, seed 0, 8 PANTHER-annotated tier-1 species): same PANTHER family for
  O2 orthologs 86.3 % (same subfamily 72.6 %), paralogs 78.8 % (subfamily 25.7 %), cross-root pairs 0.0 %.
  Reverse view: of 20,000 human–mouse pairs sharing a PANTHER subfamily, 75.0 % are O2 orthologs, 9.2 % paralogs,
  15.8 % in different root HOGs (Ancestra family splits → recall limitation; O1 values were similar).

## 7. Limitations
- Accuracy not benchmarked on this species set; Ancestra's QfO recall < OrthoFinder; deep-level splits remain
  (15.8 % of PANTHER-subfamily human–mouse pairs split across root HOGs). Flag deep-level gains as possible artifacts.
- Duplication levels are unreliable for near-identical or fast-evolving families (actins, globins).
- Context (a) proteomes have no atlas pockets yet (`pockets_source=local`); 5–17 % of their accessions lack a v4 model.
- Bacteria are under-sampled (33; protein budget went to eukaryotes); no cyclostome passed the v4 coverage check.
- Placement is best-hit (no HMM / phylogenetic placement); ~7.5 % of placements are ambiguous between root HOGs;
  9.5 % of Swiss-Prot proteins have no hit (orphans, short peptides) — presence-only evidence, never absence.
- 3,983 tier-1 atlas accessions unmapped (deleted rat/zebrafish models with no ≥95 % match); 1,040 gene-name-only.
- `run_ancestra_o2.sh` hard-codes /tmp paths (set at top); the node was shared (load up to 125), so runtimes are upper bounds.
- Release drift: everything is UniProt 2026_03; AFDB v4 tier-3 accessions deleted since then are recovered via UniParc.
