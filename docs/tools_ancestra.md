# Ancestra (orthology / HOGs) — what it is, what is precomputed, how we run it

Repo: `P/external/ancestra-mvp-ai` (A below; commit ffb0bbe, method frozen at `freeze-v1`, `A/results/freeze.json`).
P = <repo>. Install: `P/docs/ENV_ancestra.md`. Written 2026-10-30.
Caveat from the authors: everything in A is AI-generated and not yet human-verified (`A/README.md` top,
`A/docs/DISCLOSURES.md`, `A/docs/VERIFICATION.md`). Numbers marked **[measured here]** were run by me on this node.

## 1. What it does (A/README.md "What Ancestra does"; A/docs/core_v0.md)
Bottom-up reconstruction of ancestral gene sets along a **given rooted guide species tree**, one level at a time:
1. S0 `families` (`A/ancestra/families.py`): DIAMOND clustering → all proteins vs cluster reps → connected
   components → Leiden split of giant components. Reference-free (no OMA/Pfam DB).
2. S2 `scores` (`A/ancestra/scoring.py`): similarity only within families (DIAMOND, MMseqs2 re-align of missed pairs;
   pairs in huge families are capped — 22 families hit the cap in our tier-1 run).
3. S3 `merge` (`A/ancestra/merge.py`, `localtree.py`): at each guide-tree node, merge child-lineage ancestral genes by
   mutual best matches + in-paralog absorption + repair; ambiguous many-to-many clusters resolved with a small local
   tree (`accurate`: FAMSA+FastTree; `fast`: NJ) labelled by species overlap. Each merged gene = a HOG at that level.
4. Optional layers: gene order (`--gene-order`, no measurable gain), 3Di structure rescue (`--family-rescue 3di`,
   *lowers* OrthoBench precision, off by default). Expression/isoform layers were dropped.
Sequence-only by default — structures are **not** needed (our AFDB models are not an input).

## 2. Inputs / outputs (A/README.md "Input and output formats"; A/ancestra/pipeline.py, io.py, orthoxml.py)
**Input**
- `--proteomes DIR`: one FASTA per species (`.fa/.fasta/.faa/.fas/.pep`, `.gz` ok). Species name = file stem
  (`UP000005640_9606.fasta.gz` → `UP000005640_9606`). `sp|ACC|NAME` headers → `ACC` (UniProt accession works
  directly). IDs unique across species. **One protein per gene expected** (no isoform selection).
- `--tree FILE`: rooted Newick; internal node names = HOG levels (OrthoXML TaxRange); polytomies allowed (k-way merge).
  Build from NCBI taxonomy: `python -m ancestra.speciestree ncbi --proteomes DIR/*.fa.gz --taxdump TAXDUMP
  --leaf-naming qfo --out tree.nwk` (needs `nodes.dmp`+`names.dmp`; the atlas `ncbi/` dir lacks them → I used a fresh
  NCBI taxdump in /tmp/pe_big/ancestra/taxdump). Random tree costs 16 F points on OrthoBench → tree matters.
- `--gene-order DIR` (optional parquet per species; schema `A/ancestra/genomes.py`).

**Output** (`--out`)
| file | content |
|---|---|
| `hogs.orthoxml` | primary: OrthoXML 0.5, nested `orthologGroup` (speciation, TaxRange = level) / `paralogGroup` (duplication, `taxonId` = level just below the duplication); XSD-validated; loads in pyham |
| `pairs.tsv` | `#gene_a gene_b level` ortholog pairs (QfO flattening: deepest common group is an orthologGroup) |
| `orthobench_ogs.txt` | root-level groups (= orthogroups/gene families), one line per group |
| `manifest.json` | versions, argv, config, input md5, per-stage wall/RSS, per-node merge stats |
| `stages/` | checkpoints (input, families, scores, merge, output); reruns skip unchanged stages |

So it yields: **orthogroups** (root HOGs), **HOGs at every taxonomic level**, **pairwise orthologs** (and paralogs via
`hogs_to_pairs(with_paralogs=True)`), and a **gene "tree" per family as nested HOGs with S/D labels** (topology only:
multifurcating, **no branch lengths**; the local FastTree/NJ trees used in the merge are **not saved**). No separate
reconciled Newick gene trees, no species-tree inference in the default run (the guide tree is an input; species-tree
code with ASTRAL-Pro/concat exists in `A/ancestra/speciestree.py` but ASTRAL is not in the env — untested here).

## 3. API (demo `P/scripts/ancestra_api_demo.py` run on the example; same calls used for the tier-1 stats in §7)
```python
from ancestra.orthoxml import read_orthoxml, hogs_to_pairs, hogs_at_level, Group
ox = read_orthoxml("hogs.orthoxml.gz")          # tier-1: ~10 s; QfO 989k genes: 14 s, 0.7 GB RSS [measured here]
root_of = {p: i for i, g in enumerate(ox.groups) if isinstance(g, Group) for p in g.iter_genes()}  # orthogroup
orth, para = hogs_to_pairs(ox, with_paralogs=True)   # sets of (a,b), a<b; orth with with_level=True -> {pair: level}
hogs = hogs_at_level(ox, "Bilateria")            # list of gene sets = HOGs at a level (partition of that clade)
ox.species_of("P07900"); ox.species_tree          # species of a protein; guide tree (Node) from <taxonomy>
# gene tree with labels: recurse Group.children; kind "ortholog"=speciation, "paralog"=duplication, .level = taxon
```
Relation of two proteins: ortholog iff `(a,b)` in `orth`; paralog iff in `para`; otherwise different orthogroups
(or singleton). Gains/losses/duplications per branch: **pyham** (installed): `pyham.Ham(tree.nwk, hogs.orthoxml,
use_internal_name=True)`, `compare_genomes_vertically(extant, ancestral).get_lost()/get_gained()/get_duplicated()`.
CLI: `ancestra run`, `ancestra eval`; `A/scripts/orthoxml_to_pairs.py`, `A/scripts/verify_orthoxml.py`,
`A/scripts/verify_ids.py`. There is no query server / per-protein lookup API beyond these functions.

## 4. Precomputed predictions in the repo (A/predictions/README.md)
| run | species | proteins | IDs | files |
|---|---|---|---|---|
| `predictions/orthobench/{accurate,fast}` | 12 OrthoBench animals | 251,378 | **Ensembl** protein IDs (Ensembl 99) | hogs.orthoxml.gz, pairs.tsv.gz, orthobench_ogs.txt.gz |
| `predictions/qfo2022/{accurate,fast}` | 78 QfO 2022_02 ref. proteomes (`A/inputs/qfo2022/inputs.json`) | 988,778 | **UniProt accessions** | hogs.orthoxml.gz only (pairs regenerable) |
Plus guide trees and gene-order tables in `A/inputs/`. QfO accurate: 33,465 root groups, 741,738 genes (75%) in
multi-gene groups, max group 6,026 genes.

**Overlap with our atlas [measured here]** (tier-1 pocket parquet: 453,269 accessions with ≥1 kept pocket, 46 taxids):
- QfO run covers **19 / 46** tier-1 taxids (human, mouse, rat, zebrafish, fly, worm, yeast, S. pombe, C. albicans,
  Arabidopsis, rice, maize, Dicty, P. falciparum, E. coli, M. tuberculosis, P. aeruginosa, H. pylori, M. jannaschii).
  Missing: soybean and all neglected-disease/pathogen proteomes (Trypanosoma, Leishmania, Schistosoma, filarial
  nematodes, fungi of mycetoma/chromoblastomycosis, many bacteria…).
- Accessions: 218,764 / 453,269 atlas proteins (48%) are in the QfO run; within the 19 shared taxa 218,764 / 231,369
  (94.6%; zebrafish only 63% — QfO used an updated zebrafish proteome).
- Verdict: usable as a cross-check, **not** as our orthology (half the atlas missing, 2022 IDs, 78-species tree that
  does not match ours). The OrthoBench run (Ensembl IDs, 12 animals) is irrelevant for us.

## 5. Benchmarks and accuracy (A/results/orthobench_final/summary.md, A/results/qfo_final/summary.md, README)
- **Open OrthoBench** (70 RefOGs, 12 animals; dev/test split `A/splits/`): test F **76.2** accurate / 74.8 fast vs
  OrthoFinder 3.2 69.7, SonicParanoid2 70.4, OMA HOGs 60.4; CIs wide (36 test RefOGs), non-inferiority to
  OrthoFinder **not shown**; dev numbers optimistic (heavy reuse).
- **QfO 2022** (78 sp., official code run locally): SwissTree TPR/PPV 0.721/0.926, TreeFam-A 0.583/0.921, VGNC
  0.980/0.982, GO 0.475, EC 0.899, STD-Euk 0.113. **Recall clearly below OrthoFinder** (SwissTree 0.825, TreeFam-A
  0.737), similar to OMA/FastOMA; precision similar. Suspected causes: deep-level splits (Bilateria/Eukaryota), S0
  family splits (family completeness 0.72 vs ceiling 0.82), vertebrate over-merging.
- Other caveats: output depends slightly on Newick child order (DISCLOSURES #11); no human verification.

## 6. Runtime and scaling
Reported (A/results/scaling/summary.md, A/results/scaling_large/summary.md; shared nodes, 8–12 threads):
| set | proteins | preset | wall | CPU-h | RSS |
|---|---|---|---|---|---|
| QfO 64 sp | 816k | accurate | 2.65 h (8 thr) | 16.3 | 25 GB |
| QfO 78 sp | 989k | accurate / fast | 2.83 / 1.88 h (8 thr) | 18.5 / 11.6 | 29 GB |
| UniProt 250 ref. proteomes | 1.08 M | accurate / fast | 1.91 / 1.12 h (12 thr) | 17.7 / 8.6 | 23 GB |
| UniProt 500 | 2.26 M | fast | 3.3 h (12 thr) | 23 | 69 GB |
| UniProt 1,000 | 4.75 M | fast | 10.6 h (12 thr) | 64 | 221 GB |
CPU-h exponent ~1.16 (QfO subsets) → 1.35 (≥125 proteomes); S0 graph step ~quadratic at the top end; S2 memory is
the first limit (extrapolation to all AFDB, ~200 M proteins: ~10k CPU-h and tens of TB RAM → needs re-engineering).
Accuracy was only benchmarked at 12 and 78 species.

**Measured here (2026-10-30, 64-core node shared with other agents, load average 30–93 during the run):**
| run | wall | CPU-h | peak RSS |
|---|---|---|---|
| unit tests `pytest -m "not slow"` | 61 s | — | — (147 passed, 6 skipped) |
| `example/run_example.sh 4` (6 toy sp., 1,828 proteins, both presets) | 48 s | — | <0.3 GB; F1 0.9530/0.9545 = README |
| **46 tier-1 proteomes, 557,396 proteins, accurate, 48 threads** | **20.3 min** (S0 4.4, S2 6.8, merge 8.7, output 0.3 min) | 6.6 | 10.8 GB |
| same, `fast` preset re-using S0/S2 checkpoints | 65 s (merge 44 s) | <0.1 | 3.6 GB |
| `verify_orthoxml.py` on the tier-1 output | 38 s | — | all checks PASS (XSD, geneRefs, LCA, pairs, pyham) |
→ a fresh fast run ≈ 12 min. Input/command:
```bash
export PATH=P/envs/ancestra/bin:$PATH NUMBA_CACHE_DIR=/tmp/pe_big/ancestra/numba TMPDIR=/tmp/pe_big/ancestra/tmp
python -m ancestra.speciestree ncbi --proteomes /tmp/pe_big/ancestra/tier1_fa/*.fasta.gz \
    --taxdump /tmp/pe_big/ancestra/taxdump --leaf-naming qfo --out tier1_guide_tree.nwk
python -m ancestra run --proteomes /tmp/pe_big/ancestra/tier1_fa --tree tier1_guide_tree.nwk \
    --out /tmp/pe_big/ancestra/run_t1 --threads 48 --preset accurate
```
Proteomes = UniProt current-release one-protein-per-gene reference proteomes for the 46 tier-1 taxids
(UP ids: `P/data/orthology/ancestra_tier1_v0/ref_proteomes.tsv`).

## 7. Tier-1 result (preliminary, `P/data/orthology/ancestra_tier1_v0/`)
Saved: `{accurate,fast}/hogs.orthoxml.gz` (7.7 MB each), `manifest.json`, `run.log`, `accurate/verify.json`,
`tier1_guide_tree.nwk`, `ref_proteomes.tsv`. Full run dirs (pairs.tsv 194 MB, checkpoints 2.9 GB) stay in
`/tmp/pe_big/ancestra/run_t1{,_fast}` (node-local /tmp: **not persistent**; regenerate with the commands above).
- 84,845 S0 families; 36,736 root groups (orthogroups); 226,941 orthologGroups + 71,726 paralogGroups;
  451,629 / 557,396 proteins (81%) in multi-gene HOGs; 6.54 M ortholog pairs.
- Per-species fraction in HOGs: median 0.83; lowest P. falciparum 0.43, M. jannaschii 0.44, Dictyostelium 0.46,
  S. aureus 0.61, rice 0.62 (isolated lineages / orphan-rich genomes → many singletons).
- Sanity: human–mouse 23,233 ortholog pairs (16,647 1:1). HSP90AA1 (P07900)–yeast HSP82 (P02829) ortholog,
  HSP90AA1–HSP90AB1 paralog, EGFR–ERBB2 paralog, EGFR–mouse Egfr ortholog, human–mouse HBB ortholog: all correct.
- **ID coverage vs atlas:** 421,492 / 453,269 (93.0%) atlas tier-1 accessions are in these proteomes, 81.4% in a
  multi-gene HOG. Misses concentrate in zebrafish (54% present), maize (70%), rat (75%), fly (90%): AFDB v6 models
  include accessions that are not the gene-centric canonical entry (redundant TrEMBL / older entries). Not fixed yet.

## 8. Running on a new set of proteomes (recipe)
1. One FASTA per species, one protein per gene, file stem = species id (use `UPxxxx_<taxid>` so `--leaf-naming qfo`
   works). Ensure every atlas accession of interest is present (see §7 gap) — either add missing accessions
   (they will appear as species-specific in-paralogs → inflates duplications) or map them to their canonical
   gene (UniProt gene name / ≥95% identity) after the run. Proposed: map, do not add.
2. Build the guide tree from NCBI taxdump (`ancestra.speciestree ncbi`); check leaves = stems (`A/scripts/verify_ids.py`).
3. `ancestra run ... --threads 48 [--preset fast]`; work dir in /tmp (home quota!).
4. `A/scripts/verify_orthoxml.py OUT`; copy `hogs.orthoxml.gz` + manifest to P.
5. Changing only the preset/merge config reruns only merge+output (stage cache; hard-link `stages/{input,families,scores}`).

## 9. Limitations relevant to us
- Recall below OrthoFinder/SonicParanoid on QfO; deep-level (Eukaryota/LUCA) HOGs are fragmented → "gains" at deep
  nodes may be orthology splits, not real gains. Treat deep-level calls with care; cross-check with PANTHER/eggNOG.
- Needs complete proteomes and a species tree; **Swiss-Prot (tier 2) is not proteomes**: 408,716 tier-2 pocket
  proteins span 10,839 taxa (median 2 proteins per taxon) → running Ancestra on it would make absence ≠ loss and
  give a ~11k-leaf tree (untested).
- One protein per gene; isoforms/fragments not handled (AFDB human fragments are fine: same accession).
- Gene trees are nested HOGs without branch lengths; for pocket-change analyses needing branch lengths, build
  per-HOG FAMSA + FastTree trees (both in the env) and use HOG labels as constraints/reconciliation.
- Guide tree = NCBI taxonomy with polytomies (e.g. Eukaryota node with 5 children, Enterobacteriaceae with 4).
- Scaling: fine to ~1–2 M proteins / few hundred proteomes on one node; tier 3 (all AFDB, 214 M) is out of reach as is.
- Gene order layer, structure rescue: no benefit shown; skip.
- Unverified: QfO scoring was not rerun here; accuracy on our 46-species set is unknown (no benchmark); the
  UniProt-xref coverage numbers below come from 100-protein samples.

## Implications for the evolution project
- **Tier 1 (46 proteomes) is solved cheaply:** full accurate run = 20 min wall / 6.6 CPU-h / 11 GB on 48 threads of a
  loaded node, output validated. Re-running on the final AFDB v6 accession set is trivial; do it whenever inputs change.
- Fix the ID gap first (7% of atlas accessions, mostly zebrafish/maize/rat, not in gene-centric proteomes): map them
  to canonical genes, or they drop out of orthology-based analyses (biasing gains/losses in those species).
- **Densify the tree:** 46 species are sparse (4 vertebrates, 3 of them mammals; one archaeon). Budget allows adding
  complete UniProt reference proteomes as "unstructured context" species (no pockets needed for HOG inference;
  AFDB tier-3 v4 models exist for most): 250 proteomes/1.08 M proteins took 1.9 h accurate on 12 threads in the
  authors' runs → **tier 1 + ~150–200 extra reference proteomes fits in ~1–2 h on 48–64 threads**; this makes gains/
  losses "persistent across several species" testable and supports the tier-3 sampling design.
- **Tier 2 (Swiss-Prot, ~437k proteins, ~11k taxa):** do not run Ancestra on it. Place each Swiss-Prot protein into
  the tier-1(+context) HOG hierarchy by homology (DIAMOND vs HOG members → root HOG, then deepest compatible level by
  its taxonomy; OMAmer-like; minutes on 64 CPUs) and treat it as presence-only evidence (absence uninformative).
- **External orthology as cross-check/fallback** (via UniProt cross-references, no big downloads; sample of 100):
  PANTHER family/subfamily 80% (tier 1) / 92% (tier 2), eggNOG 53% / 56%, OrthoDB 77% / 39%, OMA 48% / 8%.
  PANTHER (families + gene trees with dup/spec nodes) is the best independent check of Ancestra orthogroups and of
  "misassigned orthology" artifacts; eggNOG OGs at taxonomic levels for tier 3. Fetch via UniProt REST
  `uniprotkb/accessions?...&fields=xref_panther,xref_eggnog,xref_orthodb,xref_oma` (or idmapping) — not done in bulk.
- **Tier 3 (all AFDB):** Ancestra does not scale (≥10k CPU-h, TB-scale RAM extrapolated). Use the tier-1/context HOGs
  as a backbone and place tier-3 proteins by homology (as for tier 2) or rely on eggNOG/PANTHER assignments.
- Ancestral states: use HOG nesting (S/D nodes + levels) + pyham for gene-level gain/loss/duplication; pocket presence
  is then a character on the HOG tree (Dollo/parsimony), with deep-level splits flagged as possible artifacts.
