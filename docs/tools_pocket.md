# Pocket tools brief: BURROW (pocket finder) + pocketcompare (fingerprints, AFDB atlas)

Written 2026-10-30 for the pocket-evolution project. Paths: `B` = `external/burrow`, `PC` = `external/pocketcompare`,
`PCA` = `external/pocketcompare-agents` (research record), `ATLAS` = `/cluster/scratch/gmagnusson/pocketcompare_afdb`
(read-only). "measured" = measured by me in this container (commands in `docs/ENV.md`); everything else is a claim
from the cited file, produced by AI agents and **not human-verified** (`B/README.md`, `PC/README.md` preambles).

## 1. What the tools are
- **BURROW 1.0.0** (`burrow-pockets`, import `pocketdetect`, MIT): grid (1 Å) solvent-shell points → ray-cast
  buriedness + chemistry features → LightGBM point score → multi-threshold candidate pockets → LightGBM pocket ranker
  → NMS (10 Å, Jaccard 0.75) and output floor. No neural nets; reference, Numba-CPU and GPU backends byte-identical
  (`B/README.md`, `B/docs/FOR_FURTHER_DEVELOPMENT.md` §2). Only protein heavy atoms are used; ligands, waters,
  ions, cofactors are ignored (`B/README.md` Quickstart).
- **pocketcompare 1.0.0** (MIT): every BURROW pocket → one 384-d unit float32 vector; similarity = dot product
  (`PC/README.md`, `PC/docs/METHOD.md`). Frozen method `pocketcompare-final:A4b:n1_linear|raw`, config
  `PC/configs/final_method.yaml` (sha256 b3ea8fe9…), all model files sha256-verified on load (`PC/docs/METHOD.md`).
- **Atlas**: frozen pocketcompare run over AFDB in 3 tiers (`ATLAS/POCKETOME.md` §2, `PC/afdb/README.md`).

## 2. CLI and Python APIs we need
**pocketcompare CLI** (`PC/src/pocketcompare/cli.py:470-525`):
- `pocketcompare describe IN... -o out.parquet [--af-mode] [--backend ref|fast] [--device cpu|cuda:N]
  [--workers N] [--unit merged|member] [--min-points 10] [--chains A,B] [--n1]` — every pocket of each structure.
  Use `--af-mode` for every AF model (`PC/README.md`). Also writes `out.structures.parquet` (status, timings).
- From residues: `--residues A:43,A:45` (4 Å residue shell ∩ bur10 ≥ 0.5); from ligand: `--ligand lig.sdf|.mol|.pdb`
  (4 Å ligand shell). Add `--locate burrow` to describe the BURROW pocket that best covers the site (D mode =
  how atlas pockets were made; compare like with like) (`PC/docs/METHOD.md` §2). Output id `<stem>:site`.
- `pocketcompare compare A B [--top k] [--matrix S.npy] [-o pairs.tsv]` (tables or structure files).
- `pocketcompare index FP.parquet... -o DIR [--kind flat|ivfpq] [--nlist] [--pq-m 64]`;
  `pocketcompare search Q DIR -k 10 [--rerank 100] [--nprobe] [--exclude-self] -o hits.tsv`.
- Atlas-scale index: `python PC/afdb/build_index.py --out-root ATLAS --tiers 1,2[,3] --kind flat|ivfpq -o DIR
  [--min-score] [--save-vectors] [--compact-ids] [--threads]` (`PC/afdb/build_index.py:61-73`); output dir is
  readable by `pocketcompare search`. IVF-PQ recall@10 ≈ 0.89 alone, ≈ 0.99 with exact re-rank of top-100 (needs
  `--save-vectors`) (`ATLAS/POCKETOME.md` §4).

**pocketcompare Python** (`PC/src/pocketcompare/method.py`):
- `method.describe_all(path|Structure|Featurized, af_mode=True, chains=None, unit="merged", min_points=10,
  return_n1=False)` → `(V (n,384) float32, info DataFrame)` (`method.py:220`).
- `method.describe(path, site=ligand_file | (n,3) xyz | ["A:43",...], af_mode=..., locate="shell"|"burrow"|
  "burrow_flat", return_info=True)` → `(1,384)` (`method.py:200`).
- `method.similarity(a, b)` = inner product; `similarity.csls_r(E_q, E_bg, k=10)` for optional CSLS hub correction
  (`PC/src/pocketcompare/similarity.py:50`).
- `fpindex.read_table(path)` → `(meta DataFrame, float32 V, method record)`; `read_tables(paths)` refuses mixed
  methods (`PC/src/pocketcompare/fpindex.py:83-109`). `fpindex.build_index(V, meta, out_dir, kind=...)`.
- Low level: `pockets.featurize_structure(path, chains, af_mode)` (runs BURROW once, keeps grid points, features,
  pockets) and `pockets.merge_pockets`, `shape_site(fs, idx, "sph8")` (`PC/src/pocketcompare/pockets.py:206,381,550`).
- Atlas helpers (`PC/afdb/afdb_common.py`): `open_pockets(OUT_ROOT[, "structures"])` (unified schema dataset),
  `fill_source(table)` (derives the missing `source` column of tier-1/2 files).

**BURROW CLI** (`burrow predict -h`; `B/README.md`): `burrow predict IN... --out DIR [--af-mode]
[--preset default|accurate|fast] [--ranker all|small|geometry] [--min-score X|none] [--output-rule v2|region]
[--model JSON | --model-dir D] [--af-ranker JSON] [--backend auto|fast|reference] [--device] [--workers]
[--no-points]`. Output per input: `DIR/<stem>.csv` (`rank,score,center_x,y,z,residues`), `DIR/<stem>.points.csv`
(`rank,x,y,z`), `DIR/timing.csv`. **BURROW Python**: `pocketdetect.burrow.export.pocket_representations(path,
ranker, preset, min_score=None, output_rule, af_mode)` → dict with `pockets` (76 descriptors), `points` (84 point
features), `pocket_points` memberships (`B/docs/FOR_FURTHER_DEVELOPMENT.md` §3f).

## 3. Formats
- **describe output** (parquet): `pocket_id` = `<stem>:m<merged_index>` (0-based merged index, **not** BURROW rank;
  non-contiguous because merged pockets < 10 points and, in the atlas, filtered pockets are absent;
  `PC/src/pocketcompare/cli.py:150`), `merged_index, burrow_index, burrow_rank, members, burrow_score, n_points,
  volume_A3, center_xyz` (sph8 centre), `site_n_points, site_volume_A3, sph8_fallback, residues`
  (`chain:resnum` author numbering, space-separated), `ok`, `plddt_{mean,min,frac_lt70}` (AF mode only),
  `fingerprint` (list<float>[384]). Schema metadata key `pocketcompare` = method record (`PC/README.md`).
  Site modes add `mode` (L/D), `definition`, `coverage`, `merged_coverage`, `merged_precision`, `d_match` (measured).
- **Atlas pocket table** `ATLAS/pockets/tier<T>/<unit:06d>.parquet`: above + `unit_id, tier, accession, fragment,
  model_version, taxid, organism, source`; `fingerprint` = `fixed_size_list<halffloat>[384]` (fp16; max error
  6e-5) (`ATLAS/POCKETOME.md` §3). Structure table `ATLAS/structures/tier<T>/…` has denominators (`n_res`, pLDDT
  stats, `n_pockets_all`, `n_pockets_kept`, `status`). Measured: tier-1 pocket files lack `source`; tier-1 structure
  files have `source` but lack `source_file` (POCKETOME/README list it); tier-3 files have both.
- **Atlas filters**: ok, `burrow_score ≥ 0.03`, `plddt_mean ≥ 70` (lining residues), after AF-mode pruning of
  pLDDT < 50 runs (`PC/afdb/README.md` "Pocket filters").
- **Coverage now** (measured, DuckDB over all files, 2026-10-30 ~10:00): tier 1 = 555,846 structures / 1,788,750
  pockets / 46 taxa (complete); tier 2 = 436,886 / 1,442,678 / 10,839 taxa (complete); tier 3 = 24.2 M structures /
  72.0 M pockets / 16,227 taxa from 7,621 of 80,304 units (~9.5 %, still running). Tier-3 `sum(n_pockets_kept)` was
  72,016,177 vs 71,976,247 pocket rows at query time (files being written concurrently; re-check when frozen).

## 4. The fingerprint and what similarity means
- Pipeline (`PC/docs/METHOD.md` §1-8): BURROW (default preset, `all` ranker, pinned model `burrow_all_default_s0`,
  full list, region-rule subset) → merge touching sub-pockets (points ≤ 1.5 Å or residue Jaccard ≥ 0.3), drop < 10
  points → **sph8 site** = all buried grid points (bur10 ≥ 0.5) within 8 Å of the bur16⁴-weighted centroid →
  N1 = A1 typed volumetric fingerprint (514 d: pharmacophore-typed cavity histograms, 28 typed pair-distance
  distributions, USRCAT-like moments, global shape, log volume) + env block (113 d, lining chemistry fractions) →
  z-score (fixed stats), clip ±5 → 3 linear heads 627→128, L2 each, concat/√3.
- Score = cosine = mean of three heads' cosines. Trained contrastively on BioLiP ligand-sharing pairs (0.5), same
  site ligand-shell↔detected (0.25), AF-dimer interface pairs (0.25) (`PC/docs/METHOD.md` §9). So "similar" =
  "likely to bind similar ligands", fold- and sequence-free. **Not** an alignment; gives no residue correspondence
  (`ATLAS/POCKETOME.md` §1).
- Score scale (human AFDB): random different-protein pairs median 0.11, p99 0.62, p99.9 0.77; nearest neighbour
  median 0.82; 1st/10th/50th NN 0.82/0.76/0.72. No validated "same ligand" cutoff; use ranks, per-query z-scores or
  matched backgrounds (`PC/README.md` "What a fingerprint is"; `ATLAS/POCKETOME.md` §5).
- Invariant to rigid motion, deterministic; ref/fast/GPU identical (`PC/docs/METHOD.md` §10). Measured: ref vs fast
  max|Δ| = 0.0 on 233 pockets.

## 5. Validated accuracy (held-out, pre-registered; `PC/docs/RESULTS.md`, `PCA/REPORT.md` §4)
| test | result | source |
|---|---|---|
| TOUGH-M1-test, detected pockets (21,410 pairs) | ROC-AUC 0.699 vs B1 (BURROW descriptors) 0.546, Δ +0.154 [+0.050,+0.289] | RESULTS.md C1-D |
| TOUGH-M1-test, ligand shell | 0.698 vs 0.568 (B1) / 0.568 (ESM-2); borderline (Holm p 0.0495) | C1a/C1b |
| size-stratified TOUGH-M1 | 0.719 (L) / 0.717 (D) | RESULTS.md §4.3 |
| temporal PDB-T (987 sites after 2025-11-12) | 0.752 vs B1 0.671, Δ +0.081 [+0.028,+0.132] | RESULTS.md R1 |
| ProSPECCTs P5.2+P6.2+P7 | 0.740 vs B1 0.718, not significant; P6 (cofactor) 0.622 | RESULTS.md C2 |
| BLiPS-test AF vs holo | AF D 0.723 vs holo 0.728; paired Δ within ±0.03 | RESULTS.md §4.6 |
| unseen ligand chemotypes (BLiPS) | 0.685 vs B1 0.770 (weak spot) | RESULTS.md §4.6 |
| Vertex (homologous positives) | 0.822; Foldseek TM 0.910, ESM-2 0.837 | RESULTS.md §4.8 |
| blind retrieval, 13,876 AF pockets | AUC 0.719, **hit@100 0.078**, BURROW site recall 0.700 | RESULTS.md §4.10 |
| PPI interface orthologue retrieval | MRR 0.683 vs ESM-2 0.954 | RESULTS.md §4.9 |
Caveats: leakage control removed 76 % of TOUGH-M1 → not comparable to literature; ligand heavy-atom count alone
scores 0.95 on TOUGH-M1 L-mode (size leak); TOUGH-M1-test has ~102-112 units, MDE ≈ 0.12 (`PCA/REPORT.md` §2, §10).
**BURROW** (`B/README.md` headline table): LIGYSIS top-N+2 DCC12 72.5 (vs fpocket+PRANK 60.4, P2Rank 52.6), DCA4
54.2; AF set (2,526 AFDB v6 models) 64.5 vs P2Rank-AF (+13.1, mostly ions; non-ion +1.9 n.s.); CryptoBench
(apo/cryptic) and dimer interfaces inconclusive. AF mode removes hallucinated pockets (5.0 % → 0.2 %) but recall is
equivalent to a plain pLDDT < 70 filter (`B/docs/AF_MODE.md`). Much of BURROW's margin is on ions (38 % of sites).

## 6. Runtime (measured, 1 core, warm Numba, thread caps = 1; 53 human AFDB v6 F1 models, 4.4 pockets/model)
| step | median / mean per structure |
|---|---|
| `pocketcompare describe --af-mode --backend ref` | 0.256 / 0.396 s (BURROW 0.222 s; describe 0.023 s) |
| same, `--backend fast` | 0.142 / 0.208 s; 8 workers: 0.072 s/structure wall (≈ 14 structures/s) |
| CLI start-up | ≈ 2-3 s per call; first call ever ≈ 15-17 s (Numba compile) |
| `burrow predict` HMGCR (888 aa), timing.csv | default AF 632 ms; fast 405; accurate 1,013; no-AF 821; reference backend 1,702; `--af-ranker` 1,657 |
| flat FAISS index of tier 1 (1.79 M pockets) | 28.6 s build, 2.7 GB, 4.0 GB RSS; 13 queries k=10 in 1.3 s incl. load |
Published: 0.42/0.21/0.13 s (ref/fast/GPU) per structure; 380 structures/s per 4-GPU node (`PC/README.md` Runtime).

## 7. Settings that matter for pocket-finder stability checks
- What can be varied **without** changing the fingerprint definition: nothing in BURROW — pocketcompare always runs
  its pinned BURROW model/preset (`PC/src/pocketcompare/pockets.py:206-230`; it ignores `burrow` CLI options). Its
  own knobs: `af_mode`, `unit` (`merged`/`member`), `min_points`, `chains`, site + `locate`. Fingerprints made with
  other settings are a different representation: compare only within one setting.
- BURROW knobs for detection stability (`burrow predict`): `--af-mode` on/off; `--preset fast|accurate` (grid 1.5 Å /
  sub-lattice changes); `--model models/burrow_all_default_s{1,2}.json` (independent training seeds, `B/models/`);
  `--ranker small|geometry`; `--min-score none` / `--output-rule region` (longer lists, same top pockets);
  `--af-ranker` (re-ranks; demotes pockets near disorder, −3 pts recall there, `B/docs/AF_MODE.md`).
- Measured on HMGCR (P04035; n = 1 protein, illustrative): best residue-Jaccard of the default top-5 pockets in
  other runs: repeat/reference backend 1.0 (deterministic); `--min-score none`/`region` 1.0; seed s1 0.30-1.0;
  seed s2 0.30-1.0; no AF mode 0.47-0.96; fast preset 0.31-0.83; accurate 0.37-0.82; geometry ranker 0.11-0.36.
  Rank orders also change (e.g. default #4 is #1 without AF mode).
- Measured on 53 models: of 177 atlas-filtered AF-mode pockets only 55 % have a non-AF-mode counterpart with residue
  Jaccard ≥ 0.5; among those, fingerprint cosine median 0.886 (p10 0.69). Same-site AF vs non-AF fingerprints
  therefore look like "different but similar" pockets.
- Residue-shell vs detected (`--locate burrow`) vectors of the same site differ: CA2 test 0.644 vs 1.000 (measured);
  published median 0.87, p10 0.54-0.65 (`PC/README.md` Limitations).
- Recommended matching across settings/species: residue overlap (Jaccard ≥ 0.3-0.5) after mapping numbering, not
  centre distance or pocket_id (`B/docs/FOR_FURTHER_DEVELOPMENT.md` §5; `ATLAS/POCKETOME.md` §6.6).

## 8. Residue correspondence / alignment
None of the three repos implements pocket alignment or residue correspondence (searched `align|correspond|superpos|
kabsch` over `PC/src`, `PC/afdb`, `PC/scripts`, `B/src`, `PCA/src`, `PCA/scripts`: only Fr-TM-align scores shipped
with TOUGH-M1 are read, `PCA/scripts/final_eval.py:330`). `PCA/REPORT.md` §11 lists "re-ranking with an alignment
step" and A3 pseudocentre alignment as untested follow-ups. We must build correspondence ourselves (e.g. residue
mapping via the orthogroup MSA / Foldseek-TMalign of the two models, then Jaccard of mapped lining residues).
Atlas `residues` are UniProt numbering for fragment 1 (`ATLAS/POCKETOME.md` §3), which makes MSA mapping direct.

## 9. Known weaknesses (beyond §5)
- Fixed 8 Å sph8 site: large/elongated pockets described only at their buried centre; sites < 8 Å apart share most
  of their description; ligand size deliberately not encoded (`PC/README.md` Limitations).
- Detection inherited from BURROW: ~70 % of known sites detected in AF models; one site can be split into 2-3
  pockets; one pocket can cover two ligands (`ATLAS/POCKETOME.md` §6.5-6.6; `B/docs/FOR_FURTHER_DEVELOPMENT.md` §5).
- Cofactors/metals stripped: heme/NAD sites described as empty cavities; isolated metal sites (< 10 points) dropped;
  ProSPECCTs cofactor subset P6 = 0.622 (`PC/README.md`; `PC/docs/RESULTS.md`).
- Hubness: some pockets near many others (skew 1.62 → 0.55 with CSLS; CSLS hurt retrieval) (`ATLAS/POCKETOME.md` §5).
- Homology dominates neighbour lists (75 % of kinase non-self top-10 were kinases) (`ATLAS/POCKETOME.md` §6.1).
  Measured: HMGCR m1 top hits are HMGCR orthologues (Dictyostelium P34136 0.968).
- Monomers only: no interface pockets in the atlas; PPI-modulator pockets not recognised across families.
- Shallow/flexible sites have low hit rates for every method; lipid/additive sites frequent (`B/docs/FOR_FURTHER_DEVELOPMENT.md` §5).
- Tier 3 = Foldcomp-decompressed v4: only 60 % of pockets have a Jaccard ≥ 0.5 counterpart vs the original model;
  matched cosine median 0.95 (p5 0.76) (`ATLAS/POCKETOME.md` §6.10).

## 10. Benchmark datasets shipped or referenced
Shipped in `B/benchmarks/` (25 MB, manifests + frozen site tables, no structures; `B/benchmarks/README.md`):
LIGYSIS sites (CC-BY-4.0), CryptoBench (MIT), PocketWatch (MIT), **af_set** (2,526 AFDB v6 models with LIGYSIS
sites transferred onto the model; `af_set/sites.tsv.gz`, residues in UniProt numbering; built by BURROW, MIT),
af_set_allsites, dimer (400 assemblies), COACH420/HOLO4K manifests (P2Rank datasets, no licence file). BURROW
predictions for all test sets: `B/predictions/*.csv.gz` (MIT). Licences: `B/DATA.md`, `B/benchmarks/THIRD_PARTY_LICENSES.md`.
Shipped in `PCA/data/splits/` (12 MB): BLiPS sites/pairs (from BioLiP2, "freely available", no formal licence),
PDB-T temporal sites (`pdbt_sites.parquet`), PPIPS/PPI-L tables, homology graph nodes/components; raw data
(`PCA/data/raw/`) is **not** shipped (`PCA/DATA.md`). Referenced only (fetch ourselves; `PCA/DATA.md`):
- TOUGH-M1: Dataverse doi:10.7910/DVN/L7H7JJ (CC0) + OSF https://osf.io/6ngbs/ (lists/baseline scores; licence unstated).
- ProSPECCTs: http://www.ewit.ccb.tu-dortmund.de/ag-koch/prospeccts/ (POST form) — internal use only, **no redistribution**.
- Vertex: figshare doi:10.1021/acs.jcim.6b00118.s002 (CC BY-NC 4.0). TOUGH-C1: OSF https://osf.io/enz69/ (unstated).
- KLIFS API (unstated), LIGYSIS Zenodo 10.5281/zenodo.14645504 (CC-BY-4.0), Pasteur PPI pocketome Zenodo
  10.5281/zenodo.10805580 (CC-BY-4.0), BioLiP2 https://zhanggroup.org/BioLiP/ (no formal licence).
No benchmark raw data is mounted here (only `ATLAS` and `/cluster/scratch/gmagnusson/dimers` exist under
gmagnusson's scratch; the research repo paths `/cluster/project/beltrao/gmagnus/...` are not mounted). For our
project the most directly usable sets are `B/benchmarks/af_set` (ligand sites on AFDB v6 models = positives for
"pocket present") and PDB-T/BLiPS (ligand-sharing pairs).

## 11. Unverified / flagged
- GPU backend not tested (torch not installed: home quota, `docs/ENV.md`). Tier-3 rows not re-derived (no foldcomp).
- Published accuracy numbers were not re-computed here (no benchmark data mounted).
- Stability numbers in §7 are from 1 protein (BURROW settings) and 53 proteins (AF vs non-AF): indicative only.

## Implications for the evolution project
1. Fingerprint similarity ≠ residue correspondence. To say "the same pocket in orthologues" we need our own mapping:
   lining residues (UniProt numbering) mapped through the orthogroup MSA, then Jaccard; use cosine only as a second,
   shape/chemistry axis (e.g. for "function shift" = conserved position, lower cosine).
2. Homology dominates cosine neighbours; random cross-protein pairs reach p99 0.62 / p99.9 0.77. Convergence claims
   need homology-removed pairs (Foldseek/MMseqs2), size-matched backgrounds and a size-only baseline (POCKETOME §6.8).
3. "Lost pocket" is weak evidence: BURROW recall ≈ 70 %, atlas drops score < 0.03 and pLDDT < 70 pockets, AF mode
   prunes pLDDT < 50 runs. Before calling a loss, re-run that model with `--min-score none`, without AF mode, with
   seeds s1/s2 and the fast/accurate presets, and inspect pLDDT; check `n_pockets_all` vs `n_pockets_kept`.
4. BURROW output is unstable across settings at the level of individual sub-pockets (HMGCR top-5 Jaccard 0.3-1.0
   across seeds/presets); define stability at site level (merge pockets with Jaccard ≥ 0.3, require a match in ≥ k
   of N settings), not by pocket_id or rank.
5. pocket_id `:m<k>` is a 0-based merged index with gaps; never join across models or settings on it.
6. Mixed model sources: tiers 1-2 = AFDB v6 files, tier 3 = Foldcomp-decompressed v4 (only 60 % pockets reproduce
   at Jaccard ≥ 0.5 vs the original). Human/model organisms (v6) vs most other species (v4) is a confound; for case
   studies re-describe all members from one source (e.g. original v4 or v6 files) with one setting.
7. Tier 3 is ~9.5 % done and ordered by tree round-robin: absence of an orthologue's pocket may just mean its unit
   has not run — always join with `ATLAS/state/done` / structure tables for denominators.
8. Each accession appears in one tier only; human (tier 1) vs other primates (tier 3) also differs by model version.
9. Multi-fragment human proteins (> 2,700 aa) duplicate sites with fragment-local numbering: deduplicate after
   offsetting (`_struct_ref_seq`).
10. Cofactor/metal sites are described as empty cavities and small metal sites vanish; "pocket chemistry shift"
    claims for heme/NAD/metal sites need care. Ions also inflate BURROW pocket counts.
11. The sph8 description saturates for big pockets: gains/losses of sub-pockets in large clefts can be invisible to
    the fingerprint; use BURROW points/residues for geometry.
12. Atlas schema drift: tier-1/2 pocket files lack `source`; tier-1 structure files lack `source_file`. Use
    `union_by_name=true` in DuckDB or `afdb_common.open_pockets` + `fill_source`.
13. Home quota (~50 GB, nearly full) cannot hold indices or GPU wheels; put big outputs under /tmp/pe_big/work.
14. Atlas is on Euler scratch (15-day purge, `PC/afdb/README.md`); snapshot the metadata we depend on early.
