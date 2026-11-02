# Review round 1: data, algorithms, compute (DESIGN v0)

Reviewer: data/algorithms/compute agent, 2026-10-30 ~10:30-11:15 CET. I measured feasibility on this node; I did not
edit DESIGN.md. "[m]" = measured here (commands were short probes, not committed scripts); "[e]" = estimate.

## Measured facts that change the plan
- **Flat kNN is ~5-9x slower than §9 claims.** FAISS 1.14.3 CPU, IndexFlatIP over tier-1 (1.79 M x 384): 20k queries,
  k=51, 32 threads = 36 s [m] -> tiers 1-2 all-vs-all (3.2 M^2) ~ 2.9 h on 32 threads, ~1.5 h on 64 [e], not 20 min.
  faiss in `.venv` has **0 GPUs** (CPU build) and torch is absent. A 4090 does the same GEMM in minutes [e].
  Loading the 1.79 M tier-1 fingerprints from parquet takes 18 s [m].
- **PAE and original v4 models are not fetchable as tested.** AFDB `/files/...predicted_aligned_error_v6.json` and the
  GCS bucket `public-datasets-deepmind-alphafold-v4` return HTTP 403 to curl from this node [m]. The AFDB *API*
  (`/api/prediction/<acc>`) answers and lists `paeDocUrl` and an AF2 **`msaUrl` (a3m)** [m]. The local v6 tars hold
  only `model_v6.cif.gz/pdb.gz` (no PAE) [m]. The local Foldcomp DB (1.1 TB) holds **all** v4 entries, including human
  (`AF-P00918-F1-model_v4` is in the lookup) [m], but no PAE.
- **Context species exist, but not where the experimental evidence is.** Of 3,703 near-complete tier-3 taxa, 2,141
  are named non-metagenome species (753 Archaea, 1,219 Bacteria, 169 Eukaryota) [m]. The eukaryotes are protists and
  early fungi (Perkinsus, Naegleria, Acanthamoeba, Giardia, chytrids, Lingula); there are no near-complete vertebrates,
  insects or land plants. UniProt reference proteomes overlap atlas v4 entries 95 % (Naegleria) and 91 % (Lingula) [m].
  Model length = current UniProt length for 14,888/14,888 Naegleria proteins [m].
- **Annotation sources are reachable and cheap.** UniProt REST stream (accession, PANTHER, Pfam, ft_binding, ft_act_site,
  ft_mutagen) runs at ~1.5k entries/s (15.6k in 10 s) [m]: ~12 min for tiers 1-2. PANTHER xref coverage in Naegleria is
  only 58 %, Pfam 62 % [m]. TED (Zenodo 13908086, CC BY 4.0): `ted_365m.domain_summary.cath.globularity.taxid.tsv.gz`
  is 19.9 GB and covers every AFDB v4 entry with CATH-labelled domain boundaries [m]. InterPro `protein2ipr.dat.gz` is
  13.3 GB [m]. BioLiP2 `BioLiP.txt.gz` is 79 MB [m]. SIFTS `uniprot_segments_observed.tsv.gz` is 11 MB [m]. The M-CSA API
  answers [m]. Pfam-A.regions, M-CSA and UniProt 2026_03 are already in `/tmp/pe_big/raw` [m].
- **Tools are present.** foldseek 10 (reads tar and Foldcomp input directly, `--gpu`), mmseqs 18, FAMSA 2.4.1,
  DIAMOND 2.2.8 and FastTree are in `envs/ancestra/bin` [m]. FoldMason and US-align are not installed. Foldseek
  `createdb` on a 1,773-model tar took 0.5 s, and all-vs-all took 5 s [m].
- **O1 orthogroups are mostly small.** 36,736 root groups: median 4 genes, p99 106, max 4,223. 8 % of genes are in
  groups larger than 500 [m].
- GPU 3 is at 100 % (24 GB used, the AF3 benchmark), and GPUs 0-2 are idle [m]. `/tmp` has 1.7 TB free [m].

## Numbered changes (within scope)
1. **Move the kNN to the GPU, and make it homology-filtered.** Install `faiss-gpu-cu12` (or cupy) into a /tmp venv
   (~0.5 GB, outside the home quota). Search k=1000 (GPU limit 2048), then drop H0-H2 pairs with a join on root HOG,
   and keep the top-50 remaining pairs plus the top-50 raw pairs. Why: homologs fill the raw top-50 (75 % of kinase
   top-10 hits are kinases), so F_conv computed from raw top-k is biased low and depends on family size. Cost
   [e]: ~10-20 min on one 4090. The output (3.2 B ids+scores, ~26 GB) stays on /tmp and only the filtered pairs are
   persisted. CPU fallback: ~1.5 h on 64 threads, run in the background.
2. **Do homology labelling by joins on precomputed tables; run alignment only on the residual pairs.** H1/H2: HOG join
   (free). H3: TED CATH superfamily at the pocket region (residues ∩ domain ≥ 50 %), plus Pfam clan (Pfam-A.regions,
   already downloaded), plus PANTHER (UniProt xref). Streaming TED filtered to tier-1/2 accessions takes ~15 min download
   + ~20 min DuckDB [e]. TED covers tier 3 unchanged, so it scales. MMseqs2 and Foldseek run **pair-restricted**
   (`tsv2db` of the unlabelled significant pairs -> `mmseqs align` / `foldseek structurealign --alignment-type 1`), not
   all-vs-all. For ≤ 10 M pairs this is < 1 h on 64 cores [e; measure on 100k pairs first]. The Foldseek DB for tiers 1-2
   is built straight from the v6 tars (~10-20 min [e]). For tier 3, use `--input-format 5` on the Foldcomp DB.
   Do **not** require PANTHER/Pfam to be present for an "N" label, because coverage is 58-62 % in non-model
   proteomes. Require instead: TED labels present on both sides, plus a Foldseek negative on the full chains (an E-value
   or TM-score cutoff). Otherwise label the pair "unknown".
3. **Make the structural source of each case explicit, and drop "original AFDB files" for tier 3.** The original v4
   files are not downloadable from this node, but the local Foldcomp DB has every v4 entry, including human and
   Swiss-Prot. Re-describe **all** members of a case from Foldcomp v4. That is one source, local, and the same method
   record as tier 3. Use v6 (local tars) only as a second-source stability check for tiers 1-2. Needs `pip install
   foldcomp` in a /tmp venv (minutes).
4. **Treat PAE as at risk; test it now (5 min).** Use python `requests` against the API's `paeDocUrl` for 10 accessions.
   If PAE is blocked, §5.3 uses pLDDT + disorder + the co-fold PAE. Record that PAE was unavailable; do not silently
   skip it. Similarly, test AFDB `msaUrl` downloads: if they work, AF2 a3m files are a no-rate-limit MSA source for
   non-human orthologs in co-folding (record the source; the MSA differs from AF3 jackhmmer).
5. **Add on-demand context pockets for the clades where cases will live.** Near-complete context exists only for
   prokaryotes and protists. For vertebrate, insect and plant persistence (filter 1), run the atlas worker
   (`PC/afdb/worker.py`, same method record) locally on Foldcomp entries of ~20-40 chosen reference proteomes. These are
   the tier-3 units the atlas would produce later; dedupe by accession when the atlas catches up. Cost: ~100
   structures/s on 64 cores (from 14/s on 8 workers) -> 400k structures ≈ 70 min [e]. Choose the proteomes
   deterministically from the manifest + a clade list, as for O2. O2 itself stays as designed: 557k + ~1 M proteins,
   ~1.3 h, ~30-40 GB RSS [e from the exponent 1.35].
6. **Build MSAs from model sequences, assert identity, and keep them off home.** Take sequences from the structure files
   (v6 `_entity_poly`, or Foldcomp) or check UniProt FASTA against the model by length + md5. Mismatches go to U and are
   counted. Run FAMSA per root HOG for O2 (~1.5 M sequences): ~30-60 min on 64 cores [e]. MSAs go on /tmp. Persist only
   a parquet of `site -> MSA columns` (~100 MB).
7. **Residue mapping: keep MSA as primary; replace the "Foldseek/US-align pairwise to a reference" option with
   FoldMason.** A star alignment to one reference is not transitive, so column Jaccard between two non-reference
   members becomes noisy. FoldMason (static binary, minutes per thousand structures) gives a true structural MSA from
   the same files. Decide on §8.1 by identity bin (< 30 %, 30-50 %, > 50 %) and on a dedicated multi-domain stratum
   (domain shuffling within root HOGs). Use the benchmark already shipped in `B/benchmarks/af_set/sites.tsv.gz`
   (6.7k LIGYSIS sites on AFDB v6 models, UniProt numbering) plus UniProt ECO:0000269 sites. Cost: FoldMason for all ~29k
   groups ≥ 3 members ≈ 2-6 h [e]. Run it only for low-identity groups and for case orthogroups unless the benchmark shows
   a clear gain.
8. **Compute SOG pairs sparsely.** All-pairs site Jaccard is quadratic per orthogroup: a 4,223-gene group has ~12k sites
   and ~7e7 pairs, and O2 groups will be larger. Use a column -> sites inverted index and keep only pairs that share
   ≥ 1 footprint column. Then Leiden. Total time ~minutes [e]. For tier 3, compute sites per atlas unit (per-file
   streaming, embarrassingly parallel), and compute SOG per HOG.
9. **Define tier-3 scaling now (Principle 0 needs it).** (a) The Part A query set is fixed: tiers 1-2 plus a seeded,
   clade-stratified sample of tier-3 pockets. The database grows to 0.8 B through IVF-PQ: ~64-75 GB RAM, ~0.6 TB of fp16
   vectors on /tmp for exact re-rank. Exact all-vs-all over 0.8 B (~2.5e20 FLOP) is infeasible. (b) Place tier-3
   proteins into backbone HOGs with DIAMOND `--fast`, chunked per atlas unit and idempotent. Measure on one unit (4k
   entries) before quoting a cost, which I estimate at node-days for 214 M. (c) TED/Foldseek labels already cover v4.
10. **Persist what /tmp cannot regenerate cheaply.** Copy orthoxml, sites, SOG/state tables and filtered pair tables
    (compressed parquet) to `results/`. The home quota is at ~47/50 GB, so budget < 3 GB of new results. Everything
    else must be regenerable by one command from the manifest.
11. **Experimental data: fetch subsets, record licences.** BioLiP2 (79 MB; no formal licence, so cite it and do not
    redistribute raw data): map PDB numbering to UniProt with SIFTS segments (11 MB) and the local PDB mirror. UniProt
    features: CC BY 4.0, keep ECO:0000269/7744 only. M-CSA: CC BY 4.0. ChEMBL (CC BY-SA 3.0): query the REST API for
    case proteins only, or use the SQLite dump (~5 GB compressed) on /tmp if more than ~200 targets. BindingDB:
    case-level only. PDB mirror: CC0. Add all of these to DATA.md with versions and checksums.

## Failure points (ranked)
(1) PAE and original-model access (change 4). (2) No eukaryotic context pockets near the model organisms (change 5).
(3) kNN runtime and top-k bias (change 1). (4) MSA misalignment in multi-domain or low-identity HOGs, which creates
false gains and losses (change 7, plus §8.1 strata). (5) Sequence/numbering drift between the model and the FASTA
(change 6). (6) /tmp loss on reallocation (change 10). (7) CPU contention: O2, MSAs, the kNN fallback and the worker
fill-in all compete for 64 cores, so schedule them as below.

## Ordering that de-risks the critical path
T+0 (in parallel, ~1 h): O2 Ancestra (already downloading proteomes); UniProt annotation streams; TED/BioLiP2/SIFTS
download; faiss-gpu install; PAE/msaUrl probe; foldcomp install. T+1 h: run sites + MSA + SOG on **O1 only**
together with the §8.1 benchmark. This fixes τ and the MSA-vs-FoldMason choice before O2 lands, and it can be
re-run unchanged on O2. T+1-2.5 h: GPU kNN + homology joins on GPUs 0-2, while the CPUs run the context-pocket fill-in
after O2 finishes. T+3 h: a vertical slice on 3-5 known families (e.g. carbonic anhydrase, HMGCR, kinases) through
states, reconstruction, the artifact filters and the web app view. This is the first gate before scaling to all SOGs.
Co-folding needs MSAs for non-human/yeast orthologs: request them as soon as case candidates exist.
