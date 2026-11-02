# Co-folding with AlphaFold3 (inference only)

Status: AF3 is installed and runs on the 4x RTX 4090 node with the default
Triton flash-attention kernel. Positive controls reproduce the crystal ligand
poses. See the benchmark below.

## 1. Install

```bash
scripts/cofold/setup_af3_env.sh        # idempotent, about 2 min with a warm uv cache
```

What it does:

* clones `google-deepmind/alphafold3` into `envs/af3_src`, pinned to commit
  `8f8ac383afba995f5953182d39a560ffb8bbbfd3` (2026-10-29). It needs jax 0.10.2
  with jax[cuda12], and the CUDA wheels come from PyPI, so no system CUDA is needed.
* runs `uv venv --python 3.12 /tmp/af3_env` and then `uv sync --frozen`. This
  builds the C++ extension with the system gcc 13.3, and uses cmake and ninja
  pulled in by scikit-build-core from PyPI. No conda is needed.
* runs `build_data`, which builds the CCD and chemical-component pickles into
  `envs/af3_src/src/alphafold3/constants/converters/` (about 25 s).
* decompresses the weights: `zstd -d /cluster/project/beltrao/shared/alphafold3/models/af3.bin.zst -o data/af3_weights/af3.bin` (1.1 GB).
* installs `gemmi` for analysis.

**Why the venv is on /tmp.** The home/project filesystem holding the repo has
a quota of about 50 GB, and it was at 47 GB. The first `uv sync` into
`envs/af3` failed with `Quota exceeded (os error 122)`. The venv needs 5.8 GB,
mostly CUDA wheels. It is therefore built on node-local `/tmp/af3_env`
(1.7 TB NVMe), which vanishes with the job, and you must rerun
`setup_af3_env.sh` on a new allocation. All scripts read the location from
`AF3_ENV` (default `/tmp/af3_env`). The source tree, CCD data and weights stay
in the workspace (about 1.6 GB).

RTX 4090 notes (Ada, sm_89):

* `--flash_attention_implementation=triton` (the default) works. The `xla`
  fallback is only required for pre-Ampere GPUs, and we did not need it.
* `XLA_FLAGS=--xla_gpu_enable_triton_gemm=false` is set. AF3 docs require it
  as the compile-time workaround.
* The custom-kernel-fusion workaround for 7.x GPUs is not needed.
* There is a persistent JAX compile cache at `/tmp/af3_jax_cache`
  (`--jax_compilation_cache_dir`). It is shared by all jobs on the node.

## 2. Build inputs: `scripts/cofold/make_af3_input.py`

```bash
PY=/tmp/af3_env/bin/python
$PY scripts/cofold/make_af3_input.py --acc P24941 --ligand CCD:ATP --ligand CCD:MG \
    --seeds 1,2,3 --name cdk2_atp --out /tmp/af3_inputs/cdk2_atp.json
$PY scripts/cofold/make_af3_input.py --acc P00918 --ligand CCD:ZN \
    --ligand 'SMILES:CC(=O)Nc1nnc(s1)S(N)(=O)=O' --out /tmp/af3_inputs/ca2_azm.json
```

* **Sequence to fold.** The script takes the sequence from `--seq` or
  `--seq-fasta` if given. Otherwise it fetches the current UniProt sequence
  (REST, with retries). With `--trust-cache-seq` it folds the cached sequence
  as-is.
* **Cache lookup.** It looks for `<acc>_data.json.gz` in the human
  (25.05, 20,293 files) and yeast (25.06, 6,039 files) caches. The cached
  `unpairedMsa`, `pairedMsa` and `templates` (mmCIF text plus index maps, at
  most 4) are reused **only if the cached query sequence is identical**.
* **On mismatch the script exits with code 2** and prints the substitutions
  or edits. It never silently folds another sequence.
  * `--reuse-msa-for-variant` is an explicit, recorded approximation for
    same-length variants such as ancestral states. It replaces the query row
    of the cached MSAs, keeps the templates, and appends `_variantmsa` to the
    job name.
* **No cache: exit code 3.** The script writes `<out>.NEEDS_MSA`.
  `--msa-free` writes a single-sequence input instead, which is not recommended
  for real cases.
* **Other options.**
  * `--no-templates` drops the templates, which is useful to test for template
    leakage on positive controls.
  * `--copies N` builds homo-oligomers.
  * Ligands are given as `CCD:<code>` (`+` for multi-component) or
    `SMILES:<smiles>`, which is validated with RDKit.
  * Seeds are fixed, with a default of `[1,2,3]`.
* **Sidecar file.** A `<out>.meta.json` records the provenance: sequence
  source, cache file and MSA mode.
* **Spot check.** In 60 random cached entries (40 human, 20 yeast), the cached
  sequence was identical to current UniProt in 60 of 60.
* **Templates and the PDB mirror.** The cached JSONs already contain templates,
  so the PDB mirror is not needed for cached proteins. A new template search
  would need `pdb_seqres.fasta`, which the mirror at
  `/cluster/scratch/jjaenes/pdb/` lacks; the mirror is sharded `xx/<id>.cif.gz`.
  It would also need a flat `mmcif_files/` directory. Both can be derived from
  the mirror.

### Proteins without cached MSAs: `scripts/cofold/colabfold_msa.py`

**Implemented** (the ColabFold MMseqs2 server route; the three alternatives — a local
AF3 jackhmmer pipeline over ~630 GB of databases, local MMseqs2 over the 1–2 TB
ColabFold DBs, and `--reuse-msa-for-variant` on a cached ortholog — are all rejected
here: the first two blow the 500 GB storage budget, and the third is barred across the
derived/ancestral boundary by `config/prereg.yaml`.) This is the route for non-human/non-yeast orthologs,
ancestral reconstructions and engineered variants, and — per `config/prereg.yaml`
(`cofold.msa`) — the route used for **all** members of a case unless every member
matches a cached AF3 sequence exactly.

```bash
PY=/tmp/af3_env/bin/python
# MSA only (cached by sequence SHA-256 under /tmp/pe_big/msas/<sha>/)
$PY scripts/cofold/colabfold_msa.py --acc P24941
# MSA + AF3 input JSON in one go
$PY scripts/cofold/colabfold_msa.py --acc P24941 --ligand CCD:ATP --ligand CCD:MG \
    --name cdk2_atp_mg_cfmsa --out /tmp/af3_inputs/cdk2_atp_mg_cfmsa.json
# batch (multi-FASTA, one record per sequence to fold)
$PY scripts/cofold/colabfold_msa.py --fasta members.fa
# housekeeping
$PY scripts/cofold/colabfold_msa.py --export-compact --record-data-md
```

* **Protocol.** `POST {host}/ticket/msa` with `q=<fasta>`, `mode=env`
  (UniRef30 + ColabFoldDB envdb), poll `GET /ticket/<id>`, then
  `GET /result/download/<id>` (a tar.gz). The returned `uniref.a3m` and
  `bfd.mgnify30.metaeuk30.smag30.a3m` are concatenated into one unpaired A3M with the
  query as the first row, duplicate rows and the per-file query rows dropped,
  lowercase insertion columns kept (AF3 parses A3M).
* **AF3 JSON.** `unpairedMsa` = that A3M, `pairedMsa` = `""`, `templates` = `[]`.
  For a monomer + ligand this is exactly equivalent to the cached input minus the
  templates — see the `pairedMsa` check in §2.2 below.
* **Politeness.** One submission per 20 s node-wide (a lock file in the cache dir),
  polling backs off 5 → 30 s, and `RATELIMIT`/`MAINTENANCE` replies back off 60 s+.
  Everything is cached by sequence hash, so re-runs cost nothing. Measured:
  95 s and 207 s for two 298-aa queries (server queue 10–50 at the time).
* **Storage.** Full results (raw tar + merged a3m.gz + provenance) go to
  `/tmp/pe_big/msas/`. `--export-compact` copies the gzipped a3m plus its provenance
  JSON into `data/msas/`, stopping at a hard 50 MB budget (home quota). A 298-aa
  query with 20 k sequences is 3.0 MB gzipped, so the budget holds ~16 proteins;
  beyond that the compact copy is skipped and named in the log.
* **Provenance.** `--record-data-md` appends one row per MSA to `DATA.md` with
  server, mode, ticket id, UTC timestamp, the returned database file names (the
  server exposes no `/version` endpoint, so these names *are* its DB version
  strings), sequence count and a depth proxy.

### 2.2 Does the MSA source matter? (CDK2 test)

Three CDK2 + ATP + Mg runs, all `--no-templates`, 3 seeds × 5 samples, scored against
1HCK ATP. Median ± SD over the 15 samples; "top" is the top-ranked sample.

| MSA source | n seqs | lig pLDDT | P–L PAE min | ligand RMSD (Å) | centroid (Å) | P–L ipTM |
|---|---|---|---|---|---|---|
| cached AF3 jackhmmer, unpaired + paired | 18 071 + 50 000 | 91.4 ± 0.4 | 0.84 ± 0.01 | 1.00 ± 0.15 | 0.20 ± 0.06 | 0.95 |
| cached AF3 jackhmmer, unpaired only (`pairedMsa=""`) | 18 071 | 91.4 ± 0.4 | 0.84 ± 0.01 | 1.00 ± 0.15 | 0.20 ± 0.06 | 0.95 |
| **ColabFold MMseqs2 (UniRef30 + envdb), unpaired** | 20 013 | 90.7 ± 0.3 | 0.85 ± 0.01 | 1.03 ± 0.13 | 0.29 ± 0.09 | 0.95 |
| *(reference)* CDK2 + caffeine decoy, cached MSA | | 74.6 | 1.61 | n/a | 4.11 | 0.85 |

Two conclusions.

1. **`pairedMsa` is inert for a monomer + ligand.** Emptying the 50 000-row paired MSA
   reproduced the templated-free run bit-for-bit (top ranking score identical to 16
   significant figures). AF3 only uses the paired block to pair rows *across protein
   chains*. So an AF3 JSON built from a ColabFold unpaired MSA with `pairedMsa=""` is
   directly comparable to a cached-MSA job, and the comparison above is apples to apples.
2. **The MSA source does not matter at the scale we care about.** ColabFold costs
   0.7 pLDDT, 0.01 PAE and 0.09 Å of centroid accuracy relative to jackhmmer — about
   two seed-to-seed standard deviations, and 20–30× smaller than the cognate-vs-decoy
   gap (16.8 pLDDT, 0.77 PAE, 3.9 Å) that the readouts have to resolve. Using ColabFold
   MSAs for *all* members of a case (prereg `cofold.msa`) therefore costs almost nothing
   and removes the MSA-provenance confound entirely (round-1 review change 11).

Caveat: one protein, a well-populated family (both MSAs ~20 k sequences). For a sparse
family the two pipelines could diverge, so `msa_source.tsv` records the per-member
sequence count and depth proxy for every case, and any case whose members differ by
more than ~3× in depth should be flagged.

## 3. Run inference: `scripts/cofold/run_af3.sh`

```bash
scripts/cofold/run_af3.sh -i /tmp/af3_inputs/cdk2_atp.json -j cdk2_atp -g 0
scripts/cofold/run_af3.sh -i /tmp/af3_inputs/ -j batch1 -g 1     # directory of *.json
```

* **Behaviour.** The script runs `run_alphafold.py --norun_data_pipeline`
  pinned to GPU `-g` via `CUDA_VISIBLE_DEVICES`. It uses
  `--compress_large_output_files`, so the CIF and confidence files are `.zst`.
* **Outputs (home quota rule).**
  * Full outputs go to `/tmp/pe_big/cofold/<job>/<name>/` (override with
    `AF3_OUT_ROOT`). They are regenerable and node-local, and include all
    samples, PAE matrices and the input copy.
  * Only small summaries are copied to `results/cofold/<job>/`: the logs, GPU
    memory trace, `*_ranking_scores.csv`, `*summary_confidences.json` and the
    top-ranked `*_model.cif.zst` (about 60 kB).
  * Copy any per-sample model you need into the repo explicitly.
  * Inputs live in `/tmp/pe_big/cofold/inputs/`.
* **Samples.** The default is 5 diffusion samples per seed (`-s`).
* **Memory modes (`-m`).**
  * `prealloc` (default): 95% of GPU memory is preallocated, which is fastest.
  * `grow`: no preallocation. This is used for the benchmark so that
    `nvidia-smi` shows real usage.
  * `unified`: `TF_FORCE_UNIFIED_MEMORY`, spilling to host RAM, for inputs too
    large for 24 GB.
* **Run records.**
  * `run_af3.log` ends with `# rc= wall_s= peak_gpu_mem_mib=`.
  * `gpu_mem.csv` holds 1 Hz memory samples.
  * The `<name>_data.json` that AF3 writes (a full copy of the input with
    MSAs, about 30 MB) is gzipped.
* **Parallelism.** Run one process per GPU. A directory input amortises model
  build and compile over all jobs in it.

## 4. Benchmark

Run with `scripts/cofold/benchmark_af3.sh`, which uses 4 GPU queues in
parallel, then analyse with
`/tmp/af3_env/bin/python scripts/cofold/analyze_af3.py`.

* **Analysis inputs and outputs.** The analysis reads
  `/tmp/pe_big/cofold/bench` and writes
  `results/cofold/bench/{runtime,placement_samples,placement_summary}.tsv`.
* **Memory-limit runs.** `mem_*` were run ad hoc with `run_af3.sh -m prealloc`
  and `-m unified`, on inputs built by `make_af3_input.py --copies 2|3`.
* **Settings.** All runs used 5 diffusion samples per seed and 10 recycles
  (the AF3 defaults), on one RTX 4090 each, with 4 jobs in parallel on the
  node.

### 4.1 Runtime and GPU memory (RTX 4090, 24 GB)

"Per seed" means one seed with 5 diffusion samples, timed with the compile
cache warm. The first seed in a new bucket adds the XLA compile, about 5 to
35 s. "Job wall" is the whole `run_af3.sh` call: Python start-up, model build
and weight loading (about 15 to 20 s), featurisation of all seeds, inference
and writing. Peak memory in `grow` mode is the XLA pool high-water mark. Note
that `grow` mode is still capped at 75% of the GPU (18 GB), which is why the
2048 bucket failed there.

Tokens = residues + ligand heavy atoms, padded up to the AF3 bucket.

| Job (system) | Residues | Tokens (bucket) | First seed (s) | Per extra seed (s) | Job wall, 3 seeds (s) | Peak GPU (GiB) |
|---|---|---|---|---|---|---|
| CA2 + Zn + AZM | 260 | 274 (384) | 21.8 | 16.8 | 77 | 4.4 |
| CDK2 + ATP + Mg | 298 | 330 (384) | 20.7 (cold compile: 47.4) | 16.7 | 80 | 4.4 |
| yeast ADH1 + NAD + 2 Zn | 348 | 394 (512) | 54.9 (cold) | 25.1 | 132 | 4.4 |
| PYGL + PLP | 847 | 863 (1024) | 109.3 (cold) | 78.8 | 312 | 8.4 |
| CDK2 x4 + ATP | 4x298 | 1223 (1280) | 160.7 | 124.5 | 329 (2 seeds) | 16.4 |
| CDK2 x5 + ATP | 5x298 | 1521 (1536) | 213.2 | 179.8 | 442 (2 seeds) | 16.4 |
| PYGL x2 + 2 PLP, `prealloc` | 2x847 | 1726 (2048) | 316 (cold) | about 285 (est.) | 352 (1 seed) | fits in 0.95 x 24 GB |
| PYGL x2, `grow` (75% cap) | 2x847 | 1726 (2048) | OOM | | | |
| PYGL x3 + 3 PLP, `prealloc` | 3x847 | 2589 (3072) | **OOM** (needs a 23.7 GiB alloc) | | | |
| PYGL x3, `unified` | 3x847 | 2589 (3072) | 1087 | | 1141 (1 seed) | 24 GB + host spill |

Sources: `results/cofold/bench/runtime.tsv` and `run_af3.log` files.

**Memory limit.**
* Inputs up to the **2048-token bucket fit on one 24 GB RTX 4090** with the
  default `prealloc` mode.
* The 2560/3072 buckets need `-m unified`, which spills to host RAM. This
  works, but at 3072 it is about 3.8x the 2048 time, roughly 1.5x the A100
  number.
* Unified memory works until host RAM runs out (the node has 375 GB). Above
  2048 tokens, `pair_transition_shard_spec` could also be tuned (see
  docs/performance.md); this was not tested.
* For a single chain plus small ligands, the practical maximum is about 2,000
  residues on-GPU.

**Speed vs A100.** The 4090 is about 1.3x slower than the published A100-80GB
compile-free numbers: 79 s vs 62 s at 1024 tokens.

### 4.2 GPU-hours budget formula

Per-seed time (5 samples, warm cache) scales as about the bucket size to the
power 1.8:

    t_seed(B) ≈ 79 s × (B / 1024)^1.8        # B = AF3 bucket ≥ tokens
    t_job     ≈ 35 s + n_seeds × t_seed(B)   # + ~30 s once per new bucket per node (compile, cached after)
    GPU-h     = Σ_jobs t_job / 3600 ;  wall-h ≈ GPU-h / 4  (4 GPUs, one job stream each)

Rules of thumb with 3 seeds × 5 samples:

| Bucket | Typical system | t_seed | t_job | Jobs per GPU-hour |
|---|---|---|---|---|
| 384 | ≤ about 350 aa + ligand | 17 s | about 1.4 min | about 43 |
| 512 | ≤ 500 aa | 25 s | about 1.8 min | about 33 |
| 768 | ≤ 750 aa | about 47 s | about 3 min | about 20 |
| 1024 | ≤ 1000 aa | 79 s | about 4.5 min | about 13 |
| 1536 | ≤ 1500 aa | 180 s | about 9.6 min | about 6 |
| 2048 | ≤ 2000 aa | about 285 s | about 15 min | about 4 |

The 768-bucket row is interpolated; the other rows are measured.

**Example budget.**
* Suppose 30 cases × 2 states (ancestral and derived) × (1 cognate ligand +
  2 decoys) × 1 ortholog pair, for 180 jobs at about 450 aa (bucket 512).
* That is about 180 × 1.8 min ≈ 5.4 GPU-h, or about 1.4 h wall on 4 GPUs.
* With a ±templates ablation it doubles to about 11 GPU-h.

**Reducing per-job overhead.** Batch many JSONs in one directory per GPU
(`run_af3.sh -i <dir>`). This loads the model once and reuses compiled
buckets, saving the roughly 20 s start-up per job.

### 4.3 Ligand placement sanity: positive controls and decoys

Method:
* Each sample is superposed onto the holo PDB chain using sequence-aligned
  CA atoms.
* The **site** is the set of holo residues within 5 Å of the crystal ligand.
* **centroid** is the distance between the predicted and crystal ligand
  centroids.
* **RMSD** is the heavy-atom RMSD over atoms with matching names, without
  symmetry correction, for CCD ligands only.
* **recall** is the fraction of site residues within 5 Å of the predicted
  ligand.
* Decoys are scored against the cognate ligand's site.
* The table shows the top-ranked of 15 samples; "in-site" is the fraction of
  the 15 samples with recall ≥ 0.5.

Source: `results/cofold/bench/placement_summary.tsv`, with per-sample values
in `placement_samples.tsv`.

| Job | Reference | P–L ipTM | P–L PAE min | Ligand pLDDT | Centroid (Å) | RMSD (Å) | Recall | In-site |
|---|---|---|---|---|---|---|---|---|
| **CDK2 + ATP + Mg** (pos) | 1HCK ATP | 0.95 | 0.82 | 92.2 | 0.12 | 0.74 | 1.00 | 1.00 |
| CDK2 + ATP + Mg, **no templates** | 1HCK ATP | 0.96 | 0.84 | 92.1 | 0.10 | 1.00 | 1.00 | 1.00 |
| **CA2 + Zn + AZM (CCD)** (pos) | 3HS4 AZM | 0.98 | 0.78 | 97.4 | 0.14 | 0.28 | 1.00 | 1.00 |
| CA2 + Zn + AZM (SMILES) | 3HS4 AZM | 0.97 | 0.79 | 96.7 | 0.12 | n/a | 1.00 | 1.00 |
| **PYGL (847 aa) + PLP** (pos, 1024 bucket) | 1FA9 PLP | 0.99 | 0.89 | 97.3 | 0.49 | 1.22 | 1.00 | 1.00 |
| **yeast ADH1 + NAD + 2 Zn** (yeast) | 4W6Z 8ID (NAD analogue) | 0.97 | 0.83 | 97.2 | 0.21 | 0.41 | 0.97 | 1.00 |
| CDK2 + AZM (decoy) | 1HCK ATP site | 0.86 | 1.54 | 74.4 | 4.71 | n/a | 0.54 | 1.00 |
| CDK2 + caffeine (decoy) | 1HCK ATP site | 0.85 | 1.61 | 74.6 | 4.11 | n/a | 0.50 | 1.00 |
| CA2 + Zn + ATP (decoy) | 3HS4 AZM site | 0.80 | 3.48 | 58.2 | 4.63 | n/a | 1.00 | 1.00 |
| CA2 + Zn + caffeine (decoy) | 3HS4 AZM site | 0.72 | 3.50 | 59.4 | 4.48 | n/a | 0.46 | 0.73 |
| ADH1 + 2 Zn + caffeine (decoy) | 4W6Z cofactor site | 0.77 | 4.02 | 72.3 | 7.66 | n/a | 0.29 | 0.00 |

Interpretation:
* All four positive controls are reproduced to within about 1 Å, in all 15
  samples, for both human and yeast proteins. This holds even without
  templates (CDK2), so the result is not just template copying, although the
  MSA and the training set still include these PDB complexes.
* **Decoys are not rejected by placement.** AF3 docks every decoy into, or
  next to, the known pocket, and gives it a high P–L ipTM of 0.72 to 0.86.
  The discriminating signals are:
  * ligand pLDDT: 92 to 97 for true ligands vs 58 to 75 for decoys;
  * P–L PAE min: below 0.9 for true ligands vs 1.5 to 4 for decoys;
  * a smaller ipTM gap of about 0.1 to 0.25.
* For ancestral-vs-derived comparisons, report ligand pLDDT, P–L PAE min,
  ipTM and pose consistency across seeds. Always run matched decoys on both
  states, and judge differences relative to the decoy baseline, not against an
  absolute ipTM threshold.

## 5. Caveats

* **Environment.** The venv lives on node-local `/tmp` because of the home
  quota, so rerun `setup_af3_env.sh` (about 2 min) on each new allocation.
  The compile cache and full outputs in `/tmp/pe_big` are also lost with the
  job.
* **Memory measurement.** In `grow` mode, XLA still caps the pool at 75% of
  GPU memory. Use `prealloc` (the default) for production and `unified` for
  buckets of 2560 or more.
* **Memory figures.** The peak GPU numbers are XLA pool high-water marks
  sampled at 1 Hz, not exact live-tensor peaks.
* **Template and training leakage.** The cached JSONs include up to 4 PDB
  templates, which for these controls are holo structures of the same protein.
  Positive controls are also in AF3's training set. Expect real cases, such as
  ancestral sequences or non-human orthologs, to be less confident. Use
  `--no-templates` as an ablation.
* **Ancestral or variant sequences.** These do not match the cache.
  `--reuse-msa-for-variant` is an approximation that keeps the extant
  ortholog's MSA and templates, so it is biased toward the extant state.
  State this whenever it is used, or compute a proper MSA (ColabFold server).
* **Covalent ligands.** PLP is covalent (a Schiff base to Lys680) in the
  crystal but was folded non-covalently. Covalent ligands need
  `bondedAtomPairs`, which is not exposed in `make_af3_input.py` yet.
* **Metric limitations.**
  * Ligand RMSD uses atom-name matching without symmetry correction.
  * For SMILES ligands only the centroid and recall are computed.
  * The 4W6Z reference ligand is 8-iodo-NAD (8ID), so NAD's 44 atoms are
    compared by name to it.
* **AF3 terms.** The weights are under the AF3 model terms (non-commercial).
  Outputs carry `TERMS_OF_USE.md`, which is kept in the summaries.


## 6. Negative controls on the positive controls

See `docs/COFOLD_CONTROLS.md` for the site-ablation, ligand-free refold and decoy
"exchange rates" measured on the four positive-control systems, and
`results/cofold/controls/{summary,exchange_rate,apo_burrow}.tsv` for the tables.
Scripts: `scripts/cofold/site_ablation.py` (site definition and Gly/Phe variants),
`scripts/cofold/run_controls.py` (build + GPU queue + scoring),
`scripts/cofold/cofold_readouts.py` (the preregistered readouts, including PoseBusters),
`scripts/cofold/apo_describe.py` (ligand-free refold + `pocketcompare describe`).

## 7. Running a case: `scripts/cofold/case_cofold.py`

```bash
PY=/tmp/af3_env/bin/python
$PY scripts/cofold/case_cofold.py --spec cases/cofold/<case>.yaml --build --run --unblind
```

It implements `config/prereg.yaml` (`cofold`) literally: one MSA pipeline for all
members (ColabFold unless every member is cached-exact), `--no-templates`, blinded job
names, per-seed decision statistic that must hold in all three seeds, and the
positive-control ligand-RMSD gate that labels a case `uninterpretable` rather than
dropping it. See the module docstring for the YAML schema and
`cases/cofold/cdk2_toy.yaml` for a worked example.

## 8. Calibration panel and the frozen threshold

`scripts/cofold/calibration_panel.py` builds a BioLiP2-derived panel of holo positives
with cached AF3 MSAs, by the deterministic rule R1-R7 in its docstring (also dumped to
`results/cofold/calibration/selection_rule.json`). Per decision D23 the negatives are
*both* property-matched decoys *and* the Gly site ablation of each positive, because a
case compares protein variants, not only ligands.

```bash
PY=/tmp/af3_env/bin/python
$PY scripts/cofold/calibration_panel.py --select --per-class 6
$PY scripts/cofold/calibration_panel.py --build --run --gpus 0,1,2
$PY scripts/cofold/calibration_panel.py --score      # writes threshold.json
```

Decoys are matched **within ligand class** on MW, logP, HBD, HBA and rotatable bonds
against the full in-class BioLiP ligand pool, minus anything BioLiP records as binding
that accession. Cross-class matching was tried first and rejected: it paired Ca sites
with organic molecules, which is a gross-mismatch control rather than a decoy and
flatters the readouts. The panel decoys are therefore hard near-misses (SAM/SAH,
GDP/ADP/CTP, Ca/Mg/Mn, ACO/FAD) and the panel AUC sits below the caffeine-decoy numbers
in §4.3. That lower number is the honest operating point for case scoring.

`results/cofold/calibration/threshold.json` is written **before any case is folded** and
is not to be retuned afterwards. It carries the per-readout AUCs with bootstrap CIs, the
same AUCs split by negative kind (decoy vs ablation), the D23 composite thresholds with
their sensitivity and specificity, and the measured per-job runtime.

## 9. Corrections to earlier numbers, and the blinding status of the Tier 1 cases

**Tier 1 cases are UNBLINDED.** The four pool cases in §7 (citrate synthase, MetE, choline
kinase, RNA polymerase) folded blind — inputs carried HMAC-hashed names and the state
labels lived only in `blind_key.json` — but `build()` minted a fresh salt on every call,
which destroyed the original key. Recovering it meant matching each AF3 output back to its
member by sequence identity, and that happened before the final scoring pass, so the
labels were visible while the readouts were produced. The readout code never receives a
state label, so the measurements are label-independent, but the preregistered blinding
claim cannot be made for these four. It is recorded as `blinding.status: UNBLINDED` in
each `case.yaml` and in each `results/cofold/cases/<id>/decision.json`. The salt is now
persisted; the deep dives in §10 were scored without reconstructing any mapping and are
blinded.

**Centroid distances, corrected.** An earlier commit message quoted "21–34 Å" for MetE.
That range came from its *decoy* jobs. The correct figures, per-job medians over 3 seeds:

| case | cognate-job medians | decoy-job medians |
|---|---|---|
| citrate synthase | 19.3–23.5 Å | 19.7–24.6 Å |
| MetE | 28.9–34.1 Å | 23.8–29.5 Å |
| choline kinase | 12.2–18.7 Å | 11.7–18.1 Å |

Over all cognate samples of the three cases the spread is 12.0–49.0 Å, so the summary
statement is "cognate medians 12–34 Å, individual samples 12–49 Å". The conclusion is
unchanged: the cognate ligand never reaches the atlas footprint, on either side.

**Sensitivity is visible in the toy case.** `cdk2_toy`'s true positive (human CDK2 + ATP)
sits at ligand pLDDT 90.6–90.8 against the frozen threshold of 92.0, with centroid 1.8 Å
and site recall 1.00. A real binder, placed correctly, missed by 1.2 pLDDT points. That is
the frozen sensitivity of 0.42 in one picture, and it is why a "does not bind" call is
weak evidence. The threshold was not adjusted.

## 10. Known-answer deep dives

Two 2x2 designs (protein x ligand), 3 seeds x 5 samples, no templates, ColabFold MSAs,
blinded with a persisted salt. The success criterion was committed in `21b5c9e` *before*
either dive finished folding. Results: `results/cofold/deepdives/<case>/result.json`.

### 10.1 Steroid receptor specificity — RESOLVED

MR LBD and GR LBD against aldosterone (AS4) and cortisol (HCY). Known answer: MR binds
both with high affinity; GR prefers cortisol and binds aldosterone poorly.

| | aldosterone | cortisol |
|---|---|---|
| **MR** | 98.2 | 98.5 |
| **GR** | **94.5** | 97.7 |

Ligand pLDDT, mean over 5 samples, seed 1 (seeds 2 and 3 agree to within 0.3). Gate passed
(MR/2AA2 aldosterone RMSD 0.31 Å); all four cells bind, with centroid 0.5–1.5 Å and site
recall 1.00. The interaction (GR.HCY − GR.AS4) − (MR.HCY − MR.AS4) is +2.96, +2.68, +2.80
across the three seeds: the expected sign in all three, for ligand pLDDT and for min PAE.
**Co-folding does resolve this known specificity difference.**

Caveat on size: the discriminating signal is ~3 pLDDT points, against a cognate-vs-decoy
gap of 18–42 points in the controls. It is reproducible (seed SD ~0.1) but small, so this
says co-folding can rank two real ligands at one site, not that it would detect such a
shift blind.

### 10.2 LDH vs MDH substrate specificity — biology right, gate failed

PfLDH and human MDH1 against pyruvate (PYR) and oxaloacetate (OAA), both with NADH.

| | pyruvate | oxaloacetate |
|---|---|---|
| **LDH** | 96.8 | 86.7 |
| **MDH** | 84.5 | **93.3** (only cell that binds) |

Each enzyme prefers its own substrate, and the interaction is +18.9, +9.9, +17.5 pLDDT
with the expected sign in all three seeds. But the positive-control gate FAILED
(hMDH1/7RM9 malonate RMSD 2.66 Å > 2.0 Å), so by the pre-committed rule the dive is
**uninterpretable** and is reported as such. The verdict was not overridden.

The likely cause of the gate failure is a known limitation, not a bad model: malonate is a
small symmetric dicarboxylate and our ligand RMSD matches atoms by name with no symmetry
correction (§5), so a 180° flip inflates it. The correct fix is a symmetry-corrected RMSD,
not a looser gate; until then this dive stays uninterpretable.

**A measurement bug this dive exposed.** `pick_ligand()` took the largest non-ion ligand,
so with NADH (44 atoms) present it measured the *cofactor* rather than the substrate, and
the first scoring pass showed all four cells at pLDDT 93–97 with centroid ~12 Å — NADH's
distance to the substrate site. `collect()` now focuses on the ligand that DIFFERS between
the cognate and decoy sets, treating anything shared as context. The table above is the
corrected substrate measurement; the pre-committed criterion was not changed.
