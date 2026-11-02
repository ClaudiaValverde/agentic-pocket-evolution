# ENV — Python environment for burrow + pocketcompare (built 2026-10-30)

Venv: `P/.venv` (P = <repo>), CPython 3.12.3 (`/usr/bin/python3.12`; both
tools are pinned to and tested with 3.12.3). Full pin list: `P/requirements.lock` (`uv pip freeze`, 39 packages).
Key versions: pocketcompare 1.0.0 and burrow-pockets 1.0.0 (both editable from `P/external/`), numpy 2.5.4,
scipy 1.18.1, pandas 3.1.0, pyarrow 26.0.0, numba 0.68.0, lightgbm 4.7.0, gemmi 0.7.5, faiss-cpu 1.14.3,
duckdb 2.0.0, scikit-learn 1.9.1, biopython 1.88, requests 2.34.2, matplotlib 3.11.2, pyyaml 6.0.3, pytest 9.1.1.

## Home quota: read this before installing anything
The first attempt (with `torch==2.6.0 triton==3.2.0` for the GPU backend) failed with
`Quota exceeded (os error 122)`: `<home>` (and `/home`, the same NFS export) was at ~47 GB
(Euler home quota is presumably 50 GB; `lquota`/`quota` are not available in the container), and uv's
default cache (`/home/.cache/uv`) plus the CUDA wheels (~7 GB) filled it. That cache was created by this install
and was deleted again. Rules:
- always `export UV_CACHE_DIR=/tmp/pe_big/uvcache` (PI convention; /tmp NVMe has 1.7 TB) and use `--link-mode=copy`;
  put extracted models and test outputs under /tmp/pe_big/work;
- the torch/triton GPU extra was **not** installed (CPU fast backend only). If the GPU backend is needed
  (`--device cuda:N`), install it into a separate venv under /tmp or project storage, not into home;
- large outputs (indices: a flat tier-1 index is 2.7 GB) must not go under P on home either.

## Build commands (exact)
```bash
cd <repo>
export UV_CACHE_DIR=/tmp/pe_big/uvcache
uv venv --python /usr/bin/python3.12 .venv
uv pip install --python .venv/bin/python --link-mode=copy \
    -e external/burrow -e "external/pocketcompare[search,test]" \
    pyarrow duckdb pandas numpy scipy scikit-learn biopython gemmi requests matplotlib pyyaml
uv pip freeze --python .venv/bin/python > requirements.lock
```
burrow finds its models in `external/burrow/models/` (editable install); pocketcompare uses its own pinned copy
`external/pocketcompare/models/burrow/burrow_all_default_s0.json` (sha256 f33963e7…, identical to burrow's file).
The first call compiles Numba kernels (~15 s; the cache is written next to the sources under `external/`).
Set `OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1` for timing / many workers.

## Verification (all run 2026-10-30, CPU, this container; scratch in /tmp/pe_big/work (moved there from /tmp/pe_env after the runs, per PI storage note))
Helper scripts written for this: `scripts/extract_afdb_member.py` (seeks single members of the AFDB v6 tars
via the atlas tar index, no full extraction), `scripts/verify_atlas_fingerprint.py` (fresh describe vs stored
atlas rows, applying the atlas filters score ≥ 0.03 and mean lining pLDDT ≥ 70),
`scripts/burrow_settings_jaccard.py` (pocket matching across BURROW settings by residue Jaccard).

```bash
P=<repo>; mkdir -p /tmp/pe_big/work && cd /tmp/pe_big/work
export OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1
# 1. one real AFDB v6 model (HMGCR, P04035, 888 aa) -> describe -> compare with tier-1 atlas rows
$P/.venv/bin/python $P/scripts/extract_afdb_member.py --tar UP000005640_9606_HUMAN_v6.tar --acc P04035 P07900 P00918 --out afdb
$P/.venv/bin/pocketcompare describe afdb/AF-P04035-F1-model_v6.cif.gz --af-mode -o hmgcr_ref.parquet
$P/.venv/bin/python $P/scripts/verify_atlas_fingerprint.py hmgcr_ref.parquet --tier 1
#  -> 13 pockets, 10 pass the atlas filters, the same 10 pocket_ids are in the atlas; cosine 1.0000000 for all,
#     max|Δ| 5.8e-5 (the atlas stores float16), residue lists identical. Exit 0.
# 2. 50 random human F1 models (+ the 3 above) -> same check, reference and fast backends
$P/.venv/bin/pocketcompare describe rand50 --af-mode --backend ref  -o r50_ref.parquet     # 0.401 s/structure
$P/.venv/bin/pocketcompare describe rand50 --af-mode --backend fast -o r50_fast.parquet    # 0.213 s/structure (warm)
$P/.venv/bin/python $P/scripts/verify_atlas_fingerprint.py r50_ref.parquet --tier 1
#  -> 233 pockets, 177 kept = 177 atlas rows, all cosines 1.0000000, all residue lists equal;
#     ref vs fast max|Δ| = 0.0 (bit-identical)
# 3. tier 2 (Swiss-Prot tar): B9JAQ3 B9JAR2 B9JAX6 -> 9/9 rows, cosine 1.0000000
$P/.venv/bin/python $P/scripts/extract_afdb_member.py --tar swissprot_cif_v6.tar --acc B9JAQ3 B9JAR2 B9JAX6 --out t2
$P/.venv/bin/pocketcompare describe t2 --af-mode --backend fast -o t2.parquet
$P/.venv/bin/python $P/scripts/verify_atlas_fingerprint.py t2.parquet --tier 2
# 4. site modes, compare, index, search
$P/.venv/bin/pocketcompare describe afdb/AF-P00918-F1-model_v6.cif.gz --af-mode --residues A:95,A:97,A:120,A:199,A:200 -o ca2_res_shell.parquet
$P/.venv/bin/pocketcompare describe afdb/AF-P00918-F1-model_v6.cif.gz --af-mode --residues A:95,A:97,A:120,A:199,A:200 --locate burrow -o ca2_res_burrow.parquet
$P/.venv/bin/pocketcompare compare ca2_res_burrow.parquet ca2_all.parquet --top 3     # m1 1.0000 (shell mode: 0.644)
E=$P/external/pocketcompare/examples
$P/.venv/bin/pocketcompare describe $E/pdb/1a28_A.pdb --ligand $E/pdb/1a28_STR.pdb -o site.parquet
$P/.venv/bin/pocketcompare compare site.parquet pr_af.parquet --top 2                  # m2 0.9863 = README value
cd $P/external/pocketcompare && $P/.venv/bin/python afdb/build_index.py \
   --out-root /cluster/scratch/gmagnusson/pocketcompare_afdb --tiers 1 --kind flat -o /tmp/pe_big/work/idx_t1_flat --threads 16
#  -> 1,788,750 pockets, 28.6 s, 2.7 GB index.faiss + 47 MB meta.parquet, peak RSS 4.0 GB
$P/.venv/bin/pocketcompare search /tmp/pe_big/work/hmgcr_ref.parquet /tmp/pe_big/work/idx_t1_flat -k 10 --threads 16 -o hits.tsv  # 1.3 s for 13 queries incl. load
# 5. BURROW alternative settings (pocket-finder stability), HMGCR
for args in "--af-mode" "" "--af-mode --preset fast" "--af-mode --preset accurate" "--af-mode --min-score none" \
            "--af-mode --output-rule region" "--af-mode --ranker small" "--af-mode --ranker geometry" \
            "--af-mode --af-ranker $P/external/burrow/models/afr_final_afpdb_pl.json" "--af-mode --backend reference" \
            "--af-mode --model $P/external/burrow/models/burrow_all_default_s1.json"; do
  $P/.venv/bin/burrow predict afdb/AF-P04035-F1-model_v6.cif.gz --out bw/<name> $args; done
$P/.venv/bin/python $P/scripts/burrow_settings_jaccard.py --stem AF-P04035-F1-model_v6 --top 6 --run def=bw/af_default --run noaf=bw/noaf ...
$P/.venv/bin/pocketcompare describe rand50 --backend fast -o r50_noaf.parquet                 # fingerprints without AF mode
# 6. test suites
(cd $P/external/pocketcompare && $P/.venv/bin/python -m pytest -q -p no:cacheprovider)
(cd $P/external/burrow && $P/.venv/bin/python scripts/smoke.py)
```
All BURROW settings ran (rc 0); results are in `docs/tools_pocket.md` §7. Test-suite results: see the end of
this file.

## Not verified
- GPU backend (`--device cuda:N`): torch not installed (quota). 4 × RTX 4090 are visible (`nvidia-smi -L`).
- Tier-3 rows (Foldcomp-decompressed v4 models): the `foldcomp` package is not installed and the original
  decompressed PDBs are not kept; POCKETOME.md §8 says re-describing the decompressed entry reproduces them.

## Test-suite results (2026-10-30)
- `pocketcompare` pytest: **39 passed, 2 skipped** in 808 s (the skips need torch/GPU; incl. the golden test
  reproducing research fingerprints to ≤ 1e-6). Note: run without thread caps it used ~39 cores; set
  `OMP_NUM_THREADS`/`NUMBA_NUM_THREADS` on the shared node.
- `burrow` `scripts/smoke.py`: **PASSED** in 5 s (default, fast preset and AF mode match `examples/expected/`).

## Part A env and tests (added 2026-10-30 14:42 CET)
- Part A GPU stages (`scripts/partA/02_knn.py`, `04_null.py`, `05_sig.py`, `08_topk_baselines.py`,
  `13_directional.py`) run in a separate venv on node-local disk, because faiss-gpu + torch do not fit the home quota:
  `UV_CACHE_DIR=/tmp/pe_big/uvcache uv venv -p 3.12 /tmp/pe_big/envs/faiss && uv pip install -p /tmp/pe_big/envs/faiss/bin/python faiss-gpu-cu12 numpy pandas pyarrow scipy duckdb rdkit pyyaml matplotlib scikit-learn biopython pytest`
  then `uv pip install -p /tmp/pe_big/envs/faiss/bin/python torch --index-url https://download.pytorch.org/whl/cu124`
  (faiss 1.14.1, torch 2.6.0+cu124). Use `CUDA_VISIBLE_DEVICES=3`.
- Tests: `uv run --project . pytest tests/test_partA.py tests/test_descent.py` passes in the project env. The 3 tests
  that need torch (GPD tail calibration x2, CovModel block form) are **skipped with the message "torch not installed
  in this env (GPU path optional; docs/ENV.md)"**. All 26 run in the faiss env:
  `/tmp/pe_big/envs/faiss/bin/python -m pytest tests/test_partA.py tests/test_descent.py`.
  torch was not added to the project env (about 2.5 GB with CUDA libraries, home quota).

## Part B orthology RBH check (scripts/partB/rbh_check.py)
- DIAMOND 2.2.8 (`envs/ancestra/bin/diamond`, already in the Ancestra env; blastp --sensitive, --max-target-seqs 1); per-proteome DBs in /tmp/pe_big/partB/diamond_db (regenerable). MMseqs2 18.8cc5c is also present in that env (not used).

## Cross-OG site matching (scripts/partB/crossog.py, exploratory D31)
- Foldseek 10.941cd33 (`envs/ancestra/bin/foldseek version`; already in the Ancestra env, nothing installed; binary sha256 ce5f08d8f1e15598…), `structurealign` on exact pairs (tsv2db prefilter DB, equivalent to `--exhaustive-search` for the listed pairs), default 3Di+AA, -e 10. DB: `/tmp/pe_big/partA/foldseek/afdb_v6` (Part A, AFDB v6 tier 1/2, 1,102,932 entries). Rescue reuses `pocketevo.rescue` (BURROW from .venv). CPU only.
