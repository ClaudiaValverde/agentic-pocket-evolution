#!/usr/bin/env bash
# Run AlphaFold3 inference only (no data pipeline) on one JSON or a directory of JSONs.
#
# Usage: run_af3.sh -i <input.json|input_dir> -j <job> [-g GPU] [-s NUM_SAMPLES] [-m MODE] [-- extra run_alphafold.py flags]
#   -i  AF3 input JSON or directory of *.json (from make_af3_input.py)
#   -j  job name -> full outputs in $AF3_OUT_ROOT/<job>/ (default /tmp/pe_big/cofold, bulky,
#       regenerable); small summaries (logs, ranking/summary confidences, top model CIF)
#       are copied to $P/results/cofold/<job>/
#   -g  GPU index (default 0); pinned via CUDA_VISIBLE_DEVICES
#   -s  diffusion samples per seed (default 5, AF3 default)
#   -m  memory mode: prealloc (default, fastest, XLA_CLIENT_MEM_FRACTION=0.95),
#       grow (no preallocation; used for benchmarking peak GPU memory),
#       unified (spill to host RAM for inputs that do not fit in 24 GB; slower)
# Env overrides: AF3_ENV (venv, default /tmp/af3_env), AF3_SRC, AF3_MODEL_DIR,
#                AF3_JAX_CACHE (persistent compile cache, default /tmp/af3_jax_cache),
#                AF3_KEEP_DATA_JSON=1 (keep uncompressed <name>_data.json copy of input)
set -euo pipefail
P="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"  # repo root
AF3_ENV=${AF3_ENV:-/tmp/af3_env}
AF3_SRC=${AF3_SRC:-$P/envs/af3_src}
AF3_MODEL_DIR=${AF3_MODEL_DIR:-$P/data/af3_weights}
AF3_JAX_CACHE=${AF3_JAX_CACHE:-/tmp/af3_jax_cache}
AF3_OUT_ROOT=${AF3_OUT_ROOT:-/tmp/pe_big/cofold}

INPUT= ; JOB= ; GPU=0 ; NSAMP=5 ; MODE=prealloc
while getopts "i:j:g:s:m:h" o; do
  case $o in
    i) INPUT=$OPTARG;; j) JOB=$OPTARG;; g) GPU=$OPTARG;; s) NSAMP=$OPTARG;; m) MODE=$OPTARG;;
    h|*) sed -n 2,16p "$0"; exit 1;;
  esac
done
shift $((OPTIND-1))
[[ -n $INPUT && -n $JOB ]] || { sed -n 2,16p "$0"; exit 1; }
[[ -x $AF3_ENV/bin/python ]] || { echo "AF3 env missing at $AF3_ENV; run scripts/cofold/setup_af3_env.sh" >&2; exit 1; }
[[ -f $AF3_MODEL_DIR/af3.bin ]] || { echo "weights missing: $AF3_MODEL_DIR/af3.bin" >&2; exit 1; }

OUT=$AF3_OUT_ROOT/$JOB
# Home quota is tight (storage alert, 2026-10-30): by default only the compact readout
# files go into the repo -- the log tail, the job-level ranking CSV and summary
# confidences. Models and per-sample files stay under $AF3_OUT_ROOT. Set
# AF3_SUMMARY_MODELS=1 to also copy the top-ranked model, and AF3_SUMMARY_ROOT to redirect.
AF3_SUMMARY_ROOT=${AF3_SUMMARY_ROOT:-$P/results/cofold}
AF3_SUMMARY_MODELS=${AF3_SUMMARY_MODELS:-0}
SUMM=$AF3_SUMMARY_ROOT/$JOB
mkdir -p "$OUT" "$SUMM" "$AF3_JAX_CACHE"
if [[ -d $INPUT ]]; then INFLAG=(--input_dir="$INPUT"); else INFLAG=(--json_path="$INPUT"); fi

export CUDA_VISIBLE_DEVICES=$GPU
# Compile-time workaround from AF3 docs/performance.md (set in the official Dockerfile).
export XLA_FLAGS="--xla_gpu_enable_triton_gemm=false"
case $MODE in
  prealloc) export XLA_PYTHON_CLIENT_PREALLOCATE=true XLA_CLIENT_MEM_FRACTION=0.95;;
  grow)     export XLA_PYTHON_CLIENT_PREALLOCATE=false; unset XLA_CLIENT_MEM_FRACTION || true;;
  unified)  export XLA_PYTHON_CLIENT_PREALLOCATE=false TF_FORCE_UNIFIED_MEMORY=true XLA_CLIENT_MEM_FRACTION=3.2;;
  *) echo "bad mode $MODE" >&2; exit 1;;
esac

LOG=$OUT/run_af3.log
GPULOG=$OUT/gpu_mem.csv
# Sample GPU memory every second (this GPU only) for the benchmark record.
( while true; do nvidia-smi --id="$GPU" --query-gpu=timestamp,memory.used --format=csv,noheader,nounits; sleep 1; done ) > "$GPULOG" 2>/dev/null &
SAMPLER=$!
trap 'kill $SAMPLER 2>/dev/null || true' EXIT

echo "# $(date -Is) host=$(hostname) gpu=$GPU mode=$MODE samples=$NSAMP input=$INPUT" | tee "$LOG"
T0=$(date +%s.%N)
set +e
"$AF3_ENV/bin/python" "$AF3_SRC/run_alphafold.py" \
  "${INFLAG[@]}" \
  --output_dir="$OUT" \
  --model_dir="$AF3_MODEL_DIR" \
  --norun_data_pipeline \
  --flash_attention_implementation=triton \
  --num_diffusion_samples="$NSAMP" \
  --jax_compilation_cache_dir="$AF3_JAX_CACHE" \
  --compress_large_output_files \
  --force_output_dir \
  "$@" 2>&1 | tee -a "$LOG"
RC=${PIPESTATUS[0]}
set -e
T1=$(date +%s.%N)
kill $SAMPLER 2>/dev/null || true
PEAK=$(awk -F', ' 'BEGIN{m=0}{if($2>m)m=$2}END{print m}' "$GPULOG")
WALL=$(awk -v a="$T0" -v b="$T1" 'BEGIN{printf "%.1f", b-a}')
echo "# rc=$RC wall_s=$WALL peak_gpu_mem_mib=$PEAK" | tee -a "$LOG"
# The <name>_data.json AF3 writes is a full copy of the input incl. MSAs (tens of MB): compress it.
if [[ ${AF3_KEEP_DATA_JSON:-0} != 1 ]]; then
  find "$OUT" -name '*_data.json' -size +1M -exec gzip -f {} \;
fi
# Copy small summaries into the repo (home quota is ~50 GB; full outputs stay in $OUT).
# log tail only (the full AF3 log is ~300 kB and stays in $OUT)
{ grep -E "^# |took [0-9.]+ seconds|Fold job .* done|Error|Traceback" "$LOG" | tail -200; } > "$SUMM/run_af3.log.tail"
( cd "$OUT" && find . -maxdepth 2 \( -name '*_ranking_scores.csv' -o -name '*_summary_confidences.json' \
    -o -name 'TERMS_OF_USE.md' \) -print0 | xargs -0 -r cp --parents -t "$SUMM" )
if [[ $AF3_SUMMARY_MODELS == 1 ]]; then
  ( cd "$OUT" && find . -maxdepth 2 -name '*_model.cif*' -print0 | xargs -0 -r cp --parents -t "$SUMM" )
fi
echo "# full outputs: $OUT ; summaries: $SUMM" | tee -a "$LOG" "$SUMM/run_af3.log.tail"
exit $RC
