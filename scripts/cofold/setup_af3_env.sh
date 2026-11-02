#!/usr/bin/env bash
# (Re)create the AlphaFold3 environment. Idempotent; ~2 min with a warm uv cache.
# The venv (~5.8 GB, CUDA wheels) lives on node-local /tmp because the home/project
# quota (50 GB, nearly full) cannot hold it. /tmp vanishes with the job: rerun this.
set -euo pipefail
P="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"  # repo root
AF3_ENV=${AF3_ENV:-/tmp/af3_env}
SRC=$P/envs/af3_src
AF3_COMMIT=${AF3_COMMIT:-8f8ac383afba995f5953182d39a560ffb8bbbfd3}
export UV_CACHE_DIR=${UV_CACHE_DIR:-/tmp/pe_big/uvcache} UV_LINK_MODE=copy UV_PROJECT_ENVIRONMENT=$AF3_ENV
if [[ ! -d $SRC/.git ]]; then
  git clone https://github.com/google-deepmind/alphafold3.git "$SRC"
fi
git -C "$SRC" fetch --depth 1 origin "$AF3_COMMIT" 2>/dev/null || true
git -C "$SRC" checkout -q "$AF3_COMMIT"
[[ -x $AF3_ENV/bin/python ]] || uv venv --python 3.12 "$AF3_ENV"
(cd "$SRC" && uv sync --frozen)            # builds the C++ extension (cmake/ninja from PyPI, system gcc 13)
uv pip install --python "$AF3_ENV/bin/python" gemmi==0.7.5   # analysis only
# CCD / chemical component pickles (written into $SRC/src/alphafold3/constants/converters/)
[[ -f $SRC/src/alphafold3/constants/converters/chemical_component_sets.msgpack ]] || (cd "$SRC" && "$AF3_ENV/bin/build_data")
# Weights: decompress once into the workspace
W=$P/data/af3_weights/af3.bin
[[ -f $W ]] || { mkdir -p "$(dirname "$W")"; zstd -d /cluster/project/beltrao/shared/alphafold3/models/af3.bin.zst -o "$W"; }
"$AF3_ENV/bin/python" -c "import jax, alphafold3.cpp; print('jax', jax.__version__, jax.devices())"
