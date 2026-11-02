#!/usr/bin/env bash
# Rebuild all annotation tables for a snapshot manifest. Usage: scripts/annotations/run_all.sh <manifest.parquet>
set -euo pipefail
cd "$(dirname "$0")/../.."
M=${1:-data/snapshot/20261030T085136Z/manifest.parquet}
PY=".venv/bin/python -I"
export UV_CACHE_DIR=/tmp/pe_big/uvcache
$PY scripts/annotations/accessions.py --manifest "$M" --out data/annotations/accessions.parquet
$PY scripts/annotations/uniprot.py           # needs network; ~70 min; writes sequences to /tmp/pe_big/work
for s in mcsa pfam chembl ted biolip; do $PY scripts/annotations/$s.py & done; wait
.venv/bin/python -m pytest tests/test_annotations.py -q
