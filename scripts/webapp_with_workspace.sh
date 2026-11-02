#!/usr/bin/env bash
# Start the research browser with the bulk results of a full project workspace (on the cluster).
#
#   bash scripts/webapp_with_workspace.sh /path/to/full/workspace [port]
#
# The workspace is a checkout of this project in which the pipeline has been run, so that it holds the bulk
# outputs this repository does not ship: results/webapp_index/<snap>/, results/sogs/<snap>/{O1_r2,O2_r1}/web/,
# cases/*/cofold/models/ and data/annotations/uniprot_proteins.parquet. Structures are still streamed from the
# AFDB and atlas mounts (see webapp/README.md); nothing is copied.
set -euo pipefail
P="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"  # repo root
WS="$(cd "${1:?usage: $0 WORKSPACE [port]}" && pwd)"
PORT="${2:-8765}"
SNAP="${PE_WEBAPP_SNAPSHOT:-20261030T085136Z}"

export PE_WEBAPP_INDEX_DIRS="$WS/results/webapp_index/$SNAP:/tmp/pe_big/webapp_index/$SNAP"
export PE_WEBAPP_SOG_DIRS="$WS/results/sogs/$SNAP/O1_r2/web:$WS/results/sogs/$SNAP/O2_r1/web"
export PE_WEBAPP_CASE_DIRS="$WS/cases:$P/webapp/fixtures/cases"

# Gene-name search reads data/annotations/uniprot_proteins.parquet of this repo; link the workspace copy if needed.
ANN="$P/data/annotations/uniprot_proteins.parquet"
if [[ ! -e $ANN && -e $WS/data/annotations/uniprot_proteins.parquet ]]; then
  ln -s "$WS/data/annotations/uniprot_proteins.parquet" "$ANN"
fi

exec uv run --project "$P" python -m webapp --port "$PORT"
