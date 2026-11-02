#!/usr/bin/env bash
# Part B on an O1 SOG release: events (pre-rescue) -> rescue (k=1) -> events (post-rescue), tau 0.2 and 0.3.
#   scripts/partB/run_o1_pipeline.sh O1_r3
set -euo pipefail
REL=${1:-O1_r3}
S=20261030T085136Z
P=$(cd "$(dirname "$0")/../.." && pwd); cd "$P"
TH=results/partB/sim/measured/d11_primary_known_gt10/thresholds.json
BULK=/tmp/pe_big/sogs/$S/$REL
META=results/sogs/$S/$REL/sogs_all_tau.parquet
for t in 0.2 0.3; do
  OUT=results/partB/events/$S/$REL/tau$t/primary
  PE_WORKERS=16 .venv/bin/python scripts/partB/run_events.py --sogs results/sogs/$S/$REL/web --orthology O1 \
    --snapshot $S --states-glob "$BULK/states_t$t/shard=*.parquet" --state-col state --sog-meta $META --tau $t \
    --thresholds $TH --prefilter --stratum primary --out $OUT
  .venv/bin/python scripts/partB/run_rescue.py --tau $t --events $OUT --bulk $BULK --states-dir $BULK/states_t$t \
    --workers 16
  .venv/bin/python scripts/partB/run_events.py --sogs $OUT/rescue --orthology O1 --snapshot $S \
    --guide data/orthology/ancestra_tier1_v0/tier1_guide_tree.nwk --sog-meta $META --tau $t --thresholds $TH \
    --stratum primary --workers 16 --out $OUT/rescue/events
done
