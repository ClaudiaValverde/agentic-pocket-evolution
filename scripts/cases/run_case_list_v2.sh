#!/usr/bin/env bash
# D28 co-fold pool: stage O1 and O2 case lists (filters 1-3 + literature line), then combine, dedup and select.
#   [SNAPSHOT=<id>] scripts/cases/run_case_list_v2.sh [O2_CANDIDATES]
set -euo pipefail
S=${SNAPSHOT:-20261030T085136Z}
P=$(cd "$(dirname "$0")/../.." && pwd); cd "$P"
O1C=results/partB/$S/O1_r3
O2C=${1:-results/partB/$S/O2_r1/candidates.parquet}
O2FL=$(dirname "$O2C")/false_loss_controls.parquet
STG=/tmp/pe_big/cases/stage
OUT=results/cases/$S
rm -rf "$STG"; mkdir -p "$STG/O1" "$STG/O2"
SPEC1=$O1C/specificity_candidates.parquet
.venv/bin/python scripts/cases/build_cases.py --snapshot $S --orthology O1 --sog-release O1_r3 --pool eligible --version stage \
  --candidates $O1C/candidates.parquet --false-loss $O1C/false_loss_controls.parquet --false-loss-kind bad_model \
  --n-controls 5 --events-base results/partB/events/$S/O1_r3/tau0.2/primary --specificity "$SPEC1" \
  --out "$STG/O1" --cases-dir "$STG/O1" --verify-ids
.venv/bin/python scripts/cases/build_cases.py --snapshot $S --orthology O2 --sog-release O2_r1 --pool eligible --version stage \
  --candidates "$O2C" --false-loss "$O2FL" --false-loss-kind bad_model --n-controls 5 --convergence "" \
  --specificity "$(dirname "$O2C")/specificity_candidates.parquet" \
  --dup-gains "$(dirname "$O2C")/duplication_gains_exploratory.parquet" --n-dup-gains 2 \
  --events-base results/partB/events/$S/O2_r1/tau0.2/primary \
  --orthoxml data/orthology/ancestra_O2/accurate/hogs.orthoxml.gz \
  --states-glob "/tmp/pe_big/partB/$S/O2_r1/states_resc_t0.2/shard=*.parquet" \
  --rescue-hits results/sogs/$S/O2_r1_rescue/rescue_A.parquet --out "$STG/O2" --cases-dir "$STG/O2" --verify-ids
mkdir -p $OUT/stage
cp "$STG/O1/case_list_stage.parquet" $OUT/stage/case_list_stage_O1.parquet
cp "$STG/O1/case_list_stage.parquet.json" $OUT/stage/case_list_stage_O1.parquet.json
cp "$STG/O2/case_list_stage_O2.parquet" $OUT/stage/case_list_stage_O2.parquet
cp "$STG/O2/case_list_stage_O2.parquet.json" $OUT/stage/case_list_stage_O2.parquet.json
.venv/bin/python scripts/cases/combine_pool.py --snapshot $S --lists $OUT/stage/case_list_stage_O1.parquet \
  $OUT/stage/case_list_stage_O2.parquet --stage-dirs "$STG/O1" "$STG/O2" --out $OUT/case_list_v2.parquet \
  --candidates O1=$O1C/candidates.parquet O2=$O2C
