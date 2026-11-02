#!/usr/bin/env bash
# Part B entry point (DESIGN §5): events -> rescue overlay -> events -> nulls -> candidates -> artifact cascade.
#
#   scripts/partB/run_partB.sh [--snapshot ID] [--orthology O1|O2] [--sog-run RUN] [--run-id LABEL]
#                              [--stages LIST] [--workers N] [--null-seeds N] [--results-root DIR] [--o1-run RUN]
#
# Defaults reproduce the committed results: snapshot = $PE_SNAPSHOT or config/prereg.yaml meta.snapshot_id;
# O1 -> sog-run O1_r3 (tau 0.2 primary, 0.3 sensitivity); O2 -> sog-run O2_r1 (tau 0.2; rescue overlays
# results/sogs/<snap>/<run>_rescue[2]/rescue_A.parquet, produced by scripts/rescue_sogs.py).
# Stages (comma list, default all):
#   O1: pre,rescue,post,nulls,pred11,controls,candidates,rbh,artifacts
#       (pred11 = events on the pre-D11 release --pre-d11-run, default O1_r2, both strata: source of the
#        bad-model false-loss controls)
#   O2: prep,events,nulls,agreement,report,rbh,artifacts
#       (agreement + report read the O1 run --o1-run of the SAME snapshot, default O1_r3: on a new snapshot run
#        them only after Part B O1 of that snapshot; report works without the agreement and says so)
#   (specificity: scripts/partB/run_specificity.py [--orthology O2 --release O2_r1 --run-id ID]; O2 needs tier-3
#    fingerprints from scripts/partB/extract_site_fp.py; rbh = DIAMOND reciprocal-best-hit orthology check)
# Cuts: results/partB/sim/measured/d11_primary_known_gt10/thresholds.json (run_sim.py on measured noise; the
# prereg hash is checked at load). Inputs: results/sogs/<snap>/<run>/ (web, sogs_all_tau.parquet) and the bulk
# /tmp/pe_big/sogs/<snap>/<run>/ (states_t*, msa, sogs).
set -euo pipefail
P=$(cd "$(dirname "$0")/../.." && pwd); cd "$P"
PY=.venv/bin/python
S=${PE_SNAPSHOT:-$($PY -c "import sys; sys.path.insert(0,'src'); from pocketevo import recon; print(recon.default_snapshot())")}
ORTH=O1; RUN=""; RUNID=""; STAGES=all; W=16; NSEEDS=20; RR=results/partB; PRE=O1_r2; O1RUN=O1_r3
while [ $# -gt 0 ]; do
  case $1 in
    --snapshot) S=$2; shift 2;; --orthology) ORTH=$2; shift 2;; --sog-run) RUN=$2; shift 2;;
    --run-id) RUNID=$2; shift 2;; --stages) STAGES=$2; shift 2;; --workers) W=$2; shift 2;;
    --null-seeds) NSEEDS=$2; shift 2;; --results-root) RR=$2; shift 2;; --pre-d11-run) PRE=$2; shift 2;;
    --o1-run) O1RUN=$2; shift 2;;
    -h|--help) sed -n 2,19p "$0"; exit 0;; *) echo "unknown option $1" >&2; exit 2;;
  esac
done
[ -z "$RUN" ] && { [ "$ORTH" = O1 ] && RUN=O1_r3 || RUN=O2_r1; }
RUNID=${RUNID:-$RUN}
export PE_SNAPSHOT=$S PE_WORKERS=$W
TH=results/partB/sim/measured/d11_primary_known_gt10/thresholds.json
WEB=results/sogs/$S/$RUN/web; META=results/sogs/$S/$RUN/sogs_all_tau.parquet; BULK=/tmp/pe_big/sogs/$S/$RUN
COMMON="--snapshot $S --sog-run $RUN --run-id $RUNID"
want() { [ "$STAGES" = all ] || [[ ",$STAGES," == *",$1,"* ]]; }
log() { echo "[$(date '+%H:%M:%S')] $*"; }

if [ "$ORTH" = O1 ]; then
  OX=data/orthology/ancestra_tier1_v0/accurate/hogs.orthoxml.gz; GT=data/orthology/ancestra_tier1_v0/tier1_guide_tree.nwk
  for t in 0.2 0.3; do
    OUT=$RR/events/$S/$RUN/tau$t/primary
    if want pre; then log "events (pre-rescue) tau $t"
      $PY scripts/partB/run_events.py --sogs $WEB --orthology O1 --orthoxml $OX $COMMON --states-glob "$BULK/states_t$t/shard=*.parquet" \
        --state-col state --sog-meta $META --tau $t --thresholds $TH --prefilter --stratum primary --workers $W --out $OUT; fi
    if want rescue; then log "rescue (k=1) tau $t"
      $PY scripts/partB/run_rescue.py $COMMON --tau $t --events $OUT --bulk $BULK --states-dir $BULK/states_t$t --workers $W; fi
    if want post; then log "events (post-rescue) tau $t"
      $PY scripts/partB/run_events.py --sogs $OUT/rescue --orthology O1 --orthoxml $OX $COMMON --guide $GT --sog-meta $META \
        --tau $t --thresholds $TH --stratum primary --workers $W --out $OUT/rescue/events; fi
  done
  if want nulls; then log "null-input funnel ($NSEEDS seeds)"
    $PY scripts/partB/run_nulls.py --sogs $WEB --orthology O1 --orthoxml $OX $COMMON --states-glob "$BULK/states_t0.2/shard=*.parquet" \
      --sog-meta $META --tau 0.2 --thresholds $TH --prefilter --skip-m3 --n-seeds $NSEEDS --workers $W \
      --out $RR/nulls/$S/$RUN/tau0.2/primary; fi
  if want pred11; then for sc in state state_s70; do log "pre-D11 events ($PRE, $sc)"
    $PY scripts/partB/run_events.py --sogs results/sogs/$S/$PRE/web --orthology O1 --orthoxml $OX --snapshot $S --sog-run $PRE \
      --states-glob "/tmp/pe_big/sogs/$S/$PRE/states_t0.2/shard=*.parquet" --state-col $sc --sog-meta results/sogs/$S/$PRE/sogs_all_tau.parquet \
      --tau 0.2 --thresholds $TH --prefilter --stratum pre_d11_$sc --workers $W --out $RR/events/$S/${PRE}_preD11/tau0.2/$sc; done; fi
  if want controls; then log "false-loss controls"
    $PY scripts/partB/false_loss_controls.py $COMMON --results-root $RR --pre-d11-run $PRE; fi
  if want candidates; then log "candidates"
    $PY scripts/partB/write_candidates.py $COMMON --results-root $RR --thresholds $TH; fi
  if want rbh; then log "RBH orthology check"; $PY scripts/partB/rbh_check.py --snapshot $S --runs $RUN:O1 --threads $W; fi
  if want artifacts; then log "artifact cascade"
    $PY scripts/partB/run_artifacts.py $COMMON --results-root $RR; fi
else
  OX=data/orthology/ancestra_O2/accurate/hogs.orthoxml.gz; GT=data/orthology/ancestra_O2/o2_guide_tree.nwk
  ST=/tmp/pe_big/partB/$S/$RUN/states_resc_t0.2
  if want prep; then log "O2 states + rescue overlays"; $PY scripts/partB/prep_o2_states.py $COMMON; fi
  OUT=$RR/events/$S/$RUN/tau0.2/primary
  if want events; then log "events O2"
    $PY scripts/partB/run_events.py --sogs $WEB --orthology O2 --orthoxml $OX --guide $GT $COMMON --states-glob "$ST/shard=*.parquet" \
      --sog-meta $META --tau 0.2 --thresholds $TH --prefilter --stratum primary --workers $W --out $OUT; fi
  if want nulls; then log "null-input funnel O2 ($NSEEDS seeds)"
    $PY scripts/partB/run_nulls.py --sogs $WEB --orthology O2 --orthoxml $OX --guide $GT $COMMON --states-glob "$ST/shard=*.parquet" \
      --sog-meta $META --tau 0.2 --thresholds $TH --prefilter --skip-m3 --n-seeds $NSEEDS --workers $W \
      --out $RR/nulls/$S/$RUN/tau0.2/primary; fi
  if want agreement; then log "O1 -> O2 gain agreement"; $PY scripts/partB/o1_o2_agreement.py --snapshot $S --sog-run $O1RUN --o2-run $RUN; fi
  if want report; then log "O2 report"; $PY scripts/partB/write_o2_report.py $COMMON --o1-run $O1RUN; fi
  if want rbh; then log "RBH orthology check"; $PY scripts/partB/rbh_check.py --snapshot $S --runs $RUN:O2 --threads $W; fi
  if want artifacts; then log "artifact cascade O2"; $PY scripts/partB/run_artifacts.py $COMMON --orthology O2 --results-root $RR; fi
fi
log done
