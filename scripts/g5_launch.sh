#!/usr/bin/env bash
# G5 (PI decision D27): wait until --at, take a new atlas snapshot, run the D27 stages in two lanes through
# scripts/run_pipeline.py, then compare snapshot 2 with snapshot 1 (scripts/g5_compare.py) and commit SUMMARIES only.
#
#   nohup setsid bash scripts/g5_launch.sh --at "2026-10-31 03:00" > /tmp/pe_big/logs/g5/launcher.log 2>&1 &
#   bash scripts/g5_launch.sh --now --snapshot 20261030T085136Z --run-id dry1 --dry-run --no-commit   # test
#   G5_FORCE_QUOTA_FAIL=1 bash scripts/g5_launch.sh --commit-only --snapshot <S>    # test the abort path
#   (commit-path tests only in a throwaway clone/branch under /tmp, never on the real branch)
#
# Storage (PI storage alert, home quota nearly full): ALL G5 outputs live under /tmp/pe_big.
#   - snapshot manifest   -> /tmp/pe_big/snapshot/<S>/                 (symlink data/snapshot/<S>)
#   - results/<area>/<S>  -> /tmp/pe_big/g5_results/<area>/<S>/        (symlinks, so every entry point runs unchanged)
#   - logs + live status  -> /tmp/pe_big/logs/g5/<S>/ (G5_STATUS.md refreshed after every stage)
#   - committed at the end: results/pipeline/G5_STATUS.md + results/pipeline/G5/<S>/ (small summary JSON/CSV files),
#     only after a write probe on the home file system succeeds (3x the commit size + 50 MB); otherwise the commit is
#     skipped cleanly and the reason logged to /tmp/pe_big/logs/g5/<S>/commit.log.
# Lanes (D27):
#   shared : snapshot -> sites (tiers 1-3 + L) -> web index
#   lane A : sogs_O2 -> bench_O2 -> [wait for sogs_O1] -> split_O2 -> rescue_O2 (BURROW outputs reused) -> partB_O2
#            (O2-only Part B stages; the O1 -> O2 agreement is NOT in this lane)
#   lane B : sogs_O1 -> bench_O1 -> partB_O1 -> partA (GPU 3; reproducibility) -> cases (needs Part B O1 + Part A)
#   end    : partB_O1O2 (O1 -> O2 agreement + O2 report, needs both lanes) -> og_index -> comparison -> commit
# Frozen / species-level inputs (rescue k D9, TimeTree, observation calibrations) come from pocketevo.frozen, never
# from results/<area>/<new snapshot>; `run_pipeline.py --dry-run` checks every input (tests/test_run_pipeline_dryrun.py).
# Part A writes results/partA/<S>/<run>/numbers_partA.json only (never results/numbers.json["partA"]).
# Storage guard: a stage is skipped (logged as BLOCKED) if du(/tmp/pe_big) + its estimate > --limit-gb (450), or if a
# 20 MB write probe on the home file system fails.
set -uo pipefail
{   # whole script in one compound command: bash parses it fully before running, so edits to this file
    # while the launcher waits cannot change the running instance
AT="2026-10-31 03:00"; RID=g5; LIMIT_GB=450; NOW=0; SNAP=""; DRY=""; COMMIT=1; COMMIT_ONLY=0
while [ $# -gt 0 ]; do
  case "$1" in
    --at) AT="$2"; shift 2;;
    --run-id) RID="$2"; shift 2;;
    --limit-gb) LIMIT_GB="$2"; shift 2;;
    --now) NOW=1; shift;;
    --snapshot) SNAP="$2"; shift 2;;
    --dry-run) DRY="--dry-run"; shift;;
    --no-commit) COMMIT=0; shift;;
    --commit-only) COMMIT_ONLY=1; NOW=1; shift;;   # test: only the summary commit (needs --snapshot)
    *) echo "unknown option $1"; exit 2;;
  esac
done
P=$(cd "$(dirname "$0")/.." && pwd); cd "$P"
PY="$P/.venv/bin/python"
export UV_CACHE_DIR=/tmp/pe_big/uvcache
S1=20261030T085136Z
LOGROOT=/tmp/pe_big/logs/g5
SNAPROOT=/tmp/pe_big/snapshot
RESROOT=/tmp/pe_big/g5_results
mkdir -p "$LOGROOT" "$SNAPROOT" "$RESROOT"
if [ "$NOW" = 0 ]; then
  echo $$ > "$LOGROOT/launcher.pid"   # only the scheduled instance owns the pid file
  target=$(date -d "$AT" +%s)
  echo "[$(date)] waiting until $AT (epoch $target), pid $$"
  while [ "$(date +%s)" -lt "$target" ]; do sleep 30; done
fi
echo "[$(date)] start"

home_ok() {   # home_ok <MB>: write probe on the home file system (the quota is not queryable here)
  local f="$P/.g5_quota_probe.$$"
  if dd if=/dev/zero of="$f" bs=1M count="$1" conv=fsync status=none 2>/dev/null; then rm -f "$f"; return 0; fi
  rm -f "$f"; return 1
}

if [ -z "$SNAP" ]; then
  mkdir -p "$LOGROOT/_snapshot"
  home_ok 20 || { echo "[$(date)] home write probe failed before the snapshot: abort"; exit 1; }
  t0=$(date +%s)
  $PY scripts/atlas_snapshot.py --atlas-root /cluster/scratch/gmagnusson/pocketcompare_afdb --out "$SNAPROOT" \
      --threads 32 > "$LOGROOT/_snapshot/snapshot.log" 2>&1
  rc=$?
  SNAP=$(ls -1 "$SNAPROOT" | sort | tail -1)
  if [ $rc -ne 0 ] || [ -z "$SNAP" ] || [ "$SNAP" = "$S1" ]; then echo "[$(date)] snapshot failed (rc $rc)"; exit 1; fi
  ln -sfn "$SNAPROOT/$SNAP" "data/snapshot/$SNAP"
  mkdir -p "$LOGROOT/$SNAP"; mv "$LOGROOT/_snapshot/snapshot.log" "$LOGROOT/$SNAP/"
  echo "{\"stage\": \"snapshot\", \"rc\": 0, \"wall_s\": $(( $(date +%s) - t0 ))}" >> "$LOGROOT/$SNAP/shared.jsonl"
fi
LOGD="$LOGROOT/$SNAP"; mkdir -p "$LOGD"
# results/<area>/<S> -> /tmp (only for a new snapshot; the committed snapshot's dirs are never touched)
if [ "$SNAP" != "$S1" ]; then
  for d in sogs partA partB partB/events partB/nulls partB/m3 cases webapp_index; do
    mkdir -p "$RESROOT/$d/$SNAP" "results/$d"
    [ -e "results/$d/$SNAP" ] || ln -s "$RESROOT/$d/$SNAP" "results/$d/$SNAP"
    grep -qxF "/results/$d/$SNAP" .git/info/exclude 2>/dev/null || echo "/results/$d/$SNAP" >> .git/info/exclude
  done
  grep -qxF "/data/snapshot/$SNAP" .git/info/exclude 2>/dev/null || echo "/data/snapshot/$SNAP" >> .git/info/exclude
fi
echo "[$(date)] snapshot $SNAP, logs $LOGD, results under $RESROOT"

status() { $PY scripts/g5_compare.py --s2 "$SNAP" --run-id "$RID" --status-only --out "$LOGD/G5_STATUS.md" > /dev/null 2>&1; }

guard() {   # guard <need_GB> <lane> <stage>
  local used; used=$(du -s --block-size=1G /tmp/pe_big 2>/dev/null | cut -f1)
  if [ $(( used + $1 )) -gt "$LIMIT_GB" ]; then
    echo "{\"stage\": \"$3\", \"rc\": null, \"blocked\": \"storage guard: /tmp/pe_big ${used} GB + ${1} GB > ${LIMIT_GB} GB\", \"wall_s\": 0}" >> "$LOGD/$2.jsonl"
    return 1
  fi
  if ! home_ok 20; then
    echo "{\"stage\": \"$3\", \"rc\": null, \"blocked\": \"home quota: 20 MB write probe failed\", \"wall_s\": 0}" >> "$LOGD/$2.jsonl"
    return 1
  fi
}

run() {     # run <lane> <stage> <need_GB> <workers>
  local lane=$1 stage=$2 need=$3 w=$4
  guard "$need" "$lane" "$stage" || { status; return 1; }
  echo "[$(date)] $lane: $stage"
  $PY scripts/run_pipeline.py --snapshot "$SNAP" --run-id "$RID" --stages "$stage" --workers "$w" \
      --snapshot-root "$P/data/snapshot" --msa-seed /tmp/pe_big/sogs --rescue-work /tmp/pe_big/rescue/$S1/O2_r1 \
      --log "$LOGD/$lane.jsonl" $DRY >> "$LOGD/$lane.log" 2>&1
  local rc=$?
  status
  echo "[$(date)] $lane: $stage rc=$rc"
  return $rc
}

commit_summaries() {
  local C="$LOGD/commit.log" stage="$LOGD/commit_stage" dst="results/pipeline/G5/$SNAP"
  rm -rf "$stage"; mkdir -p "$stage"
  cp "$LOGD/G5_STATUS.md" "$stage/G5_STATUS.md" 2>/dev/null
  cp "$SNAPROOT/$SNAP/manifest.json" "$stage/snapshot_manifest.json" 2>/dev/null
  for f in "sogs/$SNAP/O1_$RID/funnel.json" "sogs/$SNAP/O2_$RID/funnel.json" "sogs/$SNAP/O1_$RID/tau_benchmark.csv" \
           "sogs/$SNAP/O2_$RID/tau_benchmark.csv" "sogs/$SNAP/O2_${RID}_split/split_o2_summary.json" \
           "sogs/$SNAP/O2_${RID}_rescue/rescue_summary.json" "partA/$SNAP/$RID/numbers_partA.json" \
           "partB/$SNAP/O1_$RID/artifact_summary.json" "partB/$SNAP/O2_$RID/artifact_summary.json" \
           "partB/$SNAP/O1_$RID/o1_o2_gain_agreement.json"; do
    [ -f "$RESROOT/$f" ] && cp "$RESROOT/$f" "$stage/$(echo "$f" | cut -d/ -f3)__$(basename "$f")"
  done
  local kb; kb=$(du -sk "$stage" | cut -f1)
  local need_mb=$(( (kb * 3) / 1024 + 50 ))
  echo "[$(date)] commit: $kb KB of summaries; home write probe ${need_mb} MB" >> "$C"
  if [ "${G5_FORCE_QUOTA_FAIL:-0}" = 1 ] || ! home_ok "$need_mb"; then
    echo "[$(date)] ABORT commit: home write probe of ${need_mb} MB failed (quota); summaries stay in $stage" >> "$C"
    return 1
  fi
  mkdir -p "$dst" && cp "$stage"/* "$dst"/ && cp "$stage/G5_STATUS.md" results/pipeline/G5_STATUS.md || {
    echo "[$(date)] ABORT commit: copy into the repo failed" >> "$C"; return 1; }
  git add -- results/pipeline/G5_STATUS.md "$dst" >> "$C" 2>&1
  local bad; bad=$(git diff --cached --name-only -- results/pipeline/G5_STATUS.md "$dst" | grep -v '^results/pipeline/G5' || true)
  if [ -n "$bad" ]; then
    echo "[$(date)] ABORT commit: staged paths outside results/pipeline/G5*: $bad" >> "$C"
    git reset -q -- results/pipeline/G5_STATUS.md "$dst"; return 1
  fi
  if \
     git commit -q -m "G5 ($SNAP): rerun summaries (D27), snapshot 2 vs 1

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>" -- results/pipeline/G5_STATUS.md "$dst" >> "$C" 2>&1; then
    echo "[$(date)] committed $(git rev-parse HEAD)" >> "$C"
  else
    echo "[$(date)] ABORT commit: git add/commit failed (see above); working tree left as is" >> "$C"
    return 1
  fi
}

if [ "$COMMIT_ONLY" = 1 ]; then
  status
  commit_summaries && echo "[$(date)] summaries committed" || echo "[$(date)] commit skipped (see $LOGD/commit.log)"
  exit 0
fi

# shared stages (estimates: tier-3 site metadata ~0.71 B sites x 86 B = 61 GB; web index ~5 GB)
run shared sites 70 32 && run shared webindex 6 32 || { echo "shared stages failed"; status; [ -n "$DRY" ] || exit 1; }

lane_a() {
  run laneA sogs_O2 10 32 && run laneA bench_O2 1 8 || return 1
  while [ ! -e "$LOGD/O1_sogs.done" ]; do sleep 30; done
  run laneA split_O2 1 8
  run laneA rescue_O2 25 32 && run laneA partB_O2 6 16
}
lane_b() {
  run laneB sogs_O1 5 24; touch "$LOGD/O1_sogs.done"
  run laneB bench_O1 1 8
  run laneB partB_O1 6 16; local pb=$?
  run laneB partA 45 16 && [ $pb = 0 ] && run laneB cases 1 8
}
lane_a > "$LOGD/laneA.driver.log" 2>&1 &
PA=$!
lane_b > "$LOGD/laneB.driver.log" 2>&1 &
PB=$!
wait $PA; wait $PB
run final partB_O1O2 1 16
run final og_index 1 16
$PY scripts/g5_compare.py --s2 "$SNAP" --run-id "$RID" --note "Run finished $(date)." --out "$LOGD/G5_STATUS.md" \
  && echo "[$(date)] status: $LOGD/G5_STATUS.md"
if [ "$COMMIT" = 1 ] && [ -z "$DRY" ]; then
  commit_summaries && echo "[$(date)] summaries committed" || echo "[$(date)] commit skipped (see $LOGD/commit.log)"
fi
exit 0
}
