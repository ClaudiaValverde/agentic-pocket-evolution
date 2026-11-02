#!/usr/bin/env bash
# Run Ancestra O2 (tier-1 + context proteomes listed in $OUT/o2_proteomes.tsv, written by
# select_context_species.py). Work dir in /tmp (home quota); compact outputs copied to $OUT.
# Usage: scripts/orthology/run_ancestra_o2.sh [threads=40] [preset=accurate]
set -euo pipefail
P="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"  # repo root
THREADS=${1:-40}
PRESET=${2:-accurate}
OUT=$P/data/orthology/ancestra_O2
W=/tmp/pe_big/ancestra/o2
RUN=$W/run_O2_$PRESET
IN=$W/input_O2
TAXDUMP=/tmp/pe_big/ancestra/taxdump
export PATH=$P/envs/ancestra/bin:$PATH NUMBA_CACHE_DIR=/tmp/pe_big/ancestra/numba TMPDIR=/tmp/pe_big/ancestra/tmp
mkdir -p "$IN" "$TMPDIR" "$OUT"

# 1. input dir = symlinks named <UP>_<taxid>.fasta.gz (file stem = species id, --leaf-naming qfo)
rm -f "$IN"/*.fasta.gz
tail -n +2 "$OUT/o2_proteomes.tsv" | while IFS=$'\t' read -r set up taxid species file md5 n; do
  ln -s "$file" "$IN/${up}_${taxid}.fasta.gz"
done
echo "input proteomes: $(ls "$IN" | wc -l)"

# 2. guide tree from the NCBI taxdump
python -m ancestra.speciestree ncbi --proteomes "$IN"/*.fasta.gz --taxdump "$TAXDUMP" --leaf-naming qfo \
    --out "$W/o2_guide_tree.nwk"
cp "$W/o2_guide_tree.nwk" "$OUT/o2_guide_tree.nwk"

# 3. run (stage checkpoints in $RUN/stages; reruns skip unchanged stages)
start=$(date +%s)
/usr/bin/time -v python -m ancestra run --proteomes "$IN" --tree "$W/o2_guide_tree.nwk" --out "$RUN" \
    --threads "$THREADS" --preset "$PRESET" > "$RUN.log" 2> "$RUN.time" || { tail -30 "$RUN.time"; exit 1; }
end=$(date +%s)
echo "wall_s=$((end-start))" | tee "$RUN.wall"

# 4. verify + copy compact outputs
python $P/external/ancestra-mvp-ai/scripts/verify_orthoxml.py "$RUN" > "$RUN.verify" 2>&1 || true
mkdir -p "$OUT/$PRESET"
gzip -c "$RUN/hogs.orthoxml" > "$OUT/$PRESET/hogs.orthoxml.gz"
cp "$RUN/manifest.json" "$OUT/$PRESET/manifest.json"
cp "$RUN.log" "$OUT/$PRESET/run.log"; cp "$RUN.time" "$OUT/$PRESET/time.txt"; cp "$RUN.wall" "$OUT/$PRESET/wall.txt"
cp "$RUN.verify" "$OUT/$PRESET/verify.txt"
echo done
