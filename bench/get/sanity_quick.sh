#!/usr/bin/env bash
# Quick zebrac sanity for GET RSS/speed checks during development.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
source "$ROOT/bench/shared/tools.sh"
source "$ROOT/bench/shared/zebrac_runner.sh"

export ZEBRAC_DURATION_MS="${ZEBRAC_DURATION_MS:-800}"
export ZEBRAC_MIN_SAMPLES="${ZEBRAC_MIN_SAMPLES:-3}"
export ZEBRAC_WARMUP="${ZEBRAC_WARMUP:-1}"

DATA="$ROOT/bench/shared/data"
BED="$ROOT/bench/get/data/bed"
REG="$ROOT/bench/get/data/regions"
OUT="$ROOT/bench/get/results/sanity"
META="$OUT/meta.jsonl"
mkdir -p "$OUT"

print_result() {
    local tool="$1"
    local raw="$2"
    python3 - "$tool" "$raw" <<'PY'
import json, sys
tool, path = sys.argv[1], sys.argv[2]
r = json.load(open(path))["results"][0]
wall = r["wall_time"]["mean"] / 1e6
rss = r["peak_rss"]["mean"] / (1024 * 1024)
print(f"{tool:10} wall={wall:7.2f}ms rss={rss:7.1f}MB")
PY
}

run_one() {
    local tag="$1"
    local tool="$2"
    local script="$3"
    local raw="$OUT/${tag}_${tool}.json"
    zebrac_clear_commands
    zebrac_add_command get sanity "$tag" "$tool" "$tool" "" "" "$raw" "$(shell_command "$script")"
    zebrac_run_current_group "$raw" "$META"
    print_result "$tool" "$raw"
}

REG1="$(cat "$REG/Transcriptome_1kbp_mid.txt")"

echo "=== Transcriptome BED 100 ==="
run_one bed100 z-fasta "$ZFASTA get $DATA/REAL_Transcriptome.fa --bed $BED/Transcriptome_100.bed > /dev/null"
run_one bed100 noodles "$NOODLES get $DATA/REAL_Transcriptome.fa --bed $BED/Transcriptome_100.bed > /dev/null"
run_one bed100 rust-bio "$RUSTBIO get $DATA/REAL_Transcriptome.fa --bed $BED/Transcriptome_100.bed > /dev/null"

echo "=== Transcriptome positional 1kbp ==="
run_one pos1k z-fasta "$ZFASTA get $DATA/REAL_Transcriptome.fa \"$REG1\" > /dev/null"
run_one pos1k noodles "$NOODLES get $DATA/REAL_Transcriptome.fa \"$REG1\" > /dev/null"
run_one pos1k rust-bio "$RUSTBIO get $DATA/REAL_Transcriptome.fa \"$REG1\" > /dev/null"

echo "=== Transcriptome multi N=10 (regions file) ==="
MULTI="$REG/Transcriptome_regions_N10.txt"
if [[ -f "$MULTI" ]]; then
    run_one multi10 z-fasta "$ZFASTA get $DATA/REAL_Transcriptome.fa \$(tr '\n' ' ' < $MULTI) > /dev/null"
    run_one multi10 noodles "$NOODLES get $DATA/REAL_Transcriptome.fa \$(tr '\n' ' ' < $MULTI) > /dev/null"
fi

echo "=== Correctness: BED 10 rows vs samtools ==="
$ZFASTA get "$DATA/REAL_Transcriptome.fa" --bed "$BED/Transcriptome_10.bed" > "$OUT/zf_bed10.fa"
samtools faidx "$DATA/REAL_Transcriptome.fa" -r "$BED/Transcriptome_10.regions.txt" > "$OUT/st_bed10.fa"
diff -q "$OUT/zf_bed10.fa" "$OUT/st_bed10.fa"
echo "byte-identical OK"
