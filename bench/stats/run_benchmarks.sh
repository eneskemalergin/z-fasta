#!/bin/bash
# z-fasta STATS Benchmark Runner
# Measures statistics computation performance: runtime, memory, throughput.
# Uses hyperfine for precise timing.
#
# Usage: ./run_benchmarks.sh [--runs N] [--skip-scaling] [--skip-real]
#
# Outputs:
#   results/stats_<timestamp>/         — per-file benchmarks
#   results/indexonly_<timestamp>/      — --index-only mode benchmarks
#   results/scaling_<timestamp>/       — file-size scaling
#   results/memory_<timestamp>.csv     — memory (RSS) per tool/mode
#
# Comparison tools:
#   - seqkit stats (Go bioinformatics toolkit)
#   - samtools (fasta stats not directly comparable, but index timing is)
#
# Note: samtools and fastahack do not have a "stats" equivalent.
#       seqkit stats is the primary comparison tool.
#       We also compare z-fasta --index-only vs full scan internally.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH_ROOT="$(dirname "$SCRIPT_DIR")"
PROJECT_ROOT="$(dirname "$BENCH_ROOT")"
RESULTS_DIR="$SCRIPT_DIR/results"
DATA_DIR="$BENCH_ROOT/shared/data"
INDEX_DATA="$BENCH_ROOT/index/data"

# ── Tools ──────────────────────────────────────────────────────────
ZFASTA="$PROJECT_ROOT/zig-out/bin/z-fasta"
SEQKIT="$PROJECT_ROOT/tools/seqkit"
SEQTK="$PROJECT_ROOT/tools/seqtk/seqtk"

# ── Defaults ───────────────────────────────────────────────────────
RUNS=5
WARMUP=1
SKIP_SCALING=false
SKIP_REAL=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --runs)          RUNS="$2"; shift 2 ;;
        --warmup)        WARMUP="$2"; shift 2 ;;
        --skip-scaling)  SKIP_SCALING=true; shift ;;
        --skip-real)     SKIP_REAL=true; shift ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

# ── Preflight ──────────────────────────────────────────────────────
mkdir -p "$RESULTS_DIR"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

command -v hyperfine &>/dev/null || { echo "Error: hyperfine not found (apt install hyperfine)"; exit 1; }
[[ -x "$ZFASTA" ]]   || { echo "Error: z-fasta not found at $ZFASTA. Run: ./zig-0.16.0/zig build -Doptimize=ReleaseFast"; exit 1; }

HAS_SEQKIT=false; [[ -x "$SEQKIT" ]] && HAS_SEQKIT=true
HAS_SEQTK=false;  [[ -x "$SEQTK" ]]  && HAS_SEQTK=true

# Cache-clearing command
CACHE_CLEAR=""
BENCH_MODE="warm"
if sudo -n sh -c 'sync; echo 3 > /proc/sys/vm/drop_caches' 2>/dev/null; then
    CACHE_CLEAR="sudo sh -c 'sync; echo 3 > /proc/sys/vm/drop_caches'"
    WARMUP=0
    BENCH_MODE="cold (cache cleared before each run)"
    echo "Cache clearing: enabled → forcing --warmup 0"
else
    echo "Cache clearing: DISABLED → warm-cache mode"
fi

echo ""
echo "════════════════════════════════════════════════════════════════"
echo "  z-fasta STATS Benchmark Runner"
echo "════════════════════════════════════════════════════════════════"
echo "  Mode: $BENCH_MODE"
echo "  Runs: $RUNS | Warmup: $WARMUP"
echo "  z-fasta: $ZFASTA"
echo "  seqkit:  $HAS_SEQKIT ($SEQKIT)"
echo "  seqtk:   $HAS_SEQTK ($SEQTK)"
echo ""
echo "  NOTE: samtools/fastahack have no stats command."
echo "  Composition  : z-fasta stats vs seqtk comp"
echo "  Assembly stats: z-fasta stats vs seqkit stats -a (N50, GC%, etc.)"
echo "  Index-only   : z-fasta stats --index-only (no equivalent in seqtk/seqkit)"
echo ""

# ── Helper: ensure index exists ───────────────────────────────────
ensure_indexes() {
    local file="$1"
    [[ -f "${file}.fai" ]] || samtools faidx "$file" 2>/dev/null || true
    [[ -f "${file}.zfi" ]] || "$ZFASTA" index "$file" 2>/dev/null || true
}

# ── Helper: benchmark stats on a single file ─────────────────────
bench_stats() {
    local file="$1"
    local json_out="$2"
    local label_prefix="${3:-}"

    local prepare_cmd=""
    [[ -n "$CACHE_CLEAR" ]] && prepare_cmd="$CACHE_CLEAR"

    local cmds=()
    local names=()

    # z-fasta stats (full scan)
    names+=("${label_prefix}z-fasta-full")
    cmds+=("$ZFASTA stats '$file' > /dev/null")

    # z-fasta stats --index-only
    names+=("${label_prefix}z-fasta-indexonly")
    cmds+=("$ZFASTA stats --index-only '$file' > /dev/null")

    # seqkit stats -a (assembly stats: N50, GC%, total N, etc.)
    if $HAS_SEQKIT; then
        names+=("${label_prefix}seqkit-stats-a")
        cmds+=("$SEQKIT stats -a '$file' > /dev/null 2>&1")
    fi

    # seqtk comp (per-sequence base composition; full scan)
    if $HAS_SEQTK; then
        names+=("${label_prefix}seqtk-comp")
        cmds+=("$SEQTK comp '$file' > /dev/null 2>&1")
    fi

    local hf_args=(
        hyperfine
        --warmup "$WARMUP"
        --runs "$RUNS"
        --export-json "$json_out"
    )
    [[ -n "$prepare_cmd" ]] && hf_args+=(--prepare "$prepare_cmd")

    for i in "${!names[@]}"; do
        hf_args+=(-n "${names[$i]}" "${cmds[$i]}")
    done

    "${hf_args[@]}"
}

# ── Helper: measure memory ────────────────────────────────────────
measure_memory_stats() {
    local file="$1"
    local mode="$2"  # "full" or "indexonly"

    local tools_and_cmds=()
    if [[ "$mode" == "full" ]]; then
        tools_and_cmds+=("z-fasta-full:$ZFASTA stats '$file' > /dev/null")
    else
        tools_and_cmds+=("z-fasta-indexonly:$ZFASTA stats --index-only '$file' > /dev/null")
    fi

    if $HAS_SEQKIT && [[ "$mode" == "full" ]]; then
        tools_and_cmds+=("seqkit-stats-a:$SEQKIT stats -a '$file' > /dev/null 2>&1")
    fi
    if $HAS_SEQTK && [[ "$mode" == "full" ]]; then
        tools_and_cmds+=("seqtk-comp:$SEQTK comp '$file' > /dev/null 2>&1")
    fi

    for entry in "${tools_and_cmds[@]}"; do
        local tool="${entry%%:*}"
        local cmd="${entry#*:}"

        local tmp_time
        tmp_time=$(mktemp)
        /usr/bin/time -f "%e %M %F %R" -o "$tmp_time" bash -c "$cmd" 2>/dev/null || true

        local time_s mem_kb major_faults minor_faults
        time_s=$(awk '{print $1}' "$tmp_time")
        mem_kb=$(awk '{print $2}' "$tmp_time")
        major_faults=$(awk '{print $3}' "$tmp_time")
        minor_faults=$(awk '{print $4}' "$tmp_time")
        rm -f "$tmp_time"

        echo "$tool,$time_s,$mem_kb,$major_faults,$minor_faults"
    done
}

# ══════════════════════════════════════════════════════════════════════
#  1. Test File Stats
# ══════════════════════════════════════════════════════════════════════

echo "──────────────────────────────────────────────────"
echo " [1] Test File Stats Performance"
echo "──────────────────────────────────────────────────"

STATS_DIR="$RESULTS_DIR/stats_${TIMESTAMP}"
mkdir -p "$STATS_DIR"

for f in "$PROJECT_ROOT"/tests/data/simple.fasta \
         "$PROJECT_ROOT"/tests/data/proteome.fasta \
         "$PROJECT_ROOT"/tests/data/edge_cases.fasta \
         "$PROJECT_ROOT"/tests/data/mixed_widths.fasta; do
    if [[ -f "$f" ]]; then
        ensure_indexes "$f"
        fname=$(basename "$f" .fasta)
        echo "  $fname..."
        bench_stats "$f" "$STATS_DIR/${fname}.json" "${fname}_"
    fi
done
echo ""

# ══════════════════════════════════════════════════════════════════════
#  2. Index-Only vs Full Scan Comparison
# ══════════════════════════════════════════════════════════════════════

echo "──────────────────────────────────────────────────"
echo " [2] Index-Only vs Full Scan"
echo "──────────────────────────────────────────────────"

INDEXONLY_DIR="$RESULTS_DIR/indexonly_${TIMESTAMP}"
mkdir -p "$INDEXONLY_DIR"

# Use scaling files for meaningful comparison
for mb in 10 50 100; do
    f="$INDEX_DATA/size_${mb}mb.fasta"
    if [[ -f "$f" ]]; then
        ensure_indexes "$f"
        echo "  ${mb}MB: index-only vs full scan..."
        bench_stats "$f" "$INDEXONLY_DIR/${mb}mb.json" "${mb}mb_"
    fi
done
echo ""

# ══════════════════════════════════════════════════════════════════════
#  3. File-Size Scaling
# ══════════════════════════════════════════════════════════════════════

if ! $SKIP_SCALING; then
    echo "──────────────────────────────────────────────────"
    echo " [3] File-Size Scaling"
    echo "──────────────────────────────────────────────────"

    SCALE_DIR="$RESULTS_DIR/scale_size_${TIMESTAMP}"
    mkdir -p "$SCALE_DIR"

    SIZE_MBS=(1 5 10 25 50 100 250 500 1000)

    for mb in "${SIZE_MBS[@]}"; do
        f="$INDEX_DATA/size_${mb}mb.fasta"
        if [[ -f "$f" ]]; then
            ensure_indexes "$f"
            echo "  ${mb}MB..."
            bench_stats "$f" "$SCALE_DIR/${mb}mb.json" "${mb}mb_"
        else
            echo "  SKIP: ${mb}MB file not found (run bench/index/run_benchmarks.sh first)"
        fi
    done
    echo ""
fi

# ══════════════════════════════════════════════════════════════════════
#  4. Real Dataset Stats
# ══════════════════════════════════════════════════════════════════════

if ! $SKIP_REAL; then
    echo "──────────────────────────────────────────────────"
    echo " [4] Real Dataset Stats"
    echo "──────────────────────────────────────────────────"

    declare -A DATASETS
    [[ -f "$DATA_DIR/REAL_Genome.fa" ]]        && DATASETS["Genome"]="$DATA_DIR/REAL_Genome.fa"
    [[ -f "$DATA_DIR/REAL_Transcriptome.fa" ]] && DATASETS["Transcriptome"]="$DATA_DIR/REAL_Transcriptome.fa"
    [[ -f "$DATA_DIR/REAL_Proteome.fasta" ]]   && DATASETS["Proteome"]="$DATA_DIR/REAL_Proteome.fasta"

    if [[ ${#DATASETS[@]} -eq 0 ]]; then
        echo "  No REAL_* datasets found. Run bench/shared/download_data.sh first."
    else
        REAL_DIR="$RESULTS_DIR/real_${TIMESTAMP}"
        mkdir -p "$REAL_DIR"

        for name in "${!DATASETS[@]}"; do
            file="${DATASETS[$name]}"
            ensure_indexes "$file"
            echo "  ${name} ($(du -h "$file" | cut -f1))..."
            bench_stats "$file" "$REAL_DIR/${name}.json" "${name}_"
        done
    fi
    echo ""
fi

# ══════════════════════════════════════════════════════════════════════
#  5. Memory Measurement
# ══════════════════════════════════════════════════════════════════════

echo "──────────────────────────────────────────────────"
echo " [5] Memory Usage"
echo "──────────────────────────────────────────────────"

MEM_CSV="$RESULTS_DIR/memory_${TIMESTAMP}.csv"
echo "tool,time_s,mem_kb,major_faults,minor_faults,file_size_mb" > "$MEM_CSV"

for mb in 10 50 100; do
    f="$INDEX_DATA/size_${mb}mb.fasta"
    if [[ -f "$f" ]]; then
        ensure_indexes "$f"
        echo "  ${mb}MB full scan..."
        while IFS= read -r line; do
            echo "$line,$mb" >> "$MEM_CSV"
        done < <(measure_memory_stats "$f" "full")

        echo "  ${mb}MB index-only..."
        while IFS= read -r line; do
            echo "$line,$mb" >> "$MEM_CSV"
        done < <(measure_memory_stats "$f" "indexonly")
    fi
done

echo "  Memory data → $MEM_CSV"
echo ""

# ══════════════════════════════════════════════════════════════════════
#  6. Throughput Calculation
# ══════════════════════════════════════════════════════════════════════

echo "──────────────────────────────────────────────────"
echo " [6] Throughput Summary"
echo "──────────────────────────────────────────────────"

THROUGHPUT_CSV="$RESULTS_DIR/throughput_${TIMESTAMP}.csv"
echo "tool,file_size_mb,time_s,throughput_mbs" > "$THROUGHPUT_CSV"

for mb in 10 50 100 250 500 1000; do
    f="$INDEX_DATA/size_${mb}mb.fasta"
    if [[ -f "$f" ]]; then
        ensure_indexes "$f"
        # Time z-fasta stats full scan
        TIME_S=$( { /usr/bin/time -f "%e" "$ZFASTA" stats "$f" > /dev/null; } 2>&1 ) || true
        if [[ -n "$TIME_S" ]] && (( $(echo "$TIME_S > 0" | bc -l 2>/dev/null || echo 0) )); then
            THROUGHPUT=$(echo "scale=1; $mb / $TIME_S" | bc -l 2>/dev/null || echo "0")
            echo "z-fasta-full,$mb,$TIME_S,$THROUGHPUT" >> "$THROUGHPUT_CSV"
            echo "  z-fasta full ${mb}MB: ${TIME_S}s → ${THROUGHPUT} MB/s"
        fi

        if $HAS_SEQKIT; then
            TIME_S=$( { /usr/bin/time -f "%e" "$SEQKIT" stats -a "$f" > /dev/null; } 2>&1 ) || true
            if [[ -n "$TIME_S" ]] && (( $(echo "$TIME_S > 0" | bc -l 2>/dev/null || echo 0) )); then
                THROUGHPUT=$(echo "scale=1; $mb / $TIME_S" | bc -l 2>/dev/null || echo "0")
                echo "seqkit-stats-a,$mb,$TIME_S,$THROUGHPUT" >> "$THROUGHPUT_CSV"
                echo "  seqkit -a ${mb}MB: ${TIME_S}s → ${THROUGHPUT} MB/s"
            fi
        fi

        if $HAS_SEQTK; then
            TIME_S=$( { /usr/bin/time -f "%e" "$SEQTK" comp "$f" > /dev/null; } 2>&1 ) || true
            if [[ -n "$TIME_S" ]] && (( $(echo "$TIME_S > 0" | bc -l 2>/dev/null || echo 0) )); then
                THROUGHPUT=$(echo "scale=1; $mb / $TIME_S" | bc -l 2>/dev/null || echo "0")
                echo "seqtk-comp,$mb,$TIME_S,$THROUGHPUT" >> "$THROUGHPUT_CSV"
                echo "  seqtk comp ${mb}MB: ${TIME_S}s → ${THROUGHPUT} MB/s"
            fi
        fi
    fi
done

echo "  Throughput data → $THROUGHPUT_CSV"
echo ""

echo "════════════════════════════════════════════════════════════════"
echo "  STATS benchmarks complete."
echo "  Results in: $RESULTS_DIR"
echo "  Run: .venv/bin/python bench/stats/generate_report.py"
echo "════════════════════════════════════════════════════════════════"
