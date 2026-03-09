#!/bin/bash
# z-fasta GET Benchmark Runner
# Measures region extraction performance: runtime, memory, throughput.
# Uses hyperfine for precise timing with cache clearing between runs.
#
# Usage: ./run_benchmarks.sh [--runs N] [--skip-scaling] [--skip-real]
#
# Outputs:
#   results/single_<timestamp>.json    — single-region latency
#   results/fullseq_<timestamp>.json   — full-sequence extraction
#   results/scaling_<timestamp>/       — region-size scaling
#   results/memory_<timestamp>.csv     — memory (RSS) per tool/mode

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH_ROOT="$(dirname "$SCRIPT_DIR")"
PROJECT_ROOT="$(dirname "$BENCH_ROOT")"
RESULTS_DIR="$SCRIPT_DIR/results"
DATA_DIR="$BENCH_ROOT/shared/data"
INDEX_DATA="$BENCH_ROOT/index/data"

# ── Tools ──────────────────────────────────────────────────────────
ZFASTA="$PROJECT_ROOT/zig-out/bin/z-fasta"
SAMTOOLS="samtools"
SEQKIT="$PROJECT_ROOT/tools/seqkit"
FASTAHACK="$PROJECT_ROOT/tools/fastahack-1.0.0/fastahack"

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
[[ -x "$ZFASTA" ]]   || { echo "Error: z-fasta not found at $ZFASTA. Run: zig build -Doptimize=ReleaseFast"; exit 1; }
command -v "$SAMTOOLS" &>/dev/null || { echo "Error: samtools not found"; exit 1; }

HAS_SEQKIT=false;    [[ -x "$SEQKIT" ]]    && HAS_SEQKIT=true
HAS_FASTAHACK=false; [[ -x "$FASTAHACK" ]] && HAS_FASTAHACK=true

# Cache-clearing command
CACHE_CLEAR=""
BENCH_MODE="warm"
if sudo -n sh -c 'sync; echo 3 > /proc/sys/vm/drop_caches' 2>/dev/null; then
    CACHE_CLEAR="sudo sh -c 'sync; echo 3 > /proc/sys/vm/drop_caches'"
    WARMUP=0
    BENCH_MODE="cold (cache cleared before each run)"
    echo "Cache clearing: enabled (sudo available) → forcing --warmup 0"
else
    echo "Cache clearing: DISABLED (no passwordless sudo) → warm-cache mode"
fi

echo ""
echo "════════════════════════════════════════════════════════════════"
echo "  z-fasta GET Benchmark Runner"
echo "════════════════════════════════════════════════════════════════"
echo "  Mode: $BENCH_MODE"
echo "  Runs: $RUNS | Warmup: $WARMUP"
echo "  z-fasta:   $ZFASTA"
echo "  samtools:  $(which $SAMTOOLS) ($($SAMTOOLS --version | head -1))"
echo "  seqkit:    $HAS_SEQKIT"
echo "  fastahack: $HAS_FASTAHACK"
echo ""

# ── Helper: ensure index exists for a FASTA file ──────────────────
ensure_indexes() {
    local file="$1"
    [[ -f "${file}.fai" ]] || $SAMTOOLS faidx "$file" 2>/dev/null
    [[ -f "${file}.zfi" ]] || "$ZFASTA" index "$file" 2>/dev/null
}

# ── Helper: get first sequence name and length from .fai ──────────
first_seq_info() {
    local fai="$1.fai"
    head -1 "$fai" | awk -F'\t' '{print $1, $2}'
}

# ── Helper: get longest sequence name and length from .fai ────────
longest_seq_info() {
    local fai="$1.fai"
    sort -t$'\t' -k2 -n -r "$fai" | head -1 | awk -F'\t' '{print $1, $2}'
}

# ── Helper: run a single-region hyperfine benchmark ───────────────
bench_get_region() {
    local file="$1"
    local region="$2"
    local json_out="$3"
    local label_prefix="${4:-}"

    local prepare_cmd=""
    [[ -n "$CACHE_CLEAR" ]] && prepare_cmd="$CACHE_CLEAR"

    local cmds=()
    local names=()

    # z-fasta
    names+=("${label_prefix}z-fasta")
    cmds+=("$ZFASTA get '$file' '$region' > /dev/null")

    # samtools
    names+=("${label_prefix}samtools")
    cmds+=("$SAMTOOLS faidx '$file' '$region' > /dev/null")

    # seqkit
    if $HAS_SEQKIT; then
        names+=("${label_prefix}seqkit")
        cmds+=("$SEQKIT faidx '$file' '$region' > /dev/null 2>&1")
    fi

    # fastahack
    if $HAS_FASTAHACK; then
        names+=("${label_prefix}fastahack")
        cmds+=("$FASTAHACK '$file' '$region' > /dev/null 2>&1")
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
measure_memory_get() {
    local file="$1"
    local region="$2"

    local tools_and_cmds=(
        "z-fasta:$ZFASTA get '$file' '$region' > /dev/null"
        "samtools:$SAMTOOLS faidx '$file' '$region' > /dev/null"
    )
    $HAS_SEQKIT    && tools_and_cmds+=("seqkit:$SEQKIT faidx '$file' '$region' > /dev/null 2>&1")
    $HAS_FASTAHACK && tools_and_cmds+=("fastahack:$FASTAHACK '$file' '$region' > /dev/null 2>&1")

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
#  1. Single-Region Latency (small regions on test files)
# ══════════════════════════════════════════════════════════════════════

echo "──────────────────────────────────────────────────"
echo " [1] Single-Region Extraction Latency"
echo "──────────────────────────────────────────────────"

SINGLE_DIR="$RESULTS_DIR/single_${TIMESTAMP}"
mkdir -p "$SINGLE_DIR"

# Use synthetic scaling files from index bench data if available
# Otherwise use test files
SINGLE_TEST_FILES=()
for mb in 10 50 100; do
    f="$INDEX_DATA/size_${mb}mb.fasta"
    if [[ -f "$f" ]]; then
        SINGLE_TEST_FILES+=("$f")
    fi
done

# Fallback to test data
if [[ ${#SINGLE_TEST_FILES[@]} -eq 0 ]]; then
    for f in "$PROJECT_ROOT"/tests/data/simple.fasta "$PROJECT_ROOT"/tests/data/proteome.fasta; do
        [[ -f "$f" ]] && SINGLE_TEST_FILES+=("$f")
    done
fi

for file in "${SINGLE_TEST_FILES[@]}"; do
    ensure_indexes "$file"
    read -r SEQNAME SEQLEN <<< "$(first_seq_info "$file")"
    fname=$(basename "$file" .fasta)

    # Small region: 100 bases
    REGION_SIZE=100
    if [[ "$SEQLEN" -ge "$REGION_SIZE" ]]; then
        MID=$(( SEQLEN / 2 ))
        START=$(( MID - REGION_SIZE / 2 ))
        [[ "$START" -lt 1 ]] && START=1
        END=$(( START + REGION_SIZE - 1 ))
        REGION="${SEQNAME}:${START}-${END}"
        echo "  ${fname}: ${REGION} (${REGION_SIZE} bp)"
        bench_get_region "$file" "$REGION" "$SINGLE_DIR/${fname}_100bp.json" "${fname}_100bp_"
    fi

    # Medium region: 10,000 bases
    REGION_SIZE=10000
    if [[ "$SEQLEN" -ge "$REGION_SIZE" ]]; then
        MID=$(( SEQLEN / 2 ))
        START=$(( MID - REGION_SIZE / 2 ))
        [[ "$START" -lt 1 ]] && START=1
        END=$(( START + REGION_SIZE - 1 ))
        REGION="${SEQNAME}:${START}-${END}"
        echo "  ${fname}: ${REGION} (${REGION_SIZE} bp)"
        bench_get_region "$file" "$REGION" "$SINGLE_DIR/${fname}_10kbp.json" "${fname}_10kbp_"
    fi
done

echo ""

# ══════════════════════════════════════════════════════════════════════
#  2. Full-Sequence Extraction (varying file sizes)
# ══════════════════════════════════════════════════════════════════════

echo "──────────────────────────────────────────────────"
echo " [2] Full-Sequence Extraction"
echo "──────────────────────────────────────────────────"

FULLSEQ_DIR="$RESULTS_DIR/fullseq_${TIMESTAMP}"
mkdir -p "$FULLSEQ_DIR"

for mb in 1 10 50 100; do
    f="$INDEX_DATA/size_${mb}mb.fasta"
    if [[ -f "$f" ]]; then
        ensure_indexes "$f"
        read -r SEQNAME SEQLEN <<< "$(first_seq_info "$f")"
        echo "  ${mb}MB: full $SEQNAME ($SEQLEN bp)"
        bench_get_region "$f" "$SEQNAME" "$FULLSEQ_DIR/${mb}mb_full.json" "${mb}mb_full_"
    fi
done

echo ""

# ══════════════════════════════════════════════════════════════════════
#  3. Region-Size Scaling
# ══════════════════════════════════════════════════════════════════════

if ! $SKIP_SCALING; then
    echo "──────────────────────────────────────────────────"
    echo " [3] Region-Size Scaling"
    echo "──────────────────────────────────────────────────"

    # Use the 100MB file for scaling tests
    SCALE_FILE="$INDEX_DATA/size_100mb.fasta"
    if [[ ! -f "$SCALE_FILE" ]]; then
        SCALE_FILE="$INDEX_DATA/size_50mb.fasta"
    fi

    if [[ -f "$SCALE_FILE" ]]; then
        ensure_indexes "$SCALE_FILE"
        read -r SEQNAME SEQLEN <<< "$(first_seq_info "$SCALE_FILE")"

        SCALE_DIR="$RESULTS_DIR/scale_region_${TIMESTAMP}"
        mkdir -p "$SCALE_DIR"

        REGION_SIZES=(100 1000 10000 100000 1000000)

        for rsize in "${REGION_SIZES[@]}"; do
            if [[ "$SEQLEN" -ge "$rsize" ]]; then
                MID=$(( SEQLEN / 2 ))
                START=$(( MID - rsize / 2 ))
                [[ "$START" -lt 1 ]] && START=1
                END=$(( START + rsize - 1 ))
                REGION="${SEQNAME}:${START}-${END}"
                echo "  $rsize bp region..."
                bench_get_region "$SCALE_FILE" "$REGION" "$SCALE_DIR/${rsize}bp.json" "${rsize}bp_"
            fi
        done
    else
        echo "  SKIP: no scaling data files found"
    fi
    echo ""
fi

# ══════════════════════════════════════════════════════════════════════
#  4. Real Dataset Benchmarks
# ══════════════════════════════════════════════════════════════════════

if ! $SKIP_REAL; then
    echo "──────────────────────────────────────────────────"
    echo " [4] Real Dataset Extraction"
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
            read -r SEQNAME SEQLEN <<< "$(longest_seq_info "$file")"

            # Small region: 1000 bp
            if [[ "$SEQLEN" -ge 1000 ]]; then
                MID=$(( SEQLEN / 2 ))
                START=$(( MID - 500 ))
                [[ "$START" -lt 1 ]] && START=1
                END=$(( START + 999 ))
                REGION="${SEQNAME}:${START}-${END}"
                echo "  ${name}: 1kb region from $SEQNAME"
                bench_get_region "$file" "$REGION" "$REAL_DIR/${name}_1kb.json" "${name}_1kb_"
            fi

            # Large region: 1MB
            if [[ "$SEQLEN" -ge 1000000 ]]; then
                MID=$(( SEQLEN / 2 ))
                START=$(( MID - 500000 ))
                [[ "$START" -lt 1 ]] && START=1
                END=$(( START + 999999 ))
                REGION="${SEQNAME}:${START}-${END}"
                echo "  ${name}: 1MB region from $SEQNAME"
                bench_get_region "$file" "$REGION" "$REAL_DIR/${name}_1mb.json" "${name}_1mb_"
            fi

            # Full sequence extraction
            echo "  ${name}: full $SEQNAME ($SEQLEN bp)"
            bench_get_region "$file" "$SEQNAME" "$REAL_DIR/${name}_full.json" "${name}_full_"
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
echo "tool,time_s,mem_kb,major_faults,minor_faults,region_type" > "$MEM_CSV"

# Pick the best available file for memory measurement
MEM_FILE=""
for candidate in "$INDEX_DATA/size_100mb.fasta" "$INDEX_DATA/size_50mb.fasta" \
                 "$PROJECT_ROOT/tests/data/simple.fasta"; do
    if [[ -f "$candidate" ]]; then
        MEM_FILE="$candidate"
        break
    fi
done

if [[ -n "$MEM_FILE" ]]; then
    ensure_indexes "$MEM_FILE"
    read -r SEQNAME SEQLEN <<< "$(first_seq_info "$MEM_FILE")"

    # Small region
    if [[ "$SEQLEN" -ge 100 ]]; then
        REGION="${SEQNAME}:1-100"
        echo "  Measuring: small region (100 bp)..."
        while IFS= read -r line; do
            echo "${line},small" >> "$MEM_CSV"
        done < <(measure_memory_get "$MEM_FILE" "$REGION")
    fi

    # Large region
    if [[ "$SEQLEN" -ge 1000000 ]]; then
        REGION="${SEQNAME}:1-1000000"
        echo "  Measuring: large region (1 MB)..."
        while IFS= read -r line; do
            echo "${line},large" >> "$MEM_CSV"
        done < <(measure_memory_get "$MEM_FILE" "$REGION")
    fi

    # Full sequence
    echo "  Measuring: full sequence..."
    while IFS= read -r line; do
        echo "${line},full" >> "$MEM_CSV"
    done < <(measure_memory_get "$MEM_FILE" "$SEQNAME")
fi

echo "  Memory data → $MEM_CSV"
echo ""

echo "════════════════════════════════════════════════════════════════"
echo "  GET benchmarks complete."
echo "  Results in: $RESULTS_DIR"
echo "  Run: python3 bench/get/generate_report.py"
echo "════════════════════════════════════════════════════════════════"
