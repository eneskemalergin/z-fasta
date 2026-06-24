#!/bin/bash
# z-fasta GET Benchmark Runner
# Measures region extraction performance: runtime, memory, throughput.
# Uses hyperfine for precise timing with cache clearing between runs.
#
# Usage: ./run_benchmarks.sh [--runs N] [--skip-scaling] [--skip-real] [--skip-rc] [--rc-full]
#
# Outputs:
#   results/single_<timestamp>.json    single-region latency
#   results/fullseq_<timestamp>.json   full-sequence extraction
#   results/scaling_<timestamp>/       region-size scaling
#   results/memory_<timestamp>.csv     memory (RSS) per tool/mode

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
BEDTOOLS="bedtools"
SEQKIT="$PROJECT_ROOT/tools/seqkit"
FASTAHACK="$PROJECT_ROOT/tools/fastahack-1.0.0/fastahack"
SEQTK="$PROJECT_ROOT/tools/seqtk/seqtk"
# pyfaidx: prefer venv, fall back to conda/system python env
if [[ -x "$PROJECT_ROOT/.venv/bin/faidx" ]]; then
    PYFAIDX="$PROJECT_ROOT/.venv/bin/faidx"
elif command -v faidx &>/dev/null; then
    PYFAIDX="$(command -v faidx)"
else
    PYFAIDX=""
fi

# ── Defaults ───────────────────────────────────────────────────────
RUNS=5
WARMUP=1
SKIP_SCALING=false
SKIP_REAL=false
SKIP_RC=false
RC_PROFILE="quick"

while [[ $# -gt 0 ]]; do
    case $1 in
        --runs)          RUNS="$2"; shift 2 ;;
        --warmup)        WARMUP="$2"; shift 2 ;;
        --skip-scaling)  SKIP_SCALING=true; shift ;;
        --skip-real)     SKIP_REAL=true; shift ;;
        --skip-rc)       SKIP_RC=true; shift ;;
        --rc-full)       RC_PROFILE="full"; shift ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

# ── Preflight ──────────────────────────────────────────────────────
mkdir -p "$RESULTS_DIR"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
TMP_WORK_DIR=$(mktemp -d "$SCRIPT_DIR/.bench_tmp.XXXXXX")
trap 'rm -rf "$TMP_WORK_DIR"' EXIT

command -v hyperfine &>/dev/null || { echo "Error: hyperfine not found (apt install hyperfine)"; exit 1; }
[[ -x "$ZFASTA" ]]   || { echo "Error: z-fasta not found at $ZFASTA. Run: ./zig build -Doptimize=ReleaseFast"; exit 1; }
command -v "$SAMTOOLS" &>/dev/null || { echo "Error: samtools not found"; exit 1; }

HAS_BEDTOOLS=false;  command -v "$BEDTOOLS" &>/dev/null && HAS_BEDTOOLS=true
HAS_SEQKIT=false;    [[ -x "$SEQKIT" ]]    && HAS_SEQKIT=true
HAS_FASTAHACK=false; [[ -x "$FASTAHACK" ]] && HAS_FASTAHACK=true
HAS_SEQTK=false;     [[ -x "$SEQTK" ]]     && HAS_SEQTK=true
HAS_PYFAIDX=false;   [[ -n "$PYFAIDX" && -x "$PYFAIDX" ]] && HAS_PYFAIDX=true

# Cache-clearing command
CACHE_CLEAR=""
BENCH_MODE="warm"
if sudo -n sh -c 'sync; echo 3 > /proc/sys/vm/drop_caches' 2>/dev/null; then
    CACHE_CLEAR="sudo sh -c 'sync; echo 3 > /proc/sys/vm/drop_caches'"
    WARMUP=0
    BENCH_MODE="cold (cache cleared before each run)"
    echo "Cache clearing: enabled (sudo available) -> forcing --warmup 0"
else
    echo "Cache clearing: DISABLED (no passwordless sudo) -> warm-cache mode"
fi

echo ""
echo "════════════════════════════════════════════════════════════════"
echo "  z-fasta GET Benchmark Runner"
echo "════════════════════════════════════════════════════════════════"
echo "  Mode: $BENCH_MODE"
echo "  Runs: $RUNS | Warmup: $WARMUP"
echo "  z-fasta:   $ZFASTA"
echo "  samtools:  $(which $SAMTOOLS) ($($SAMTOOLS --version | head -1))"
echo "  bedtools:  $HAS_BEDTOOLS"
echo "  seqkit:    $HAS_SEQKIT"
echo "  fastahack: $HAS_FASTAHACK"
echo "  seqtk:     $HAS_SEQTK ($SEQTK)"
echo "  pyfaidx:   $HAS_PYFAIDX ($PYFAIDX)"
echo "  rc bench:  $([[ "$SKIP_RC" == true ]] && echo skipped || echo "$RC_PROFILE")"
echo ""

# ── Helper: ensure index exists for a FASTA file ──────────────────
ensure_indexes() {
    local file="$1"
    [[ -f "${file}.fai" ]] || $SAMTOOLS faidx "$file" 2>/dev/null
    [[ -f "${file}.zfi" ]] || "$ZFASTA" index "$file" 2>/dev/null
}
# ── Helper: convert samtools region to seqtk BED file ──────────────
# seqtk subseq takes either a sequence-name list or a BED file (0-based start).
# Samtools regions are 1-based inclusive, so start_bed = start - 1.
make_seqtk_region_file() {
    local region="$1"
    local tmpfile
    tmpfile=$(mktemp /tmp/seqtk_region.XXXXXX.bed)
    if [[ "$region" == *:* ]]; then
        local seqname="${region%%:*}"
        local coords="${region#*:}"
        local start="${coords%-*}"
        local end="${coords#*-}"
        local start0=$(( start - 1 ))
        printf '%s\t%s\t%s\n' "$seqname" "$start0" "$end" > "$tmpfile"
    else
        # Full sequence: just the name
        printf '%s\n' "$region" > "$tmpfile"
    fi
    echo "$tmpfile"
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
        cmds+=("$FASTAHACK -r '$region' '$file' > /dev/null 2>&1")
    fi

    # seqtk subseq (requires a BED/name file for region specification)
    local seqtk_region_file=""
    if $HAS_SEQTK; then
        seqtk_region_file=$(make_seqtk_region_file "$region")
        names+=("${label_prefix}seqtk")
        cmds+=("$SEQTK subseq '$file' '$seqtk_region_file' > /dev/null 2>&1")
    fi

    # pyfaidx (accepts samtools-style regions directly)
    if $HAS_PYFAIDX; then
        names+=("${label_prefix}pyfaidx")
        cmds+=("$PYFAIDX '$file' '$region' > /dev/null 2>&1")
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

    # Clean up seqtk temp file after benchmark
    [[ -n "$seqtk_region_file" ]] && rm -f "$seqtk_region_file"
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
    $HAS_FASTAHACK && tools_and_cmds+=("fastahack:$FASTAHACK -r '$region' '$file' > /dev/null 2>&1")

    # seqtk: pre-create the region file once (not per /usr/bin/time run)
    local seqtk_mem_file=""
    if $HAS_SEQTK; then
        seqtk_mem_file=$(make_seqtk_region_file "$region")
        tools_and_cmds+=("seqtk:$SEQTK subseq '$file' '$seqtk_mem_file' > /dev/null 2>&1")
    fi
    $HAS_PYFAIDX && tools_and_cmds+=("pyfaidx:$PYFAIDX '$file' '$region' > /dev/null 2>&1")

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

    # Clean up seqtk temp file
    [[ -n "$seqtk_mem_file" ]] && rm -f "$seqtk_mem_file"
}

generate_bed_benchmark_files() {
    local file="$1"
    local count="$2"
    local bed_out="$3"
    local regions_out="$4"

    read -r seqname seqlen <<< "$(longest_seq_info "$file")"

    local span=100
    [[ "$seqlen" -lt "$span" ]] && span="$seqlen"

    local max_start=$(( seqlen - span ))
    [[ "$max_start" -lt 0 ]] && max_start=0

    : > "$bed_out"
    : > "$regions_out"

    for ((i = 0; i < count; i++)); do
        local start0=0
        if [[ "$count" -gt 1 && "$max_start" -gt 0 ]]; then
            start0=$(( (i * max_start) / (count - 1) ))
        fi
        local end0=$(( start0 + span ))
        local strand='+'
        (( i % 3 == 0 )) && strand='-'

        printf '%s\t%d\t%d\tbed_%d\t0\t%s\n' "$seqname" "$start0" "$end0" "$i" "$strand" >> "$bed_out"
        printf '%s:%d-%d\n' "$seqname" $(( start0 + 1 )) "$end0" >> "$regions_out"
    done
}

bench_get_bed() {
    local file="$1"
    local bed_file="$2"
    local regions_file="$3"
    local chunk_size="$4"
    local json_out="$5"
    local label_prefix="$6"
    local honor_strand="$7"

    local prepare_cmd=""
    [[ -n "$CACHE_CLEAR" ]] && prepare_cmd="$CACHE_CLEAR"

    local cmds=()
    local names=()

    local zf_cmd="$ZFASTA get '$file' --bed '$bed_file' --chunk-size '$chunk_size'"
    if [[ "$honor_strand" == true ]]; then
        zf_cmd+=" --strand-aware"
    fi
    zf_cmd+=" > /dev/null"
    names+=("${label_prefix}z-fasta")
    cmds+=("$zf_cmd")

    names+=("${label_prefix}samtools")
    cmds+=("$SAMTOOLS faidx -r '$regions_file' '$file' > /dev/null")

    if $HAS_BEDTOOLS; then
        local bt_cmd="$BEDTOOLS getfasta -fi '$file' -bed '$bed_file'"
        if [[ "$honor_strand" == true ]]; then
            bt_cmd+=" -s"
        fi
        bt_cmd+=" > /dev/null"
        names+=("${label_prefix}bedtools")
        cmds+=("$bt_cmd")
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
#  5. BED Batch Extraction
# ══════════════════════════════════════════════════════════════════════

echo "──────────────────────────────────────────────────"
echo " [5] BED Batch Extraction"
echo "──────────────────────────────────────────────────"

BED_FILE=""
for candidate in \
    "$INDEX_DATA/size_100mb.fasta" \
    "$INDEX_DATA/size_50mb.fasta" \
    "$DATA_DIR/REAL_Genome.fa" \
    "$PROJECT_ROOT/tests/data/simple.fasta"; do
    if [[ -f "$candidate" ]]; then
        BED_FILE="$candidate"
        break
    fi
done

if [[ -n "$BED_FILE" ]]; then
    ensure_indexes "$BED_FILE"
    BED_DIR="$RESULTS_DIR/bed_${TIMESTAMP}"
    mkdir -p "$BED_DIR"

    BED_COUNTS=(100 1000 10000 100000)
    BED_CHUNKS=(31 97 257 4096)

    echo "  File: $BED_FILE"
    for i in "${!BED_COUNTS[@]}"; do
        count="${BED_COUNTS[$i]}"
        chunk="${BED_CHUNKS[$i]}"
        bed_path="$TMP_WORK_DIR/${count}regions.bed"
        regions_path="$TMP_WORK_DIR/${count}regions.txt"

        generate_bed_benchmark_files "$BED_FILE" "$count" "$bed_path" "$regions_path"

        echo "  ${count} BED regions (default, chunk-size=${chunk})..."
        bench_get_bed "$BED_FILE" "$bed_path" "$regions_path" "$chunk" \
            "$BED_DIR/${count}regions_default.json" "${count}regions_default_" false

        if $HAS_BEDTOOLS; then
            echo "  ${count} BED regions (stranded, chunk-size=${chunk})..."
            bench_get_bed "$BED_FILE" "$bed_path" "$regions_path" "$chunk" \
                "$BED_DIR/${count}regions_stranded.json" "${count}regions_stranded_" true
        fi
    done
else
    echo "  SKIP: no suitable FASTA found for BED benchmarks"
fi

echo ""

# ══════════════════════════════════════════════════════════════════════
#  6. Memory Measurement
# ══════════════════════════════════════════════════════════════════════

echo "──────────────────────────────────────────────────"
echo " [6] Memory Usage"
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

echo "  Memory data -> $MEM_CSV"
echo ""

# ══════════════════════════════════════════════════════════════════════
#  7. Multi-Region Extraction (v0.2.4)
# ══════════════════════════════════════════════════════════════════════

echo "──────────────────────────────────────────────────"
echo " [7] Multi-Region Extraction (v0.2.4)"
echo "──────────────────────────────────────────────────"

# Helper: make a seqtk BED file from multiple samtools-style regions
make_seqtk_multi_bed() {
    local tmpfile
    tmpfile=$(mktemp /tmp/seqtk_multi.XXXXXX.bed)
    for region in "$@"; do
        if [[ "$region" == *:* ]]; then
            local seqname="${region%%:*}"
            local coords="${region#*:}"
            local start="${coords%-*}"
            local end="${coords#*-}"
            printf '%s\t%s\t%s\n' "$seqname" "$(( start - 1 ))" "$end"
        else
            printf '%s\n' "$region"
        fi
    done > "$tmpfile"
    echo "$tmpfile"
}

# Helper: bench_get_multi <file> <region_count_label> <json_out> <label_prefix> region1 region2 ...
bench_get_multi() {
    local file="$1"
    local json_out="$2"
    local label_prefix="$3"
    shift 3
    local regions=("$@")

    local prepare_cmd=""
    [[ -n "$CACHE_CLEAR" ]] && prepare_cmd="$CACHE_CLEAR"

    local cmds=()
    local names=()

    # z-fasta (native multi-region)
    names+=("${label_prefix}z-fasta")
    local zf_cmd="$ZFASTA get '$file'"
    for r in "${regions[@]}"; do zf_cmd+=" '$r'"; done
    zf_cmd+=" > /dev/null"
    cmds+=("$zf_cmd")

    # samtools (accepts multiple region args natively)
    names+=("${label_prefix}samtools")
    local st_cmd="$SAMTOOLS faidx '$file'"
    for r in "${regions[@]}"; do st_cmd+=" '$r'"; done
    st_cmd+=" > /dev/null"
    cmds+=("$st_cmd")

    # seqtk (BED file with all regions)
    if $HAS_SEQTK; then
        local seqtk_bed
        seqtk_bed=$(make_seqtk_multi_bed "${regions[@]}")
        names+=("${label_prefix}seqtk")
        cmds+=("$SEQTK subseq '$file' '$seqtk_bed' > /dev/null 2>&1")
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

    [[ $HAS_SEQTK == true ]] && rm -f "$seqtk_bed" 2>/dev/null || true
}

# Pick a file with enough sequences
MULTI_FILE=""
MULTI_SEQ_COUNT=0
for candidate in \
    "$DATA_DIR/REAL_Proteome.fasta" \
    "$DATA_DIR/REAL_Transcriptome.fa" \
    "$INDEX_DATA/size_100mb.fasta" \
    "$INDEX_DATA/size_50mb.fasta" \
    "$PROJECT_ROOT/tests/data/edge_cases.fasta" \
    "$PROJECT_ROOT/tests/data/simple.fasta"; do
    if [[ -f "${candidate}.fai" ]] || { [[ -f "$candidate" ]] && "$SAMTOOLS" faidx "$candidate" 2>/dev/null; }; then
        CNT=$(wc -l < "${candidate}.fai" 2>/dev/null || echo 0)
        if [[ "$CNT" -gt "$MULTI_SEQ_COUNT" ]]; then
            MULTI_FILE="$candidate"
            MULTI_SEQ_COUNT=$CNT
        fi
    fi
done

if [[ -n "$MULTI_FILE" && "$MULTI_SEQ_COUNT" -gt 0 ]]; then
    ensure_indexes "$MULTI_FILE"
    MULTI_DIR="$RESULTS_DIR/multi_${TIMESTAMP}"
    mkdir -p "$MULTI_DIR"

    echo "  File: $MULTI_FILE ($MULTI_SEQ_COUNT sequences)"

    # Read all sequence names and lengths from .fai
    mapfile -t ALL_NAMES < <(awk -F'\t' '{print $1}' "${MULTI_FILE}.fai")
    mapfile -t ALL_LENS  < <(awk -F'\t' '{print $2}' "${MULTI_FILE}.fai")

    # For each region count, pick evenly spaced sequences + 100bp sub-regions
    for RCOUNT in 1 10 50 100; do
        if [[ "$RCOUNT" -gt "$MULTI_SEQ_COUNT" ]]; then
            echo "  SKIP: $RCOUNT regions (only $MULTI_SEQ_COUNT sequences available)"
            continue
        fi

        REGIONS=()
        STEP=$(( MULTI_SEQ_COUNT / RCOUNT ))
        [[ "$STEP" -lt 1 ]] && STEP=1
        idx2=0
        while [[ ${#REGIONS[@]} -lt $RCOUNT && $idx2 -lt $MULTI_SEQ_COUNT ]]; do
            NAME="${ALL_NAMES[$idx2]}"
            LEN="${ALL_LENS[$idx2]}"
            if [[ "$LEN" -ge 100 ]]; then
                MID=$(( LEN / 2 ))
                START=$(( MID - 50 ))
                [[ "$START" -lt 1 ]] && START=1
                END=$(( START + 99 ))
                REGIONS+=("${NAME}:${START}-${END}")
            else
                REGIONS+=("$NAME")
            fi
            idx2=$(( idx2 + STEP ))
        done

        echo "  $RCOUNT region(s)..."
        bench_get_multi "$MULTI_FILE" \
            "$MULTI_DIR/${RCOUNT}regions.json" \
            "${RCOUNT}regions_" \
            "${REGIONS[@]}"
    done
else
    echo "  SKIP: no suitable multi-sequence file found"
fi

if ! $SKIP_RC; then
    echo ""
    echo "──────────────────────────────────────────────────"
    echo " [8] RC Orientation Review"
    echo "──────────────────────────────────────────────────"

    rc_run_hyperfine_markdown() {
        local markdown_out="$1"
        local json_out="$2"
        shift 2
        hyperfine --warmup "$WARMUP" --runs "$RUNS" --export-markdown "$markdown_out" --export-json "$json_out" "$@"
    }

    rc_make_region() {
        local start="$1"
        local size="$2"
        local end=$(( start + size - 1 ))
        printf '%s:%d-%d' "$RC_SEQ_NAME" "$start" "$end"
    }

    rc_build_multi_command() {
        local count="$1"
        local with_rc="$2"
        local args=()
        local size=1000
        local start=1000
        local max_start=$(( RC_SEQ_LEN - size - 1 ))
        local usable_span=$(( max_start - start ))
        local step

        if (( usable_span <= 0 )); then
            echo "RC review fixture is too short for multi-region benchmarking" >&2
            exit 1
        fi

        step=$(( usable_span / (count + 1) ))
        if (( step < size )); then
            step=$size
        fi

        for ((i = 0; i < count; i += 1)); do
            args+=("$(rc_make_region $(( start + (i * step) )) "$size")")
        done

        printf '%q get %q' "$ZFASTA" "$RC_FASTA"
        for region in "${args[@]}"; do
            printf ' %q' "$region"
        done
        if [[ "$with_rc" == "true" ]]; then
            printf ' --rc'
        fi
        printf ' > /dev/null'
    }

    rc_samtools_multi_command() {
        local count="$1"
        local size=1000
        local start=1000
        local max_start=$(( RC_SEQ_LEN - size - 1 ))
        local usable_span=$(( max_start - start ))
        local step=$(( usable_span / (count + 1) ))

        if (( step < size )); then
            step=$size
        fi

        printf '%q faidx -i --mark-strand no %q' "$SAMTOOLS" "$RC_FASTA"
        for ((i = 0; i < count; i += 1)); do
            printf ' %q' "$(rc_make_region $(( start + (i * step) )) "$size")"
        done
        printf ' > /dev/null'
    }

    rc_bedtools_seqtk_command() {
        local bed_file="$1"
        local stranded="$2"

        if [[ "$stranded" == "true" ]]; then
            printf '%q getfasta -fi %q -bed %q -s | %q seq -r > /dev/null' "$BEDTOOLS" "$RC_FASTA" "$bed_file" "$SEQTK"
        else
            printf '%q getfasta -fi %q -bed %q | %q seq -r > /dev/null' "$BEDTOOLS" "$RC_FASTA" "$bed_file" "$SEQTK"
        fi
    }

    rc_measure_rss() {
        local label="$1"
        local command="$2"
        local tmpfile="$TMP_WORK_DIR/rc_${label}.time"

        bash -lc "/usr/bin/time -v $command" >/dev/null 2> "$tmpfile"
        awk -v label="$label" '
            /Elapsed \(wall clock\) time/ { elapsed = $NF }
            /Maximum resident set size/ { rss = $NF }
            END { printf "%s\t%s\t%s\n", label, elapsed, rss }
        ' "$tmpfile"
    }

    RC_DIR="$RESULTS_DIR/rc_review_${TIMESTAMP}"
    RC_JSON_DIR="$RC_DIR/json"
    RC_WORK_DIR="$TMP_WORK_DIR/rc_review"
    mkdir -p "$RC_DIR" "$RC_JSON_DIR" "$RC_WORK_DIR"

    RC_FASTA="$RC_WORK_DIR/rc_review.fa"
    RC_BED="$RC_WORK_DIR/rc_review.bed"
    RC_BED_STRANDED="$RC_WORK_DIR/rc_review_stranded.bed"
    RC_SEQ_NAME="chrSynthetic"

    if [[ "$RC_PROFILE" == "full" ]]; then
        RC_SEQ_LEN=8000000
        RC_BED_COUNT=10000
        RC_LARGE_REGION_SIZE=1000000
        RC_MULTI_COUNT_A=10
        RC_MULTI_COUNT_B=100
    else
        RC_SEQ_LEN=2000000
        RC_BED_COUNT=2000
        RC_LARGE_REGION_SIZE=250000
        RC_MULTI_COUNT_A=10
        RC_MULTI_COUNT_B=50
    fi

    echo "  Preparing synthetic RC fixture..."
    {
        echo ">${RC_SEQ_NAME} reverse-complement review fixture"
        awk -v total="$RC_SEQ_LEN" -v width=60 '
            BEGIN {
                pattern = "ACGTNRYWSKMBDHVacgtnrywskmbdhv"
                seq = ""
                while (length(seq) < total) {
                    seq = seq pattern
                }
                seq = substr(seq, 1, total)
                while (length(seq) > 0) {
                    print substr(seq, 1, width)
                    seq = substr(seq, width + 1)
                }
            }
        '
    } > "$RC_FASTA"
    "$ZFASTA" index "$RC_FASTA" >/dev/null

    : > "$RC_BED"
    : > "$RC_BED_STRANDED"
    for ((i = 0; i < RC_BED_COUNT; i += 1)); do
        s=$(( 1000 + (i * 400) ))
        e=$(( s + 120 ))
        printf '%s\t%d\t%d\n' "$RC_SEQ_NAME" "$s" "$e" >> "$RC_BED"
        if ((( i % 2 ) == 0 )); then
            strand='+'
        else
            strand='-'
        fi
        printf '%s\t%d\t%d\tregion_%05d\t0\t%s\n' "$RC_SEQ_NAME" "$s" "$e" "$i" "$strand" >> "$RC_BED_STRANDED"
    done

    RC_REGION_SMALL="$(rc_make_region 1000 100)"
    RC_REGION_MEDIUM="$(rc_make_region 200000 10000)"
    RC_REGION_LARGE="$(rc_make_region 200000 "$RC_LARGE_REGION_SIZE")"
    RC_REGION_FULL="$RC_SEQ_NAME"
    RC_MULTI_A_FORWARD="$(rc_build_multi_command "$RC_MULTI_COUNT_A" false)"
    RC_MULTI_A_RC="$(rc_build_multi_command "$RC_MULTI_COUNT_A" true)"
    RC_MULTI_B_FORWARD="$(rc_build_multi_command "$RC_MULTI_COUNT_B" false)"
    RC_MULTI_B_RC="$(rc_build_multi_command "$RC_MULTI_COUNT_B" true)"
    RC_SAMTOOLS_MULTI_A_RC=""
    RC_SAMTOOLS_MULTI_B_RC=""
    RC_BEDTOOLS_BED_RC=""
    RC_BEDTOOLS_BED_STRANDED_RC=""

    if [[ "$HAS_SEQTK" == true && "$HAS_BEDTOOLS" == true ]]; then
        RC_BEDTOOLS_BED_RC="$(rc_bedtools_seqtk_command "$RC_BED" false)"
        RC_BEDTOOLS_BED_STRANDED_RC="$(rc_bedtools_seqtk_command "$RC_BED_STRANDED" true)"
    fi
    RC_HAS_SAMTOOLS=true
    if ! command -v "$SAMTOOLS" >/dev/null 2>&1; then
        RC_HAS_SAMTOOLS=false
    else
        RC_SAMTOOLS_MULTI_A_RC="$(rc_samtools_multi_command "$RC_MULTI_COUNT_A")"
        RC_SAMTOOLS_MULTI_B_RC="$(rc_samtools_multi_command "$RC_MULTI_COUNT_B")"
    fi

    echo "  profile: $RC_PROFILE (runs=$RUNS, warmup=$WARMUP, seq_len=$RC_SEQ_LEN, bed_count=$RC_BED_COUNT)"
    echo "  [8.1] orientation small"
    RC_ORIENTATION_SMALL_ARGS=(
        -n forward "$ZFASTA get '$RC_FASTA' '$RC_REGION_SMALL' > /dev/null"
        -n zf-rc "$ZFASTA get '$RC_FASTA' '$RC_REGION_SMALL' --rc > /dev/null"
        -n reverse-only "$ZFASTA get '$RC_FASTA' '$RC_REGION_SMALL' --reverse-only > /dev/null"
        -n complement-only "$ZFASTA get '$RC_FASTA' '$RC_REGION_SMALL' --complement-only > /dev/null"
        -n rc-annotate "$ZFASTA get '$RC_FASTA' '$RC_REGION_SMALL' --rc --annotate-rc > /dev/null"
    )
    if [[ "$RC_HAS_SAMTOOLS" == true ]]; then
        RC_ORIENTATION_SMALL_ARGS+=(
            -n samtools-rc "$SAMTOOLS faidx -i --mark-strand no '$RC_FASTA' '$RC_REGION_SMALL' > /dev/null"
        )
    fi
    rc_run_hyperfine_markdown \
        "$RC_DIR/orientation_small.md" \
        "$RC_JSON_DIR/orientation_small.json" \
        "${RC_ORIENTATION_SMALL_ARGS[@]}"

    echo "  [8.2] orientation medium"
    RC_ORIENTATION_MEDIUM_ARGS=(
        -n forward "$ZFASTA get '$RC_FASTA' '$RC_REGION_MEDIUM' > /dev/null"
        -n zf-rc "$ZFASTA get '$RC_FASTA' '$RC_REGION_MEDIUM' --rc > /dev/null"
        -n reverse-only "$ZFASTA get '$RC_FASTA' '$RC_REGION_MEDIUM' --reverse-only > /dev/null"
        -n complement-only "$ZFASTA get '$RC_FASTA' '$RC_REGION_MEDIUM' --complement-only > /dev/null"
        -n rc-annotate "$ZFASTA get '$RC_FASTA' '$RC_REGION_MEDIUM' --rc --annotate-rc > /dev/null"
    )
    if [[ "$RC_HAS_SAMTOOLS" == true ]]; then
        RC_ORIENTATION_MEDIUM_ARGS+=(
            -n samtools-rc "$SAMTOOLS faidx -i --mark-strand no '$RC_FASTA' '$RC_REGION_MEDIUM' > /dev/null"
        )
    fi
    rc_run_hyperfine_markdown \
        "$RC_DIR/orientation_medium.md" \
        "$RC_JSON_DIR/orientation_medium.json" \
        "${RC_ORIENTATION_MEDIUM_ARGS[@]}"

    echo "  [8.3] orientation large"
    RC_ORIENTATION_LARGE_ARGS=(
        -n forward "$ZFASTA get '$RC_FASTA' '$RC_REGION_LARGE' > /dev/null"
        -n zf-rc "$ZFASTA get '$RC_FASTA' '$RC_REGION_LARGE' --rc > /dev/null"
        -n reverse-only "$ZFASTA get '$RC_FASTA' '$RC_REGION_LARGE' --reverse-only > /dev/null"
        -n complement-only "$ZFASTA get '$RC_FASTA' '$RC_REGION_LARGE' --complement-only > /dev/null"
        -n rc-annotate "$ZFASTA get '$RC_FASTA' '$RC_REGION_LARGE' --rc --annotate-rc > /dev/null"
    )
    if [[ "$RC_HAS_SAMTOOLS" == true ]]; then
        RC_ORIENTATION_LARGE_ARGS+=(
            -n samtools-rc "$SAMTOOLS faidx -i --mark-strand no '$RC_FASTA' '$RC_REGION_LARGE' > /dev/null"
        )
    fi
    rc_run_hyperfine_markdown \
        "$RC_DIR/orientation_large.md" \
        "$RC_JSON_DIR/orientation_large.json" \
        "${RC_ORIENTATION_LARGE_ARGS[@]}"

    echo "  [8.4] full-sequence forward vs rc"
    RC_FULL_ARGS=(
        -n forward "$ZFASTA get '$RC_FASTA' '$RC_REGION_FULL' > /dev/null"
        -n zf-rc "$ZFASTA get '$RC_FASTA' '$RC_REGION_FULL' --rc > /dev/null"
        -n rc-annotate "$ZFASTA get '$RC_FASTA' '$RC_REGION_FULL' --rc --annotate-rc > /dev/null"
    )
    if [[ "$RC_HAS_SAMTOOLS" == true ]]; then
        RC_FULL_ARGS+=(
            -n samtools-rc "$SAMTOOLS faidx -i --mark-strand no '$RC_FASTA' '$RC_REGION_FULL' > /dev/null"
        )
    fi
    rc_run_hyperfine_markdown \
        "$RC_DIR/full_sequence.md" \
        "$RC_JSON_DIR/full_sequence.json" \
        "${RC_FULL_ARGS[@]}"

    echo "  [8.5] multi-region no-flag regression vs rc"
    RC_MULTI_ARGS=(
        -n multi${RC_MULTI_COUNT_A}-forward "$RC_MULTI_A_FORWARD"
        -n multi${RC_MULTI_COUNT_A}-zf-rc "$RC_MULTI_A_RC"
        -n multi${RC_MULTI_COUNT_B}-forward "$RC_MULTI_B_FORWARD"
        -n multi${RC_MULTI_COUNT_B}-zf-rc "$RC_MULTI_B_RC"
    )
    if [[ "$RC_HAS_SAMTOOLS" == true ]]; then
        RC_MULTI_ARGS+=(
            -n multi${RC_MULTI_COUNT_A}-samtools-rc "$RC_SAMTOOLS_MULTI_A_RC"
            -n multi${RC_MULTI_COUNT_B}-samtools-rc "$RC_SAMTOOLS_MULTI_B_RC"
        )
    fi
    rc_run_hyperfine_markdown \
        "$RC_DIR/multi_region.md" \
        "$RC_JSON_DIR/multi_region.json" \
        "${RC_MULTI_ARGS[@]}"

    echo "  [8.6] bed batch orientation overhead"
    RC_BED_ARGS=(
        -n bed-forward "$ZFASTA get '$RC_FASTA' --bed '$RC_BED' > /dev/null"
        -n bed-zf-rc "$ZFASTA get '$RC_FASTA' --bed '$RC_BED' --rc > /dev/null"
        -n bed-reverse-only "$ZFASTA get '$RC_FASTA' --bed '$RC_BED' --reverse-only > /dev/null"
        -n bed-complement-only "$ZFASTA get '$RC_FASTA' --bed '$RC_BED' --complement-only > /dev/null"
        -n bed-honor-strand-zf-rc "$ZFASTA get '$RC_FASTA' --bed '$RC_BED_STRANDED' --honor-strand --rc > /dev/null"
    )
    if [[ -n "$RC_BEDTOOLS_BED_RC" ]]; then
        RC_BED_ARGS+=(
            -n bed-bedtools-seqtk-rc "$RC_BEDTOOLS_BED_RC"
            -n bed-honor-strand-bedtools-seqtk-rc "$RC_BEDTOOLS_BED_STRANDED_RC"
        )
    fi
    rc_run_hyperfine_markdown \
        "$RC_DIR/bed_batch.md" \
        "$RC_JSON_DIR/bed_batch.json" \
        "${RC_BED_ARGS[@]}"

    {
        printf 'label\telapsed\tmaxrss_kb\n'
        rc_measure_rss forward "$ZFASTA get '$RC_FASTA' '$RC_REGION_LARGE' > /dev/null"
        rc_measure_rss zf_rc "$ZFASTA get '$RC_FASTA' '$RC_REGION_LARGE' --rc > /dev/null"
        rc_measure_rss rc_annotate "$ZFASTA get '$RC_FASTA' '$RC_REGION_LARGE' --rc --annotate-rc > /dev/null"
        rc_measure_rss multi_b_forward "$RC_MULTI_B_FORWARD"
        rc_measure_rss multi_b_zf_rc "$RC_MULTI_B_RC"
        rc_measure_rss bed_forward "$ZFASTA get '$RC_FASTA' --bed '$RC_BED' > /dev/null"
        rc_measure_rss bed_honor_strand_zf_rc "$ZFASTA get '$RC_FASTA' --bed '$RC_BED_STRANDED' --honor-strand --rc > /dev/null"
        if [[ "$RC_HAS_SAMTOOLS" == true ]]; then
            rc_measure_rss samtools_rc "$SAMTOOLS faidx -i --mark-strand no '$RC_FASTA' '$RC_REGION_LARGE' > /dev/null"
            rc_measure_rss multi_b_samtools_rc "$RC_SAMTOOLS_MULTI_B_RC"
        fi
        if [[ -n "$RC_BEDTOOLS_BED_STRANDED_RC" ]]; then
            rc_measure_rss bed_honor_strand_bedtools_seqtk_rc "$RC_BEDTOOLS_BED_STRANDED_RC"
        fi
    } > "$RC_DIR/rss.tsv"

    echo "  RC outputs -> $RC_DIR"
fi

echo ""
echo "════════════════════════════════════════════════════════════════"
echo "  All GET benchmark slices complete."
echo "  Results in: $RESULTS_DIR"
echo "  Run: .venv/bin/python bench/get/generate_report.py"
echo "════════════════════════════════════════════════════════════════"

MANIFEST_PATH="$RESULTS_DIR/run_${TIMESTAMP}.json"
PRODUCED=(single fullseq bed memory multi)
if ! $SKIP_SCALING; then
    PRODUCED+=(scale_region)
fi
if ! $SKIP_REAL; then
    PRODUCED+=(real)
fi
if ! $SKIP_RC; then
    PRODUCED+=(rc_review)
fi

{
    echo "{"
    printf '  "timestamp": "%s",\n' "$TIMESTAMP"
    printf '  "runs": %s,\n' "$RUNS"
    printf '  "warmup": %s,\n' "$WARMUP"
    printf '  "bench_mode": "%s",\n' "$BENCH_MODE"
    printf '  "skip_scaling": %s,\n' "$([[ "$SKIP_SCALING" == true ]] && echo true || echo false)"
    printf '  "skip_real": %s,\n' "$([[ "$SKIP_REAL" == true ]] && echo true || echo false)"
    printf '  "skip_rc": %s,\n' "$([[ "$SKIP_RC" == true ]] && echo true || echo false)"
    printf '  "rc_profile": "%s",\n' "$RC_PROFILE"
    printf '  "produced": ['
    for i in "${!PRODUCED[@]}"; do
        if [[ "$i" -gt 0 ]]; then
            printf ', '
        fi
        printf '"%s"' "${PRODUCED[$i]}"
    done
    printf ']\n'
    echo "}"
} > "$MANIFEST_PATH"
