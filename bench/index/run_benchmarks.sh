#!/bin/bash
# z-fasta Benchmark Runner
# Uses hyperfine for precise timing with cache clearing between runs.
# Outputs raw JSON data — all analysis done by generate_report.py.
#
# Usage: ./run_benchmarks.sh [--runs N] [--skip-scaling] [--skip-real]
#
# Outputs:
#   results/perf_<timestamp>.json       — real dataset benchmarks
#   results/scaling_<timestamp>.json    — file-size scaling
#   results/seqscale_<timestamp>.json   — sequence-count scaling
#   results/memory_<timestamp>.csv      — memory (RSS) per tool/mode

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH_ROOT="$(dirname "$SCRIPT_DIR")"
PROJECT_ROOT="$(dirname "$BENCH_ROOT")"
RESULTS_DIR="$SCRIPT_DIR/results"
SCALING_DIR="$SCRIPT_DIR/data"
DATA_DIR="$BENCH_ROOT/shared/data"

# ── Tools ──────────────────────────────────────────────────────────
ZFASTA="$PROJECT_ROOT/zig-out/bin/z-fasta"
SAMTOOLS="samtools"
SEQKIT="$PROJECT_ROOT/tools/seqkit"
FASTAHACK="$PROJECT_ROOT/tools/fastahack-1.0.0/fastahack"
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
mkdir -p "$RESULTS_DIR" "$SCALING_DIR"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

command -v hyperfine &>/dev/null || { echo "Error: hyperfine not found (apt install hyperfine)"; exit 1; }
[[ -x "$ZFASTA" ]]   || { echo "Error: z-fasta not found at $ZFASTA. Run: ./zig build -Doptimize=ReleaseFast"; exit 1; }
command -v "$SAMTOOLS" &>/dev/null || { echo "Error: samtools not found"; exit 1; }

HAS_SEQKIT=false;    [[ -x "$SEQKIT" ]]    && HAS_SEQKIT=true
HAS_FASTAHACK=false; [[ -x "$FASTAHACK" ]] && HAS_FASTAHACK=true
HAS_PYFAIDX=false;   [[ -n "$PYFAIDX" && -x "$PYFAIDX" ]] && HAS_PYFAIDX=true

# Cache-clearing command (needs sudo for drop_caches)
CACHE_CLEAR=""
BENCH_MODE="warm"
if sudo -n sh -c 'sync; echo 3 > /proc/sys/vm/drop_caches' 2>/dev/null; then
    CACHE_CLEAR="sudo sh -c 'sync; echo 3 > /proc/sys/vm/drop_caches'"
    # Warmup is contradictory with cache clearing: the prepare command
    # drops the page cache before EVERY run (including warmup), so warmup
    # never actually warms anything — it just wastes time.
    WARMUP=0
    BENCH_MODE="cold (cache cleared before each run)"
    echo "Cache clearing: enabled (sudo available) → forcing --warmup 0"
else
    echo "Cache clearing: DISABLED (no passwordless sudo) → warm-cache mode"
    echo "  Tip: echo '$(whoami) ALL=(ALL) NOPASSWD: /usr/bin/sh -c sync*' | sudo tee /etc/sudoers.d/benchcache"
fi

echo ""
echo "════════════════════════════════════════════════════════════════"
echo "  z-fasta Benchmark Runner"
echo "════════════════════════════════════════════════════════════════"
echo "  Mode: $BENCH_MODE"
echo "  Runs: $RUNS | Warmup: $WARMUP"
echo "  z-fasta:   $ZFASTA ($(ls -lh "$ZFASTA" | awk '{print $5}'))"
echo "  samtools:  $(which $SAMTOOLS) ($($SAMTOOLS --version | head -1))"
echo "  seqkit:    $HAS_SEQKIT ($SEQKIT)"
echo "  fastahack: $HAS_FASTAHACK ($FASTAHACK)"
echo "  pyfaidx:   $HAS_PYFAIDX ($PYFAIDX)"
echo ""

# ── Helper: build hyperfine command for a single file ──────────────
# Runs all tools against one FASTA file, exports to JSON.
bench_file() {
    local file="$1"
    local json_out="$2"
    local label_prefix="${3:-}"  # e.g. "100MB_" for scaling labels

    local prepare_cmd="rm -f '${file}.fai' '${file}.zfi'"
    [[ -n "$CACHE_CLEAR" ]] && prepare_cmd="$prepare_cmd; $CACHE_CLEAR"

    local cmds=()
    local names=()

    # z-fasta (default: mmap + dedup)
    names+=("${label_prefix}z-fasta-default")
    cmds+=("$ZFASTA index --emit-fai '$file' > /dev/null")

    # z-fasta --no-dedup
    names+=("${label_prefix}z-fasta-nodedup")
    cmds+=("$ZFASTA index --emit-fai --no-dedup '$file' > /dev/null")

    # z-fasta --low-mem
    names+=("${label_prefix}z-fasta-lowmem")
    cmds+=("$ZFASTA index --low-mem '$file' > /dev/null")

    # samtools
    names+=("${label_prefix}samtools")
    cmds+=("$SAMTOOLS faidx '$file'")

    # seqkit
    if $HAS_SEQKIT; then
        names+=("${label_prefix}seqkit")
        cmds+=("$SEQKIT faidx '$file' > /dev/null 2>&1")
    fi

    # fastahack (indexes via first access: -i flag)
    if $HAS_FASTAHACK; then
        names+=("${label_prefix}fastahack")
        cmds+=("$FASTAHACK -i '$file' > /dev/null 2>&1")
    fi

    # pyfaidx (creates .fai on first access; --no-output suppresses sequence extraction)
    if $HAS_PYFAIDX; then
        names+=("${label_prefix}pyfaidx")
        cmds+=("$PYFAIDX '$file' --no-output 2>&1")
    fi

    # Build the hyperfine invocation
    local hf_args=(
        hyperfine
        --warmup "$WARMUP"
        --runs "$RUNS"
        --prepare "$prepare_cmd"
        --export-json "$json_out"
    )

    for i in "${!names[@]}"; do
        hf_args+=(-n "${names[$i]}" "${cmds[$i]}")
    done

    "${hf_args[@]}"
}

# ── Helper: measure memory + page faults with /usr/bin/time ────────
# Format: %e=elapsed %M=MaxRSS(KB) %F=major-faults %R=minor-faults
# Major faults = blocking disk reads (should be ~0 for warm cache).
# Minor faults = non-blocking page mappings (high for mmap tools).
measure_memory() {
    local file="$1"
    local csv_out="$2"

    local tools_and_cmds=(
        "z-fasta-default:$ZFASTA index --emit-fai '$file' > /dev/null"
        "z-fasta-nodedup:$ZFASTA index --emit-fai --no-dedup '$file' > /dev/null"
        "z-fasta-lowmem:$ZFASTA index --low-mem '$file' > /dev/null"
        "samtools:$SAMTOOLS faidx '$file'"
    )
    $HAS_SEQKIT    && tools_and_cmds+=("seqkit:$SEQKIT faidx '$file' > /dev/null 2>&1")
    $HAS_FASTAHACK && tools_and_cmds+=("fastahack:$FASTAHACK -i '$file' > /dev/null 2>&1")
    $HAS_PYFAIDX   && tools_and_cmds+=("pyfaidx:$PYFAIDX '$file' --no-output 2>&1")

    for entry in "${tools_and_cmds[@]}"; do
        local tool="${entry%%:*}"
        local cmd="${entry#*:}"

        rm -f "${file}.fai" "${file}.zfi" 2>/dev/null || true

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
#  1. Real Dataset Performance
# ══════════════════════════════════════════════════════════════════════

if ! $SKIP_REAL; then
    echo "──────────────────────────────────────────────────"
    echo " [1] Real Dataset Performance"
    echo "──────────────────────────────────────────────────"

    declare -A DATASETS
    [[ -f "$DATA_DIR/REAL_Genome.fa" ]]          && DATASETS["Genome"]="$DATA_DIR/REAL_Genome.fa"
    [[ -f "$DATA_DIR/REAL_Transcriptome.fa" ]]   && DATASETS["Transcriptome"]="$DATA_DIR/REAL_Transcriptome.fa"
    [[ -f "$DATA_DIR/REAL_Proteome.fasta" ]]     && DATASETS["Proteome"]="$DATA_DIR/REAL_Proteome.fasta"

    if [[ ${#DATASETS[@]} -eq 0 ]]; then
        echo "  No REAL_* datasets found in $DATA_DIR"
        echo "  Running download_data.sh to fetch real datasets..."
        echo ""
        bash "$BENCH_ROOT/shared/download_data.sh"
        echo ""
        # Re-check after download
        [[ -f "$DATA_DIR/REAL_Genome.fa" ]]          && DATASETS["Genome"]="$DATA_DIR/REAL_Genome.fa"
        [[ -f "$DATA_DIR/REAL_Transcriptome.fa" ]]   && DATASETS["Transcriptome"]="$DATA_DIR/REAL_Transcriptome.fa"
        [[ -f "$DATA_DIR/REAL_Proteome.fasta" ]]     && DATASETS["Proteome"]="$DATA_DIR/REAL_Proteome.fasta"
    fi

    if [[ ${#DATASETS[@]} -eq 0 ]]; then
        echo "  ERROR: download_data.sh failed to fetch any datasets. Skipping real benchmarks."
    else

    PERF_DIR="$RESULTS_DIR/perf_${TIMESTAMP}"
    mkdir -p "$PERF_DIR"

    for name in "${!DATASETS[@]}"; do
        file="${DATASETS[$name]}"
        echo ""
        echo "  Benchmarking: $name ($(du -h "$file" | cut -f1))"
        bench_file "$file" "$PERF_DIR/${name}.json"
    done

    # Memory measurement on largest dataset
    echo ""
    echo "  Measuring memory + page faults..."
    MEM_CSV="$RESULTS_DIR/memory_${TIMESTAMP}.csv"
    echo "tool,time_s,mem_kb,major_faults,minor_faults" > "$MEM_CSV"

    # Pick the largest dataset for memory measurement
    LARGEST_FILE=""
    LARGEST_SIZE=0
    for name in "${!DATASETS[@]}"; do
        file="${DATASETS[$name]}"
        sz=$(stat --printf='%s' "$file" 2>/dev/null || stat -f '%z' "$file" 2>/dev/null || echo 0)
        if (( sz > LARGEST_SIZE )); then
            LARGEST_SIZE=$sz
            LARGEST_FILE="$file"
        fi
    done

    if [[ -n "$LARGEST_FILE" ]]; then
        measure_memory "$LARGEST_FILE" "$MEM_CSV" >> "$MEM_CSV"
        echo "  Memory data written to $MEM_CSV"
    fi
    echo ""
    fi  # end: datasets available
fi

# ══════════════════════════════════════════════════════════════════════
#  2. File-Size Scaling
# ══════════════════════════════════════════════════════════════════════

if ! $SKIP_SCALING; then
    echo "──────────────────────────────────────────────────"
    echo " [2] File-Size Scaling"
    echo "──────────────────────────────────────────────────"

    SIZE_MBS=(1 5 10 25 50 100 250 500 1000)

    for mb in "${SIZE_MBS[@]}"; do
        f="$SCALING_DIR/size_${mb}mb.fasta"
        if [[ ! -f "$f" ]]; then
            echo "  Generating ${mb}MB file..."
            python3 -c "
mb=$mb; total=mb*1024*1024; line='ACGTACGTAC'*8
with open('$f','w') as fh:
    for i in range(1,101):
        fh.write(f'>seq{i}\n')
        r=total//100-50
        while r>0:
            c=min(80,r); fh.write(line[:c]+'\n'); r-=c
"
        fi
    done

    SCALE_SIZE_DIR="$RESULTS_DIR/scale_size_${TIMESTAMP}"
    mkdir -p "$SCALE_SIZE_DIR"

    for mb in "${SIZE_MBS[@]}"; do
        f="$SCALING_DIR/size_${mb}mb.fasta"
        echo "  ${mb}MB..."
        bench_file "$f" "$SCALE_SIZE_DIR/${mb}mb.json"
    done
    echo ""

    # ── 3. Sequence-Count Scaling ──────────────────────────────────
    echo "──────────────────────────────────────────────────"
    echo " [3] Sequence-Count Scaling"
    echo "──────────────────────────────────────────────────"

    SEQ_COUNTS=(10 100 1000 10000 100000)

    for count in "${SEQ_COUNTS[@]}"; do
        f="$SCALING_DIR/seqs_${count}.fasta"
        if [[ ! -f "$f" ]]; then
            echo "  Generating ${count} sequences..."
            python3 -c "
count=$count; total=50*1024*1024; seq_len=max(80,total//count-20); line='ACGTACGTAC'*8
with open('$f','w') as fh:
    for i in range(1,count+1):
        fh.write(f'>seq{i}\n')
        r=seq_len
        while r>0:
            c=min(80,r); fh.write(line[:c]+'\n'); r-=c
"
        fi
    done

    SCALE_SEQ_DIR="$RESULTS_DIR/scale_seqs_${TIMESTAMP}"
    mkdir -p "$SCALE_SEQ_DIR"

    for count in "${SEQ_COUNTS[@]}"; do
        f="$SCALING_DIR/seqs_${count}.fasta"
        echo "  ${count} seqs..."
        bench_file "$f" "$SCALE_SEQ_DIR/${count}.json"
    done
    echo ""
fi

echo "════════════════════════════════════════════════════════════════"
echo "  Benchmarks complete."
echo "  Results in: $RESULTS_DIR"
echo "  Run: .venv/bin/python bench/index/generate_report.py"
echo "════════════════════════════════════════════════════════════════"
