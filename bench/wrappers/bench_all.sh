#!/bin/bash
# bench/wrappers/bench_all.sh
# Hyperfine benchmarks for Tier 2 tools across all datasets.
# Tests: index and get on Genome (3GB), Transcriptome (459MB), Proteome (14MB).
#
# Tools:
#   z-fasta   - our tool (SIMD + mmap)
#   noodles   - noodles-fasta 0.61 (library indexer + query)
#   rustbio   - rust-bio 2.3 (custom index, library fetch)
#   samtools  - system samtools 1.13 (htslib 1.13)
#
# Note: bedtools getfasta is benchmarked separately in bench/get/run_benchmarks.sh
# for BED batch extraction.
#
# Usage: bash bench/wrappers/bench_all.sh [--runs N]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

ZFASTA="$REPO_ROOT/zig-out/bin/z-fasta"
NOODLES="$REPO_ROOT/tools/noodles_wrapper/target/release/noodles_wrapper"
RUSTBIO="$REPO_ROOT/tools/rustbio_wrapper/target/release/rustbio_wrapper"
SAMTOOLS="samtools"

DATA_DIR="$REPO_ROOT/bench/shared/data"
RESULTS_DIR="$SCRIPT_DIR/results"
mkdir -p "$RESULTS_DIR"

RUNS=5

while [[ $# -gt 0 ]]; do
    case $1 in
        --runs) RUNS="$2"; shift 2 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

# Verify all tools exist
for bin in "$ZFASTA" "$NOODLES" "$RUSTBIO"; do
    [[ -x "$bin" ]] || { echo "Missing: $bin"; exit 1; }
done
command -v "$SAMTOOLS" &>/dev/null || { echo "samtools not found"; exit 1; }
command -v hyperfine &>/dev/null || { echo "hyperfine not found"; exit 1; }

TIMESTAMP=$(date +%Y%m%d_%H%M%S)

echo "================================================================"
echo " Tier 2 Benchmark Suite - $TIMESTAMP"
echo " Runs per tool: $RUNS"
echo "================================================================"
echo ""
echo "Tools:"
echo "  z-fasta   - SIMD + mmap (.zfi index)"
echo "  noodles   - noodles-fasta 0.61 (library indexer)"
echo "  rustbio   - rust-bio 2.3 (custom indexer)"
echo "  samtools  - samtools 1.13 / htslib 1.13"
echo ""

# ── Dataset definitions ──────────────────────────────────────────────
# name:fasta:region
DATASETS=(
    "Genome:$DATA_DIR/REAL_Genome.fa:1:10000-10099"
    "Transcriptome:$DATA_DIR/REAL_Transcriptome.fa:ENST00000456328.2|ENSG00000290825.1|-|OTTHUMT00000362751.1|DDX11L2-202|DDX11L2|1657|lncRNA|:1-100"
    "Proteome:$DATA_DIR/REAL_Proteome.fasta:sp|A0A0J9YXV3|GREP1_HUMAN:1-100"
)

# ── Benchmark functions ──────────────────────────────────────────────

bench_index() {
    local name="$1"
    local fasta="$2"

    echo ""
    echo "── INDEX: $name ──"
    echo "   File: $fasta ($(du -h "$fasta" | cut -f1))"

    # Prepare clean state: remove all index files
    rm -f "${fasta}.zfi" "${fasta}.fai"

    hyperfine \
        --runs "$RUNS" \
        --warmup 1 \
        --export-json "$RESULTS_DIR/index_${name}_${TIMESTAMP}.json" \
        --export-markdown "$RESULTS_DIR/index_${name}_${TIMESTAMP}.md" \
        -n "z-fasta"   "$ZFASTA index $fasta" \
        -n "noodles"   "$NOODLES index $fasta" \
        -n "rustbio"   "$RUSTBIO index $fasta" \
        -n "samtools"  "$SAMTOOLS faidx $fasta" \
        2>&1

    echo ""
}

bench_get() {
    local name="$1"
    local fasta="$2"
    local region="$3"

    echo ""
    echo "── GET: $name ──"
    echo "   Region: $region"

    # Ensure all index files exist
    [[ -f "${fasta}.fai" ]] || "$SAMTOOLS" faidx "$fasta"
    [[ -f "${fasta}.zfi" ]] || "$ZFASTA" index "$fasta" 2>/dev/null

    hyperfine \
        --runs "$RUNS" \
        --warmup 1 \
        --export-json "$RESULTS_DIR/get_${name}_${TIMESTAMP}.json" \
        --export-markdown "$RESULTS_DIR/get_${name}_${TIMESTAMP}.md" \
        -n "z-fasta"   "$ZFASTA get '$fasta' '$region'" \
        -n "noodles"   "$NOODLES get '$fasta' '$region'" \
        -n "rustbio"   "$RUSTBIO get '$fasta' '$region'" \
        -n "samtools"  "$SAMTOOLS faidx '$fasta' '$region'" \
        2>&1

    echo ""
}

# ── Verify correctness first ─────────────────────────────────────────
echo "── Correctness verification ──"
errors=0
for entry in "${DATASETS[@]}"; do
    IFS=':' read -r name fasta region <<< "$entry"
    echo -n "  $name: "

    # Ensure index exists
    [[ -f "${fasta}.fai" ]] || "$SAMTOOLS" faidx "$fasta" 2>/dev/null
    [[ -f "${fasta}.zfi" ]] || "$ZFASTA" index "$fasta" 2>/dev/null

    ref=$("$SAMTOOLS" faidx "$fasta" "$region" 2>/dev/null | grep -v '^>' | tr -d '\n')
    zf=$("$ZFASTA" get "$fasta" "$region" 2>/dev/null | grep -v '^>' | tr -d '\n')
    no=$("$NOODLES" get "$fasta" "$region" 2>/dev/null | grep -v '^>' | tr -d '\n')
    rb=$("$RUSTBIO" get "$fasta" "$region" 2>/dev/null | grep -v '^>' | tr -d '\n')

    all_match=true
    for seq in "$zf" "$no" "$rb"; do
        [[ "$seq" == "$ref" ]] || { all_match=false; break; }
    done

    if $all_match; then
        echo "OK (${#ref} bp)"
    else
        echo "MISMATCH"
        ((errors++))
    fi
done

if [[ "$errors" -gt 0 ]]; then
    echo ""
    echo "Correctness check failed. Aborting benchmarks."
    exit 1
fi
echo ""

# ── Run benchmarks ────────────────────────────────────────────────────
for entry in "${DATASETS[@]}"; do
    IFS=':' read -r name fasta region <<< "$entry"
    bench_index "$name" "$fasta"
    bench_get "$name" "$fasta" "$region"
done

# ── Summary ───────────────────────────────────────────────────────────
echo ""
echo "================================================================"
echo " Results saved to: $RESULTS_DIR/"
echo "================================================================"
echo ""

# Print markdown summaries
for f in "$RESULTS_DIR"/*_${TIMESTAMP}.md; do
    [[ -f "$f" ]] || continue
    echo "$(basename "$f" .md)"
    cat "$f"
    echo ""
done
