#!/bin/bash
# bench/wrappers/sanity_check.sh
# Quick sanity checks for Tier 2 wrappers + z-fasta against the Genome FASTA.
# Verifies: index creation, region extraction, output matches samtools.
#
# Usage: bash bench/wrappers/sanity_check.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

ZFASTA="$REPO_ROOT/zig-out/bin/z-fasta"
HTSLIB="$REPO_ROOT/tools/htslib_wrapper/htslib_wrapper"
NOODLES="$REPO_ROOT/tools/noodles_wrapper/target/release/noodles_wrapper"
RUSTBIO="$REPO_ROOT/tools/rustbio_wrapper/target/release/rustbio_wrapper"
SAMTOOLS="samtools"

GENOME="$REPO_ROOT/bench/shared/data/REAL_Genome.fa"
TEST_REGION="1:10000-10099"

RED='\033[0;31m'
GREEN='\033[0;32m'
CYAN='\033[0;36m'
BOLD='\033[1m'
NC='\033[0m'

ok()   { echo -e "${GREEN}  OK${NC} $*"; }
fail() { echo -e "${RED}  FAIL${NC} $*"; }

errors=0

echo -e "${BOLD}=== Tier 2 + z-fasta Sanity Checks ===${NC}"
echo ""

# Check prerequisites
for tool_bin in "$ZFASTA" "$HTSLIB" "$NOODLES" "$RUSTBIO"; do
    if [[ ! -x "$tool_bin" ]]; then
        fail "Binary not found: $tool_bin"
        ((errors++))
    fi
done

if [[ ! -f "$GENOME" ]]; then
    fail "Genome FASTA not found: $GENOME"
    ((errors++))
fi

if ! command -v "$SAMTOOLS" &>/dev/null; then
    fail "samtools not found in PATH"
    ((errors++))
fi

if [[ "$errors" -gt 0 ]]; then
    echo ""
    echo "Prerequisites missing. Exiting."
    exit 1
fi

echo "Genome: $GENOME ($(du -h "$GENOME" | cut -f1))"
echo "Region: $TEST_REGION"
echo ""

# Get reference output from samtools
echo -e "${CYAN}--- Reference (samtools) ---${NC}"
REF_OUTPUT=$("$SAMTOOLS" faidx "$GENOME" "$TEST_REGION")
echo "$REF_OUTPUT"
echo ""

# Collect timing results for summary table
declare -A INDEX_TIMES
declare -A GET_TIMES
declare -A MATCH_STATUS

# Test each tool
for name_tool in "z-fasta:$ZFASTA" "htslib:$HTSLIB" "noodles:$NOODLES" "rustbio:$RUSTBIO"; do
    name="${name_tool%%:*}"
    tool="${name_tool##*:}"

    echo -e "${CYAN}--- $name ---${NC}"

    # Test index
    INDEX_START=$(date +%s%N)
    if "$tool" index "$GENOME" 2>&1; then
        INDEX_END=$(date +%s%N)
        INDEX_MS=$(( (INDEX_END - INDEX_START) / 1000000 ))
        INDEX_TIMES[$name]="$INDEX_MS"
        ok "index completed in ${INDEX_MS}ms"
    else
        fail "index failed"
        INDEX_TIMES[$name]="FAIL"
        ((errors++))
        echo ""
        continue
    fi

    # Test get
    GET_START=$(date +%s%N)
    TOOL_OUTPUT=$("$tool" get "$GENOME" "$TEST_REGION" 2>&1)
    GET_EXIT=$?
    GET_END=$(date +%s%N)
    GET_MS=$(( (GET_END - GET_START) / 1000000 ))

    if [[ "$GET_EXIT" -ne 0 ]]; then
        fail "get failed (exit=$GET_EXIT)"
        GET_TIMES[$name]="FAIL"
        MATCH_STATUS[$name]="FAIL"
        ((errors++))
    else
        GET_TIMES[$name]="$GET_MS"
        # Compare with samtools output (sequence only, ignore header)
        TOOL_SEQ=$(echo "$TOOL_OUTPUT" | grep -v '^>' | tr -d '\n')
        REF_SEQ=$(echo "$REF_OUTPUT" | grep -v '^>' | tr -d '\n')

        if [[ "$TOOL_SEQ" == "$REF_SEQ" ]]; then
            ok "get matches samtools (${#TOOL_SEQ} bp, ${GET_MS}ms)"
            MATCH_STATUS[$name]="PASS"
        else
            fail "get output differs from samtools"
            echo "    Expected length: ${#REF_SEQ}"
            echo "    Got length:      ${#TOOL_SEQ}"
            MATCH_STATUS[$name]="DIFF"
            ((errors++))
        fi
    fi
    echo ""
done

# Summary table
echo -e "${BOLD}=== Timing Summary ===${NC}"
echo ""
printf "%-12s  %-12s  %-12s  %-8s\n" "Tool" "Index (ms)" "Get (ms)" "Match"
printf "%-12s  %-12s  %-12s  %-8s\n" "------------" "------------" "------------" "--------"

for name in "z-fasta" "htslib" "noodles" "rustbio"; do
    idx="${INDEX_TIMES[$name]:-N/A}"
    getv="${GET_TIMES[$name]:-N/A}"
    match="${MATCH_STATUS[$name]:-N/A}"
    if [[ "$idx" != "FAIL" && "$getv" != "FAIL" ]]; then
        printf "%-12s  %-12s  %-12s  %-8s\n" "$name" "${idx}" "${getv}" "${match}"
    else
        printf "%-12s  %-12s  %-12s  %-8s\n" "$name" "FAIL" "FAIL" "FAIL"
    fi
done

echo ""
echo -e "${BOLD}=== Summary ===${NC}"
if [[ "$errors" -eq 0 ]]; then
    echo -e "${GREEN}All checks passed.${NC}"
else
    echo -e "${RED}${errors} check(s) failed.${NC}"
    exit 1
fi
