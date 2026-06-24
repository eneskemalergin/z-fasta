#!/bin/bash
# z-fasta Edge Case / Correctness Tests
# Generates 20 edge case FASTA files, runs each tool, diffs outputs.
# Outputs raw CSV; all formatting done by generate_report.py.
#
# Usage: ./run_tests.sh
#
# Outputs:
#   results/tests_<timestamp>.csv  per-case exit codes + match status

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH_ROOT="$(dirname "$SCRIPT_DIR")"
PROJECT_ROOT="$(dirname "$BENCH_ROOT")"
RESULTS_DIR="$SCRIPT_DIR/results"
EDGE_DIR="$SCRIPT_DIR/edge_cases"

# ── Tools ──────────────────────────────────────────────────────────
ZFASTA="$PROJECT_ROOT/zig-out/bin/z-fasta"
SAMTOOLS="samtools"
SEQKIT="$PROJECT_ROOT/tools/seqkit"
FASTAHACK="$PROJECT_ROOT/tools/fastahack-1.0.0/fastahack"

HAS_SEQKIT=false;    [[ -x "$SEQKIT" ]]    && HAS_SEQKIT=true
HAS_FASTAHACK=false; [[ -x "$FASTAHACK" ]] && HAS_FASTAHACK=true

mkdir -p "$RESULTS_DIR" "$EDGE_DIR"

TIMESTAMP=$(date +%Y%m%d_%H%M%S)
CSV="$RESULTS_DIR/tests_${TIMESTAMP}.csv"

echo "════════════════════════════════════════════════════════════════"
echo "  z-fasta Edge Case Tests"
echo "════════════════════════════════════════════════════════════════"
echo ""

# ══════════════════════════════════════════════════════════════════════
#  Generate edge case files
# ══════════════════════════════════════════════════════════════════════
echo "Generating edge case files..."

# 1. Empty file
> "$EDGE_DIR/empty.fasta"

# 2. Single byte (just >)
echo -n ">" > "$EDGE_DIR/single_byte.fasta"

# 3. Header only, no newline
echo -n ">seq1" > "$EDGE_DIR/header_no_newline.fasta"

# 4. Header with newline but no sequence
echo ">seq1" > "$EDGE_DIR/header_only.fasta"

# 5. Very long header (10KB)
printf ">%s description\nACGT\n" "$(head -c 10000 < /dev/zero | tr '\0' 'A')" \
    > "$EDGE_DIR/long_header.fasta"

# 6. Very long single line sequence (1MB, no wrapping)
printf ">longseq\n%s\n" "$(head -c 1000000 < /dev/zero | tr '\0' 'A')" \
    > "$EDGE_DIR/no_wrap_1mb.fasta"

# 7. Sequence with special characters
printf ">special\nACGT*-NXYZ\n" > "$EDGE_DIR/special_chars.fasta"

# 8. CRLF line endings (Windows)
printf ">seq1\r\nACGT\r\nGGGG\r\n" > "$EDGE_DIR/crlf.fasta"

# 9. Mixed line endings
printf ">seq1\nACGT\r\nGGGG\n" > "$EDGE_DIR/mixed_endings.fasta"

# 10. No final newline
printf ">seq1\nACGT" > "$EDGE_DIR/no_final_newline.fasta"

# 11. Binary garbage in sequence
printf ">binary\nACGT\x00\x01\x02GGGG\n" > "$EDGE_DIR/binary_data.fasta"

# 12. Unicode in header
printf ">seq_émoji_🧬\nACGT\n" > "$EDGE_DIR/unicode_header.fasta"

# 13. Tab in sequence name
printf ">seq1\tmore\ttabs\nACGT\n" > "$EDGE_DIR/tab_in_name.fasta"

# 14. Multiple > on same line (not headers)
printf ">seq1 description with > symbols > here\nACGT\n" \
    > "$EDGE_DIR/gt_in_description.fasta"

# 15. > in middle of sequence
printf ">seq1\nACGT>GGGG\n" > "$EDGE_DIR/gt_in_sequence.fasta"

# 16. Space before >
printf " >seq1\nACGT\n" > "$EDGE_DIR/space_before_gt.fasta"

# 17. All bases are N
printf ">all_n\nNNNNNNNNNN\n" > "$EDGE_DIR/all_n.fasta"

# 18. Lowercase bases
printf ">lower\nacgtacgt\n" > "$EDGE_DIR/lowercase.fasta"

# 19. Pipe in name (proteome style)
printf ">sp|P12345|PROT_HUMAN description\nMKWVTFISLL\n" \
    > "$EDGE_DIR/pipe_name.fasta"

# 20. Duplicate sequence names
printf ">dup\nAAAA\n>dup\nCCCC\n>dup\nGGGG\n" > "$EDGE_DIR/triple_dup.fasta"

echo "  Generated $(ls "$EDGE_DIR"/*.fasta 2>/dev/null | wc -l) edge case files"
echo ""

# ══════════════════════════════════════════════════════════════════════
#  Run tests
# ══════════════════════════════════════════════════════════════════════
echo "test_case,zfasta_exit,samtools_exit,seqkit_exit,fastahack_exit,output_match" \
    > "$CSV"

FAIL_LOG="$RESULTS_DIR/failures.log"
> "$FAIL_LOG"  # truncate

pass=0
fail=0
total=0

for file in "$EDGE_DIR"/*.fasta; do
    [[ -f "$file" ]] || continue
    name=$(basename "$file" .fasta)
    total=$((total + 1))

    # ── Clean indexes ──────────────────────────────────────────────
    rm -f "${file}.fai" "${file}.zfi" 2>/dev/null || true

    # ── z-fasta ────────────────────────────────────────────────────
    zf_exit=0
    "$ZFASTA" index --emit-fai "$file" > /tmp/zf_edge_test.fai 2>/dev/null || zf_exit=$?

    # ── samtools ───────────────────────────────────────────────────
    rm -f "${file}.fai" 2>/dev/null || true
    sam_exit=0
    $SAMTOOLS faidx "$file" 2>/dev/null || sam_exit=$?

    # ── seqkit ─────────────────────────────────────────────────────
    seq_exit=0
    if $HAS_SEQKIT; then
        rm -f "${file}.fai" 2>/dev/null || true
        "$SEQKIT" faidx "$file" > /dev/null 2>&1 || seq_exit=$?
    fi

    # ── fastahack ──────────────────────────────────────────────────
    fh_exit=0
    if $HAS_FASTAHACK; then
        rm -f "${file}.fai" 2>/dev/null || true
        "$FASTAHACK" -i "$file" > /dev/null 2>&1 || fh_exit=$?
    fi

    # ── Compare z-fasta vs samtools output ─────────────────────────
    # Re-run samtools to get a clean .fai for comparison
    rm -f "${file}.fai" "${file}.zfi" 2>/dev/null || true
    $SAMTOOLS faidx "$file" 2>/dev/null || true

    if [[ $zf_exit -ne 0 && $sam_exit -ne 0 ]]; then
        # Both rejected; that's agreement
        match="MATCH"
    elif [[ $zf_exit -ne 0 || $sam_exit -ne 0 ]]; then
        # One succeeded, one failed
        match="DIFF"
        echo "--- FAILURE: $name (exit codes differ: zf=$zf_exit sam=$sam_exit) ---" >> "$FAIL_LOG"
    elif [[ -f "${file}.fai" ]]; then
        # Both succeeded; diff the .fai content
        if diff -q /tmp/zf_edge_test.fai "${file}.fai" &>/dev/null; then
            match="MATCH"
        else
            match="DIFF"
            echo "--- FAILURE: $name (output differs) ---" >> "$FAIL_LOG"
            diff -u "${file}.fai" /tmp/zf_edge_test.fai >> "$FAIL_LOG" || true
        fi
    else
        # Both succeeded but no samtools .fai (empty input?)
        if [[ ! -s /tmp/zf_edge_test.fai ]]; then
            match="MATCH"
        else
            match="DIFF"
            echo "--- FAILURE: $name (z-fasta produced output, samtools did not) ---" >> "$FAIL_LOG"
        fi
    fi

    echo "$name,$zf_exit,$sam_exit,$seq_exit,$fh_exit,$match" >> "$CSV"

    if [[ "$match" == "MATCH" ]]; then
        sym="✓"; pass=$((pass + 1))
    else
        sym="✗"; fail=$((fail + 1))
    fi
    echo "  $sym  $name  (zf=$zf_exit sam=$sam_exit seq=$seq_exit fh=$fh_exit) $match"

    # Clean up
    rm -f "${file}.fai" "${file}.zfi" /tmp/zf_edge_test.fai 2>/dev/null || true
done

echo ""
echo "════════════════════════════════════════════════════════════════"
echo "  Results: $pass / $total match   ($fail diff)"
echo "  CSV: $CSV"
if [[ -s "$FAIL_LOG" ]]; then
    echo "  Failures: $FAIL_LOG"
else
    rm -f "$FAIL_LOG"
fi
echo "════════════════════════════════════════════════════════════════"
