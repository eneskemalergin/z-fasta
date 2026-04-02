#!/usr/bin/env bash
# verify_multi_get.sh — Verify z-fasta multi-region get output is byte-identical
# to samtools faidx with the same arguments.
#
# Tests:
#   - 2 regions on the same sequence
#   - 2 regions on different sequences
#   - Reversed CLI order (seq2 before seq1 by file position)
#   - Overlapping regions
#   - Duplicate regions (same region twice)
#   - Mixed full-sequence + sub-region
#   - 10 regions (stays under 16 → direct path)
#   - 20 regions (>= 16 → offset-sort path)
#   - All sequences from a file
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
TMPDIR_LOCAL="$SCRIPT_DIR/.verify_multi_tmp"
mkdir -p "$TMPDIR_LOCAL"
trap 'rm -rf "$TMPDIR_LOCAL"' EXIT

ZFASTA="${ZFASTA:-$PROJECT_DIR/zig-out/bin/z-fasta}"
SAMTOOLS="${SAMTOOLS:-samtools}"
PASS=0
FAIL=0

pass() { PASS=$((PASS + 1)); printf "  PASS: %s\n" "$1"; }
fail() {
    FAIL=$((FAIL + 1)); printf "  FAIL: %s\n" "$1"
    echo "    expected (samtools):"
    head -5 "$TMPDIR_LOCAL/st.tmp" 2>/dev/null | sed 's/^/      /'
    echo "    got (z-fasta):"
    head -5 "$TMPDIR_LOCAL/zf.tmp" 2>/dev/null | sed 's/^/      /'
}

# verify_multi <desc> <fasta> region1 region2 ...
verify_multi() {
    local desc="$1"
    local fasta="$2"
    shift 2
    local regions=("$@")

    "$SAMTOOLS" faidx "$fasta" "${regions[@]}" > "$TMPDIR_LOCAL/st.tmp" 2>/dev/null \
        || { fail "$desc (samtools err)"; return; }
    "$ZFASTA" get "$fasta" "${regions[@]}" > "$TMPDIR_LOCAL/zf.tmp" 2>/dev/null \
        || { fail "$desc (z-fasta err)"; return; }
    diff -q "$TMPDIR_LOCAL/st.tmp" "$TMPDIR_LOCAL/zf.tmp" > /dev/null 2>&1 \
        && pass "$desc" || fail "$desc"
}

echo "z-fasta multi-region get verification against samtools faidx"
echo "============================================================="

cd "$PROJECT_DIR"

# ── Helpers ──────────────────────────────────────────────────────────
ensure_index() {
    local f="$1"
    [[ -f "${f}.fai" ]] || "$SAMTOOLS" faidx "$f"
    [[ -f "${f}.zfi" ]] || "$ZFASTA" index "$f"
}

SIMPLE="tests/data/simple.fasta"
PROTEOME="tests/data/proteome.fasta"
EDGE="tests/data/edge_cases.fasta"

ensure_index "$SIMPLE"
ensure_index "$PROTEOME"
ensure_index "$EDGE"

echo ""
echo "=== simple.fasta (seq1=24bp, seq2=12bp) ==="

# 2 regions, same sequence (direct path, < 16)
verify_multi "two sub-regions same seq" "$SIMPLE" \
    "seq1:1-10" "seq1:13-24"

# 2 regions, different sequences, CLI order matches file order
verify_multi "two seqs, CLI order = file order" "$SIMPLE" \
    "seq1:1-12" "seq2:1-6"

# 2 regions, reversed CLI order (seq2 before seq1 by file position)
verify_multi "two seqs, reversed CLI order" "$SIMPLE" \
    "seq2:1-6" "seq1:1-12"

# Overlapping regions on same sequence
verify_multi "overlapping regions same seq" "$SIMPLE" \
    "seq1:1-15" "seq1:10-24"

# Duplicate identical region
verify_multi "duplicate region" "$SIMPLE" \
    "seq1:1-12" "seq1:1-12"

# Mixed full-sequence + sub-region
verify_multi "mixed full-seq and sub-region" "$SIMPLE" \
    "seq1" "seq2:3-10"

# Full sequences for both
verify_multi "both full sequences" "$SIMPLE" \
    "seq1" "seq2"

# Same region 3 times
verify_multi "triple duplicate" "$SIMPLE" \
    "seq1:1-5" "seq1:1-5" "seq1:1-5"

# 3 regions: full + 2 sub-regions
verify_multi "full + 2 sub-regions" "$SIMPLE" \
    "seq1" "seq2:1-6" "seq1:10-20"

# All sequences in CLI reverse order
verify_multi "all seqs reversed" "$SIMPLE" \
    "seq2" "seq1"

echo ""
echo "=== edge_cases.fasta ==="

# Long sequence name
LONGNAME="AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"

verify_multi "long header name + normal name" "$EDGE" \
    "single_line" "$LONGNAME"

verify_multi "long header name reversed" "$EDGE" \
    "$LONGNAME" "single_line"

verify_multi "lowercase + nonstandard" "$EDGE" \
    "lowercase" "nonstandard"

# Duplicate + lowercase
verify_multi "same sequence twice (different seqs)" "$EDGE" \
    "single_line:1-5" "single_line:3-8"

echo ""
echo "=== proteome.fasta (pipe-delimited names) ==="

SEQ1="sp|P12345|PROT_HUMAN"
SEQ2="sp|Q98765|ANOT_MOUSE"

verify_multi "two pipe-delimited names" "$PROTEOME" \
    "$SEQ1" "$SEQ2"

verify_multi "pipe names reversed" "$PROTEOME" \
    "$SEQ2" "$SEQ1"

verify_multi "pipe name sub-region + full" "$PROTEOME" \
    "${SEQ1}:1-10" "$SEQ2"

verify_multi "duplicate pipe-delimited name" "$PROTEOME" \
    "$SEQ1" "$SEQ1"

echo ""
echo "=== Offset-sort path (>= 16 regions) ==="

# Build 20 regions from simple.fasta by alternating seq1/seq2 sub-ranges
# This exercises the sort path (>= 16 regions)
REGIONS_20=()
for i in 1 2 3 4 5 6 7 8 9 10; do
    REGIONS_20+=("seq1:${i}-$((i))")
    REGIONS_20+=("seq2:$((i <= 12 ? i : 12))-$((i <= 12 ? i : 12))")
done

echo "  Testing 20 regions (sort path)..."
verify_multi "20 regions triggers sort path" "$SIMPLE" "${REGIONS_20[@]}"

# 20 regions in reverse file order to test that sort path preserves CLI order
REGIONS_20_REV=()
for i in 10 9 8 7 6 5 4 3 2 1; do
    REGIONS_20_REV+=("seq2:$((i <= 12 ? i : 12))-$((i <= 12 ? i : 12))")
    REGIONS_20_REV+=("seq1:${i}-$((i))")
done

verify_multi "20 regions reversed (sort path preserves CLI order)" "$SIMPLE" "${REGIONS_20_REV[@]}"

echo ""
echo "=== All sequences from edge_cases.fasta (6 sequences) ==="

ALL_EDGE=()
while IFS=$'\t' read -r name _rest; do
    ALL_EDGE+=("$name")
done < "${EDGE}.fai"

if [[ ${#ALL_EDGE[@]} -gt 0 ]]; then
    verify_multi "all edge_cases sequences" "$EDGE" "${ALL_EDGE[@]}"
    # Reversed
    REVERSED_EDGE=()
    for (( i=${#ALL_EDGE[@]}-1; i>=0; i-- )); do
        REVERSED_EDGE+=("${ALL_EDGE[$i]}")
    done
    verify_multi "all edge_cases reversed" "$EDGE" "${REVERSED_EDGE[@]}"
fi

echo ""
echo "============================================================="
echo "Results: $PASS passed, $FAIL failed"
echo "============================================================="
[ "$FAIL" -gt 0 ] && { echo "VERIFICATION FAILED"; exit 1; }
echo "ALL PASSED"
