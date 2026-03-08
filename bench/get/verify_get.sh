#!/usr/bin/env bash
# verify_get.sh — Verify z-fasta get output is byte-identical to samtools faidx
# Targets: non-REAL test files in tests/data/ only
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BENCH_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
PROJECT_DIR="$(cd "$BENCH_ROOT/.." && pwd)"
TMPDIR_LOCAL="$SCRIPT_DIR/.verify_tmp"
mkdir -p "$TMPDIR_LOCAL"
trap 'rm -rf "$TMPDIR_LOCAL"' EXIT

ZFASTA="${ZFASTA:-$PROJECT_DIR/zig-out/bin/z-fasta}"
SAMTOOLS="${SAMTOOLS:-samtools}"
PASS=0
FAIL=0

pass() { PASS=$((PASS + 1)); printf "  PASS: %s\n" "$1"; }
fail() {
    FAIL=$((FAIL + 1)); printf "  FAIL: %s\n" "$1"
    echo "    samtools:"
    head -3 "$TMPDIR_LOCAL/st.tmp" 2>/dev/null | sed 's/^/      /'
    echo "    z-fasta:"
    head -3 "$TMPDIR_LOCAL/zf.tmp" 2>/dev/null | sed 's/^/      /'
}

verify() {
    local fasta="$1" region="$2" desc="$3"
    "$SAMTOOLS" faidx "$fasta" "$region" > "$TMPDIR_LOCAL/st.tmp" 2>/dev/null || { fail "$desc (samtools err)"; return; }
    "$ZFASTA" get "$fasta" "$region" > "$TMPDIR_LOCAL/zf.tmp" 2>/dev/null || { fail "$desc (z-fasta err)"; return; }
    diff -q "$TMPDIR_LOCAL/st.tmp" "$TMPDIR_LOCAL/zf.tmp" > /dev/null 2>&1 && pass "$desc" || fail "$desc"
}

test_sequence() {
    local fasta="$1" name="$2" length="$3"
    local label="${fasta##*/}:${name}"

    # 1. Full sequence
    verify "$fasta" "$name" "$label full"
    # 2. First 10 bases
    if [ "$length" -ge 10 ]; then
        verify "$fasta" "${name}:1-10" "$label :1-10"
    else
        verify "$fasta" "${name}:1-${length}" "$label :1-${length}"
    fi
    # 3. Last 10 bases
    if [ "$length" -ge 10 ]; then
        local ls=$((length - 9))
        verify "$fasta" "${name}:${ls}-${length}" "$label :${ls}-${length}"
    fi
    # 4. Middle 101 bases
    if [ "$length" -ge 101 ]; then
        local m=$((length / 2)) ms=$((length / 2 - 50)) me=$((length / 2 + 50))
        [ "$ms" -lt 1 ] && ms=1
        [ "$me" -gt "$length" ] && me=$length
        verify "$fasta" "${name}:${ms}-${me}" "$label mid-span"
    fi
    # 5. Single base
    verify "$fasta" "${name}:1-1" "$label :1-1"
    # 6. Last base
    verify "$fasta" "${name}:${length}-${length}" "$label :${length}-${length}"
    # 7. Full via region syntax
    verify "$fasta" "${name}:1-${length}" "$label :1-${length}"
    # 8. Clamp: END > seq_len
    local over=$((length + 100))
    verify "$fasta" "${name}:1-${over}" "$label :1-${over} (clamp)"
}

test_file() {
    local fasta="$1"
    [ -f "${fasta}.fai" ] || { echo "SKIP $fasta (no .fai)"; return; }
    echo ""
    echo "=== ${fasta##*/} ==="
    while IFS=$'\t' read -r name length _offset _lb _lby _rest; do
        [ "$length" -gt 0 ] && test_sequence "$fasta" "$name" "$length"
    done < "${fasta}.fai"
}

echo "z-fasta get verification against samtools faidx"
echo "================================================"

cd "$PROJECT_DIR"
for f in \
    tests/data/simple.fasta \
    tests/data/proteome.fasta \
    tests/data/single.fasta \
    tests/data/edge_cases.fasta \
    tests/data/mixed_widths.fasta; do
    [ -f "$f" ] && test_file "$f"
done

echo ""
echo "================================================"
echo "Results: $PASS passed, $FAIL failed"
echo "================================================"
[ "$FAIL" -gt 0 ] && { echo "VERIFICATION FAILED"; exit 1; }
echo "ALL PASSED"
