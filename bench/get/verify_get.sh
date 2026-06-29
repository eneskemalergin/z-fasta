#!/usr/bin/env bash
# verify_get.sh - Verify z-fasta get output is byte-identical to samtools faidx
# Targets: non-REAL test files in tests/data/ only
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BENCH_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
PROJECT_DIR="$(cd "$BENCH_ROOT/.." && pwd)"
TMPDIR_LOCAL="$SCRIPT_DIR/.verify_tmp"
MESSY_DIR="$PROJECT_DIR/bench/index/messy_variants"
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

verify_expected() {
    local fasta="$1" region="$2" desc="$3" expected="$4"
    printf "%s" "$expected" > "$TMPDIR_LOCAL/expected.tmp"
    "$ZFASTA" get "$fasta" "$region" > "$TMPDIR_LOCAL/zf.tmp" 2>/dev/null || { fail "$desc (z-fasta err)"; return; }
    diff -q "$TMPDIR_LOCAL/expected.tmp" "$TMPDIR_LOCAL/zf.tmp" > /dev/null 2>&1 && pass "$desc" || fail "$desc"
}

write_expected_region() {
    local fasta="$1" name="$2" start="$3" end="$4" out="$5"
    python - "$fasta" "$name" "$start" "$end" "$out" <<'PY'
from pathlib import Path
import sys

fasta, name, start, end, out = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
seqs = {}
current = None
for raw in Path(fasta).read_text().splitlines():
    if raw.startswith(">"):
        current = raw[1:].split()[0]
        seqs[current] = []
    elif current is not None:
        seqs[current].append("".join(ch for ch in raw if not ch.isspace()))

seq = "".join(seqs[name])[start - 1:end]
header = f">{name}:{start}-{end}"
Path(out).write_text(header + "\n" + seq + "\n")
PY
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

test_mixed_width_fixture() {
    local fasta="$TMPDIR_LOCAL/mixed_widths.fasta"
    cp tests/data/mixed_widths.fasta "$fasta"
    "$SAMTOOLS" faidx "$fasta" >/dev/null 2>&1 || { echo "SKIP mixed_widths temp fixture (samtools faidx failed)"; return; }
    "$ZFASTA" index "$fasta" >/dev/null 2>&1 || { fail "mixed_widths temp fixture (z-fasta index err)"; return; }
    test_file "$fasta"
}

test_non_uniform_fixture() {
    local source_fasta="$1" name="$2" start="$3" end="$4" desc="$5"
    local fasta="$TMPDIR_LOCAL/${name}_${start}_${end}.fasta"
    cp "$source_fasta" "$fasta"
    "$ZFASTA" index "$fasta" >/dev/null 2>&1 || { fail "$desc (z-fasta index err)"; return; }
    write_expected_region "$fasta" "$name" "$start" "$end" "$TMPDIR_LOCAL/expected.tmp"
    "$ZFASTA" get "$fasta" "${name}:${start}-${end}" > "$TMPDIR_LOCAL/zf.tmp" 2>/dev/null || { fail "$desc (z-fasta err)"; return; }
    diff -q "$TMPDIR_LOCAL/expected.tmp" "$TMPDIR_LOCAL/zf.tmp" > /dev/null 2>&1 && pass "$desc" || fail "$desc"
}

echo "z-fasta get verification against samtools faidx"
echo "================================================"

cd "$PROJECT_DIR"
for f in \
    tests/data/simple.fasta \
    tests/data/proteome.fasta \
    tests/data/single.fasta \
    tests/data/edge_cases.fasta; do
    [ -f "$f" ] && test_file "$f"
done
test_mixed_width_fixture
test_non_uniform_fixture "$MESSY_DIR/mixed_line_widths.fasta" mixed_line_widths 3 24 "non-uniform mixed-width region"
test_non_uniform_fixture "$MESSY_DIR/trailing_whitespace.fasta" trailing_whitespace 7 20 "non-uniform trailing-whitespace region"
test_non_uniform_fixture "$MESSY_DIR/blank_lines.fasta" blank_lines 5 16 "non-uniform blank-line region"
test_non_uniform_fixture "$MESSY_DIR/mixed_crlf_lf.fasta" mixed_crlf_lf 3 18 "non-uniform mixed-CRLF region"

echo ""
echo "================================================"
echo "Results: $PASS passed, $FAIL failed"
echo "================================================"
[ "$FAIL" -gt 0 ] && { echo "VERIFICATION FAILED"; exit 1; }
echo "ALL PASSED"
