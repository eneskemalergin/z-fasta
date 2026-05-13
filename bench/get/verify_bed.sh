#!/usr/bin/env bash
# verify_bed.sh — Verify z-fasta --bed output against bedtools getfasta.
#
# Coverage:
#   - chunked BED extraction for small / medium / large / x-large inputs
#   - default BED extraction and strand-aware extraction (-s / --honor-strand)
#   - comments / track lines / duplicate intervals / overlapping intervals
#   - stdin BED path for z-fasta
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
TMPDIR_LOCAL="$SCRIPT_DIR/.verify_bed_tmp"
mkdir -p "$TMPDIR_LOCAL"
trap 'rm -rf "$TMPDIR_LOCAL"' EXIT

ZFASTA="${ZFASTA:-$PROJECT_DIR/zig-out/bin/z-fasta}"
BEDTOOLS="${BEDTOOLS:-bedtools}"
FASTA="${BED_FASTA:-tests/data/simple.fasta}"

SMALL_COUNT="${BED_VERIFY_SMALL:-10}"
MEDIUM_COUNT="${BED_VERIFY_MEDIUM:-1000}"
LARGE_COUNT="${BED_VERIFY_LARGE:-10000}"
XLARGE_COUNT="${BED_VERIFY_XLARGE:-100000}"

SMALL_CHUNK="${BED_VERIFY_SMALL_CHUNK:-3}"
MEDIUM_CHUNK="${BED_VERIFY_MEDIUM_CHUNK:-97}"
LARGE_CHUNK="${BED_VERIFY_LARGE_CHUNK:-257}"
XLARGE_CHUNK="${BED_VERIFY_XLARGE_CHUNK:-4096}"

PASS=0
FAIL=0

pass() { PASS=$((PASS + 1)); printf "  PASS: %s\n" "$1"; }
fail() {
    FAIL=$((FAIL + 1))
    printf "  FAIL: %s\n" "$1"
    echo "    expected (bedtools normalized):"
    head -5 "$TMPDIR_LOCAL/bedtools.norm" 2>/dev/null | sed 's/^/      /'
    echo "    got (z-fasta normalized):"
    head -5 "$TMPDIR_LOCAL/zf.norm" 2>/dev/null | sed 's/^/      /'
}

require_inputs() {
    if [[ ! -x "$ZFASTA" ]]; then
        echo "SKIP: z-fasta binary not found at $ZFASTA"
        echo "      Build it first with: ./zig build -Doptimize=ReleaseFast"
        exit 0
    fi

    if ! command -v "$BEDTOOLS" >/dev/null 2>&1; then
        echo "SKIP: bedtools not found in PATH"
        echo "      Install bedtools to run BED verification against getfasta"
        exit 0
    fi

    if [[ ! -f "$PROJECT_DIR/$FASTA" ]]; then
        echo "SKIP: FASTA file not found at $PROJECT_DIR/$FASTA"
        exit 0
    fi
}

generate_bed() {
    local out="$1" count="$2"

    : > "$out"
    printf '# synthetic BED for %s rows\n' "$count" >> "$out"
    printf 'track name=verify_bed_%s\n' "$count" >> "$out"

    # Fixed edge cases first.
    cat >> "$out" <<'EOF'
seq1	0	24	full_seq1	0	+
seq1	0	5	dup_a	0	+
seq1	0	5	dup_b	0	+
seq1	2	8	overlap_minus	0	-
seq2	0	12	full_seq2	0	-
seq2	0	4	short_plus	0	+
EOF

    local current=6
    while (( current < count )); do
        local idx="$current"
        local chrom len start span end strand name

        if (( idx % 2 == 0 )); then
            chrom="seq1"
            len=24
        else
            chrom="seq2"
            len=12
        fi

        if (( len <= 4 )); then
            span=1
        else
            span=$(( (idx % 4) + 1 ))
        fi
        start=$(( idx % (len - span + 1) ))
        end=$(( start + span ))
        if (( idx % 3 == 0 )); then
            strand='-'
        else
            strand='+'
        fi
        name="${chrom}_${start}_${end}_${strand}_${idx}"

        printf '%s\t%d\t%d\t%s\t0\t%s\n' "$chrom" "$start" "$end" "$name" "$strand" >> "$out"
        current=$((current + 1))
    done
}

normalize_fasta_against_bed() {
    local fasta_in="$1" bed_in="$2" honor_strand="$3" fasta_out="$4"

    awk -v bed_file="$bed_in" -v honor_strand="$honor_strand" '
        function flush_record(    header, seq_copy) {
            if (record_idx == 0) return
            if (!(record_idx in chrom)) {
                printf("record count exceeds BED rows at FASTA record %d\n", record_idx) > "/dev/stderr"
                exit 1
            }
            header = chrom[record_idx] ":" start1[record_idx] "-" end1[record_idx]
            if (honor_strand == "1" && strand[record_idx] == "-") {
                header = header ":rc"
            }
            print ">" header
            seq_copy = seq
            while (length(seq_copy) > 60) {
                print substr(seq_copy, 1, 60)
                seq_copy = substr(seq_copy, 61)
            }
            print seq_copy
        }
        BEGIN {
            FS = "\t"
            while ((getline < bed_file) > 0) {
                sub(/\r$/, "", $0)
                if ($0 == "" || $0 ~ /^#/ || $0 ~ /^track/ || $0 ~ /^browser/) continue
                bed_count += 1
                chrom[bed_count] = $1
                start1[bed_count] = $2 + 1
                end1[bed_count] = $3
                strand[bed_count] = (NF >= 6 ? $6 : ".")
            }
            close(bed_file)
            record_idx = 0
            seq = ""
        }
        /^>/ {
            flush_record()
            record_idx += 1
            seq = ""
            next
        }
        {
            gsub(/[[:space:]]/, "", $0)
            seq = seq $0
        }
        END {
            flush_record()
            if (record_idx != bed_count) {
                printf("record count mismatch: fasta=%d bed=%d\n", record_idx, bed_count) > "/dev/stderr"
                exit 1
            }
        }
    ' "$fasta_in" > "$fasta_out"
}

verify_bed_case() {
    local desc="$1" fasta="$2" bed="$3" chunk_size="$4" honor_strand="$5" use_stdin="$6"
    local bedtools_args=()
    local zf_args=(get "$fasta")

    if [[ "$honor_strand" == "1" ]]; then
        bedtools_args+=(-s)
        zf_args+=(--honor-strand)
    fi
    zf_args+=(--chunk-size "$chunk_size")

    if [[ "$use_stdin" == "1" ]]; then
        zf_args+=(--bed -)
        "$BEDTOOLS" getfasta -fi "$fasta" -bed "$bed" "${bedtools_args[@]}" > "$TMPDIR_LOCAL/bedtools.raw" 2>/dev/null \
            || { fail "$desc (bedtools err)"; return; }
        cat "$bed" | "$ZFASTA" "${zf_args[@]}" > "$TMPDIR_LOCAL/zf.raw" 2>/dev/null \
            || { fail "$desc (z-fasta err)"; return; }
    else
        zf_args+=(--bed "$bed")
        "$BEDTOOLS" getfasta -fi "$fasta" -bed "$bed" "${bedtools_args[@]}" > "$TMPDIR_LOCAL/bedtools.raw" 2>/dev/null \
            || { fail "$desc (bedtools err)"; return; }
        "$ZFASTA" "${zf_args[@]}" > "$TMPDIR_LOCAL/zf.raw" 2>/dev/null \
            || { fail "$desc (z-fasta err)"; return; }
    fi

    normalize_fasta_against_bed "$TMPDIR_LOCAL/bedtools.raw" "$bed" "$honor_strand" "$TMPDIR_LOCAL/bedtools.norm" \
        || { fail "$desc (bedtools normalize err)"; return; }
    normalize_fasta_against_bed "$TMPDIR_LOCAL/zf.raw" "$bed" "$honor_strand" "$TMPDIR_LOCAL/zf.norm" \
        || { fail "$desc (z-fasta normalize err)"; return; }

    if diff -q "$TMPDIR_LOCAL/bedtools.norm" "$TMPDIR_LOCAL/zf.norm" >/dev/null 2>&1; then
        pass "$desc"
    else
        fail "$desc"
    fi
}

run_size_suite() {
    local label="$1" count="$2" chunk_size="$3"
    local bed_file="$TMPDIR_LOCAL/${label}.bed"

    generate_bed "$bed_file" "$count"
    echo ""
    echo "=== ${label} (${count} BED rows, chunk-size=${chunk_size}) ==="
    verify_bed_case "${label} default BED file" "$FASTA" "$bed_file" "$chunk_size" 0 0
    verify_bed_case "${label} stranded BED file" "$FASTA" "$bed_file" "$chunk_size" 1 0
}

echo "z-fasta BED verification against bedtools getfasta"
echo "================================================"

cd "$PROJECT_DIR"
require_inputs

run_size_suite "small" "$SMALL_COUNT" "$SMALL_CHUNK"
verify_bed_case "small default via stdin" "$FASTA" "$TMPDIR_LOCAL/small.bed" "$SMALL_CHUNK" 0 1
verify_bed_case "small stranded via stdin" "$FASTA" "$TMPDIR_LOCAL/small.bed" "$SMALL_CHUNK" 1 1
run_size_suite "medium" "$MEDIUM_COUNT" "$MEDIUM_CHUNK"
run_size_suite "large" "$LARGE_COUNT" "$LARGE_CHUNK"
run_size_suite "x-large" "$XLARGE_COUNT" "$XLARGE_CHUNK"

echo ""
echo "================================================"
echo "Results: $PASS passed, $FAIL failed"
echo "================================================"
[ "$FAIL" -gt 0 ] && { echo "VERIFICATION FAILED"; exit 1; }
echo "ALL PASSED"