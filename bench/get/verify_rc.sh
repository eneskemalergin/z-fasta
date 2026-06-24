#!/usr/bin/env bash
# verify_rc.sh - Verify z-fasta RC / reverse / complement extraction behavior.
#
# Coverage:
#   - --rc for single regions against samtools faidx -i --mark-strand no
#   - --rc for names-file extraction against samtools faidx -i --mark-strand no
#   - synthetic IUPAC-heavy fixture, including lowercase u
#   - synthetic chromosome-like full-sequence fixture to cover long wrapped whole-sequence RC
#   - multi-region --rc output matches per-region concatenation, including sort-path cases
#   - manual checks for --annotate-rc, --reverse-only, and --complement-only
#   - BED composition with --strand-aware --rc against bedtools getfasta -s plus seqtk seq -r
#   - protein rejection for complement-based transforms
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
TMPDIR_LOCAL="$SCRIPT_DIR/.verify_rc_tmp"
mkdir -p "$TMPDIR_LOCAL"
trap 'rm -rf "$TMPDIR_LOCAL"' EXIT

ZFASTA="${ZFASTA:-$PROJECT_DIR/zig-out/bin/z-fasta}"
SAMTOOLS="${SAMTOOLS:-samtools}"
BEDTOOLS="${BEDTOOLS:-bedtools}"
SEQTK="${SEQTK:-$PROJECT_DIR/tools/seqtk/seqtk}"
PASS=0
FAIL=0

pass() {
    PASS=$((PASS + 1))
    printf "  PASS: %s\n" "$1"
}

fail_diff() {
    local desc="$1" expected_label="$2" expected_file="$3" got_label="$4" got_file="$5"
    FAIL=$((FAIL + 1))
    printf "  FAIL: %s\n" "$desc"
    echo "    ${expected_label}:"
    head -5 "$expected_file" 2>/dev/null | sed 's/^/      /'
    echo "    ${got_label}:"
    head -5 "$got_file" 2>/dev/null | sed 's/^/      /'
}

fail_msg() {
    local desc="$1" msg="$2"
    FAIL=$((FAIL + 1))
    printf "  FAIL: %s\n" "$desc"
    echo "      $msg"
}

require_inputs() {
    if [[ ! -x "$ZFASTA" ]]; then
        echo "SKIP: z-fasta binary not found at $ZFASTA"
        echo "      Build it first with: ./zig build -Doptimize=ReleaseFast"
        exit 0
    fi

    if ! command -v "$SAMTOOLS" >/dev/null 2>&1; then
        echo "SKIP: samtools not found in PATH"
        echo "      Install samtools to run reverse-complement verification"
        exit 0
    fi

    if ! "$SAMTOOLS" faidx --help 2>&1 | grep -q -- '--reverse-complement'; then
        echo "SKIP: installed samtools does not support faidx --reverse-complement"
        echo "      Need samtools with: faidx -i / --reverse-complement"
        exit 0
    fi

    if ! command -v "$BEDTOOLS" >/dev/null 2>&1; then
        echo "SKIP: bedtools not found in PATH"
        echo "      Install bedtools to run RC BED composition verification"
        exit 0
    fi

    if [[ ! -x "$SEQTK" ]]; then
        echo "SKIP: seqtk binary not found at $SEQTK"
        echo "      Build or point SEQTK to a usable seqtk binary"
        exit 0
    fi
}

write_iupac_fixture() {
    cat > "$TMPDIR_LOCAL/iupac.fasta" <<'EOF'
>iupac_all synthetic IUPAC coverage
ACGTURYSWKMBDHVNacgturyswkmbdhvnu
EOF

    "$ZFASTA" index "$TMPDIR_LOCAL/iupac.fasta" >/dev/null 2>&1 \
        || { fail_msg "synthetic IUPAC fixture setup" "z-fasta index failed for generated IUPAC FASTA"; exit 1; }
}

write_chrom_fixture() {
    local pattern="ACGTNRYWSKMBDHVacgtnrywskmbdhv"

    {
        echo ">chrSynthetic chromosome-like full-sequence RC coverage"
        awk -v pattern="$pattern" -v total=16384 '
            BEGIN {
                seq = ""
                while (length(seq) < total) {
                    seq = seq pattern
                }
                seq = substr(seq, 1, total)
                while (length(seq) > 0) {
                    print substr(seq, 1, 71)
                    seq = substr(seq, 72)
                }
            }
        '
    } > "$TMPDIR_LOCAL/chrom_like.fasta"

    "$ZFASTA" index "$TMPDIR_LOCAL/chrom_like.fasta" >/dev/null 2>&1 \
        || { fail_msg "synthetic chromosome fixture setup" "z-fasta index failed for generated chromosome-like FASTA"; exit 1; }
}

verify_rc_region() {
    local fasta="$1" region="$2" desc="$3"

    "$SAMTOOLS" faidx -i --mark-strand no "$fasta" "$region" > "$TMPDIR_LOCAL/st.tmp" 2>/dev/null \
        || { fail_msg "$desc (samtools err)" "samtools faidx failed"; return; }
    "$ZFASTA" get "$fasta" "$region" --rc > "$TMPDIR_LOCAL/zf.tmp" 2>/dev/null \
        || { fail_msg "$desc (z-fasta err)" "z-fasta get --rc failed"; return; }

    if diff -q "$TMPDIR_LOCAL/st.tmp" "$TMPDIR_LOCAL/zf.tmp" >/dev/null 2>&1; then
        pass "$desc"
    else
        fail_diff "$desc" "samtools" "$TMPDIR_LOCAL/st.tmp" "z-fasta" "$TMPDIR_LOCAL/zf.tmp"
    fi
}

verify_rc_names() {
    local fasta="$1" names_file="$2" desc="$3"

    "$SAMTOOLS" faidx -i --mark-strand no -r "$names_file" "$fasta" > "$TMPDIR_LOCAL/st.tmp" 2>/dev/null \
        || { fail_msg "$desc (samtools err)" "samtools faidx -r failed"; return; }
    "$ZFASTA" get "$fasta" --names "$names_file" --rc > "$TMPDIR_LOCAL/zf.tmp" 2>/dev/null \
        || { fail_msg "$desc (z-fasta err)" "z-fasta get --names --rc failed"; return; }

    if diff -q "$TMPDIR_LOCAL/st.tmp" "$TMPDIR_LOCAL/zf.tmp" >/dev/null 2>&1; then
        pass "$desc"
    else
        fail_diff "$desc" "samtools" "$TMPDIR_LOCAL/st.tmp" "z-fasta" "$TMPDIR_LOCAL/zf.tmp"
    fi
}

verify_exact_output() {
    local desc="$1" expected_file="$2"
    shift 2

    "$ZFASTA" "$@" > "$TMPDIR_LOCAL/zf.tmp" 2> "$TMPDIR_LOCAL/zf.err" \
        || { fail_msg "$desc (z-fasta err)" "z-fasta command failed"; return; }

    if diff -q "$expected_file" "$TMPDIR_LOCAL/zf.tmp" >/dev/null 2>&1; then
        pass "$desc"
    else
        fail_diff "$desc" "expected" "$expected_file" "z-fasta" "$TMPDIR_LOCAL/zf.tmp"
    fi
}

verify_multi_rc_concat() {
    local desc="$1" fasta="$2"
    shift 2
    local regions=("$@")

    : > "$TMPDIR_LOCAL/expected.tmp"
    for region in "${regions[@]}"; do
        "$ZFASTA" get "$fasta" "$region" --rc >> "$TMPDIR_LOCAL/expected.tmp" 2> "$TMPDIR_LOCAL/zf.err" \
            || { fail_msg "$desc (per-region z-fasta err)" "z-fasta single-region --rc failed"; return; }
    done

    "$ZFASTA" get "$fasta" "${regions[@]}" --rc > "$TMPDIR_LOCAL/zf.tmp" 2> "$TMPDIR_LOCAL/zf.err" \
        || { fail_msg "$desc (multi-region z-fasta err)" "z-fasta multi-region --rc failed"; return; }

    if diff -q "$TMPDIR_LOCAL/expected.tmp" "$TMPDIR_LOCAL/zf.tmp" >/dev/null 2>&1; then
        pass "$desc"
    else
        fail_diff "$desc" "per-region concatenated" "$TMPDIR_LOCAL/expected.tmp" "multi-region" "$TMPDIR_LOCAL/zf.tmp"
    fi
}

normalize_fasta_against_bed() {
    local fasta_in="$1" bed_in="$2" fasta_out="$3"

    awk -v bed_file="$bed_in" '
        function flush_record(    header, seq_copy) {
            if (record_idx == 0) return
            if (!(record_idx in chrom)) {
                printf("record count exceeds BED rows at FASTA record %d\n", record_idx) > "/dev/stderr"
                exit 1
            }
            header = chrom[record_idx] ":" start1[record_idx] "-" end1[record_idx]
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

verify_bed_rc_composed() {
    local desc="$1" fasta="$2" bed="$3" chunk_size="$4" use_stdin="$5"
    local zf_args=(get "$fasta" --strand-aware --rc --chunk-size "$chunk_size")

    "$BEDTOOLS" getfasta -fi "$fasta" -bed "$bed" -s > "$TMPDIR_LOCAL/bedtools.raw" 2>/dev/null \
        || { fail_msg "$desc (bedtools err)" "bedtools getfasta -s failed"; return; }
    "$SEQTK" seq -r "$TMPDIR_LOCAL/bedtools.raw" > "$TMPDIR_LOCAL/bedtools.transformed" 2>/dev/null \
        || { fail_msg "$desc (seqtk err)" "seqtk seq -r failed"; return; }

    if [[ "$use_stdin" == "1" ]]; then
        zf_args+=(--bed -)
        cat "$bed" | "$ZFASTA" "${zf_args[@]}" > "$TMPDIR_LOCAL/zf.raw" 2>/dev/null \
            || { fail_msg "$desc (z-fasta err)" "z-fasta BED --strand-aware --rc failed via stdin"; return; }
    else
        zf_args+=(--bed "$bed")
        "$ZFASTA" "${zf_args[@]}" > "$TMPDIR_LOCAL/zf.raw" 2>/dev/null \
            || { fail_msg "$desc (z-fasta err)" "z-fasta BED --strand-aware --rc failed"; return; }
    fi

    normalize_fasta_against_bed "$TMPDIR_LOCAL/bedtools.transformed" "$bed" "$TMPDIR_LOCAL/bedtools.norm" \
        || { fail_msg "$desc (bedtools normalize err)" "failed to normalize bedtools+seqtk output"; return; }
    normalize_fasta_against_bed "$TMPDIR_LOCAL/zf.raw" "$bed" "$TMPDIR_LOCAL/zf.norm" \
        || { fail_msg "$desc (z-fasta normalize err)" "failed to normalize z-fasta BED RC output"; return; }

    if diff -q "$TMPDIR_LOCAL/bedtools.norm" "$TMPDIR_LOCAL/zf.norm" >/dev/null 2>&1; then
        pass "$desc"
    else
        fail_diff "$desc" "bedtools+seqtk normalized" "$TMPDIR_LOCAL/bedtools.norm" "z-fasta normalized" "$TMPDIR_LOCAL/zf.norm"
    fi
}

verify_protein_error() {
    local desc="$1"

    if "$ZFASTA" get tests/data/proteome.fasta 'sp|P12345|PROT_HUMAN:1-10' --rc > "$TMPDIR_LOCAL/zf.tmp" 2> "$TMPDIR_LOCAL/zf.err"; then
        fail_msg "$desc" "expected z-fasta get --rc on protein to fail"
        return
    fi

    if grep -q 'reverse complement is not defined for protein sequences' "$TMPDIR_LOCAL/zf.err"; then
        pass "$desc"
    else
        fail_diff "$desc" "expected stderr" <(printf '%s\n' 'reverse complement is not defined for protein sequences') "z-fasta stderr" "$TMPDIR_LOCAL/zf.err"
    fi
}

echo "z-fasta reverse/complement verification"
echo "====================================="

cd "$PROJECT_DIR"
require_inputs

./zig build -Doptimize=ReleaseFast >/dev/null

verify_rc_region tests/data/simple.fasta seq1:1-5 "simple seq1:1-5 --rc vs samtools"
verify_rc_region tests/data/simple.fasta seq1:10-15 "simple seq1:10-15 line-boundary --rc vs samtools"
verify_rc_region tests/data/simple.fasta seq2:1-12 "simple seq2 full --rc vs samtools"
verify_rc_region tests/data/single.fasta single_sequence:1-4 "single-sequence span --rc vs samtools"
verify_rc_region tests/data/edge_cases.fasta lowercase:1-12 "lowercase preservation --rc vs samtools"

write_iupac_fixture
verify_rc_region "$TMPDIR_LOCAL/iupac.fasta" iupac_all:1-33 "synthetic IUPAC full-span --rc vs samtools"

write_chrom_fixture
verify_rc_region "$TMPDIR_LOCAL/chrom_like.fasta" chrSynthetic "synthetic chromosome-like full-sequence --rc vs samtools"

cat > "$TMPDIR_LOCAL/expected_iupac_complement.fa" <<'EOF'
>iupac_all:1-33 (complement)
TGCAAYRSWMKVHDBNtgcaayrswmkvhdbna
EOF
verify_exact_output "synthetic IUPAC --complement-only exact output" "$TMPDIR_LOCAL/expected_iupac_complement.fa" \
    get "$TMPDIR_LOCAL/iupac.fasta" iupac_all:1-33 --complement-only --annotate-rc

verify_multi_rc_concat "multi-region --rc same seq concatenation" tests/data/simple.fasta \
    seq1:1-5 seq1:10-15 seq1:20-24
verify_multi_rc_concat "multi-region --rc mixed seq order concatenation" tests/data/simple.fasta \
    seq2:1-4 seq1:1-5 seq2:5-12 seq1:10-15

REGIONS_20_RC=()
for i in 10 9 8 7 6 5 4 3 2 1; do
    REGIONS_20_RC+=("seq2:${i}-${i}")
    REGIONS_20_RC+=("seq1:${i}-${i}")
done
verify_multi_rc_concat "multi-region --rc 20-region sort-path concatenation" tests/data/simple.fasta "${REGIONS_20_RC[@]}"

cat > "$TMPDIR_LOCAL/names.txt" <<'EOF'
seq2
seq1
EOF
verify_rc_names tests/data/simple.fasta "$TMPDIR_LOCAL/names.txt" "names-file --rc vs samtools"

cat > "$TMPDIR_LOCAL/expected_rc_annotated.fa" <<'EOF'
>seq1:1-5 (reverse complement)
TACGT
EOF
verify_exact_output "annotated --rc header" "$TMPDIR_LOCAL/expected_rc_annotated.fa" \
    get tests/data/simple.fasta seq1:1-5 --rc --annotate-rc

cat > "$TMPDIR_LOCAL/expected_reverse_only.fa" <<'EOF'
>seq1:1-5 (reverse)
ATGCA
EOF
verify_exact_output "annotated --reverse-only" "$TMPDIR_LOCAL/expected_reverse_only.fa" \
    get tests/data/simple.fasta seq1:1-5 --reverse-only --annotate-rc

cat > "$TMPDIR_LOCAL/expected_complement_only.fa" <<'EOF'
>seq1:1-5 (complement)
TGCAT
EOF
verify_exact_output "annotated --complement-only" "$TMPDIR_LOCAL/expected_complement_only.fa" \
    get tests/data/simple.fasta seq1:1-5 --complement-only --annotate-rc

cat > "$TMPDIR_LOCAL/bed.bed" <<'EOF'
seq1	0	5	plus	0	+
seq1	0	5	minus	0	-
seq1	2	8	overlap_minus	0	-
seq2	0	4	short_plus	0	+
EOF
cat > "$TMPDIR_LOCAL/expected_bed_rc.fa" <<'EOF'
>seq1:1-5
TACGT
>seq1:1-5
ACGTA
>seq1:3-8
GTACGT
>seq2:1-4
CCCC
EOF
verify_exact_output "BED --strand-aware --rc composition" "$TMPDIR_LOCAL/expected_bed_rc.fa" \
    get tests/data/simple.fasta --bed "$TMPDIR_LOCAL/bed.bed" --strand-aware --rc

verify_protein_error "protein rc rejection"
verify_bed_rc_composed "BED --strand-aware --rc vs bedtools+seqtk" tests/data/simple.fasta "$TMPDIR_LOCAL/bed.bed" 2 0
verify_bed_rc_composed "BED --strand-aware --rc via stdin vs bedtools+seqtk" tests/data/simple.fasta "$TMPDIR_LOCAL/bed.bed" 2 1

echo ""
echo "====================================="
echo "Results: $PASS passed, $FAIL failed"
echo "====================================="
[ "$FAIL" -gt 0 ] && { echo "VERIFICATION FAILED"; exit 1; }
echo "ALL PASSED"
