#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
RESULTS_DIR="$SCRIPT_DIR/results"
ZFASTA="${ZFASTA:-$PROJECT_ROOT/zig-out/bin/z-fasta}"
RUNS=3
WARMUP=1
PROFILE="quick"
TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
OUT_DIR="$RESULTS_DIR/rc_review_$TIMESTAMP"
TMPDIR_LOCAL="$(mktemp -d "$SCRIPT_DIR/.rc_review_tmp.XXXXXX")"
trap 'rm -rf "$TMPDIR_LOCAL"' EXIT

while [[ $# -gt 0 ]]; do
    case "$1" in
        --runs)
            RUNS="$2"
            shift 2
            ;;
        --warmup)
            WARMUP="$2"
            shift 2
            ;;
        --quick)
            PROFILE="quick"
            shift
            ;;
        --full)
            PROFILE="full"
            shift
            ;;
        *)
            echo "Unknown option: $1" >&2
            exit 1
            ;;
    esac
done

command -v hyperfine >/dev/null 2>&1 || {
    echo "Error: hyperfine not found" >&2
    exit 1
}

command -v /usr/bin/time >/dev/null 2>&1 || {
    echo "Error: /usr/bin/time not found" >&2
    exit 1
}

[[ -x "$ZFASTA" ]] || {
    echo "Error: z-fasta not found at $ZFASTA. Run: ./zig build -Doptimize=ReleaseFast" >&2
    exit 1
}

mkdir -p "$OUT_DIR"
mkdir -p "$OUT_DIR/json"

FASTA="$TMPDIR_LOCAL/rc_review.fa"
BED="$TMPDIR_LOCAL/rc_review.bed"
BED_STRANDED="$TMPDIR_LOCAL/rc_review_stranded.bed"
SEQ_NAME="chrSynthetic"

if [[ "$PROFILE" == "full" ]]; then
    SEQ_LEN=8000000
    BED_COUNT=10000
    LARGE_REGION_SIZE=1000000
    MULTI_COUNT_A=10
    MULTI_COUNT_B=100
else
    SEQ_LEN=2000000
    BED_COUNT=2000
    LARGE_REGION_SIZE=250000
    MULTI_COUNT_A=10
    MULTI_COUNT_B=50
fi

write_fixture() {
    {
        echo ">$SEQ_NAME reverse-complement review fixture"
        awk -v total="$SEQ_LEN" -v width=60 '
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
    } > "$FASTA"

    "$ZFASTA" index "$FASTA" >/dev/null
}

write_bed() {
    : > "$BED"
    : > "$BED_STRANDED"

    local count="$BED_COUNT"
    local width=120
    local step=400
    local start=1000
    local strand

    for ((i = 0; i < count; i += 1)); do
        local s=$(( start + (i * step) ))
        local e=$(( s + width ))
        printf '%s\t%d\t%d\n' "$SEQ_NAME" "$s" "$e" >> "$BED"
        if ((( i % 2 ) == 0 )); then
            strand='+'
        else
            strand='-'
        fi
        printf '%s\t%d\t%d\tregion_%05d\t0\t%s\n' "$SEQ_NAME" "$s" "$e" "$i" "$strand" >> "$BED_STRANDED"
    done
}

make_region() {
    local start="$1"
    local size="$2"
    local end=$(( start + size - 1 ))
    printf '%s:%d-%d' "$SEQ_NAME" "$start" "$end"
}

build_multi_command() {
    local count="$1"
    local with_rc="$2"
    local args=()
    local size=1000
    local start=1000
    local max_start=$(( SEQ_LEN - size - 1 ))
    local usable_span=$(( max_start - start ))
    local step

    if (( usable_span <= 0 )); then
        echo "Synthetic RC review fixture is too short for multi-region benchmarking" >&2
        exit 1
    fi

    step=$(( usable_span / (count + 1) ))
    if (( step < size )); then
        step=$size
    fi

    for ((i = 0; i < count; i += 1)); do
        args+=("$(make_region $(( start + (i * step) )) "$size")")
    done

    printf '%q get %q' "$ZFASTA" "$FASTA"
    for region in "${args[@]}"; do
        printf ' %q' "$region"
    done
    if [[ "$with_rc" == "true" ]]; then
        printf ' --rc'
    fi
    printf ' > /dev/null'
}

run_hyperfine_markdown() {
    local markdown_out="$1"
    local json_out="$2"
    shift 2
    hyperfine --warmup "$WARMUP" --runs "$RUNS" --export-markdown "$markdown_out" --export-json "$json_out" "$@"
}

measure_rss() {
    local label="$1"
    local command="$2"
    local tmpfile="$TMPDIR_LOCAL/${label}.time"

    bash -lc "/usr/bin/time -v $command" >/dev/null 2> "$tmpfile"
    awk -v label="$label" '
        /Elapsed \(wall clock\) time/ { elapsed = $NF }
        /Maximum resident set size/ { rss = $NF }
        END { printf "%s\t%s\t%s\n", label, elapsed, rss }
    ' "$tmpfile"
}

write_fixture
write_bed

REGION_SMALL="$(make_region 1000 100)"
REGION_MEDIUM="$(make_region 200000 10000)"
REGION_LARGE="$(make_region 200000 "$LARGE_REGION_SIZE")"
REGION_FULL="$SEQ_NAME"
MULTI_A_FORWARD="$(build_multi_command "$MULTI_COUNT_A" false)"
MULTI_A_RC="$(build_multi_command "$MULTI_COUNT_A" true)"
MULTI_B_FORWARD="$(build_multi_command "$MULTI_COUNT_B" false)"
MULTI_B_RC="$(build_multi_command "$MULTI_COUNT_B" true)"

echo "Running focused RC review benchmarks"
echo "  output: $OUT_DIR"
echo "  profile: $PROFILE (runs=$RUNS, warmup=$WARMUP, seq_len=$SEQ_LEN, bed_count=$BED_COUNT)"

echo "  [1/6] orientation small"
run_hyperfine_markdown \
    "$OUT_DIR/orientation_small.md" \
    "$OUT_DIR/json/orientation_small.json" \
    -n forward "$ZFASTA get '$FASTA' '$REGION_SMALL' > /dev/null" \
    -n rc "$ZFASTA get '$FASTA' '$REGION_SMALL' --rc > /dev/null" \
    -n reverse-only "$ZFASTA get '$FASTA' '$REGION_SMALL' --reverse-only > /dev/null" \
    -n complement-only "$ZFASTA get '$FASTA' '$REGION_SMALL' --complement-only > /dev/null" \
    -n rc-annotate "$ZFASTA get '$FASTA' '$REGION_SMALL' --rc --annotate-rc > /dev/null"

echo "  [2/6] orientation medium"
run_hyperfine_markdown \
    "$OUT_DIR/orientation_medium.md" \
    "$OUT_DIR/json/orientation_medium.json" \
    -n forward "$ZFASTA get '$FASTA' '$REGION_MEDIUM' > /dev/null" \
    -n rc "$ZFASTA get '$FASTA' '$REGION_MEDIUM' --rc > /dev/null" \
    -n reverse-only "$ZFASTA get '$FASTA' '$REGION_MEDIUM' --reverse-only > /dev/null" \
    -n complement-only "$ZFASTA get '$FASTA' '$REGION_MEDIUM' --complement-only > /dev/null" \
    -n rc-annotate "$ZFASTA get '$FASTA' '$REGION_MEDIUM' --rc --annotate-rc > /dev/null"

echo "  [3/6] orientation large"
run_hyperfine_markdown \
    "$OUT_DIR/orientation_large.md" \
    "$OUT_DIR/json/orientation_large.json" \
    -n forward "$ZFASTA get '$FASTA' '$REGION_LARGE' > /dev/null" \
    -n rc "$ZFASTA get '$FASTA' '$REGION_LARGE' --rc > /dev/null" \
    -n reverse-only "$ZFASTA get '$FASTA' '$REGION_LARGE' --reverse-only > /dev/null" \
    -n complement-only "$ZFASTA get '$FASTA' '$REGION_LARGE' --complement-only > /dev/null" \
    -n rc-annotate "$ZFASTA get '$FASTA' '$REGION_LARGE' --rc --annotate-rc > /dev/null"

echo "  [4/6] full-sequence forward vs rc"
run_hyperfine_markdown \
    "$OUT_DIR/full_sequence.md" \
    "$OUT_DIR/json/full_sequence.json" \
    -n forward "$ZFASTA get '$FASTA' '$REGION_FULL' > /dev/null" \
    -n rc "$ZFASTA get '$FASTA' '$REGION_FULL' --rc > /dev/null" \
    -n rc-annotate "$ZFASTA get '$FASTA' '$REGION_FULL' --rc --annotate-rc > /dev/null"

echo "  [5/6] multi-region no-flag regression vs rc"
run_hyperfine_markdown \
    "$OUT_DIR/multi_region.md" \
    "$OUT_DIR/json/multi_region.json" \
    -n multi${MULTI_COUNT_A}-forward "$MULTI_A_FORWARD" \
    -n multi${MULTI_COUNT_A}-rc "$MULTI_A_RC" \
    -n multi${MULTI_COUNT_B}-forward "$MULTI_B_FORWARD" \
    -n multi${MULTI_COUNT_B}-rc "$MULTI_B_RC"

echo "  [6/6] bed batch orientation overhead"
run_hyperfine_markdown \
    "$OUT_DIR/bed_batch.md" \
    "$OUT_DIR/json/bed_batch.json" \
    -n bed-forward "$ZFASTA get '$FASTA' --bed '$BED' > /dev/null" \
    -n bed-rc "$ZFASTA get '$FASTA' --bed '$BED' --rc > /dev/null" \
    -n bed-reverse-only "$ZFASTA get '$FASTA' --bed '$BED' --reverse-only > /dev/null" \
    -n bed-complement-only "$ZFASTA get '$FASTA' --bed '$BED' --complement-only > /dev/null" \
    -n bed-honor-strand-rc "$ZFASTA get '$FASTA' --bed '$BED_STRANDED' --honor-strand --rc > /dev/null"

{
    printf 'label\telapsed\tmaxrss_kb\n'
    measure_rss forward "$ZFASTA get '$FASTA' '$REGION_LARGE' > /dev/null"
    measure_rss rc "$ZFASTA get '$FASTA' '$REGION_LARGE' --rc > /dev/null"
    measure_rss rc_annotate "$ZFASTA get '$FASTA' '$REGION_LARGE' --rc --annotate-rc > /dev/null"
    measure_rss multi_b_forward "$MULTI_B_FORWARD"
    measure_rss multi_b_rc "$MULTI_B_RC"
    measure_rss bed_forward "$ZFASTA get '$FASTA' --bed '$BED' > /dev/null"
    measure_rss bed_honor_strand_rc "$ZFASTA get '$FASTA' --bed '$BED_STRANDED' --honor-strand --rc > /dev/null"
} > "$OUT_DIR/rss.tsv"

cat <<EOF

Focused RC review complete.

Outputs:
  $OUT_DIR/orientation_small.md
  $OUT_DIR/orientation_medium.md
  $OUT_DIR/orientation_large.md
  $OUT_DIR/full_sequence.md
  $OUT_DIR/multi_region.md
  $OUT_DIR/bed_batch.md
  $OUT_DIR/rss.tsv
EOF