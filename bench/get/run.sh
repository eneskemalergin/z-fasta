#!/usr/bin/env bash
# GET benchmark runner: correctness (run_tests), then zebrac perf and report.
#
# Usage:
#   bash bench/get/run.sh [options]
#
# Defaults: correctness first (418 checks), then perf (--runs 5, --warmup 5, --duration 5000).
# A full perf pass (all sections, default zebrac settings) typically takes several hours.
# Multi-region seqtk/fastahack reference loops are omitted by default (O(N) subprocesses;
# cost grows linearly and duplicates positional seqtk reference data). Set
# GET_MULTI_REFERENCE=1 to re-enable them for a one-off reference sweep.
# Use --skip-* for partial re-runs; resume with GET_RUN_TIMESTAMP=<ts> and GET_MULTI_COUNTS.
# Use --skip-tests (alias --skip-verify) to run perf/report only. Correctness failure aborts before perf or report.
# --skip-messy skips messy zebrac perf only (never skips messy cases inside run_tests).
# --focused runs the complete development benchmark without regenerating the checked report.
#
# Canonical flags (index-aligned): --skip-tests --skip-benchmarks --skip-messy --skip-report
# Deprecated aliases (one release cycle): --skip-verify --skip-perf
#
#   bash bench/get/run.sh
#   bash bench/get/run.sh --skip-tests --skip-report   # perf only
#   bash bench/get/run.sh --skip-verify --skip-perf    # same via aliases
#   bash bench/get/run.sh --skip-messy --skip-bed
#   bash bench/get/run.sh --focused --skip-tests --runs 5 --warmup 5
#   bash bench/get/run.sh --regenerate-fixtures
#   bash bench/get/run.sh --skip-benchmarks --skip-report   # correctness only
#
# ── Workload fixtures (generated; bench/get/data/ is gitignored) ───
#   regions/{Dataset}_{label}.txt        positional slices on REAL_* (incl. 1kbp_mid for RC)
#   regions/{Dataset}_regions_N{n}.txt   multi-region lists
#   bed/{Dataset}_{n}.bed                BED batch (+ .regions.txt for samtools)
#   Created on demand by generate_get_fixtures in this script.
#
#   -h|--help  print this header

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH_ROOT="$(dirname "$SCRIPT_DIR")"
PROJECT_ROOT="$(dirname "$BENCH_ROOT")"
RESULTS_DIR="$SCRIPT_DIR/results"
DATA_DIR="$BENCH_ROOT/shared/data"
MESSY_DIR="$BENCH_ROOT/shared/cache/messy_perf"
MESSY_PERF_JSON="$SCRIPT_DIR/messy_perf.json"
FIXTURE_DIR="$SCRIPT_DIR/data"
REGIONS_DIR="$FIXTURE_DIR/regions"
BED_DIR="$FIXTURE_DIR/bed"
FOCUSED_DIR="$FIXTURE_DIR/focused"

source "$BENCH_ROOT/shared/runner_common.sh"

RUNS=5
WARMUP=5
ZEBRAC_DURATION_MS="${ZEBRAC_DURATION_MS:-5000}"
# Messy GET uses sub-millisecond timings on small slices; large proteome fixtures use
# higher sampling. Override via GET_MESSY_* or --messy-runs etc.
MESSY_RUNS="${GET_MESSY_RUNS:-50}"
MESSY_WARMUP="${GET_MESSY_WARMUP:-10}"
MESSY_DURATION_MS="${GET_MESSY_DURATION_MS:-30000}"
DO_TESTS=true
DO_BENCHMARKS=true
DO_POS=true
DO_MULTI=true
DO_BED=true
DO_RC=true
DO_MESSY=true
DO_REPORT=true
REGENERATE_FIXTURES=false
FOCUSED=false

BENCH_VIEW_DIR=""
TEST_DATA_DIR=""
declare -A ZFI_VIEWS=()
declare -A FAI_VIEWS=()

# Multi-region perf sweep (log-spaced). 15/16 boundary lives in run_tests, not timed here.
MULTI_COUNTS=(1 10 100 1000)
BED_COUNTS=(10 100 1000 10000)
POS_LABELS=(100bp 1kbp_mid 10kbp full_seq)

# Mirror getter.parseRegion: only the rightmost colon can start START-END or START-.
region_span_bases() {
    local region="$1"
    local fasta="${2:-}"
    python3 - "$region" "$fasta" <<'PY'
import sys

region = sys.argv[1].strip()
fasta = sys.argv[2] or None


def fai_length(name: str) -> int | None:
    if not fasta:
        return None
    for line in open(fasta + ".fai", encoding="ascii", errors="replace"):
        parts = line.split("\t")
        if parts[0] == name:
            return int(parts[1])
    return None


def parse_suffix(suffix: str) -> tuple[int, int | None] | None:
    dash = suffix.find("-")
    if dash < 0:
        return None
    start_s, end_s = suffix[:dash], suffix[dash + 1 :]
    if not start_s.isdigit():
        return None
    start = int(start_s)
    if end_s == "":
        return start, None
    if not end_s.isdigit():
        return None
    return start, int(end_s)


colon_pos = region.rfind(":")
if colon_pos >= 0:
    parsed = parse_suffix(region[colon_pos + 1 :])
    if parsed is not None:
        start, end = parsed
        if end is None:
            length = fai_length(region[:colon_pos])
            print(length - start + 1 if length is not None else 1)
        else:
            print(end - start + 1)
        raise SystemExit(0)

length = fai_length(region)
if length is not None:
    print(length)
    raise SystemExit(0)
raise SystemExit(f"region_span_bases: {region!r} not found in index")
PY
}

region_bases() {
    region_span_bases "$1" "$2"
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --runs) RUNS="$2"; shift 2 ;;
        --warmup) WARMUP="$2"; shift 2 ;;
        --duration) ZEBRAC_DURATION_MS="$2"; shift 2 ;;
        --messy-runs) MESSY_RUNS="$2"; shift 2 ;;
        --messy-warmup) MESSY_WARMUP="$2"; shift 2 ;;
        --messy-duration) MESSY_DURATION_MS="$2"; shift 2 ;;
        --skip-tests|--skip-verify) DO_TESTS=false; shift ;;
        --skip-benchmarks|--skip-perf) DO_BENCHMARKS=false; shift ;;
        --skip-pos) DO_POS=false; shift ;;
        --skip-multi) DO_MULTI=false; shift ;;
        --skip-bed) DO_BED=false; shift ;;
        --skip-rc) DO_RC=false; shift ;;
        --skip-messy) DO_MESSY=false; shift ;;
        --skip-report) DO_REPORT=false; shift ;;
        --regenerate-fixtures) REGENERATE_FIXTURES=true; shift ;;
        --focused) FOCUSED=true; DO_MESSY=false; DO_REPORT=false; shift ;;
        -h|--help)
            sed -n '2,/^set -euo pipefail$/p' "$0" | head -n -1
            exit 0
            ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

if $FOCUSED && { ! $DO_POS || ! $DO_MULTI || ! $DO_BED || ! $DO_RC; }; then
    echo "error: --focused cannot be combined with focused section skip flags" >&2
    exit 1
fi

get_add_command() {
    local section="$1" workload="$2" tool="$3" family="$4"
    local json_out="$5" script="$6" input_bytes="$7" output_bases="${8:-}"
    zebrac_add_command "get" "$section" "$workload" "$tool" "$family" \
        "$input_bytes" "$output_bases" "$json_out" "$(shell_command "$script")"
}

cleanup_benchmark_views() {
    if [[ -n "$BENCH_VIEW_DIR" && -d "$BENCH_VIEW_DIR" ]]; then
        rm -rf -- "$BENCH_VIEW_DIR"
    fi
}

prepare_benchmark_views() {
    bench_require_tool z-fasta
    bench_require_tool samtools

    BENCH_VIEW_DIR="$SCRIPT_DIR/.bench_views_${TIMESTAMP}_$$"
    mkdir -p "$BENCH_VIEW_DIR/zfi" "$BENCH_VIEW_DIR/fai"

    local ds source name zfi_view fai_view
    for ds in Genome Transcriptome Proteome; do
        source="${REAL_DATASETS[$ds]:-}"
        [[ -n "$source" ]] || continue
        name="${source##*/}"
        zfi_view="$BENCH_VIEW_DIR/zfi/$name"
        fai_view="$BENCH_VIEW_DIR/fai/$name"

        ln -s "$source" "$zfi_view"
        ln -s "$source" "$fai_view"

        "$ZFASTA" index "$zfi_view" > /dev/null
        "$SAMTOOLS" faidx "$fai_view" > /dev/null 2>&1

        [[ -f "${zfi_view}.zfi" && ! -e "${zfi_view}.fai" ]] || {
            echo "error: invalid private .zfi benchmark view for $ds" >&2
            return 1
        }
        [[ -f "${fai_view}.fai" && ! -e "${fai_view}.zfi" ]] || {
            echo "error: invalid private .fai benchmark view for $ds" >&2
            return 1
        }

        ZFI_VIEWS[$ds]="$zfi_view"
        FAI_VIEWS[$ds]="$fai_view"
    done
}

bed_to_regions() {
    awk -F'\t' '!/^[[:space:]]*$|^#|^track|^browser/{sub(/\r$/,""); printf "%s:%d-%d\n",$1,$2+1,$3}' "$1"
}

multi_positional_get_core() {
    local tool_q="$1"
    local fa_q="$2"
    local regions_path="$3"
    local extra="${4:-}"
    local qr
    qr="$(quote_arg "$regions_path")"
    echo "mapfile -t _regs < ${qr}; ${tool_q} get ${fa_q} \"\${_regs[@]}\"${extra}"
}

# Keep large request sets out of zebrac's argv. Bash reads them at sample time.
multi_positional_get_cmd() {
    echo "$(multi_positional_get_core "$1" "$2" "$3") > /dev/null"
}


export_manifest_tool_versions() {
    export_manifest_core_versions
    export_manifest_required_tool_versions bedtools seqtk
    export_manifest_optional_tool_versions fastahack noodles rustbio
}

write_run_manifest() {
    local manifest="$1" timestamp="$2" metadata="$3"
    python3 - "$manifest" "$timestamp" "$metadata" "$RUNS" "$WARMUP" "$ZEBRAC_DURATION_MS" \
        "$DO_TESTS" "$DO_POS" "$DO_MULTI" "$DO_BED" "$DO_RC" "$DO_MESSY" \
        "${SECTION_POS:-}" "${SECTION_MULTI:-}" "${SECTION_BED:-}" "${SECTION_RC:-}" "${SECTION_MESSY:-}" \
        "$FOCUSED" "${SECTION_FOCUSED:-}" <<'PY'
import json, os, sys
from pathlib import Path

manifest, ts, metadata, runs, warmup, duration = sys.argv[1:7]
do_verify = sys.argv[7] == "true"
section_vals = sys.argv[13:18]
focused = sys.argv[18] == "true"
focused_section = sys.argv[19]

sections = {}
keys = ("perf_pos", "perf_multi", "perf_bed", "perf_rc", "messy")
for key, val in zip(keys, section_vals):
    if val:
        sections[key] = val
if focused and focused_section:
    sections = {"focused": focused_section}

skip_pos = "perf_pos" not in sections
skip_multi = "perf_multi" not in sections
skip_bed = "perf_bed" not in sections
skip_rc = "perf_rc" not in sections
skip_messy = "messy" not in sections

tools = {}
for name in ("samtools", "bedtools", "seqtk", "fastahack", "noodles", "rustbio"):
    val = os.environ.get(f"BENCH_VER_{name.upper()}")
    if val:
        tools[name] = val

out = {
    "schema_version": "get-run.v1",
    "timestamp": ts,
    "runner": "zebrac",
    "mode": "warm",
    "profile": "focused" if focused else "publication",
    "zebrac": os.environ.get("BENCH_VER_ZEBRAC", ""),
    "z_fasta": os.environ.get("BENCH_VER_ZFASTA", ""),
    "runs": int(runs),
    "warmup": int(warmup),
    "duration_ms": int(duration),
    "verify_skipped": not do_verify,
    "verify_pass": os.environ.get("BENCH_VERIFY_PASS"),
    "index_preload": True,
    "sidecars": "private single-format views",
    "skip_pos": skip_pos,
    "skip_multi": skip_multi,
    "skip_bed": skip_bed,
    "skip_rc": skip_rc,
    "skip_messy": skip_messy,
    "skip_names": not focused,
    "metadata": metadata,
    "tools": tools,
    "sections": sections,
}
for key, env_key in (
    ("messy_runs", "GET_MESSY_RUNS"),
    ("messy_warmup", "GET_MESSY_WARMUP"),
    ("messy_duration_ms", "GET_MESSY_DURATION_MS"),
):
    val = os.environ.get(env_key)
    if val:
        out[key] = int(val)
manifest_path = Path(manifest)
if manifest_path.is_file():
    try:
        prev = json.loads(manifest_path.read_text())
    except (json.JSONDecodeError, OSError):
        prev = {}
    for key in ("messy_runs", "messy_warmup", "messy_duration_ms"):
        if key not in out and key in prev:
            out[key] = prev[key]
manifest_path.write_text(json.dumps(out, indent=2) + "\n")
PY
}

generate_get_fixtures() {
    local mode="${1:-all}"
    mkdir -p "$REGIONS_DIR" "$BED_DIR" "$FOCUSED_DIR"
    python3 - "$REGIONS_DIR" "$BED_DIR" "$mode" "$REGENERATE_FIXTURES" \
        "${FAI_VIEWS[Genome]:-${REAL_DATASETS[Genome]:-}}" \
        "${FAI_VIEWS[Transcriptome]:-${REAL_DATASETS[Transcriptome]:-}}" \
        "${FAI_VIEWS[Proteome]:-${REAL_DATASETS[Proteome]:-}}" \
        "${BED_COUNTS[*]}" \
        "${MULTI_COUNTS[*]}" \
        "$FOCUSED_DIR" <<'PY'
import random
import sys
from pathlib import Path

regions_dir = Path(sys.argv[1])
bed_dir = Path(sys.argv[2])
mode = sys.argv[3]
regenerate = sys.argv[4] == "true"
genome_fa = sys.argv[5] or None
transcriptome_fa = sys.argv[6] or None
proteome_fa = sys.argv[7] or None
bed_counts = [int(x) for x in sys.argv[8].split()] if sys.argv[8] else []
multi_counts = [int(x) for x in sys.argv[9].split()] if len(sys.argv) > 9 and sys.argv[9] else []
focused_dir = Path(sys.argv[10])

MULTI_COUNTS = tuple(multi_counts) if multi_counts else (1, 10, 100, 1000)
POS_LABELS = (("100bp", 100), ("1kbp_mid", 1000), ("10kbp", 10000))
SPAN = 1000
MULTI_SEED = 42
FULL_SEQ_GENOME_MAX = 1000
FULL_SEQ_PROTEOME_MAX = 5000
FULL_SEQ_TRANSCRIPTOME_MAX = 20000
DATASET_SEEDS = {
    "Genome": 1665052211,
    "Transcriptome": 3571123842,
    "Proteome": 3098860823,
}


def full_seq_region(dataset: str, rows: list[tuple[str, int]]) -> str | None:
    """Name-only or capped whole-entry fetch; avoids pathological entries (e.g. titin)."""
    if not rows:
        return None
    if dataset == "Genome":
        short = [r for r in rows if 100 <= r[1] <= FULL_SEQ_GENOME_MAX]
        if not short:
            short = [r for r in rows if r[1] >= 100]
        name, length = min(short, key=lambda r: r[1])
        end = min(length, FULL_SEQ_GENOME_MAX)
        if end == length:
            return name
        return f"{name}:1-{end}"
    if dataset == "Transcriptome":
        eligible = [r for r in rows if 1000 <= r[1] <= FULL_SEQ_TRANSCRIPTOME_MAX]
        if not eligible:
            eligible = [r for r in rows if r[1] >= 1000]
        if not eligible:
            return None
        name, _length = min(eligible, key=lambda r: r[1])
        return name
    if dataset == "Proteome":
        eligible = [r for r in rows if 500 <= r[1] <= FULL_SEQ_PROTEOME_MAX]
        if not eligible:
            eligible = [r for r in rows if r[1] >= 500]
        if not eligible:
            return None
        name, _length = min(eligible, key=lambda r: r[1])
        return name
    return None


def read_fai(path: Path) -> list[tuple[str, int]]:
    rows: list[tuple[str, int]] = []
    for line in path.read_text().splitlines():
        parts = line.split("\t")
        if len(parts) >= 2 and parts[1].isdigit():
            rows.append((parts[0], int(parts[1])))
    if not rows:
        raise SystemExit(f"empty or missing fai: {path}")
    return rows


def pick_row(rows: list[tuple[str, int]], min_len: int) -> tuple[str, int]:
    candidates = [r for r in rows if r[1] >= min_len]
    if not candidates:
        raise SystemExit(f"no sequence with length >= {min_len}")
    return max(candidates, key=lambda r: r[1])


def region_span(name: str, length: int, span: int) -> str:
    if span >= length:
        return f"{name}:1-{length}"
    start = max(1, length // 2 - span // 2)
    end = start + span - 1
    if end > length:
        end = length
        start = max(1, end - span + 1)
    return f"{name}:{start}-{end}"


def write_if_needed(path: Path, text: str) -> None:
    if regenerate or not path.is_file():
        path.write_text(text)


def pos_regions(dataset: str, fasta: str) -> None:
    fai = Path(fasta + ".fai")
    rows = read_fai(fai)
    name, length = pick_row(rows, 100)
    for label, span in POS_LABELS:
        if span > length and label != "100bp":
            continue
        region = region_span(name, length, span)
        write_if_needed(regions_dir / f"{dataset}_{label}.txt", region + "\n")
    fs = full_seq_region(dataset, rows)
    if fs:
        write_if_needed(regions_dir / f"{dataset}_full_seq.txt", fs + "\n")


def multi_regions(dataset: str, fasta: str) -> None:
    fai = Path(fasta + ".fai")
    rows = [r for r in read_fai(fai) if r[1] >= SPAN]
    if not rows:
        raise SystemExit(f"{dataset}: no sequences with length >= {SPAN}")
    eligible = len(rows)
    rng = random.Random(MULTI_SEED + DATASET_SEEDS[dataset] % 10000)
    for n in MULTI_COUNTS:
        if n > eligible:
            continue
        path = regions_dir / f"{dataset}_regions_N{n}.txt"
        if not regenerate and path.is_file():
            continue
        lines = [region_span(name, length, SPAN) for name, length in rng.sample(rows, n)]
        path.write_text("\n".join(lines) + "\n")


def bed_file(dataset: str, fasta: str, out_path: Path, rows: int, seed: int) -> None:
    fai = Path(fasta + ".fai")
    eligible = [r for r in read_fai(fai) if r[1] >= SPAN]
    if not eligible:
        raise SystemExit("no sequences long enough for BED spans")
    rng = random.Random(seed)
    if regenerate or not out_path.is_file():
        with out_path.open("w", encoding="ascii") as handle:
            handle.write(f"# synthetic BED {dataset} {rows} rows span={SPAN}\n")
            for i in range(rows):
                name, length = rng.choice(eligible)
                start = rng.randint(0, max(0, length - SPAN))
                end = start + SPAN
                strand = "+" if i % 2 == 0 else "-"
                handle.write(f"{name}\t{start}\t{end}\t{name}_{i}\t0\t{strand}\n")


def focused_files(dataset: str, fasta: str) -> None:
    rows = [(name, length) for name, length in read_fai(Path(fasta + ".fai")) if length > 0]
    eligible = [(name, length) for name, length in rows if length >= SPAN]
    if not eligible:
        raise SystemExit(f"{dataset}: no focused records with length >= {SPAN}")

    rng = random.Random(DATASET_SEEDS[dataset])
    if dataset != "Genome":
        name_pool = [row for row in rows if row[1] <= (20000 if dataset == "Transcriptome" else 5000)]
        count = min(100, len(name_pool))
        names = rng.sample(name_pool, count)
        write_if_needed(
            focused_dir / f"{dataset}_names_{count}.txt",
            "".join(f"{name}\n" for name, _length in names),
        )

    longest_name, longest_len = max(eligible, key=lambda row: row[1])
    clustered = []
    start0 = max(0, longest_len // 2 - 2500)
    for i in range(100):
        start = min(start0 + i * 32, longest_len - SPAN)
        clustered.append(f"{longest_name}:{start + 1}-{start + SPAN}")
    write_if_needed(
        focused_dir / f"{dataset}_locality_clustered.txt",
        "\n".join(clustered) + "\n",
    )

    alternating_rows = sorted(eligible, key=lambda row: row[1], reverse=True)[:2]
    if len(alternating_rows) == 2:
        alternating = [region_span(*alternating_rows[i % 2], SPAN) for i in range(100)]
        write_if_needed(
            focused_dir / f"{dataset}_locality_alternating.txt",
            "\n".join(alternating) + "\n",
        )

    bed_path = focused_dir / f"{dataset}_strand_100.bed"
    if regenerate or not bed_path.is_file():
        with bed_path.open("w", encoding="ascii") as handle:
            handle.write(f"# focused BED {dataset}\n")
            for i in range(100):
                name, length = rng.choice(eligible)
                start = rng.randint(0, length - SPAN)
                strand = "+" if i % 2 == 0 else "-"
                handle.write(f"{name}\t{start}\t{start + SPAN}\t.\t0\t{strand}\n")


if mode in ("all", "pos"):
    for dataset, fasta in (
        ("Genome", genome_fa),
        ("Transcriptome", transcriptome_fa),
        ("Proteome", proteome_fa),
    ):
        if fasta:
            pos_regions(dataset, fasta)
dataset_fasta_pairs = (
    ("Genome", genome_fa),
    ("Transcriptome", transcriptome_fa),
    ("Proteome", proteome_fa),
)
if mode in ("all", "multi"):
    for dataset, fasta in dataset_fasta_pairs:
        if fasta:
            multi_regions(dataset, fasta)
if mode in ("all", "bed"):
    for dataset, fasta in dataset_fasta_pairs:
        if not fasta:
            continue
        for n in bed_counts:
            bed_file(dataset, fasta, bed_dir / f"{dataset}_{n}.bed", n, n + DATASET_SEEDS[dataset] % 1000)
if mode == "focused":
    focused_dir.mkdir(parents=True, exist_ok=True)
    for dataset, fasta in dataset_fasta_pairs:
        if fasta:
            pos_regions(dataset, fasta)
            multi_regions(dataset, fasta)
            focused_files(dataset, fasta)
PY
    if [[ "$mode" == "all" || "$mode" == "bed" ]]; then
        local ds n bed regions
        for ds in Genome Transcriptome Proteome; do
            for n in "${BED_COUNTS[@]}"; do
                bed="$BED_DIR/${ds}_${n}.bed"
                regions="$BED_DIR/${ds}_${n}.regions.txt"
                [[ -f "$bed" ]] || continue
                if $REGENERATE_FIXTURES || [[ ! -f "$regions" ]]; then
                    bed_to_regions "$bed" > "$regions"
                fi
            done
        done
    fi
}



# ══════════════════════════════════════════════════════════════════════
#  Correctness: z-fasta get vs oracles
# ══════════════════════════════════════════════════════════════════════

SKIP_INDEX=false SKIP_GET=false SKIP_MULTI=false SKIP_BED=false SKIP_RC=false SKIP_EDGE=false

# Correctness fixture paths (distinct from perf MESSY_DIR=shared/cache/messy_perf)
MESSY_TEST_DIR="$BENCH_ROOT/shared/cache/messy_fixtures"

MESSY_POS_CASES=(
    "mixed_widths:1-1,1-32,3-24,10-20"
    "trailing_whitespace:1-8,7-20,24-24"
    "blank_lines:1-4,9-12,17-20"
    "mixed_crlf:1-8,3-18,7-14"
)

PASS=0 FAIL=0

# --- primitives ---
pass() { PASS=$((PASS + 1)); printf "  PASS: %s\n" "$1"; }

fail() {
    FAIL=$((FAIL + 1)); printf "  FAIL: %s\n" "$1"
    local ef="${2:-$TMPDIR/expected.tmp}" gf="${3:-$TMPDIR/got.tmp}"
    echo "    expected:"; head -3 "$ef" 2>/dev/null | sed 's/^/      /'
    echo "    got:";      head -3 "$gf" 2>/dev/null | sed 's/^/      /'
}

diff_oracle() {
    diff -q "$1" "$2" >/dev/null 2>&1 && pass "$3" || fail "$3" "${4:-$1}" "$2"
}

# samtools 1.24 formats extracted FASTA at a fixed 60-column width. Older
# samtools releases followed the source record's wrapping, which is also what
# z-fasta and the wrappers emit. Compare headers and bases here while leaving
# the repository's own byte-level index/output checks unchanged.
diff_fasta_oracle() {
    local expected="$1" got="$2" label="$3"
    if python3 - "$expected" "$got" <<'PY'
import sys
from pathlib import Path


def canonical(path: str) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    sequence: list[str] = []
    for line in Path(path).read_text(encoding="ascii", errors="replace").splitlines():
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(sequence)))
            header = line
            sequence = []
        elif header is not None:
            sequence.append("".join(line.split()))
    if header is not None:
        records.append((header, "".join(sequence)))
    return records


raise SystemExit(0 if canonical(sys.argv[1]) == canonical(sys.argv[2]) else 1)
PY
    then
        pass "$label"
    else
        fail "$label" "$expected" "$got"
    fi
}

section_hdr() { echo ""; echo "--- [$1] $2 ---"; }

ensure_index() {
    local f="$1"
    [[ -f "${f}.fai" ]] || "$SAMTOOLS" faidx "$f" 2>/dev/null
    [[ -f "${f}.zfi" ]] || "$ZFASTA" index "$f" 2>/dev/null
}

# Messy-FASTA oracles (region/rc/rev/comp/bed). argv: MODE FASTA ... OUT
oracle() {
    python3 - "$@" <<'PY'
import sys
from pathlib import Path

RC = str.maketrans(
    "ACGTURYSWKMBDHVNacgturyswkmbdhvn",
    "TGCAAYRSWMKVHDBNtgcaayrswmkvhdbn",
)


def load_seqs(fasta: str) -> dict[str, str]:
    seqs: dict[str, list[str]] = {}
    cur: str | None = None
    for line in Path(fasta).read_text().splitlines():
        if line.startswith(">"):
            cur = line[1:].split()[0]
            seqs[cur] = []
        elif cur is not None:
            seqs[cur].append("".join(c for c in line if not c.isspace()))
    return {name: "".join(parts) for name, parts in seqs.items()}


def write_region(path: str, name: str, start: int, end: int, seq: str) -> None:
    Path(path).write_text(f">{name}:{start}-{end}\n{seq}\n")


def cmd_region(argv: list[str], transform) -> None:
    fasta, name, start, end, out = argv[2], argv[3], int(argv[4]), int(argv[5]), argv[6]
    seq = load_seqs(fasta)[name][start - 1 : end]
    if transform is not None:
        seq = transform(seq)
    write_region(out, name, start, end, seq)


def cmd_bed(argv: list[str], stranded_rc: bool = False) -> None:
    fasta, bed, out = argv[2], argv[3], argv[4]
    seqs = load_seqs(fasta)
    parts: list[str] = []
    for line in Path(bed).read_text().splitlines():
        if not line or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
            continue
        fields = line.split("\t")
        chrom, s0, e0 = fields[0], int(fields[1]), int(fields[2])
        seq = seqs[chrom][s0:e0]
        if stranded_rc and (fields[5] if len(fields) >= 6 else ".") != "-":
            seq = seq.translate(RC)[::-1]
        parts.append(f">{chrom}:{s0 + 1}-{e0}\n{seq}\n")
    Path(out).write_text("".join(parts))


COMMANDS = {
    "region": lambda a: cmd_region(a, None),
    "rc": lambda a: cmd_region(a, lambda s: s.translate(RC)[::-1]),
    "rev": lambda a: cmd_region(a, lambda s: s[::-1]),
    "comp": lambda a: cmd_region(a, lambda s: s.translate(RC)),
    "bed": lambda a: cmd_bed(a),
    "bed-stranded-rc": lambda a: cmd_bed(a, stranded_rc=True),
}

if len(sys.argv) < 2 or sys.argv[1] not in COMMANDS:
    print("usage: oracle <region|rc|rev|comp|bed|bed-stranded-rc> ...", file=sys.stderr)
    sys.exit(2)
COMMANDS[sys.argv[1]](sys.argv)
PY
}

_fixture_paths() {
    local prefix="$1" src="$2" desc="$3"
    local base="${src##*/}"; base="${base%.fasta}"
    local safe="${desc// /_}"; safe="${safe//\//_}"; safe="${safe//|/_}"
    FIXTURE_FASTA="$TMPDIR/${prefix}_${base}_${safe}.fasta"
    FIXTURE_STASH="$TMPDIR/${prefix}_${base}_${safe}.zfi"
}

_zfasta_get() {
    local fasta="$1" stdin_file="${2:-}"; shift 2
    if [[ -n "$stdin_file" ]]; then
        cat "$stdin_file" | "$ZFASTA" get "$fasta" "$@"
    else
        "$ZFASTA" get "$fasta" "$@"
    fi
}

names_to_regions() { awk '!/^[[:space:]]*$|^#/{sub(/\r$/,""); print}' "$1"; }

# Split name:start-end on the last ':' (handles sp|P12345|PROT_HUMAN:1-10).
# Sets PARSE_NAME, PARSE_START, PARSE_END; returns 0 when coords present.
parse_region_spec() {
    local region="$1"
    PARSE_NAME="" PARSE_START="" PARSE_END=""
    [[ "$region" == *:* ]] || { PARSE_NAME="$region"; return 1; }
    local coords="${region##*:}"
    PARSE_NAME="${region%:"$coords"}"
    [[ "$coords" == *-* ]] || return 1
    PARSE_START="${coords%%-*}"
    PARSE_END="${coords#*-}"
    [[ -n "$PARSE_START" && -n "$PARSE_END" ]] || return 1
    return 0
}

fasta_concat_sequence() {
    awk '/^>/{next} {gsub(/[[:space:]]/, ""); printf "%s", $0} END{print ""}' "$1" > "$2"
}

normalize_bed_output() {
    local raw_in="$1" bed_in="$2" honor_strand="$3" fasta_out="$4"
    awk -v bed_file="$bed_in" -v hs="$honor_strand" '
        function flush(    h, s2) {
            if (!idx) return
            if (!(idx in c)) { printf("record %d exceeds BED rows\n", idx) > "/dev/stderr"; exit 1 }
            h = c[idx] ":" s1[idx] "-" e[idx]
            if (hs == "1" && str[idx] == "-") h = h ":rc"
            print ">" h; s2 = seq
            while (length(s2) > 60) { print substr(s2, 1, 60); s2 = substr(s2, 61) }
            print s2; seq = ""
        }
        BEGIN {
            FS = "\t"
            while ((getline < bed_file) > 0) {
                sub(/\r$/, "", $0)
                if ($0 == "" || $0 ~ /^#/ || $0 ~ /^track/ || $0 ~ /^browser/) continue
                n++; c[n] = $1; s1[n] = $2 + 1; e[n] = $3; str[n] = NF >= 6 ? $6 : "."
            }
            close(bed_file); idx = 0; seq = ""
        }
        /^>/ { flush(); idx++; next }
        { gsub(/[[:space:]]/, "", $0); seq = seq $0 }
        END { flush(); if (idx != n) { printf("mismatch: fasta=%d bed=%d\n", idx, n) > "/dev/stderr"; exit 1 } }
    ' "$raw_in" > "$fasta_out"
}

expect_fail() {
    local label="$1"; shift
    "$@" >/dev/null 2>/dev/null && fail "$label (expected error)" || pass "$label"
}

verify_oracle_transform() {
    local mode="$1" label="$2" fasta="$3" name="$4" start="$5" end="$6"
    shift 6
    oracle "$mode" "$fasta" "$name" "$start" "$end" "$TMPDIR/expected.tmp"
    local zf_args=(get "$fasta" "${name}:${start}-${end}")
    zf_args+=("$@")
    "$ZFASTA" "${zf_args[@]}" > "$TMPDIR/got.tmp" 2>/dev/null \
        || { fail "$label (z-fasta err)"; return; }
    diff_oracle "$TMPDIR/expected.tmp" "$TMPDIR/got.tmp" "$label"
}

verify_exact_get() {
    local label="$1" expected="$2" fasta="$3"; shift 3
    "$ZFASTA" get "$fasta" "$@" > "$TMPDIR/got.tmp" 2>/dev/null \
        && diff_oracle "$expected" "$TMPDIR/got.tmp" "$label" \
        || fail "$label"
}

# --- parity: samtools ---
parity_samtools() {
    local label="$1" fasta="$2"; shift 2
    "$SAMTOOLS" faidx "$fasta" "$@" > "$TMPDIR/expected.tmp" 2>/dev/null \
        || { fail "$label (samtools err)"; return; }
    "$ZFASTA" get "$fasta" "$@" > "$TMPDIR/got.tmp" 2>/dev/null \
        || { fail "$label (z-fasta err)"; return; }
    diff_fasta_oracle "$TMPDIR/expected.tmp" "$TMPDIR/got.tmp" "$label"
}

parity_samtools_region() { parity_samtools "[parity:samtools] $3" "$1" "$2"; }

parity_samtools_regions() {
    local fasta="$1" desc="$2"; shift 2
    parity_samtools "[parity:samtools] $desc" "$fasta" "$@"
}

parity_samtools_rc() {
    local label="[parity:samtools] $3" fasta="$1" region="$2"
    "$SAMTOOLS" faidx -i --mark-strand no "$fasta" "$region" > "$TMPDIR/expected.tmp" 2>/dev/null \
        || { fail "$label (samtools err)"; return; }
    "$ZFASTA" get "$fasta" "$region" --rc > "$TMPDIR/got.tmp" 2>/dev/null \
        || { fail "$label (z-fasta err)"; return; }
    diff_fasta_oracle "$TMPDIR/expected.tmp" "$TMPDIR/got.tmp" "$label"
}

verify_open_ended_region() {
    local fasta="$1" name="$2" start="$3" len="$4" desc="$5"
    local region="${name}:${start}-"
    local seq_label="[parity:seq] $desc bases vs samtools"
    local hdr_label="[extended:header] $desc open suffix (samtools) vs clamped end (z-fasta)"

    "$SAMTOOLS" faidx "$fasta" "$region" > "$TMPDIR/expected.tmp" 2>/dev/null \
        || { fail "$seq_label (samtools err)"; return; }
    "$ZFASTA" get "$fasta" "$region" > "$TMPDIR/got.tmp" 2>/dev/null \
        || { fail "$seq_label (z-fasta err)"; return; }

    fasta_concat_sequence "$TMPDIR/expected.tmp" "$TMPDIR/exp.seq"
    fasta_concat_sequence "$TMPDIR/got.tmp" "$TMPDIR/got.seq"
    diff_oracle "$TMPDIR/exp.seq" "$TMPDIR/got.seq" "$seq_label"

    local st_hdr zf_hdr
    st_hdr=$(head -1 "$TMPDIR/expected.tmp")
    zf_hdr=$(head -1 "$TMPDIR/got.tmp")
    if [[ "$st_hdr" == ">${name}:${start}-" && "$zf_hdr" == ">${name}:${start}-${len}" ]]; then
        pass "$hdr_label"
    else
        fail "$hdr_label" "$TMPDIR/expected.tmp" "$TMPDIR/got.tmp"
    fi
}

# verify_index_cross: byte-identical get with `.zfi` present vs `.fai` selected.
# Optional third arg: existing file → read GET input from stdin (e.g. BED via pipe).
verify_index_cross() {
    local desc="$1" src="$2"
    shift 2
    local stdin_file=""
    [[ $# -gt 0 && -f "$1" ]] && { stdin_file="$1"; shift; }

    _fixture_paths xidx "$src" "$desc"
    local fasta="$FIXTURE_FASTA" stash="$FIXTURE_STASH"
    local tag="[index:cross] $desc"

    cp "$src" "$fasta"
    ensure_index "$fasta"
    [[ -f "${fasta}.zfi" && -f "${fasta}.fai" ]] || { fail "[index:cross] $desc sidecars missing"; return; }

    if ! _zfasta_get "$fasta" "$stdin_file" "$@" > "$TMPDIR/zfi.out" 2>/dev/null; then
        fail "[index:zfi] $desc get failed"
        return
    fi
    mv "${fasta}.zfi" "$stash"
    if ! _zfasta_get "$fasta" "$stdin_file" "$@" > "$TMPDIR/fai.out" 2>/dev/null; then
        mv "$stash" "${fasta}.zfi"
        fail "[index:fai] $desc get via .fai failed"
        return
    fi
    diff_oracle "$TMPDIR/zfi.out" "$TMPDIR/fai.out" "$tag"
    mv "$stash" "${fasta}.zfi"
}

verify_index() {
    local mode="$1" src="$2" target="$3" desc="$4"
    local base="${src##*/}"; base="${base%.fasta}"
    _fixture_paths "idx_${mode}" "$src" "$desc"
    local fasta="$FIXTURE_FASTA" stash="$FIXTURE_STASH"
    local tag="idx $mode $base $desc"
    local zf_get=(get "$fasta") st_get=(faidx)

    cp "$src" "$fasta"
    ensure_index "$fasta"
    [[ -f "${fasta}.zfi" && -f "${fasta}.fai" ]] || { fail "[index:zfi] $tag sidecars missing"; return; }

    if [[ "$mode" == bed ]]; then
        bed_to_regions "$target" > "$TMPDIR/idx_regions.txt"
        zf_get+=(--bed "$target")
        st_get+=(-r "$TMPDIR/idx_regions.txt" "$fasta")
    else
        zf_get+=("$target")
        st_get+=("$fasta" "$target")
    fi

    "$ZFASTA" "${zf_get[@]}" > "$TMPDIR/zfi.out" 2>/dev/null \
        || { fail "[index:zfi] $tag get failed"; return; }
    pass "[index:zfi] $tag get with .zfi present"

    "$SAMTOOLS" "${st_get[@]}" > "$TMPDIR/st.out" 2>/dev/null \
        || { fail "[parity:samtools] $tag (samtools err)"; return; }
    diff_fasta_oracle "$TMPDIR/st.out" "$TMPDIR/zfi.out" "[parity:samtools] $tag via .zfi"

    mv "${fasta}.zfi" "$stash"
    "$ZFASTA" "${zf_get[@]}" > "$TMPDIR/fai.out" 2>/dev/null || {
        mv "$stash" "${fasta}.zfi"; fail "[index:fai] $tag get via .fai failed"; return
    }
    pass "[index:fai] $tag get with .zfi absent"
    diff_fasta_oracle "$TMPDIR/st.out" "$TMPDIR/fai.out" "[parity:samtools] $tag via .fai"
    diff_oracle "$TMPDIR/zfi.out" "$TMPDIR/fai.out" "[index:cross] $tag .zfi == .fai"
    mv "$stash" "${fasta}.zfi"
}

verify_messy_zfi_required() {
    local name="$1" start="$2" end="$3"
    local fasta="$TMPDIR/messy_nozfi_${name}.fasta"
    local stash="$TMPDIR/messy_nozfi_${name}.zfi"
    [[ -f "$MESSY_TEST_DIR/${name}.fasta" ]] || { fail "[extended:messy] missing fixture $name"; return; }
    prepare_messy_fixture "$name" "$fasta" || { fail "[extended:messy] index failed $name"; return; }
    mv "${fasta}.zfi" "$stash"
    rm -f "${fasta}.fai"
    expect_fail "[extended:messy] $name:${start}-${end} without .zfi" \
        "$ZFASTA" get "$fasta" "${name}:${start}-${end}"
    mv "$stash" "${fasta}.zfi"
}

# --- extended: messy (Python oracle) ---
prepare_messy_fixture() {
    cp "$MESSY_TEST_DIR/${1}.fasta" "$2"
    "$ZFASTA" index "$2" >/dev/null 2>&1
}

verify_messy_region() {
    local fasta="$1" name="$2" start="$3" end="$4" desc="$5"
    local mode=region label="[extended:messy] $desc" extra=()
    [[ "${6:-}" == "--rc" ]] && { mode=rc; label="$label --rc"; extra=(--rc); }
    verify_oracle_transform "$mode" "$label" "$fasta" "$name" "$start" "$end" "${extra[@]}"
}

verify_expected_bed() {
    local label="[extended:messy] $1" fasta="$2" bed="$3" expected="$4"
    shift 4
    "$ZFASTA" get "$fasta" --bed "$bed" "$@" > "$TMPDIR/got.tmp" 2>/dev/null \
        || { fail "$label (z-fasta err)"; return; }
    diff_oracle "$expected" "$TMPDIR/got.tmp" "$label"
}

# --- parity: bedtools ---
parity_bedtools() {
    local label="[parity:bedtools] $1" bed="$2" honor_strand="$3" fasta="$4"
    local stdin="${5:-0}" strand_flag="${6:---strand-aware}"
    local bt=("$BEDTOOLS" getfasta -fi "$fasta" -bed "$bed")
    [[ "$honor_strand" == "1" ]] && bt+=(-s)
    local zf=(get "$fasta" --bed "$bed")
    [[ "$honor_strand" == "1" ]] && zf+=("$strand_flag")
    if [[ "$stdin" == "1" ]]; then
        zf[2]="--bed"
        zf[3]="-"
        "${bt[@]}" > "$TMPDIR/bt.raw" 2>/dev/null || { fail "$label (bedtools err)"; return; }
        cat "$bed" | "$ZFASTA" "${zf[@]}" > "$TMPDIR/zf.raw" 2>/dev/null || { fail "$label (z-fasta err)"; return; }
    else
        "${bt[@]}" > "$TMPDIR/bt.raw" 2>/dev/null || { fail "$label (bedtools err)"; return; }
        "$ZFASTA" "${zf[@]}" > "$TMPDIR/zf.raw" 2>/dev/null || { fail "$label (z-fasta err)"; return; }
    fi
    normalize_bed_output "$TMPDIR/bt.raw" "$bed" "$honor_strand" "$TMPDIR/bt.norm"
    normalize_bed_output "$TMPDIR/zf.raw" "$bed" "$honor_strand" "$TMPDIR/zf.norm"
    diff_oracle "$TMPDIR/bt.norm" "$TMPDIR/zf.norm" "$label"
}

parity_bed_samtools() {
    local label="[parity:samtools] $1" bed="$2" fasta="$3"
    local stdin="${4:-0}"
    local zf=(get "$fasta")
    bed_to_regions "$bed" > "$TMPDIR/st_regions.txt"
    if [[ "$stdin" == "1" ]]; then
        zf+=(--bed -)
        cat "$bed" | "$ZFASTA" "${zf[@]}" > "$TMPDIR/got.tmp" 2>/dev/null || { fail "$label (z-fasta err)"; return; }
    else
        zf+=(--bed "$bed")
        "$ZFASTA" "${zf[@]}" > "$TMPDIR/got.tmp" 2>/dev/null || { fail "$label (z-fasta err)"; return; }
    fi
    "$SAMTOOLS" faidx -r "$TMPDIR/st_regions.txt" "$fasta" > "$TMPDIR/expected.tmp" 2>/dev/null \
        || { fail "$label (samtools err)"; return; }
    diff_fasta_oracle "$TMPDIR/expected.tmp" "$TMPDIR/got.tmp" "$label"
}

parity_bed_rc_composed() {
    local label="[parity:bedtools] $1" fasta="$2" bed="$3"
    local stdin="${4:-0}"
    local zf=(get "$fasta" --strand-aware --rc)
    if ! "$BEDTOOLS" getfasta -fi "$fasta" -bed "$bed" -s 2>/dev/null \
        | "$SEQTK" seq -r > "$TMPDIR/bt.raw"; then
        fail "$label (bedtools/seqtk err)"
        return
    fi
    if [[ "$stdin" == "1" ]]; then
        cat "$bed" | "$ZFASTA" "${zf[@]}" --bed - > "$TMPDIR/zf.raw" 2>/dev/null || { fail "$label (z-fasta err)"; return; }
    else
        "$ZFASTA" "${zf[@]}" --bed "$bed" > "$TMPDIR/zf.raw" 2>/dev/null || { fail "$label (z-fasta err)"; return; }
    fi
    normalize_bed_output "$TMPDIR/bt.raw" "$bed" 0 "$TMPDIR/bt.norm"
    normalize_bed_output "$TMPDIR/zf.raw" "$bed" 0 "$TMPDIR/zf.norm"
    diff_oracle "$TMPDIR/bt.norm" "$TMPDIR/zf.norm" "$label"
}

# --- optional wrapper cross-checks ---
wrapper_vs() {
    local desc="$1" cmp="$2" bin="$3"; shift 3
    [[ -x "$bin" ]] || return 0
    if "$bin" "$@" > "$TMPDIR/wrapper.tmp" 2>/dev/null; then
        diff_fasta_oracle "$cmp" "$TMPDIR/wrapper.tmp" "$desc"
    else
        fail "$desc (wrapper err)"
    fi
}

wrapper_region() {
    local desc="$1" fasta="$2" region="$3"
    wrapper_vs "$desc noodles" "$TMPDIR/got.tmp" "$NOODLES" get "$fasta" "$region"
    wrapper_vs "$desc rustbio" "$TMPDIR/got.tmp" "$RUSTBIO" get "$fasta" "$region"
}

wrapper_regions() {
    local desc="$1" fasta="$2"; shift 2
    "$SAMTOOLS" faidx "$fasta" "$@" > "$TMPDIR/expected.tmp" 2>/dev/null || return 0
    wrapper_vs "$desc noodles" "$TMPDIR/expected.tmp" "$NOODLES" get "$fasta" "$@"
    wrapper_vs "$desc rustbio" "$TMPDIR/expected.tmp" "$RUSTBIO" get "$fasta" "$@"
}

wrapper_rc() {
    local desc="$1" fasta="$2" region="$3"
    local st
    st=$("$SAMTOOLS" faidx -i --mark-strand no "$fasta" "$region" 2>/dev/null) || return 0
    printf "%s\n" "$st" > "$TMPDIR/expected.tmp"
    wrapper_vs "$desc noodles" "$TMPDIR/expected.tmp" "$NOODLES" get "$fasta" "$region" --rc
    wrapper_vs "$desc rustbio" "$TMPDIR/expected.tmp" "$RUSTBIO" get "$fasta" "$region" --rc
}

wrapper_bed() {
    local desc="$1" bed="$2" fasta="${3:-$TEST_DATA_DIR/simple.fasta}"
    bed_to_regions "$bed" > "$TMPDIR/st_regions.txt"
    "$SAMTOOLS" faidx -r "$TMPDIR/st_regions.txt" "$fasta" > "$TMPDIR/expected.tmp" 2>/dev/null || return 0
    for pair in "noodles:$NOODLES" "rustbio:$RUSTBIO"; do
        IFS=: read -r name bin <<< "$pair"
        [[ -x "$bin" ]] || continue
        if "$bin" get "$fasta" --bed "$bed" > "$TMPDIR/${name}.raw" 2>/dev/null; then
            normalize_bed_output "$TMPDIR/expected.tmp" "$bed" 0 "$TMPDIR/st.norm"
            normalize_bed_output "$TMPDIR/${name}.raw" "$bed" 0 "$TMPDIR/${name}.norm"
            diff_oracle "$TMPDIR/st.norm" "$TMPDIR/${name}.norm" "$desc $name"
        else
            fail "$desc ($name err)"
        fi
    done
}

# --- fixtures ---
gen_names_file() { printf '# names file\nseq2\n\nseq1\nseq2\n' > "$1"; }

gen_bed_file() {
    local out="$1" count="$2"
    : > "$out"
    printf '# synthetic BED for %s\ntrack name=verify_bed_%s\n' "$count" "$count" >> "$out"
    cat >> "$out" <<'BED'
seq1	0	24	full_seq1	0	+
seq1	0	5	dup_a	0	+
seq1	0	5	dup_b	0	+
seq1	2	8	overlap_minus	0	-
seq2	0	12	full_seq2	0	-
seq2	0	4	short_plus	0	+
BED
    for ((i=6; i<count; i++)); do
        local idx=$i chrom len span start end strand
        ((idx%2==0)) && chrom=seq1 len=24 || chrom=seq2 len=12
        ((len<=4)) && span=1 || span=$((idx%4+1))
        start=$((idx%(len-span+1))); end=$((start+span))
        ((idx%3==0)) && strand=- || strand=+
        printf '%s\t%d\t%d\t%s_%d\t0\t%s\n' "$chrom" "$start" "$end" "$chrom" "$idx" "$strand" >> "$out"
    done
}

gen_iupac_fixture() {
    local out="$TMPDIR/iupac.fasta"
    if [[ ! -f "$out" ]]; then
        printf '>iupac_all\nACGTURYSWKMBDHVNacgturyswkmbdhvnu\n' > "$out"
    fi
    [[ -f "${out}.zfi" ]] || "$ZFASTA" index "$out" >/dev/null 2>&1
    printf '%s\n' "$out"
}

gen_chrom_fixture() {
    local total="${1:-16384}"
    local out="$TMPDIR/chrom_${total}.fasta"
    if [[ ! -f "$out" ]]; then
        {
            echo ">chrSynthetic"
            awk -v total="$total" 'BEGIN {
                pattern = "ACGTNRYWSKMBDHVacgtnrywskmbdhv"
                pattern_len = length(pattern)
                window = pattern
                while (length(window) < pattern_len + 71) window = window pattern
                for (offset = 0; offset < total; offset += 71) {
                    width = total - offset < 71 ? total - offset : 71
                    print substr(window, offset % pattern_len + 1, width)
                }
            }'
        } > "$out"
    fi
    [[ -f "${out}.zfi" ]] || "$ZFASTA" index "$out" >/dev/null 2>&1
    printf '%s\n' "$out"
}

gen_bed_rc_file() {
    local out="$1" stranded_out="$2"
    : > "$out"; : > "$stranded_out"
    for ((i=0; i<1200; i++)); do
        local s=$((1000+i*400)) e strand
        e=$((s+120))
        printf '%s\t%d\t%d\n' "chrSynthetic" "$s" "$e" >> "$out"
        ((i%2==0)) && strand=+ || strand=-
        printf '%s\t%d\t%d\tregion_%05d\t0\t%s\n' "chrSynthetic" "$s" "$e" "$i" "$strand" >> "$stranded_out"
    done
}

REG20=(); for i in {1..10}; do REG20+=("seq1:${i}-${i}" "seq2:${i}-${i}"); done
REG20_REV=(); for i in {10..1}; do REG20_REV+=("seq2:${i}-${i}" "seq1:${i}-${i}"); done

messy_bed_for() {
    case "$1" in
        mixed_widths) cat <<'BED'
mixed_widths	0	12	line1	0	+
mixed_widths	9	20	cross	0	+
BED
            ;;
        trailing_whitespace) cat <<'BED'
trailing_whitespace	0	8	start	0	+
trailing_whitespace	16	24	end	0	+
BED
            ;;
        blank_lines) cat <<'BED'
blank_lines	0	4	start	0	+
blank_lines	8	12	mid	0	+
BED
            ;;
        mixed_crlf) cat <<'BED'
mixed_crlf	0	8	start	0	+
mixed_crlf	6	14	cross	0	+
BED
            ;;
    esac
}

# --- sections ---
section0_index() {
    section_hdr 0 "Index path coverage (.zfi vs .fai)"
    echo "  Uniform fixtures only; cache messy_fixtures stay .zfi-only ([extended:messy])"
    local simple="$TEST_DATA_DIR/simple.fasta" proteome="$TEST_DATA_DIR/proteome.fasta"
    local edge="$TEST_DATA_DIR/edge_cases.fasta" mixed="$TEST_DATA_DIR/mixed_widths.fasta"
    local bed_small="$TMPDIR/idx_bed_small.bed" bed_medium="$TMPDIR/idx_bed_medium.bed"
    local bed_large="$TMPDIR/idx_bed_large.bed"
    local mx_bed="$TMPDIR/idx_mx.bed" names="$TMPDIR/idx_names.txt"
    local bed_rc="$TMPDIR/idx_bed_rc.bed"
    local -a p16=() reg15=() reg1024=()
    local prot1="sp|P12345|PROT_HUMAN" prot2="sp|Q98765|ANOT_MOUSE"

    gen_bed_file "$bed_small" 10
    gen_bed_file "$bed_medium" 100
    gen_bed_file "$bed_large" 1000
    gen_names_file "$names"
    for i in {1..8}; do
        p16+=("${prot1}:${i}-$((i + 2))" "${prot2}:${i}-$((i + 2))")
    done
    reg15=("${REG20[@]:0:15}")
    for _ in {1..64}; do reg1024+=("seq1:1-1"); done

    cat > "$mx_bed" <<'BED'
mixed1	54	75	mixed1_span	0	+
mixed2	74	95	mixed2_span	0	-
mixed3	57	72	mixed3_span	0	+
BED
    cat > "$bed_rc" <<'BED'
seq1	0	5	plus	0	+
seq1	0	5	minus	0	-
seq1	2	8	overlap_minus	0	-
seq2	0	4	short_plus	0	+
BED

    # Spot checks with samtools oracle on both index paths.
    verify_index positional "$simple" "seq1:1-10" "sub-region"
    verify_index positional "$simple" "seq1" "full sequence"
    verify_index positional "$TEST_DATA_DIR/single.fasta" "single_sequence:1-4" "single record"
    verify_index positional "$proteome" "sp|P12345|PROT_HUMAN:1-10" "pipe name"
    verify_index positional "$mixed" "mixed1:55-75" "mixed-width record"
    verify_index bed "$simple" "$bed_small" "small BED"

    echo "  [index:cross] single-region"
    verify_index_cross "sub-region" "$simple" "seq1:1-10"
    verify_index_cross "full sequence" "$simple" "seq1"
    verify_index_cross "open-ended region" "$simple" "seq1:10-"
    verify_index_cross "pipe name sub-region" "$proteome" "sp|P12345|PROT_HUMAN:1-10"
    verify_index_cross "mixed-width sub-region" "$mixed" "mixed1:55-75"
    verify_index_cross "duplicate name (first wins)" "$edge" "dupname"
    verify_index_cross "lowercase sub-region" "$edge" "lowercase:1-6"

    echo "  [index:cross] multi-region"
    verify_index_cross "two sub-regions" "$simple" "seq1:1-10" "seq1:13-24"
    verify_index_cross "full + sub-region" "$simple" "seq1" "seq2:3-10"
    verify_index_cross "15 regions (hash map; was scan threshold)" "$simple" "${reg15[@]}"
    verify_index_cross "16 regions (hash map path)" "$proteome" "${p16[@]}"
    verify_index_cross "20 regions (sort path)" "$simple" "${REG20[@]}"
    verify_index_cross "64 regions (hash map path)" "$simple" "${reg1024[@]}"

    echo "  [index:cross] BED and names"
    verify_index_cross "small BED" "$simple" --bed "$bed_small"
    verify_index_cross "medium BED" "$simple" --bed "$bed_medium"
    verify_index_cross "large BED" "$simple" --bed "$bed_large"
    verify_index_cross "BED via stdin" "$simple" "$bed_small" --bed -
    verify_index_cross "BED stranded" "$simple" --bed "$bed_small" --strand-aware
    verify_index_cross "BED --honor-strand alias" "$simple" --bed "$bed_small" --honor-strand
    verify_index_cross "large BED stranded" "$simple" --bed "$bed_large" --strand-aware
    verify_index_cross "mixed-width BED" "$mixed" --bed "$mx_bed"
    verify_index_cross "names file" "$simple" --names "$names"

    echo "  [index:cross] RC / complement / reverse"
    verify_index_cross "--rc" "$simple" "seq1:1-5" --rc
    verify_index_cross "--complement-only" "$simple" "seq1:1-5" --complement-only
    verify_index_cross "--reverse-only" "$simple" "seq1:1-5" --reverse-only
    verify_index_cross "multi --rc" "$simple" "seq1:1-5" "seq1:10-15" "seq1:20-24" --rc
    verify_index_cross "names --rc" "$simple" --names "$names" --rc
    verify_index_cross "BED --strand-aware --rc" "$simple" --bed "$bed_rc" --strand-aware --rc
    verify_index_cross "--rc --annotate-rc" "$simple" "seq1:1-5" --rc --annotate-rc
    verify_index_cross "--complement-only --annotate-rc" "$simple" "seq1:1-5" --complement-only --annotate-rc
    verify_index_cross "IUPAC --rc" "$(gen_iupac_fixture)" "iupac_all:1-33" --rc
    verify_index_cross "IUPAC --complement-only" "$(gen_iupac_fixture)" "iupac_all:1-33" --complement-only

    verify_messy_zfi_required mixed_widths 1 8
}

section1() {
    section_hdr 1 "Single-region extraction"

    test_seq() {
        local fasta="$1" name="$2" len="$3"
        local label="${fasta##*/}:${name}"
        parity_samtools_region "$fasta" "$name" "$label full"
        wrapper_region "$label full wrapper" "$fasta" "$name"
        if ((len>=10)); then
            parity_samtools_region "$fasta" "${name}:1-10" "$label :1-10"
            wrapper_region "$label :1-10 wrapper" "$fasta" "${name}:1-10"
            local ls=$((len-9))
            parity_samtools_region "$fasta" "${name}:${ls}-${len}" "$label :${ls}-${len}"
            wrapper_region "$label :${ls}-${len} wrapper" "$fasta" "${name}:${ls}-${len}"
        else
            parity_samtools_region "$fasta" "${name}:1-${len}" "$label :1-${len}"
        fi
        if ((len>=101)); then
            local ms=$((len/2-50)) me=$((len/2+50))
            ((ms<1)) && ms=1; ((me>len)) && me=$len
            parity_samtools_region "$fasta" "${name}:${ms}-${me}" "$label mid-span"
            wrapper_region "$label mid-span wrapper" "$fasta" "${name}:${ms}-${me}"
        fi
        parity_samtools_region "$fasta" "${name}:1-1" "$label :1-1"
        parity_samtools_region "$fasta" "${name}:${len}-${len}" "$label :${len}-${len}"
        parity_samtools_region "$fasta" "${name}:1-${len}" "$label :1-${len}"
        parity_samtools_region "$fasta" "${name}:1-$((len+100))" "$label clamp"
    }

    test_file() {
        local fasta="$1"
        [[ -f "${fasta}.fai" ]] || return
        echo "  ${fasta##*/}"
        while IFS=$'\t' read -r name length _; do
            ((length>0)) && test_seq "$fasta" "$name" "$length"
        done < "${fasta}.fai"
    }

    for f in "$TEST_DATA_DIR/simple.fasta" "$TEST_DATA_DIR/proteome.fasta" "$TEST_DATA_DIR/single.fasta" "$TEST_DATA_DIR/edge_cases.fasta"; do
        [[ -f "$f" ]] && { ensure_index "$f"; test_file "$f"; }
    done
    cp "$TEST_DATA_DIR/mixed_widths.fasta" "$TMPDIR/mixed_widths.fasta"
    ensure_index "$TMPDIR/mixed_widths.fasta"
    test_file "$TMPDIR/mixed_widths.fasta"

    verify_open_ended_region "$TEST_DATA_DIR/simple.fasta" seq1 10 24 "simple open-ended region"

    echo "  cache/messy_fixtures [extended:messy]"
    for case in "${MESSY_POS_CASES[@]}"; do
        IFS=: read -r name regions <<< "$case"
        local tgt="$TMPDIR/messy_${name}.fasta"
        [[ -f "$MESSY_TEST_DIR/${name}.fasta" ]] || { fail "[extended:messy] missing fixture $name"; continue; }
        prepare_messy_fixture "$name" "$tgt" || { fail "[extended:messy] index failed $name"; continue; }
        IFS=, read -ra spans <<< "$regions"
        for span in "${spans[@]}"; do
            IFS=- read -r start end <<< "$span"
            verify_messy_region "$tgt" "$name" "$start" "$end" "$name:${start}-${end}"
        done
    done
}

section2() {
    section_hdr 2 "Multi-region extraction"
    local simple="$TEST_DATA_DIR/simple.fasta" proteome="$TEST_DATA_DIR/proteome.fasta" edge="$TEST_DATA_DIR/edge_cases.fasta"
    ensure_index "$simple"; ensure_index "$proteome"; ensure_index "$edge"
    local LONG="AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"

    echo "  simple.fasta"
    parity_samtools_regions "$simple" "same seq, two sub-regions" "seq1:1-10" "seq1:13-24"
    parity_samtools_regions "$simple" "seqs file order" "seq1:1-12" "seq2:1-6"
    parity_samtools_regions "$simple" "seqs reversed order" "seq2:1-6" "seq1:1-12"
    parity_samtools_regions "$simple" "overlapping" "seq1:1-15" "seq1:10-24"
    parity_samtools_regions "$simple" "duplicate" "seq1:1-12" "seq1:1-12"
    parity_samtools_regions "$simple" "full + sub" "seq1" "seq2:3-10"
    parity_samtools_regions "$simple" "both full" "seq1" "seq2"
    parity_samtools_regions "$simple" "triple duplicate" "seq1:1-5" "seq1:1-5" "seq1:1-5"
    parity_samtools_regions "$simple" "full + 2 sub" "seq1" "seq2:1-6" "seq1:10-20"
    parity_samtools_regions "$simple" "all reversed" "seq2" "seq1"
    parity_samtools_regions "$simple" "20 regions sort path" "${REG20[@]}"
    parity_samtools_regions "$simple" "20 regions reversed" "${REG20_REV[@]}"
    parity_samtools_regions "$simple" "15 regions (hash map)" "${REG20[@]:0:15}"
    wrapper_regions "multi same seq" "$simple" "seq1:1-10" "seq1:13-24"
    wrapper_regions "multi file order" "$simple" "seq1:1-12" "seq2:1-6"
    wrapper_regions "multi 20 regions" "$simple" "${REG20[@]}"

    echo "  edge_cases.fasta"
    parity_samtools_regions "$edge" "long + normal" "single_line" "$LONG"
    parity_samtools_regions "$edge" "long reversed" "$LONG" "single_line"
    parity_samtools_regions "$edge" "lowercase + nonstandard" "lowercase" "nonstandard"
    parity_samtools_regions "$edge" "overlapping same seq" "single_line:1-5" "single_line:3-8"

    echo "  proteome.fasta"
    parity_samtools_regions "$proteome" "pipe names" "sp|P12345|PROT_HUMAN" "sp|Q98765|ANOT_MOUSE"
    parity_samtools_regions "$proteome" "pipe reversed" "sp|Q98765|ANOT_MOUSE" "sp|P12345|PROT_HUMAN"
    parity_samtools_regions "$proteome" "pipe sub + full" "sp|P12345|PROT_HUMAN:1-10" "sp|Q98765|ANOT_MOUSE"
    parity_samtools_regions "$proteome" "pipe duplicate" "sp|P12345|PROT_HUMAN" "sp|P12345|PROT_HUMAN"
    local prot1="sp|P12345|PROT_HUMAN" prot2="sp|Q98765|ANOT_MOUSE" P16=()
    for i in {1..8}; do P16+=("${prot1}:${i}-$((i+2))" "${prot2}:${i}-$((i+2))"); done
    parity_samtools_regions "$proteome" "16 proteome sub-regions" "${P16[@]}"

    ALL_EDGE=(); while IFS=$'\t' read -r name _; do ALL_EDGE+=("$name"); done < "${edge}.fai"
    (( ${#ALL_EDGE[@]} > 0 )) && {
        parity_samtools_regions "$edge" "all edge seqs" "${ALL_EDGE[@]}"
        REV_EDGE=(); for ((i=${#ALL_EDGE[@]}-1; i>=0; i--)); do REV_EDGE+=("${ALL_EDGE[$i]}"); done
        parity_samtools_regions "$edge" "all edge reversed" "${REV_EDGE[@]}"
    }

    REG1024=()
    for _ in {1..1024}; do REG1024+=("seq1:1-1"); done
    parity_samtools_regions "$simple" "1024 positional regions" "${REG1024[@]}"
    REG1025=("${REG1024[@]}" "seq1:1-1")
    expect_fail "1025 positional regions rejected" "$ZFASTA" get "$simple" "${REG1025[@]}"
}

section3() {
    section_hdr 3 "BED and names extraction"
    local simple="$TEST_DATA_DIR/simple.fasta"
    ensure_index "$simple"
    local BED_SMALL="$TMPDIR/small.bed" BED_MEDIUM="$TMPDIR/medium.bed"

    run_size_suite() {
        local label="$1" count="$2"
        local bed="$TMPDIR/${label}.bed"
        gen_bed_file "$bed" "$count"
        echo "  suite $label ($count rows)"
        parity_bedtools "$label default" "$bed" 0 "$simple"
        parity_bedtools "$label stranded" "$bed" 1 "$simple"
        parity_bed_samtools "$label vs samtools" "$bed" "$simple"
        parity_bedtools "$label via stdin" "$bed" 0 "$simple" 1
        parity_bed_samtools "$label stdin vs samtools" "$bed" "$simple" 1
        parity_bedtools "$label stranded stdin" "$bed" 1 "$simple" 1
    }

    run_size_suite small 10
    run_size_suite medium 100
    run_size_suite large 1000

    if command -v "$BEDTOOLS" &>/dev/null; then
        parity_bedtools "--honor-strand alias" "$BED_SMALL" 1 "$simple" 0 --honor-strand
    fi

    "$ZFASTA" get "$simple" seq1:1-10 seq2:1-6 --summary > "$TMPDIR/got.tmp" 2>"$TMPDIR/summary.err" \
        && grep -Eq '^summary: regions=[0-9]+ total_bases=[0-9]+ elapsed_s=' "$TMPDIR/summary.err" \
        && pass "[extended:summary] --summary stderr" || fail "[extended:summary] --summary stderr"

    local NAMES="$TMPDIR/names.txt"
    gen_names_file "$NAMES"
    names_to_regions "$NAMES" > "$TMPDIR/st_names.txt"
    "$ZFASTA" get "$simple" --names "$NAMES" > "$TMPDIR/got.tmp" 2>/dev/null \
        && "$SAMTOOLS" faidx -r "$TMPDIR/st_names.txt" "$simple" > "$TMPDIR/expected.tmp" 2>/dev/null \
        && diff_fasta_oracle "$TMPDIR/expected.tmp" "$TMPDIR/got.tmp" "[parity:samtools] names file" \
        || fail "[parity:samtools] names file"
    cat "$NAMES" | "$ZFASTA" get "$simple" --names - > "$TMPDIR/got.tmp" 2>/dev/null \
        && diff_fasta_oracle "$TMPDIR/expected.tmp" "$TMPDIR/got.tmp" "[parity:samtools] names stdin" \
        || fail "[parity:samtools] names stdin"

    printf '# names file\r\nseq2\r\n\r\nseq1' > "$TMPDIR/names_crlf.txt"
    printf 'seq2\nseq1\n' > "$TMPDIR/st_names_crlf.txt"
    "$ZFASTA" get "$simple" --names "$TMPDIR/names_crlf.txt" > "$TMPDIR/got.tmp" 2>/dev/null \
        && "$SAMTOOLS" faidx -r "$TMPDIR/st_names_crlf.txt" "$simple" > "$TMPDIR/expected.tmp" 2>/dev/null \
        && diff_fasta_oracle "$TMPDIR/expected.tmp" "$TMPDIR/got.tmp" "[parity:samtools] names CRLF and final line" \
        || fail "[parity:samtools] names CRLF and final line"

    names_to_regions "$NAMES" > "$TMPDIR/st_names.txt"
    "$SAMTOOLS" faidx -r "$TMPDIR/st_names.txt" "$simple" > "$TMPDIR/expected.tmp" 2>/dev/null
    "$ZFASTA" get "$simple" --names "$NAMES" --summary > "$TMPDIR/got.tmp" 2>"$TMPDIR/summary.err" \
        && grep -Eq '^summary: regions=3 total_bases=48 elapsed_s=' "$TMPDIR/summary.err" \
        && diff_fasta_oracle "$TMPDIR/expected.tmp" "$TMPDIR/got.tmp" "[extended:summary] names output and summary" \
        || fail "[extended:summary] names output and summary"

    "$ZFASTA" get "$simple" --bed "$BED_SMALL" --summary > "$TMPDIR/got.tmp" 2>"$TMPDIR/summary.err" \
        && grep -Eq '^summary: regions=10 total_bases=[0-9]+ elapsed_s=' "$TMPDIR/summary.err" \
        && pass "[extended:summary] BED output and summary" \
        || fail "[extended:summary] BED output and summary"

    cat > "$TMPDIR/annotate.bed" <<'BED'
seq1	0	5	plus	0	+
seq1	0	5	minus	0	-
BED
    printf '>seq1:1-5\nACGTA\n>seq1:1-5 (reverse complement)\nTACGT\n' > "$TMPDIR/annotate_expected.fa"
    "$ZFASTA" get "$simple" --bed "$TMPDIR/annotate.bed" --strand-aware --annotate-rc > "$TMPDIR/got.tmp" 2>/dev/null \
        && diff_oracle "$TMPDIR/annotate_expected.fa" "$TMPDIR/got.tmp" "[extended:header] BED final orientation annotation" \
        || fail "[extended:header] BED final orientation annotation"
    printf '>seq1:1-5 (reverse complement)\nTACGT\n>seq1:1-5\nACGTA\n' > "$TMPDIR/annotate_expected.fa"
    "$ZFASTA" get "$simple" --bed "$TMPDIR/annotate.bed" --strand-aware --rc --annotate-rc > "$TMPDIR/got.tmp" 2>/dev/null \
        && diff_oracle "$TMPDIR/annotate_expected.fa" "$TMPDIR/got.tmp" "[extended:header] BED composed identity annotation" \
        || fail "[extended:header] BED composed identity annotation"

    local mx="$TMPDIR/mx.fasta" mx_bed="$TMPDIR/mx.bed"
    cp "$TEST_DATA_DIR/mixed_widths.fasta" "$mx"
    ensure_index "$mx"
    cat > "$mx_bed" <<'BED'
mixed1	54	75	mixed1_span	0	+
mixed2	74	95	mixed2_span	0	-
mixed3	57	72	mixed3_span	0	+
BED
    command -v "$BEDTOOLS" &>/dev/null && {
        parity_bedtools "mixed-width default" "$mx_bed" 0 "$mx"
        parity_bedtools "mixed-width stranded" "$mx_bed" 1 "$mx"
    }
    parity_bed_samtools "mixed-width vs samtools" "$mx_bed" "$mx"

    for messy_name in mixed_widths trailing_whitespace blank_lines mixed_crlf; do
        local mf="$TMPDIR/messy_bed_${messy_name}.fasta" bed="$TMPDIR/messy_bed_${messy_name}.bed" exp="$TMPDIR/messy_bed_${messy_name}.fa"
        prepare_messy_fixture "$messy_name" "$mf" || { fail "[extended:messy] bed index $messy_name"; continue; }
        messy_bed_for "$messy_name" > "$bed"
        oracle bed "$mf" "$bed" "$exp"
        verify_expected_bed "messy BED $messy_name" "$mf" "$bed" "$exp"
    done

    local nu="$TMPDIR/nu_rc_bed.fasta" nu_bed="$TMPDIR/nu_rc_bed.bed" nu_exp="$TMPDIR/nu_rc_bed_expected.fa"
    prepare_messy_fixture mixed_widths "$nu" || fail "[extended:messy] stranded RC index"
    cat > "$nu_bed" <<'BED'
mixed_widths	9	20	plus	0	+
mixed_widths	14	28	minus	0	-
BED
    oracle bed-stranded-rc "$nu" "$nu_bed" "$nu_exp"
    verify_expected_bed "messy BED --strand-aware --rc" "$nu" "$nu_bed" "$nu_exp" --strand-aware --rc

    wrapper_bed "small BED" "$BED_SMALL"
    wrapper_bed "medium BED" "$BED_MEDIUM"
    wrapper_bed "mixed-width BED" "$mx_bed" "$mx"

    printf 'seq1:1-5\n' > "$TMPDIR/names_reg.txt"
    expect_fail "names file region syntax rejected" "$ZFASTA" get "$simple" --names "$TMPDIR/names_reg.txt"
    printf 'seq1\t5\t5\tzero_len\t0\t+\n' > "$TMPDIR/bad.bed"
    expect_fail "0-length BED rejected" "$ZFASTA" get "$simple" --bed "$TMPDIR/bad.bed"
    printf 'seq1\t5\n' > "$TMPDIR/bad.bed"
    expect_fail "short BED rejected" "$ZFASTA" get "$simple" --bed "$TMPDIR/bad.bed"
    printf 'seq1\t0\t5\tname\t0\t?\n' > "$TMPDIR/bad.bed"
    expect_fail "invalid BED strand rejected" "$ZFASTA" get "$simple" --bed "$TMPDIR/bad.bed" --strand-aware
}

section4() {
    section_hdr 4 "Reverse complement"
    local simple="$TEST_DATA_DIR/simple.fasta" single="$TEST_DATA_DIR/single.fasta"
    local edge="$TEST_DATA_DIR/edge_cases.fasta"
    local iupac chrom
    iupac=$(gen_iupac_fixture); chrom=$(gen_chrom_fixture 500000)
    ensure_index "$iupac"; ensure_index "$chrom"

    for r in seq1:1-5 seq1:10-15 seq2:1-12; do
        parity_samtools_rc "$simple" "$r" "simple $r"
        wrapper_rc "simple $r wrapper" "$simple" "$r"
    done
    parity_samtools_rc "$single" "single_sequence:1-4" "single span"
    parity_samtools_rc "$edge" "lowercase:1-12" "lowercase preservation"
    parity_samtools_rc "$iupac" "iupac_all:1-33" "IUPAC full"
    parity_samtools_rc "$chrom" "chrSynthetic" "chrom full"

    local mx="$TMPDIR/mx_rc.fasta"
    cp "$TEST_DATA_DIR/mixed_widths.fasta" "$mx"
    ensure_index "$mx"
    parity_samtools_rc "$mx" "mixed1:55-75" "mixed-width region"

    for spec in "mixed_widths:3:24" "trailing_whitespace:7:20" "blank_lines:9:12" "mixed_crlf:7:14"; do
        IFS=: read -r n s e <<< "$spec"
        local mf="$TMPDIR/messy_rc_${n}.fasta"
        prepare_messy_fixture "$n" "$mf" || { fail "[extended:messy] RC index $n"; continue; }
        verify_messy_region "$mf" "$n" "$s" "$e" "messy RC $n:${s}-${e}" --rc
    done

    verify_oracle_transform rev "[parity:seq] simple --reverse-only" "$simple" seq1 1 5 --reverse-only
    verify_oracle_transform comp "[parity:seq] simple --complement-only" "$simple" seq1 1 5 --complement-only
    printf '>seq1:1-5 (reverse)\nATGCA\n' > "$TMPDIR/rev_annot.fa"
    verify_exact_get "[extended:header] annotated --reverse-only" "$TMPDIR/rev_annot.fa" "$simple" seq1:1-5 --reverse-only --annotate-rc
    printf '>seq1:1-5 (complement)\nTGCAT\n' > "$TMPDIR/comp_annot.fa"
    verify_exact_get "[extended:header] annotated --complement-only" "$TMPDIR/comp_annot.fa" "$simple" seq1:1-5 --complement-only --annotate-rc

    printf '>iupac_all:1-33 (complement)\nTGCAAYRSWMKVHDBNtgcaayrswmkvhdbna\n' > "$TMPDIR/iupac_comp.fa"
    "$ZFASTA" get "$iupac" iupac_all:1-33 --complement-only --annotate-rc > "$TMPDIR/got.tmp" 2>/dev/null \
        && diff_oracle "$TMPDIR/iupac_comp.fa" "$TMPDIR/got.tmp" "[extended:header] IUPAC --complement-only annotated" \
        || fail "[extended:header] IUPAC --complement-only annotated"

    : > "$TMPDIR/expected.tmp"
    for r in seq1:1-5 seq1:10-15 seq1:20-24; do "$ZFASTA" get "$simple" "$r" --rc >> "$TMPDIR/expected.tmp" 2>/dev/null; done
    "$ZFASTA" get "$simple" seq1:1-5 seq1:10-15 seq1:20-24 --rc > "$TMPDIR/got.tmp" 2>/dev/null
    diff_oracle "$TMPDIR/expected.tmp" "$TMPDIR/got.tmp" "[parity:samtools] multi RC self-consistency"

    local nf="$TMPDIR/names_rc.txt"
    printf 'seq2\nseq1\n' > "$nf"
    names_to_regions "$nf" > "$TMPDIR/st_rc_names.txt"
    "$SAMTOOLS" faidx -i --mark-strand no -r "$TMPDIR/st_rc_names.txt" "$simple" > "$TMPDIR/expected.tmp" 2>/dev/null \
        && "$ZFASTA" get "$simple" --names "$nf" --rc > "$TMPDIR/got.tmp" 2>/dev/null \
        && diff_fasta_oracle "$TMPDIR/expected.tmp" "$TMPDIR/got.tmp" "[parity:samtools] names --rc" \
        || fail "[parity:samtools] names --rc"

    local bed_comp="$TMPDIR/bed_comp.bed" bed_exp="$TMPDIR/bed_comp_expected.fa"
    cat > "$bed_comp" <<'BED'
seq1	0	5	plus	0	+
seq1	0	5	minus	0	-
seq1	2	8	overlap_minus	0	-
seq2	0	4	short_plus	0	+
BED
    printf '>seq1:1-5\nTACGT\n>seq1:1-5\nACGTA\n>seq1:3-8\nGTACGT\n>seq2:1-4\nCCCC\n' > "$bed_exp"
    "$ZFASTA" get "$simple" --bed "$bed_comp" --strand-aware --rc > "$TMPDIR/got.tmp" 2>/dev/null \
        && diff_oracle "$bed_exp" "$TMPDIR/got.tmp" "[parity:seq] BED --strand-aware --rc exact" \
        || fail "[parity:seq] BED --strand-aware --rc exact"

    printf '>seq1:1-5 (reverse complement)\nTACGT\n' > "$TMPDIR/annot_rc.fa"
    "$ZFASTA" get "$simple" seq1:1-5 --rc --annotate-rc > "$TMPDIR/got.tmp" 2>/dev/null \
        && diff_oracle "$TMPDIR/annot_rc.fa" "$TMPDIR/got.tmp" "[extended:header] annotated --rc header" \
        || fail "[extended:header] annotated --rc header"

    if command -v "$BEDTOOLS" &>/dev/null && [[ -x "$SEQTK" ]]; then
        parity_bed_rc_composed "BED stranded RC" "$simple" "$bed_comp" 0
        parity_bed_rc_composed "BED stranded RC stdin" "$simple" "$bed_comp" 1
        gen_bed_rc_file "$TMPDIR/bed_rc.bed" "$TMPDIR/bed_str.bed"
        parity_bed_rc_composed "BED stranded RC chrom" "$chrom" "$TMPDIR/bed_str.bed" 0
    fi

    oracle rc "$iupac" iupac_all 1 33 "$TMPDIR/iupac_rc_py.fa"
    "$ZFASTA" get "$iupac" iupac_all:1-33 --rc > "$TMPDIR/got.tmp" 2>/dev/null \
        && diff_oracle "$TMPDIR/iupac_rc_py.fa" "$TMPDIR/got.tmp" "[parity:seq] IUPAC --rc Python oracle" \
        || fail "[parity:seq] IUPAC --rc Python oracle"

    "$ZFASTA" get "$TEST_DATA_DIR/proteome.fasta" 'sp|P12345|PROT_HUMAN:1-10' --rc > "$TMPDIR/got.tmp" 2> "$TMPDIR/got.err" \
        && fail "protein RC rejected (should fail)" \
        || grep -q 'reverse complement is not defined' "$TMPDIR/got.err" \
        && pass "protein RC rejected" || fail "protein RC rejected (wrong err)"
}

section5() {
    section_hdr 5 "Edge cases and error paths"
    local E="$TEST_DATA_DIR"

    verify_edge_case() {
        local desc="$1" fasta="$2" region="$3" exp_zf="$4" exp_st="$5"
        local exp_noo="${6:--}" exp_rb="${7:--}"
        local tag="[parity:samtools]"
        (( exp_st == 1 && exp_zf == 0 )) && tag="[extended:messy]"

        "$ZFASTA" index "$fasta" 2>/dev/null || true
        [[ -f "${fasta}.fai" ]] || "$SAMTOOLS" faidx "$fasta" 2>/dev/null || true

        local zf_exit=0 st_exit=0 zf_out st_out
        zf_out=$("$ZFASTA" get "$fasta" "$region" 2>/dev/null) || zf_exit=$?
        st_out=$("$SAMTOOLS" faidx "$fasta" "$region" 2>/dev/null) || st_exit=$?

        if (( (exp_zf == 0 && zf_exit == 0) || (exp_zf == 1 && zf_exit != 0) )); then pass "$tag $desc z-fasta exit ok"
        else fail "$tag $desc z-fasta (expected exit $exp_zf, got $zf_exit)"; fi
        if (( (exp_st == 0 && st_exit == 0) || (exp_st == 1 && st_exit != 0) )); then pass "$tag $desc samtools exit ok"
        else fail "$tag $desc samtools (expected exit $exp_st, got $st_exit)"; fi

        if (( zf_exit == 0 && st_exit == 0 )); then
            printf "%s" "$zf_out" > "$TMPDIR/zf_edge.tmp"
            printf "%s" "$st_out" > "$TMPDIR/st_edge.tmp"
            diff_fasta_oracle "$TMPDIR/st_edge.tmp" "$TMPDIR/zf_edge.tmp" "$tag $desc FASTA parity"
        fi
        if (( zf_exit == 0 && st_exit != 0 )); then
            if parse_region_spec "$region"; then
                oracle region "$fasta" "$PARSE_NAME" "$PARSE_START" "$PARSE_END" "$TMPDIR/expected.tmp"
                printf "%s\n" "$zf_out" > "$TMPDIR/zf_edge.tmp"
                diff_oracle "$TMPDIR/expected.tmp" "$TMPDIR/zf_edge.tmp" "$tag $desc Python oracle"
            else
                fail "$tag $desc Python oracle skipped (region has no coords)"
            fi
        fi
        if [[ -x "$NOODLES" && "$exp_noo" != "-" ]]; then
            local noo_exit=0; noo_out=$("$NOODLES" get "$fasta" "$region" 2>/dev/null) || noo_exit=$?
            (( (exp_noo == 0 && noo_exit == 0) || (exp_noo == 1 && noo_exit != 0) )) \
                && pass "$tag $desc noodles exit ok" || fail "$tag $desc noodles (expected $exp_noo, got $noo_exit)"
        fi
        if [[ -x "$RUSTBIO" && "$exp_rb" != "-" ]]; then
            local rb_exit=0; rb_out=$("$RUSTBIO" get "$fasta" "$region" 2>/dev/null) || rb_exit=$?
            (( (exp_rb == 0 && rb_exit == 0) || (exp_rb == 1 && rb_exit != 0) )) \
                && pass "$tag $desc rustbio exit ok" || fail "$tag $desc rustbio (expected $exp_rb, got $rb_exit)"
        fi
    }

    verify_edge_case "empty.fasta" "$E/empty.fasta" "nonexistent" 1 1 1 1
    verify_edge_case "not_fasta.txt" "$E/not_fasta.txt" "simple" 1 1 1 1
    verify_edge_case "edge_cases:single_line" "$E/edge_cases.fasta" "single_line" 0 0 0 0
    verify_edge_case "edge_cases:empty_seq" "$E/edge_cases.fasta" "empty_seq" 1 1 1 1
    verify_edge_case "edge_cases:dupname" "$E/edge_cases.fasta" "dupname" 0 0 0 0
    verify_edge_case "edge_cases:lowercase" "$E/edge_cases.fasta" "lowercase:1-6" 0 0 0 0
    verify_edge_case "edge_cases:nonstandard" "$E/edge_cases.fasta" "nonstandard" 0 0 0 0

    for case in "${MESSY_POS_CASES[@]}"; do
        IFS=: read -r name regions <<< "$case"
        IFS=, read -r first _ <<< "$regions"
        IFS=- read -r start end <<< "$first"
        local tgt="$TMPDIR/e_${name}.fasta"
        [[ -f "$MESSY_TEST_DIR/${name}.fasta" ]] || continue
        cp "$MESSY_TEST_DIR/${name}.fasta" "$tgt"
        verify_edge_case "messy:$name:$start-$end" "$tgt" "${name}:${start}-${end}" 0 1 1 1
    done

    local dup="$TMPDIR/dup_edge.fasta"
    printf '>dupname\nAAAA\n>other\nGG\n>dupname\nCCCCCC\n' > "$dup"
    "$ZFASTA" index --no-dedup "$dup" >/dev/null 2>&1 || fail "[extended:dedup] index"
    "$ZFASTA" index --emit-fai --no-dedup "$dup" > "${dup}.fai" 2>/dev/null || fail "[extended:dedup] fai index"
    printf '>dupname\nAAAA\n' > "$TMPDIR/dup_expected.fa"
    printf '>dupname:1-4\nAAAA\n' > "$TMPDIR/dup_bed_expected.fa"
    printf 'dupname\n' > "$TMPDIR/dup_names.txt"
    printf 'dupname\t0\t4\n' > "$TMPDIR/dup.bed"
    for sidecar in zfi fai; do
        if [[ "$sidecar" == fai ]]; then
            mv "${dup}.zfi" "$TMPDIR/dup.zfi"
        fi
        "$ZFASTA" get "$dup" dupname > "$TMPDIR/got.tmp" 2>/dev/null \
            && diff_oracle "$TMPDIR/dup_expected.fa" "$TMPDIR/got.tmp" "[extended:dedup] $sidecar positional first record" \
            || fail "[extended:dedup] $sidecar positional first record"
        "$ZFASTA" get "$dup" --names "$TMPDIR/dup_names.txt" > "$TMPDIR/got.tmp" 2>/dev/null \
            && diff_oracle "$TMPDIR/dup_expected.fa" "$TMPDIR/got.tmp" "[extended:dedup] $sidecar names first record" \
            || fail "[extended:dedup] $sidecar names first record"
        "$ZFASTA" get "$dup" --bed "$TMPDIR/dup.bed" > "$TMPDIR/got.tmp" 2>/dev/null \
            && diff_oracle "$TMPDIR/dup_bed_expected.fa" "$TMPDIR/got.tmp" "[extended:dedup] $sidecar BED first record" \
            || fail "[extended:dedup] $sidecar BED first record"
        if [[ "$sidecar" == fai ]]; then
            mv "$TMPDIR/dup.zfi" "${dup}.zfi"
        fi
    done

    local empty_name="$TMPDIR/empty_name.fasta"
    printf '>\nAAAA\n>\nCCCCCC\n' > "$empty_name"
    "$ZFASTA" index --no-dedup "$empty_name" >/dev/null 2>&1 || fail "[extended:empty-name] zfi index"
    "$ZFASTA" index --emit-fai --no-dedup "$empty_name" > "${empty_name}.fai" 2>/dev/null || fail "[extended:empty-name] fai index"
    printf '>\nAAAA\n' > "$TMPDIR/empty_name_expected.fa"
    "$ZFASTA" get "$empty_name" "" > "$TMPDIR/got.tmp" 2>/dev/null \
        && diff_oracle "$TMPDIR/empty_name_expected.fa" "$TMPDIR/got.tmp" "[extended:empty-name] zfi positional first record" \
        || fail "[extended:empty-name] zfi positional first record"
    mv "${empty_name}.zfi" "$TMPDIR/empty_name.zfi"
    "$ZFASTA" get "$empty_name" "" > "$TMPDIR/got.tmp" 2>/dev/null \
        && diff_oracle "$TMPDIR/empty_name_expected.fa" "$TMPDIR/got.tmp" "[extended:empty-name] fai positional first record" \
        || fail "[extended:empty-name] fai positional first record"
    "$SAMTOOLS" faidx "$empty_name" "" > "$TMPDIR/got.tmp" 2>/dev/null \
        && diff_fasta_oracle "$TMPDIR/empty_name_expected.fa" "$TMPDIR/got.tmp" "[parity:samtools] empty identifier" \
        || fail "[parity:samtools] empty identifier"
    mv "$TMPDIR/empty_name.zfi" "${empty_name}.zfi"

    printf '\n# comment\n\r\n' > "$TMPDIR/blank_names.txt"
    expect_fail "blank names input has no regions" "$ZFASTA" get "$E/simple.fasta" --names "$TMPDIR/blank_names.txt"
    printf '\t0\t4\n' > "$TMPDIR/empty_chrom.bed"
    expect_fail "empty BED chromosome rejected" "$ZFASTA" get "$E/simple.fasta" --bed "$TMPDIR/empty_chrom.bed"

    expect_fail "get on nonexistent seq" "$ZFASTA" get "$E/simple.fasta" NOSUCHSEQ
    expect_fail "get on invalid region" "$ZFASTA" get "$E/simple.fasta" seq1:invalid
    expect_fail "get with zero args" "$ZFASTA" get "$E/simple.fasta"
    expect_fail "mutually exclusive --rc and --complement-only" "$ZFASTA" get "$E/simple.fasta" seq1:1-5 --rc --complement-only
    expect_fail "mutually exclusive --rc and --reverse-only" "$ZFASTA" get "$E/simple.fasta" seq1:1-5 --rc --reverse-only
    expect_fail "mutually exclusive --complement-only and --reverse-only" "$ZFASTA" get "$E/simple.fasta" seq1:1-5 --complement-only --reverse-only
    expect_fail "--annotate-rc without transform" "$ZFASTA" get "$E/simple.fasta" seq1:1-5 --annotate-rc
    expect_fail "mixed positional and names sources" "$ZFASTA" get "$E/simple.fasta" seq1 --names "$TMPDIR/blank_names.txt"
    expect_fail "repeated BED source" "$ZFASTA" get "$E/simple.fasta" --bed "$TMPDIR/empty_chrom.bed" --bed "$TMPDIR/empty_chrom.bed"
    expect_fail "strand handling without BED" "$ZFASTA" get "$E/simple.fasta" seq1 --strand-aware
}

# --- run_tests ---
run_tests() {
    PASS=0
    FAIL=0
    local verify_tmp="$SCRIPT_DIR/.verify_tmp"
    rm -rf "$verify_tmp"
    mkdir -p "$verify_tmp"
    TMPDIR="$verify_tmp"
    TEST_DATA_DIR="$TMPDIR/tests_data"
    mkdir -p "$TEST_DATA_DIR"
    cp -a "$PROJECT_ROOT/tests/data/." "$TEST_DATA_DIR/"
    rm -f "$TEST_DATA_DIR"/*.zfi "$TEST_DATA_DIR"/*.fai
    bench_ensure_messy --fixtures


    echo "z-fasta get verification"
    echo "  z-fasta:  $ZFASTA"
    echo "  samtools: $(command -v "$SAMTOOLS" 2>/dev/null || echo missing)"
    echo "  bedtools: $(command -v "$BEDTOOLS" 2>/dev/null || echo missing)"
    echo "  skip: index=$SKIP_INDEX get=$SKIP_GET multi=$SKIP_MULTI bed=$SKIP_BED rc=$SKIP_RC edge=$SKIP_EDGE"

    cd "$PROJECT_ROOT"
    [[ -x "$ZFASTA" ]] || { echo "Error: z-fasta not found"; rm -rf "$TMPDIR"; return 1; }
    command -v "$SAMTOOLS" &>/dev/null || { echo "Error: samtools not found"; rm -rf "$TMPDIR"; return 1; }

    ! $SKIP_INDEX && section0_index
    ! $SKIP_GET   && section1
    ! $SKIP_MULTI && section2
    ! $SKIP_BED   && section3
    ! $SKIP_RC    && section4
    ! $SKIP_EDGE  && section5

    echo ""
    echo "Results: $PASS passed, $FAIL failed"
    if [[ "$FAIL" -gt 0 ]]; then
        echo "VERIFICATION FAILED"
        rm -rf "$TMPDIR"
        return 1
    fi
    echo "ALL PASSED"
    rm -rf "$TMPDIR"
    return 0
}

focused_names_bases() {
    python3 - "$1" "$2" <<'PY'
import sys

fai_path, names_path = sys.argv[1:]
lengths = {}
with open(fai_path + ".fai", encoding="ascii") as handle:
    for line in handle:
        name, length, *_rest = line.rstrip("\n").split("\t")
        lengths.setdefault(name, int(length))

total = 0
with open(names_path, encoding="ascii") as handle:
    for line in handle:
        name = line.rstrip("\n").removesuffix("\r")
        if not name or name.startswith("#"):
            continue
        total += lengths[name]
print(total)
PY
}

focused_record_cell() {
    printf '%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "${4:-}" >> "$FOCUSED_MATRIX"
}

focused_check_output() {
    local workload="$1" expected_bases="$2"
    shift 2
    local reference="" command output actual_bases
    local i=0
    for command in "$@"; do
        output="$BENCH_VIEW_DIR/focused_output_${i}.fa"
        bash -c "$command" > "$output" || {
            echo "error: focused output command failed: $workload" >&2
            return 1
        }
        if [[ -z "$reference" ]]; then
            reference="$output"
            [[ -s "$reference" ]] || {
                echo "error: focused output is empty: $workload" >&2
                return 1
            }
        elif ! cmp -s "$reference" "$output"; then
            echo "error: focused output differs: $workload" >&2
            return 1
        fi
        i=$((i + 1))
    done

    actual_bases="$(awk '!/^>/{sub(/\r$/, ""); bases += length($0)} END{print bases + 0}' "$reference")"
    [[ "$actual_bases" -eq "$expected_bases" ]] || {
        echo "error: focused output base count differs: $workload ($actual_bases != $expected_bases)" >&2
        return 1
    }
}

focused_run_group() {
    local workload="$1" file_id="$2" input_bytes="$3" output_bases="$4"
    shift 4
    [[ $# -gt 0 && $(( $# % 3 )) -eq 0 ]] || {
        echo "error: invalid focused command group: $workload" >&2
        return 1
    }
    local json="$FOCUSED_OUT/${file_id}.json"
    local tool family command
    local -a commands=()

    zebrac_clear_commands
    while [[ $# -gt 0 ]]; do
        tool="$1"
        family="$2"
        command="$3"
        shift 3
        commands+=("$command")
        get_add_command focused "$workload" "$tool" "$family" "$json" \
            "$command > /dev/null" "$input_bytes" "$output_bases"
        focused_record_cell "$workload" "$tool" run
    done
    focused_check_output "$workload" "$output_bases" "${commands[@]}"
    bench_group "$json"
    echo "  focused $workload"
}

validate_focused_matrix() {
    python3 - "$FOCUSED_MATRIX" "$METADATA_JSONL" "$FOCUSED_OUT" "$RUNS" <<'PY'
import csv
import json
import sys
from pathlib import Path

matrix_path, metadata_path, output_dir, runs = sys.argv[1:]
output_dir = Path(output_dir).resolve()
runs = int(runs)

with open(matrix_path, newline="", encoding="utf-8") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))

tools = {"z-fasta-default", "z-fasta-fai", "noodles"}
matched = {
    f"{dataset}/{feature}"
    for dataset in ("Genome", "Transcriptome", "Proteome")
    for feature in ("direct", "multi-scattered", "bed", "locality-clustered", "locality-alternating")
}
matched.update({
    "Genome/transform-rc",
    "Genome/bed-strand",
    "Transcriptome/names",
    "Transcriptome/transform-rc",
    "Transcriptome/bed-strand",
    "Proteome/names",
})
zf_only = {
    "Transcriptome/names-stdin",
    "Transcriptome/bed-stdin",
    "Transcriptome/summary",
    "Transcriptome/annotation",
    "Transcriptome/composed-transform",
}
all_skipped = {"Genome/names", "Proteome/transform-rc", "Proteome/bed-strand"}

required = {(workload, tool): "run" for workload in matched for tool in tools}
required.update({
    (workload, tool): "run" if tool != "noodles" else "skip"
    for workload in zf_only
    for tool in tools
})
required.update({(workload, tool): "skip" for workload in all_skipped for tool in tools})
required.update({
    ("Mixed/side-table", "z-fasta-default"): "run",
    ("Mixed/side-table", "z-fasta-fai"): "skip",
    ("Mixed/side-table", "noodles"): "skip",
})

declared = {(row["workload"], row["tool"]): row["state"] for row in rows}
if len(declared) != len(rows):
    raise SystemExit("duplicate focused matrix cell")
if declared != required:
    missing = sorted(required.keys() - declared.keys())
    extra = sorted(declared.keys() - required.keys())
    wrong = sorted(key for key in required.keys() & declared.keys() if required[key] != declared[key])
    raise SystemExit(f"focused declaration mismatch: missing={missing}, extra={extra}, wrong_state={wrong}")

expected = {key for key, state in required.items() if state == "run"}
skipped = {key for key, state in required.items() if state == "skip"}
if any(row["state"] == "skip" and not row["reason"] for row in rows):
    raise SystemExit("focused matrix skip lacks a reason")

metadata = []
with open(metadata_path, encoding="utf-8") as handle:
    for line in handle:
        row = json.loads(line)
        if row.get("section") != "focused":
            raise SystemExit(f"unexpected metadata section: {row.get('section')}")
        metadata.append(row)

actual = {(row["workload"], row["tool"]) for row in metadata}
if len(actual) != len(metadata):
    raise SystemExit("duplicate focused metadata cell")
if actual != expected:
    missing = sorted(expected - actual)
    extra = sorted(actual - expected)
    raise SystemExit(f"focused matrix mismatch: missing={missing}, extra={extra}")

by_json = {}
for row in metadata:
    expected_family = "noodles" if row["tool"] == "noodles" else "z-fasta"
    if row.get("suite") != "get" or row.get("tool_family") != expected_family:
        raise SystemExit(f"invalid suite or tool family: {row['workload']} {row['tool']}")
    path = Path(row["raw_json"]).resolve()
    if output_dir not in path.parents:
        raise SystemExit(f"focused JSON outside output directory: {path}")
    if not isinstance(row.get("output_bases"), int) or row["output_bases"] <= 0:
        raise SystemExit(f"invalid output_bases: {row['workload']} {row['tool']}")
    by_json.setdefault(path, {})[row["command"]] = row

seen = set()
for path, command_rows in by_json.items():
    data = json.loads(path.read_text())
    results = data.get("results", [])
    if len(results) != len(command_rows):
        raise SystemExit(f"command count mismatch: {path}")
    for result in results:
        command = result.get("command")
        if command not in command_rows:
            raise SystemExit(f"unrecorded command in {path}: {command}")
        row = command_rows[command]
        cell = (row["workload"], row["tool"])
        seen.add(cell)
        if result.get("sample_count") != runs or result.get("failed_sample_count") != 0:
            raise SystemExit(f"invalid samples or failures: {cell}")
        for metric in ("wall_time", "peak_rss", "minor_faults", "major_faults"):
            if result.get(metric, {}).get("sample_count") != runs:
                raise SystemExit(f"invalid {metric} samples: {cell}")

if seen != expected:
    raise SystemExit(f"focused result mismatch: missing={sorted(expected - seen)}")
print(f"  focused benchmark: {len(expected)} timed cells, {len(skipped)} explicit skips, zero failures")
PY
}

run_perf_focused() {
    bench_require_tool noodles
    local required
    for required in Genome Transcriptome Proteome; do
        [[ -n "${REAL_DATASETS[$required]:-}" ]] || {
            echo "error: focused GET matrix requires $required" >&2
            return 1
        }
    done
    generate_get_fixtures focused

    FOCUSED_OUT="$RESULTS_DIR/focused_${TIMESTAMP}"
    FOCUSED_MATRIX="$FOCUSED_OUT/matrix.tsv"
    SECTION_FOCUSED="focused_${TIMESTAMP}"
    mkdir -p "$FOCUSED_OUT"
    printf 'workload\ttool\tstate\treason\n' > "$FOCUSED_MATRIX"

    local qz qn
    qz="$(quote_arg "$ZFASTA")"
    qn="$(quote_arg "$NOODLES")"

    local ds fa zfi_fa fai_fa qzfi qfai nbytes spec region qrgn out_bases
    local multi_count multi_path zfi_core fai_core noodles_core bed qbed names names_count
    for ds in Genome Transcriptome Proteome; do
        fa="${REAL_DATASETS[$ds]}"
        zfi_fa="${ZFI_VIEWS[$ds]}"
        fai_fa="${FAI_VIEWS[$ds]}"
        qzfi="$(quote_arg "$zfi_fa")"
        qfai="$(quote_arg "$fai_fa")"
        nbytes="$(file_size_bytes "$fa")"

        spec="$REGIONS_DIR/${ds}_1kbp_mid.txt"
        [[ -f "$spec" ]] || spec="$REGIONS_DIR/${ds}_100bp.txt"
        region="$(tr -d '\n' < "$spec")"
        qrgn="$(quote_arg "$region")"
        out_bases="$(region_bases "$region" "$fai_fa")"
        focused_run_group "$ds/direct" "${ds}_direct" "$nbytes" "$out_bases" \
            z-fasta-default z-fasta "$qz get $qzfi $qrgn" \
            z-fasta-fai z-fasta "$qz get $qfai $qrgn" \
            noodles noodles "$qn get $qfai $qrgn"

        multi_count=100
        [[ "$ds" == "Genome" ]] && multi_count=10
        multi_path="$REGIONS_DIR/${ds}_regions_N${multi_count}.txt"
        zfi_core="$(multi_positional_get_core "$qz" "$qzfi" "$multi_path")"
        fai_core="$(multi_positional_get_core "$qz" "$qfai" "$multi_path")"
        noodles_core="$(multi_positional_get_core "$qn" "$qfai" "$multi_path")"
        focused_run_group "$ds/multi-scattered" "${ds}_multi_scattered" "$nbytes" "$((multi_count * 1000))" \
            z-fasta-default z-fasta "$zfi_core" \
            z-fasta-fai z-fasta "$fai_core" \
            noodles noodles "$noodles_core"

        if [[ "$ds" == "Genome" ]]; then
            focused_record_cell "$ds/names" z-fasta-default skip "whole chromosome output is not a focused development cell"
            focused_record_cell "$ds/names" z-fasta-fai skip "whole chromosome output is not a focused development cell"
            focused_record_cell "$ds/names" noodles skip "whole chromosome output is not a focused development cell"
        else
            names_count=100
            names="$FOCUSED_DIR/${ds}_names_${names_count}.txt"
            out_bases="$(focused_names_bases "$fai_fa" "$names")"
            focused_run_group "$ds/names" "${ds}_names" "$nbytes" "$out_bases" \
                z-fasta-default z-fasta "$qz get $qzfi --names $(quote_arg "$names")" \
                z-fasta-fai z-fasta "$qz get $qfai --names $(quote_arg "$names")" \
                noodles noodles "$qn get $qfai --names $(quote_arg "$names")"
        fi

        bed="$FOCUSED_DIR/${ds}_strand_100.bed"
        qbed="$(quote_arg "$bed")"
        focused_run_group "$ds/bed" "${ds}_bed" "$nbytes" 100000 \
            z-fasta-default z-fasta "$qz get $qzfi --bed $qbed" \
            z-fasta-fai z-fasta "$qz get $qfai --bed $qbed" \
            noodles noodles "$qn get $qfai --bed $qbed"

        for locality in clustered alternating; do
            multi_path="$FOCUSED_DIR/${ds}_locality_${locality}.txt"
            zfi_core="$(multi_positional_get_core "$qz" "$qzfi" "$multi_path")"
            fai_core="$(multi_positional_get_core "$qz" "$qfai" "$multi_path")"
            noodles_core="$(multi_positional_get_core "$qn" "$qfai" "$multi_path")"
            focused_run_group "$ds/locality-$locality" "${ds}_locality_${locality}" "$nbytes" 100000 \
                z-fasta-default z-fasta "$zfi_core" \
                z-fasta-fai z-fasta "$fai_core" \
                noodles noodles "$noodles_core"
        done

        if [[ "$ds" == "Proteome" ]]; then
            for tool in z-fasta-default z-fasta-fai noodles; do
                focused_record_cell "$ds/transform-rc" "$tool" skip "complement transforms reject protein"
                focused_record_cell "$ds/bed-strand" "$tool" skip "strand-aware complement rejects protein"
            done
        else
            multi_path="$REGIONS_DIR/${ds}_regions_N${multi_count}.txt"
            zfi_core="$(multi_positional_get_core "$qz" "$qzfi" "$multi_path" " --rc")"
            fai_core="$(multi_positional_get_core "$qz" "$qfai" "$multi_path" " --rc")"
            noodles_core="$(multi_positional_get_core "$qn" "$qfai" "$multi_path" " --rc")"
            focused_run_group "$ds/transform-rc" "${ds}_transform_rc" "$nbytes" "$((multi_count * 1000))" \
                z-fasta-default z-fasta "$zfi_core" \
                z-fasta-fai z-fasta "$fai_core" \
                noodles noodles "$noodles_core"
            focused_run_group "$ds/bed-strand" "${ds}_bed_strand" "$nbytes" 100000 \
                z-fasta-default z-fasta "$qz get $qzfi --bed $qbed --strand-aware" \
                z-fasta-fai z-fasta "$qz get $qfai --bed $qbed --strand-aware" \
                noodles noodles "$qn get $qfai --bed $qbed --rc --honor-strand"
        fi
    done

    ds=Transcriptome
    zfi_fa="${ZFI_VIEWS[$ds]}"
    fai_fa="${FAI_VIEWS[$ds]}"
    qzfi="$(quote_arg "$zfi_fa")"
    qfai="$(quote_arg "$fai_fa")"
    nbytes="$(file_size_bytes "${REAL_DATASETS[$ds]}")"
    names="$FOCUSED_DIR/${ds}_names_100.txt"
    out_bases="$(focused_names_bases "$fai_fa" "$names")"
    focused_run_group "$ds/names-stdin" "${ds}_names_stdin" "$nbytes" "$out_bases" \
        z-fasta-default z-fasta "cat $(quote_arg "$names") | $qz get $qzfi --names -" \
        z-fasta-fai z-fasta "cat $(quote_arg "$names") | $qz get $qfai --names -"
    focused_record_cell "$ds/names-stdin" noodles skip "Noodles wrapper has no names stdin mode"

    bed="$FOCUSED_DIR/${ds}_strand_100.bed"
    qbed="$(quote_arg "$bed")"
    focused_run_group "$ds/bed-stdin" "${ds}_bed_stdin" "$nbytes" 100000 \
        z-fasta-default z-fasta "cat $qbed | $qz get $qzfi --bed -" \
        z-fasta-fai z-fasta "cat $qbed | $qz get $qfai --bed -"
    focused_record_cell "$ds/bed-stdin" noodles skip "Noodles wrapper has no BED stdin mode"

    spec="$REGIONS_DIR/${ds}_1kbp_mid.txt"
    region="$(tr -d '\n' < "$spec")"
    qrgn="$(quote_arg "$region")"
    out_bases="$(region_bases "$region" "$fai_fa")"
    focused_run_group "$ds/summary" "${ds}_summary" "$nbytes" "$out_bases" \
        z-fasta-default z-fasta "$qz get $qzfi $qrgn --summary 2>/dev/null" \
        z-fasta-fai z-fasta "$qz get $qfai $qrgn --summary 2>/dev/null"
    focused_record_cell "$ds/summary" noodles skip "Summary text and timing are z-fasta-specific"
    focused_run_group "$ds/annotation" "${ds}_annotation" "$nbytes" "$out_bases" \
        z-fasta-default z-fasta "$qz get $qzfi $qrgn --rc --annotate-rc" \
        z-fasta-fai z-fasta "$qz get $qfai $qrgn --rc --annotate-rc"
    focused_record_cell "$ds/annotation" noodles skip "Transform annotations are z-fasta-specific"
    focused_run_group "$ds/composed-transform" "${ds}_composed_transform" "$nbytes" 100000 \
        z-fasta-default z-fasta "$qz get $qzfi --bed $qbed --strand-aware --rc --annotate-rc" \
        z-fasta-fai z-fasta "$qz get $qfai --bed $qbed --strand-aware --rc --annotate-rc"
    focused_record_cell "$ds/composed-transform" noodles skip "Composed transform annotations are z-fasta-specific"

    bench_ensure_messy --fixtures
    local messy_source="$MESSY_TEST_DIR/mixed_widths.fasta"
    local messy_view="$BENCH_VIEW_DIR/zfi/focused_mixed_widths.fasta"
    ln -s "$messy_source" "$messy_view"
    "$ZFASTA" index "$messy_view" > /dev/null
    focused_run_group "Mixed/side-table" "Mixed_side_table" \
        "$(file_size_bytes "$messy_source")" 22 \
        z-fasta-default z-fasta "$qz get $(quote_arg "$messy_view") mixed_widths:3-24"
    focused_record_cell "Mixed/side-table" z-fasta-fai skip "Non-uniform side tables are available only in .zfi"
    focused_record_cell "Mixed/side-table" noodles skip "Non-uniform side tables are available only in .zfi"

    validate_focused_matrix
}

run_perf_pos() {
    local out_dir="$RESULTS_DIR/perf_pos_${TIMESTAMP}"
    mkdir -p "$out_dir"
    SECTION_POS="perf_pos_${TIMESTAMP}"
    generate_get_fixtures pos

    local ds fa zfi_fa fai_fa qzfi qfai qz qs qn qr qh qt
    for ds in Genome Transcriptome Proteome; do
        fa="${REAL_DATASETS[$ds]:-}"
        [[ -n "$fa" ]] || continue
        zfi_fa="${ZFI_VIEWS[$ds]}"
        fai_fa="${FAI_VIEWS[$ds]}"
        qzfi="$(quote_arg "$zfi_fa")"
        qfai="$(quote_arg "$fai_fa")"
        qz="$(quote_arg "$ZFASTA")"
        qs="$(quote_arg "$SAMTOOLS")"
        qn="$(quote_arg "$NOODLES")"
        qr="$(quote_arg "$RUSTBIO")"
        qh="$(quote_arg "$FASTAHACK")"
        qt="$(quote_arg "$SEQTK")"
        nbytes="$(file_size_bytes "$fa")"

        for label in "${POS_LABELS[@]}"; do
            local spec="$REGIONS_DIR/${ds}_${label}.txt"
            [[ -f "$spec" ]] || continue
            local region; region="$(tr -d '\n' < "$spec")"
            local qrgn; qrgn="$(quote_arg "$region")"
            local out_bases; out_bases="$(region_bases "$region" "$fai_fa")"
            local json="$out_dir/${ds}_${label}.json"
            local workload="${ds}/${label}"

            zebrac_clear_commands
            get_add_command perf_pos "$workload" z-fasta-default z-fasta "$json" \
                "$qz get $qzfi $qrgn > /dev/null" "$nbytes" "$out_bases"
            get_add_command perf_pos "$workload" z-fasta-fai z-fasta "$json" \
                "$qz get $qfai $qrgn > /dev/null" "$nbytes" "$out_bases"
            bench_has_tool samtools && get_add_command perf_pos "$workload" samtools samtools "$json" \
                "$qs faidx $qfai $qrgn > /dev/null" "$nbytes" "$out_bases"
            bench_has_tool noodles && get_add_command perf_pos "$workload" noodles noodles "$json" \
                "$qn get $qfai $qrgn > /dev/null" "$nbytes" "$out_bases"
            bench_has_tool rustbio && get_add_command perf_pos "$workload" rustbio-custom-get rustbio "$json" \
                "$qr get $qfai $qrgn > /dev/null" "$nbytes" "$out_bases"
            bench_has_tool seqtk && get_add_command perf_pos "$workload" seqtk-reference seqtk "$json" \
                "printf '%s\\n' $qrgn | $qt subseq $qfai - > /dev/null" "$nbytes" "$out_bases"
            if [[ "$label" == "full_seq" ]] && bench_has_tool fastahack; then
                get_add_command perf_pos "$workload" fastahack fastahack "$json" \
                    "$qh -r $qrgn $qfai > /dev/null" "$nbytes" "$out_bases"
            fi
            bench_group "$json"
            echo "  perf_pos $workload"
        done
    done
}

run_perf_multi() {
    local out_dir="$RESULTS_DIR/perf_multi_${TIMESTAMP}"
    mkdir -p "$out_dir"
    SECTION_MULTI="perf_multi_${TIMESTAMP}"
    generate_get_fixtures multi

    local qz qs qn qr qt qh
    qz="$(quote_arg "$ZFASTA")"
    qs="$(quote_arg "$SAMTOOLS")"
    qn="$(quote_arg "$NOODLES")"
    qr="$(quote_arg "$RUSTBIO")"
    qt="$(quote_arg "$SEQTK")"
    qh="$(quote_arg "$FASTAHACK")"

    local -a multi_counts=("${MULTI_COUNTS[@]}")
    if [[ -n "${GET_MULTI_COUNTS:-}" ]]; then
        read -r -a multi_counts <<< "${GET_MULTI_COUNTS}"
    fi

    local ds fa zfi_fa fai_fa qzfi qfai nbytes n path json workload total_bases
    for ds in Genome Transcriptome Proteome; do
        fa="${REAL_DATASETS[$ds]:-}"
        [[ -n "$fa" ]] || continue
        zfi_fa="${ZFI_VIEWS[$ds]}"
        fai_fa="${FAI_VIEWS[$ds]}"
        qzfi="$(quote_arg "$zfi_fa")"
        qfai="$(quote_arg "$fai_fa")"
        nbytes="$(file_size_bytes "$fa")"

        for n in "${multi_counts[@]}"; do
            path="$REGIONS_DIR/${ds}_regions_N${n}.txt"
            [[ -f "$path" ]] || continue
            json="$out_dir/${ds}_N${n}.json"
            workload="${ds}/N=${n}"
            total_bases=$((n * 1000))

            local zf_cmd zf_fai_cmd nd_cmd rb_cmd st_cmd
            zf_cmd="$(multi_positional_get_cmd "$qz" "$qzfi" "$path")"
            zf_fai_cmd="$(multi_positional_get_cmd "$qz" "$qfai" "$path")"
            nd_cmd="$(multi_positional_get_cmd "$qn" "$qfai" "$path")"
            rb_cmd="$(multi_positional_get_cmd "$qr" "$qfai" "$path")"
            st_cmd="$qs faidx -r $(quote_arg "$path") $qfai > /dev/null"

            zebrac_clear_commands
            get_add_command perf_multi "$workload" z-fasta-default z-fasta "$json" "$zf_cmd" "$nbytes" "$total_bases"
            get_add_command perf_multi "$workload" z-fasta-fai z-fasta "$json" "$zf_fai_cmd" "$nbytes" "$total_bases"
            bench_has_tool samtools && get_add_command perf_multi "$workload" samtools samtools "$json" "$st_cmd" "$nbytes" "$total_bases"
            bench_has_tool noodles && get_add_command perf_multi "$workload" noodles noodles "$json" "$nd_cmd" "$nbytes" "$total_bases"
            bench_has_tool rustbio && get_add_command perf_multi "$workload" rustbio-custom-get rustbio "$json" "$rb_cmd" "$nbytes" "$total_bases"
            if [[ "${GET_MULTI_REFERENCE:-}" == "1" ]]; then
                local seqtk_loop="while IFS= read -r r; do printf '%s\\n' \"\$r\" | $qt subseq $qfai - > /dev/null; done < $(quote_arg "$path")"
                local fh_loop="while IFS= read -r r; do $qh -r \"\$r\" $qfai > /dev/null; done < $(quote_arg "$path")"
                bench_has_tool seqtk && get_add_command perf_multi "$workload" seqtk-reference seqtk "$json" "$seqtk_loop" "$nbytes" "$total_bases"
                bench_has_tool fastahack && get_add_command perf_multi "$workload" fastahack-reference fastahack "$json" "$fh_loop" "$nbytes" "$total_bases"
            fi
            bench_group "$json"
            echo "  perf_multi ${ds} ${workload}"
        done
    done
}

run_perf_bed() {
    local out_dir="$RESULTS_DIR/perf_bed_${TIMESTAMP}"
    mkdir -p "$out_dir"
    SECTION_BED="perf_bed_${TIMESTAMP}"
    generate_get_fixtures bed

    local qz qs qb qn qr
    qz="$(quote_arg "$ZFASTA")"
    qs="$(quote_arg "$SAMTOOLS")"
    qb="$(quote_arg "$BEDTOOLS")"
    qn="$(quote_arg "$NOODLES")"
    qr="$(quote_arg "$RUSTBIO")"

    local ds fa zfi_fa fai_fa qzfi qfai nbytes n bed regions json workload out_bases
    for ds in Genome Transcriptome Proteome; do
        fa="${REAL_DATASETS[$ds]:-}"
        [[ -n "$fa" ]] || continue
        zfi_fa="${ZFI_VIEWS[$ds]}"
        fai_fa="${FAI_VIEWS[$ds]}"
        qzfi="$(quote_arg "$zfi_fa")"
        qfai="$(quote_arg "$fai_fa")"
        nbytes="$(file_size_bytes "$fa")"

        for n in "${BED_COUNTS[@]}"; do
            bed="$BED_DIR/${ds}_${n}.bed"
            regions="$BED_DIR/${ds}_${n}.regions.txt"
            [[ -f "$bed" && -f "$regions" ]] || continue
            json="$out_dir/${ds}_rows_${n}.json"
            workload="${ds}/rows=${n}"
            out_bases=$((n * 1000))
            local qbed qreg
            qbed="$(quote_arg "$bed")"
            qreg="$(quote_arg "$regions")"

            zebrac_clear_commands
            get_add_command perf_bed "$workload" z-fasta-default z-fasta "$json" \
                "$qz get $qzfi --bed $qbed > /dev/null" "$nbytes" "$out_bases"
            get_add_command perf_bed "$workload" z-fasta-fai z-fasta "$json" \
                "$qz get $qfai --bed $qbed > /dev/null" "$nbytes" "$out_bases"
            bench_has_tool bedtools && get_add_command perf_bed "$workload" bedtools bedtools "$json" \
                "$qb getfasta -fi $qfai -bed $qbed > /dev/null" "$nbytes" "$out_bases"
            bench_has_tool samtools && get_add_command perf_bed "$workload" samtools samtools "$json" \
                "$qs faidx -r $qreg $qfai > /dev/null" "$nbytes" "$out_bases"
            bench_has_tool noodles && get_add_command perf_bed "$workload" noodles noodles "$json" \
                "$qn get $qfai --bed $qbed > /dev/null" "$nbytes" "$out_bases"
            bench_has_tool rustbio && get_add_command perf_bed "$workload" rustbio-custom-get rustbio "$json" \
                "$qr get $qfai --bed $qbed > /dev/null" "$nbytes" "$out_bases"
            bench_group "$json"
            echo "  perf_bed ${workload}"

            if [[ "$n" == "10" ]]; then
                json="$out_dir/${ds}_rows_${n}_stdin.json"
                zebrac_clear_commands
                get_add_command perf_bed "${workload}/stdin" z-fasta-default z-fasta "$json" \
                    "cat $qbed | $qz get $qzfi --bed - > /dev/null" "$nbytes" "$out_bases"
                get_add_command perf_bed "${workload}/stdin" z-fasta-fai z-fasta "$json" \
                    "cat $qbed | $qz get $qfai --bed - > /dev/null" "$nbytes" "$out_bases"
                bench_group "$json"
                echo "  perf_bed ${workload}/stdin"
            fi
        done
    done
}

run_perf_rc() {
    local out_dir="$RESULTS_DIR/perf_rc_${TIMESTAMP}"
    mkdir -p "$out_dir"
    SECTION_RC="perf_rc_${TIMESTAMP}"
    generate_get_fixtures pos

    local qz qs qt qn qr
    qz="$(quote_arg "$ZFASTA")"
    qs="$(quote_arg "$SAMTOOLS")"
    qt="$(quote_arg "$SEQTK")"
    qn="$(quote_arg "$NOODLES")"
    qr="$(quote_arg "$RUSTBIO")"

    local ds fa zfi_fa fai_fa region_file region workload json qzfi qfai qrgn nbytes out_bases
    # Nucleotide datasets only: Proteome omitted (--rc / complement N/A; reverse-only is trivial).
    for ds in Genome Transcriptome; do
        fa="${REAL_DATASETS[$ds]:-}"
        [[ -n "$fa" ]] || continue
        region_file="$REGIONS_DIR/${ds}_1kbp_mid.txt"
        [[ -f "$region_file" ]] || continue
        region="$(tr -d '\n' < "$region_file")"
        workload="${ds}/1kbp_mid"
        json="$out_dir/${ds}_1kbp_mid.json"
        zfi_fa="${ZFI_VIEWS[$ds]}"
        fai_fa="${FAI_VIEWS[$ds]}"
        qzfi="$(quote_arg "$zfi_fa")"
        qfai="$(quote_arg "$fai_fa")"
        qrgn="$(quote_arg "$region")"
        nbytes="$(file_size_bytes "$fa")"
        out_bases="$(region_bases "$region" "$fai_fa")"

        zebrac_clear_commands
        get_add_command perf_rc "$workload" z-fasta-default z-fasta "$json" \
            "$qz get $qzfi $qrgn > /dev/null" "$nbytes" "$out_bases"
        get_add_command perf_rc "$workload" z-fasta-rc z-fasta "$json" \
            "$qz get $qzfi $qrgn --rc > /dev/null" "$nbytes" "$out_bases"
        get_add_command perf_rc "$workload" z-fasta-fai-rc z-fasta "$json" \
            "$qz get $qfai $qrgn --rc > /dev/null" "$nbytes" "$out_bases"
        get_add_command perf_rc "$workload" z-fasta-complement-only z-fasta "$json" \
            "$qz get $qzfi $qrgn --complement-only > /dev/null" "$nbytes" "$out_bases"
        get_add_command perf_rc "$workload" z-fasta-reverse-only z-fasta "$json" \
            "$qz get $qzfi $qrgn --reverse-only > /dev/null" "$nbytes" "$out_bases"
        bench_has_tool seqtk && get_add_command perf_rc "$workload" seqtk-reference seqtk "$json" \
            "printf '%s\\n' $qrgn | $qt subseq $qfai - | $qt seq -r > /dev/null" "$nbytes" "$out_bases"
        bench_has_tool noodles && get_add_command perf_rc "$workload" noodles-rc noodles "$json" \
            "$qn get $qfai $qrgn --rc > /dev/null" "$nbytes" "$out_bases"
        bench_has_tool rustbio && get_add_command perf_rc "$workload" rustbio-custom-get-rc rustbio "$json" \
            "$qr get $qfai $qrgn --rc > /dev/null" "$nbytes" "$out_bases"
        bench_group "$json"
        echo "  perf_rc ${workload}"
    done
}

run_perf_messy() {
    local out_dir="$RESULTS_DIR/messy_${TIMESTAMP}"
    mkdir -p "$out_dir"

    bench_ensure_messy --perf

    # Drop stale JSON from prior messy perf formats (e.g. verify-region spans).
    while IFS= read -r stale; do
        rm -f "$out_dir/$stale"
    done < <(python3 - "$out_dir" "$MESSY_PERF_JSON" <<'PY'
import json
import sys
from pathlib import Path

out_dir = Path(sys.argv[1])
cfg = json.load(open(sys.argv[2]))
expected = {
    f"{variant}_{span_id}.json"
    for variant in cfg["variants"]
    for span_id in cfg["spans"]
}
for path in sorted(out_dir.glob("*.json")):
    if path.name not in expected:
        print(path.name)
PY
    )
    SECTION_MESSY="messy_${TIMESTAMP}"
    local work="$RESULTS_DIR/messy_work_${TIMESTAMP}"
    mkdir -p "$work"

    [[ -f "$MESSY_PERF_JSON" ]] || { echo "error: missing $MESSY_PERF_JSON" >&2; exit 1; }

    local saved_min saved_max saved_warm saved_dur
    saved_min="${ZEBRAC_MIN_SAMPLES:-5}"
    saved_max="${ZEBRAC_MAX_SAMPLES:-5}"
    saved_warm="${ZEBRAC_WARMUP:-1}"
    saved_dur="${ZEBRAC_DURATION_MS:-5000}"
    ZEBRAC_MIN_SAMPLES="$MESSY_RUNS"
    ZEBRAC_MAX_SAMPLES="$MESSY_RUNS"
    ZEBRAC_WARMUP="$MESSY_WARMUP"
    ZEBRAC_DURATION_MS="$MESSY_DURATION_MS"
    export GET_MESSY_RUNS="$MESSY_RUNS"
    export GET_MESSY_WARMUP="$MESSY_WARMUP"
    export GET_MESSY_DURATION_MS="$MESSY_DURATION_MS"

    echo "  messy zebrac: runs=$MESSY_RUNS warmup=$MESSY_WARMUP duration=${MESSY_DURATION_MS}ms"
    echo "  messy fixtures: $MESSY_DIR (see $MESSY_PERF_JSON)"

    local qz; qz="$(quote_arg "$ZFASTA")"
    while read -r name; do
        local src="$MESSY_DIR/${name}.fasta"
        [[ -f "$src" ]] || { echo "error: missing messy fixture $src" >&2; exit 1; }
        local messy_fasta="$work/${name}.fasta"
        local uniform_fasta="$work/${name}_uniform.fasta"
        if [[ ! -f "$messy_fasta" ]]; then
            echo "  preparing $name (copy + validate --fix + index)..."
            cp "$src" "$messy_fasta"
            # stdout is the event dump (can be huge); keep stderr so real failures surface.
            if ! "$ZFASTA" validate --fix -o "$uniform_fasta" "$messy_fasta" >/dev/null; then
                echo "error: validate --fix failed for $name" >&2
                exit 1
            fi
            "$ZFASTA" index "$messy_fasta" >/dev/null
            "$ZFASTA" index "$uniform_fasta" >/dev/null
        fi
    done < <(python3 - "$MESSY_PERF_JSON" <<'PY'
import json, sys
for name in json.load(open(sys.argv[1]))["variants"]:
    print(name)
PY
    )

    while IFS=$'\t' read -r name span_id region_full; do
        local messy_fasta="$work/${name}.fasta"
        local uniform_fasta="$work/${name}_uniform.fasta"
        local qf_messy qf_uniform messy_nbytes uniform_nbytes qrgn out_bases json base_workload
        qf_messy="$(quote_arg "$messy_fasta")"
        qf_uniform="$(quote_arg "$uniform_fasta")"
        messy_nbytes="$(file_size_bytes "$messy_fasta")"
        uniform_nbytes="$(file_size_bytes "$uniform_fasta")"
        qrgn="$(quote_arg "$region_full")"
        out_bases="$(region_span_bases "$region_full")"
        json="$out_dir/${name}_${span_id}.json"
        base_workload="${name}/${span_id}"

        zebrac_clear_commands
        get_add_command perf_messy "${base_workload}/uniform" z-fasta-uniform z-fasta "$json" \
            "$qz get $qf_uniform $qrgn > /dev/null" "$uniform_nbytes" "$out_bases"
        get_add_command perf_messy "${base_workload}/messy" z-fasta-messy z-fasta "$json" \
            "$qz get $qf_messy $qrgn > /dev/null" "$messy_nbytes" "$out_bases"
        bench_group "$json"
        echo "  perf_messy ${base_workload}/uniform ${base_workload}/messy"
    done < <(python3 - "$MESSY_PERF_JSON" <<'PY'
import json
import sys

cfg = json.load(open(sys.argv[1]))
for name in cfg["variants"]:
    for span_id, span in cfg["spans"].items():
        print(f"{name}\t{span_id}\t{span['region']}")
PY
    )
    rm -rf "$work"

    ZEBRAC_MIN_SAMPLES="$saved_min"
    ZEBRAC_MAX_SAMPLES="$saved_max"
    ZEBRAC_WARMUP="$saved_warm"
    ZEBRAC_DURATION_MS="$saved_dur"
}

run_perf() {
    bench_require_tool zebrac
    bench_require_tool z-fasta
    METADATA_JSONL="$RESULTS_DIR/metadata_${TIMESTAMP}.jsonl"
    RUN_MANIFEST="$RESULTS_DIR/run_${TIMESTAMP}.json"
    if $FOCUSED && { [[ -f "$METADATA_JSONL" ]] || [[ -d "$RESULTS_DIR/focused_${TIMESTAMP}" ]]; }; then
        echo "error: focused GET result already exists for $TIMESTAMP" >&2
        return 1
    elif [[ -n "${GET_RUN_TIMESTAMP:-}" ]] && [[ -f "$METADATA_JSONL" ]]; then
        :
    else
        : > "$METADATA_JSONL"
    fi
    ZEBRAC_MIN_SAMPLES="$RUNS"
    ZEBRAC_MAX_SAMPLES="$RUNS"
    ZEBRAC_WARMUP="$WARMUP"

    echo ""
    echo "================================================================"
    echo "  GET zebrac performance"
    echo "================================================================"
    echo "  Runs: $RUNS | Warmup: $WARMUP | Duration: ${ZEBRAC_DURATION_MS}ms"
    echo "  metadata: $METADATA_JSONL"
    echo ""

    ensure_real_data
    prepare_benchmark_views

    if $FOCUSED; then
        run_perf_focused
    else
        $DO_POS && run_perf_pos
        $DO_MULTI && run_perf_multi
        $DO_BED && run_perf_bed
        $DO_RC && run_perf_rc
        $DO_MESSY && run_perf_messy

        : "${SECTION_POS:=$(existing_section_dir perf_pos)}"
        : "${SECTION_MULTI:=$(existing_section_dir perf_multi)}"
        : "${SECTION_BED:=$(existing_section_dir perf_bed)}"
        : "${SECTION_RC:=$(existing_section_dir perf_rc)}"
        : "${SECTION_MESSY:=$(existing_section_dir messy)}"
    fi

    export_manifest_tool_versions
    write_run_manifest "$RUN_MANIFEST" "$TIMESTAMP" "$(basename "$METADATA_JSONL")"
    if $FOCUSED; then
        printf '%s\n' "$TIMESTAMP" > "$RESULTS_DIR/FOCUSED_LATEST"
    else
        printf '%s\n' "$TIMESTAMP" > "$RESULTS_DIR/LATEST"
    fi
    echo "  manifest: $RUN_MANIFEST"
}

run_report() {
    local py; py="$(report_python)"
    echo ""
    echo "================================================================"
    echo "  Generating REPORT.md"
    echo "================================================================"
    "$py" "$SCRIPT_DIR/generate_report.py"
}

TIMESTAMP="${GET_RUN_TIMESTAMP:-$(date +%Y%m%d_%H%M%S)}"
mkdir -p "$RESULTS_DIR"
trap cleanup_benchmark_views EXIT

if $DO_TESTS; then
    echo ""
    echo "================================================================"
    echo "  run_tests (correctness)"
    echo "================================================================"
    if run_tests; then
        export BENCH_VERIFY_PASS="$PASS"
    else
        export BENCH_VERIFY_PASS=0
        echo "error: run_tests failed" >&2
        exit 1
    fi
fi
$DO_BENCHMARKS && run_perf
$DO_REPORT && run_report

echo ""
echo "================================================================"
echo "  GET suite complete"
echo "================================================================"
