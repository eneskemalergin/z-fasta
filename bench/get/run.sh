#!/usr/bin/env bash
# GET benchmark runner: verify (L2) then zebrac perf (L3) and report.
#
# Usage:
#   bash bench/get/run.sh [options]
#
# Defaults: verify first (409 checks), then perf (--runs 5, --warmup 1, --duration 5000).
# A full perf pass (all sections, default zebrac settings) typically takes several hours.
# pyfaidx is omitted from timed positional runs (typically 30-100x slower than z-fasta;
# may be enabled for a final polish pass).
# Multi-region seqtk/fastahack reference loops are omitted by default (O(N) subprocesses;
# cost grows linearly and duplicates positional seqtk reference data). Set
# GET_MULTI_REFERENCE=1 to re-enable them for a one-off reference sweep.
# Use --skip-* for partial re-runs; resume with GET_RUN_TIMESTAMP=<ts> and GET_MULTI_COUNTS.
# Use --skip-verify to run perf/report only. Verify failure aborts before perf or report.
#
#   bash bench/get/run.sh
#   bash bench/get/run.sh --skip-verify --skip-report   # perf only
#   bash bench/get/run.sh --skip-messy --skip-bed
#   bash bench/get/run.sh --regenerate-fixtures
#   bash bench/get/run.sh --allow-incomplete
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
MESSY_DIR="$BENCH_ROOT/shared/messy_perf"
MESSY_PERF_JSON="$SCRIPT_DIR/messy_perf.json"
FIXTURE_DIR="$SCRIPT_DIR/data"
REGIONS_DIR="$FIXTURE_DIR/regions"
BED_DIR="$FIXTURE_DIR/bed"

source "$BENCH_ROOT/shared/zebrac_runner.sh"

RUNS=5
WARMUP=1
ZEBRAC_DURATION_MS="${ZEBRAC_DURATION_MS:-5000}"
# Messy GET uses sub-millisecond timings on small slices; large proteome fixtures use
# higher sampling. Override via GET_MESSY_* or --messy-runs etc.
MESSY_RUNS="${GET_MESSY_RUNS:-50}"
MESSY_WARMUP="${GET_MESSY_WARMUP:-10}"
MESSY_DURATION_MS="${GET_MESSY_DURATION_MS:-30000}"
DO_VERIFY=true
DO_PERF=true
DO_POS=true
DO_MULTI=true
DO_BED=true
DO_RC=true
DO_MESSY=true
DO_REPORT=true
REGENERATE_FIXTURES=false
ALLOW_INCOMPLETE=false

# Multi-region perf sweep (log-spaced). 15/16 boundary lives in verify.sh, not timed here.
MULTI_COUNTS=(1 10 100 1000)
BED_COUNTS=(10 100 1000 10000)
POS_LABELS=(100bp 1kbp_mid 10kbp full_seq)

# Mirror getter.parseRegion: scan colons from the right; suffix is START-END or START-.
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
while colon_pos >= 0:
    parsed = parse_suffix(region[colon_pos + 1 :])
    if parsed is not None:
        start, end = parsed
        if end is None:
            length = fai_length(region[:colon_pos])
            print(length - start + 1 if length is not None else 1)
        else:
            print(end - start + 1)
        raise SystemExit(0)
    colon_pos = region.rfind(":", 0, colon_pos)

length = fai_length(region)
if length is not None:
    print(length)
    raise SystemExit(0)
raise SystemExit(f"region_span_bases: {region!r} not found in index")
PY
}

region_out_bases() {
    region_span_bases "$1"
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
        --skip-verify) DO_VERIFY=false; shift ;;
        --skip-perf) DO_PERF=false; shift ;;
        --skip-pos) DO_POS=false; shift ;;
        --skip-multi) DO_MULTI=false; shift ;;
        --skip-bed) DO_BED=false; shift ;;
        --skip-rc) DO_RC=false; shift ;;
        --skip-messy) DO_MESSY=false; shift ;;
        --skip-report) DO_REPORT=false; shift ;;
        --regenerate-fixtures) REGENERATE_FIXTURES=true; shift ;;
        --allow-incomplete) ALLOW_INCOMPLETE=true; shift ;;
        -h|--help)
            sed -n '2,/^set -euo pipefail$/p' "$0" | head -n -1
            exit 0
            ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

report_python() {
    if [[ -x "$PROJECT_ROOT/.venv/bin/python" ]]; then
        echo "$PROJECT_ROOT/.venv/bin/python"
    else
        echo python3
    fi
}

get_add_command() {
    local section="$1" workload="$2" tool="$3" family="$4"
    local json_out="$5" script="$6" input_bytes="$7" output_bases="${8:-}"
    zebrac_add_command "get" "$section" "$workload" "$tool" "$family" \
        "$input_bytes" "$output_bases" "$json_out" "$(shell_command "$script")"
}

ensure_real_data() {
    declare -gA REAL_DATASETS=()
    [[ -f "$DATA_DIR/REAL_Genome.fa" ]] && REAL_DATASETS["Genome"]="$DATA_DIR/REAL_Genome.fa"
    [[ -f "$DATA_DIR/REAL_Transcriptome.fa" ]] && REAL_DATASETS["Transcriptome"]="$DATA_DIR/REAL_Transcriptome.fa"
    [[ -f "$DATA_DIR/REAL_Proteome.fasta" ]] && REAL_DATASETS["Proteome"]="$DATA_DIR/REAL_Proteome.fasta"
    if [[ ${#REAL_DATASETS[@]} -eq 0 ]]; then
        echo "  Fetching REAL_* datasets..."
        bash "$BENCH_ROOT/shared/download_data.sh"
        [[ -f "$DATA_DIR/REAL_Genome.fa" ]] && REAL_DATASETS["Genome"]="$DATA_DIR/REAL_Genome.fa"
        [[ -f "$DATA_DIR/REAL_Transcriptome.fa" ]] && REAL_DATASETS["Transcriptome"]="$DATA_DIR/REAL_Transcriptome.fa"
        [[ -f "$DATA_DIR/REAL_Proteome.fasta" ]] && REAL_DATASETS["Proteome"]="$DATA_DIR/REAL_Proteome.fasta"
    fi
    if [[ ${#REAL_DATASETS[@]} -eq 0 ]]; then
        echo "error: no REAL_* datasets under $DATA_DIR" >&2
        exit 1
    fi
}

preload_real_indexes() {
    bench_require_tool z-fasta
    echo "  Preloading REAL_* indexes (.zfi + .fai)..."
    local fa
    for fa in "${REAL_DATASETS[@]}"; do
        # Rebuild when FASTA is newer, sidecar is missing, z-fasta binary is newer than .zfi,
        # or .zfi predates the embedded name-table footer (ZFNM).
        local need_zfi=false
        if [[ ! -f "${fa}.zfi" ]] || [[ "$fa" -nt "${fa}.zfi" ]] || [[ "$ZFASTA" -nt "${fa}.zfi" ]]; then
            need_zfi=true
        elif ! tail -c 12 "${fa}.zfi" | grep -q 'ZFNM'; then
            need_zfi=true
        fi
        if $need_zfi; then
            "$ZFASTA" index "$fa" > /dev/null
        fi
        if bench_has_tool samtools; then
            if [[ ! -f "${fa}.fai" ]] || [[ "$fa" -nt "${fa}.fai" ]]; then
                samtools faidx "$fa" > /dev/null 2>&1 || true
            fi
        fi
    done
}

bed_to_regions() {
    awk -F'\t' '!/^[[:space:]]*$|^#|^track|^browser/{sub(/\r$/,""); printf "%s:%d-%d\n",$1,$2+1,$3}' "$1"
}

read_region_file() {
    local path="$1"
    local -n _out=$2
    _out=()
    while IFS= read -r line || [[ -n "$line" ]]; do
        [[ -z "$line" || "$line" == \#* ]] && continue
        _out+=("$line")
    done < "$path"
}

bench_group() {
    local json_out="$1"
    zebrac_run_current_group "$json_out" "$METADATA_JSONL"
    zebrac_clear_commands
}

# Build one timed multi-region GET command without expanding regions on the zebrac argv
# (avoids ARG_MAX when N is large). Regions are read inside bash at sample time.
multi_positional_get_cmd() {
    local tool_q="$1"
    local fa_q="$2"
    local regions_path="$3"
    local qr
    qr="$(quote_arg "$regions_path")"
    echo "mapfile -t _regs < ${qr}; ${tool_q} get ${fa_q} \"\${_regs[@]}\" > /dev/null"
}

existing_section_dir() {
    local prefix="$1"
    local dir="$RESULTS_DIR/${prefix}_${TIMESTAMP}"
    if [[ -d "$dir" ]] && compgen -G "${dir}/*.json" > /dev/null; then
        echo "${prefix}_${TIMESTAMP}"
    fi
}

export_manifest_tool_versions() {
    export BENCH_VER_ZEBRAC="$(bench_tool_version zebrac)"
    export BENCH_VER_ZFASTA="$(bench_tool_version z-fasta 2>/dev/null || echo unknown)"
    export BENCH_VER_SAMTOOLS="$(bench_tool_version samtools 2>/dev/null || echo unknown)"
    export BENCH_VER_BEDTOOLS="$(bench_tool_version bedtools 2>/dev/null || echo unknown)"
    export BENCH_VER_SEQTK="$(bench_tool_version seqtk 2>/dev/null | head -1 || echo unknown)"
    for tool in fastahack pyfaidx noodles rustbio; do
        if bench_has_tool "$tool"; then
            upper="${tool^^}"
            export "BENCH_VER_${upper}=$(bench_tool_version "$tool")"
        fi
    done
}

write_run_manifest() {
    local manifest="$1" timestamp="$2" metadata="$3"
    python3 - "$manifest" "$timestamp" "$metadata" "$RUNS" "$WARMUP" "$ZEBRAC_DURATION_MS" \
        "$DO_VERIFY" "$DO_POS" "$DO_MULTI" "$DO_BED" "$DO_RC" "$DO_MESSY" \
        "${SECTION_POS:-}" "${SECTION_MULTI:-}" "${SECTION_BED:-}" "${SECTION_RC:-}" "${SECTION_MESSY:-}" <<'PY'
import json, os, sys
from pathlib import Path

manifest, ts, metadata, runs, warmup, duration = sys.argv[1:7]
do_verify = sys.argv[7] == "true"
section_vals = sys.argv[13:18]

sections = {}
keys = ("perf_pos", "perf_multi", "perf_bed", "perf_rc", "messy")
for key, val in zip(keys, section_vals):
    if val:
        sections[key] = val

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
    "zebrac": os.environ.get("BENCH_VER_ZEBRAC", ""),
    "z_fasta": os.environ.get("BENCH_VER_ZFASTA", ""),
    "runs": int(runs),
    "warmup": int(warmup),
    "duration_ms": int(duration),
    "verify_skipped": not do_verify,
    "verify_pass": os.environ.get("BENCH_VERIFY_PASS"),
    "index_preload": True,
    "skip_pos": skip_pos,
    "skip_multi": skip_multi,
    "skip_bed": skip_bed,
    "skip_rc": skip_rc,
    "skip_messy": skip_messy,
    "skip_names": True,
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

ensure_fai() {
    local fa="$1"
    [[ -f "${fa}.fai" ]] || samtools faidx "$fa"
}

generate_get_fixtures() {
    local mode="${1:-all}"
    mkdir -p "$REGIONS_DIR" "$BED_DIR"
    local fa
    for fa in "${REAL_DATASETS[@]}"; do
        ensure_fai "$fa"
    done
    python3 - "$REGIONS_DIR" "$BED_DIR" "$mode" "$REGENERATE_FIXTURES" \
        "${REAL_DATASETS[Genome]:-}" \
        "${REAL_DATASETS[Transcriptome]:-}" \
        "${REAL_DATASETS[Proteome]:-}" \
        "${BED_COUNTS[*]}" \
        "${MULTI_COUNTS[*]}" <<'PY'
import hashlib
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

MULTI_COUNTS = tuple(multi_counts) if multi_counts else (1, 10, 100, 1000)
POS_LABELS = (("100bp", 100), ("1kbp_mid", 1000), ("10kbp", 10000))
SPAN = 1000
MULTI_SEED = 42
FULL_SEQ_GENOME_MAX = 1000
FULL_SEQ_PROTEOME_MAX = 5000
FULL_SEQ_TRANSCRIPTOME_MAX = 20000


def stable_u32(text: str) -> int:
    """Process/platform-stable u32 from text. Do not use Python hash() (salted per process)."""
    return int.from_bytes(hashlib.sha256(text.encode("utf-8")).digest()[:4], "little")


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
    rng = random.Random(MULTI_SEED + stable_u32(dataset) % 10000)
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
            bed_file(dataset, fasta, bed_dir / f"{dataset}_{n}.bed", n, n + stable_u32(dataset) % 1000)
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

run_verify() {
    echo ""
    echo "================================================================"
    echo "  verify.sh"
    echo "================================================================"
    if bash "$SCRIPT_DIR/verify.sh"; then
        export BENCH_VERIFY_PASS=409
    else
        export BENCH_VERIFY_PASS=0
        echo "error: verify.sh failed" >&2
        exit 1
    fi
}

run_perf_pos() {
    local out_dir="$RESULTS_DIR/perf_pos_${TIMESTAMP}"
    mkdir -p "$out_dir"
    SECTION_POS="perf_pos_${TIMESTAMP}"
    generate_get_fixtures pos

    local ds fa qf qz qs qn qr qh qt
    for ds in Genome Transcriptome Proteome; do
        fa="${REAL_DATASETS[$ds]:-}"
        [[ -n "$fa" ]] || continue
        qf="$(quote_arg "$fa")"
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
            local out_bases; out_bases="$(region_bases "$region" "$fa")"
            local json="$out_dir/${ds}_${label}.json"
            local json_fai="$out_dir/${ds}_${label}_zfai.json"
            local workload="${ds}/${label}"

            zebrac_clear_commands
            get_add_command perf_pos "$workload" z-fasta-default z-fasta "$json" \
                "$qz get $qf $qrgn > /dev/null" "$nbytes" "$out_bases"
            bench_has_tool samtools && get_add_command perf_pos "$workload" samtools samtools "$json" \
                "$qs faidx $qf $qrgn > /dev/null" "$nbytes" "$out_bases"
            bench_has_tool noodles && get_add_command perf_pos "$workload" noodles noodles "$json" \
                "$qn get $qf $qrgn > /dev/null" "$nbytes" "$out_bases"
            bench_has_tool rustbio && get_add_command perf_pos "$workload" rustbio-custom-get rustbio "$json" \
                "$qr get $qf $qrgn > /dev/null" "$nbytes" "$out_bases"
            bench_has_tool seqtk && get_add_command perf_pos "$workload" seqtk-reference seqtk "$json" \
                "printf '%s\\n' $qrgn | $qt subseq $qf - > /dev/null" "$nbytes" "$out_bases"
            if [[ "$label" == "full_seq" ]] && bench_has_tool fastahack; then
                get_add_command perf_pos "$workload" fastahack fastahack "$json" \
                    "$qh -r $qrgn $qf > /dev/null" "$nbytes" "$out_bases"
            fi
            bench_group "$json"

            # .fai lane: stash .zfi outside zebrac so timed get does not include mv cost.
            mv -f "${fa}.zfi" "${fa}.zfi.stash" 2>/dev/null || true
            zebrac_clear_commands
            get_add_command perf_pos "$workload" z-fasta-fai z-fasta "$json_fai" \
                "$qz get $qf $qrgn > /dev/null" "$nbytes" "$out_bases"
            bench_group "$json_fai"
            mv -f "${fa}.zfi.stash" "${fa}.zfi" 2>/dev/null || true
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

    local ds fa qf nbytes n path json workload total_bases
    for ds in Genome Transcriptome Proteome; do
        fa="${REAL_DATASETS[$ds]:-}"
        [[ -n "$fa" ]] || continue
        qf="$(quote_arg "$fa")"
        nbytes="$(file_size_bytes "$fa")"

        for n in "${multi_counts[@]}"; do
            path="$REGIONS_DIR/${ds}_regions_N${n}.txt"
            [[ -f "$path" ]] || continue
            json="$out_dir/${ds}_N${n}.json"
            workload="${ds}/N=${n}"
            total_bases=$((n * 1000))

            local zf_cmd nd_cmd rb_cmd st_cmd
            zf_cmd="$(multi_positional_get_cmd "$qz" "$qf" "$path")"
            nd_cmd="$(multi_positional_get_cmd "$qn" "$qf" "$path")"
            rb_cmd="$(multi_positional_get_cmd "$qr" "$qf" "$path")"
            st_cmd="$qs faidx -r $(quote_arg "$path") $qf > /dev/null"

            zebrac_clear_commands
            get_add_command perf_multi "$workload" z-fasta-default z-fasta "$json" "$zf_cmd" "$nbytes" "$total_bases"
            bench_has_tool samtools && get_add_command perf_multi "$workload" samtools samtools "$json" "$st_cmd" "$nbytes" "$total_bases"
            bench_has_tool noodles && get_add_command perf_multi "$workload" noodles noodles "$json" "$nd_cmd" "$nbytes" "$total_bases"
            bench_has_tool rustbio && get_add_command perf_multi "$workload" rustbio-custom-get rustbio "$json" "$rb_cmd" "$nbytes" "$total_bases"
            if [[ "${GET_MULTI_REFERENCE:-}" == "1" ]]; then
                local seqtk_loop="while IFS= read -r r; do printf '%s\\n' \"\$r\" | $qt subseq $qf - > /dev/null; done < $(quote_arg "$path")"
                local fh_loop="while IFS= read -r r; do $qh -r \"\$r\" $qf > /dev/null; done < $(quote_arg "$path")"
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

    local ds fa qf nbytes n bed regions json workload out_bases
    for ds in Genome Transcriptome Proteome; do
        fa="${REAL_DATASETS[$ds]:-}"
        [[ -n "$fa" ]] || continue
        qf="$(quote_arg "$fa")"
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
                "$qz get $qf --bed $qbed > /dev/null" "$nbytes" "$out_bases"
            if (( n <= 10000 )); then
                get_add_command perf_bed "$workload" z-fasta-chunk-all z-fasta "$json" \
                    "$qz get $qf --bed $qbed --chunk-size -1 > /dev/null" "$nbytes" "$out_bases"
            fi
            if (( n <= 100 )); then
                get_add_command perf_bed "$workload" z-fasta-chunk-1 z-fasta "$json" \
                    "$qz get $qf --bed $qbed --chunk-size 1 > /dev/null" "$nbytes" "$out_bases"
            fi
            bench_has_tool bedtools && get_add_command perf_bed "$workload" bedtools bedtools "$json" \
                "$qb getfasta -fi $qf -bed $qbed > /dev/null" "$nbytes" "$out_bases"
            bench_has_tool samtools && get_add_command perf_bed "$workload" samtools samtools "$json" \
                "$qs faidx -r $qreg $qf > /dev/null" "$nbytes" "$out_bases"
            bench_has_tool noodles && get_add_command perf_bed "$workload" noodles noodles "$json" \
                "$qn get $qf --bed $qbed > /dev/null" "$nbytes" "$out_bases"
            bench_has_tool rustbio && get_add_command perf_bed "$workload" rustbio-custom-get rustbio "$json" \
                "$qr get $qf --bed $qbed > /dev/null" "$nbytes" "$out_bases"
            bench_group "$json"
            echo "  perf_bed ${workload}"

            if [[ "$n" == "10" ]]; then
                json="$out_dir/${ds}_rows_${n}_stdin.json"
                zebrac_clear_commands
                get_add_command perf_bed "${workload}/stdin" z-fasta-default z-fasta "$json" \
                    "cat $qbed | $qz get $qf --bed - > /dev/null" "$nbytes" "$out_bases"
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

    local ds fa region_file region workload json qf qrgn nbytes out_bases
    # Nucleotide datasets only: Proteome omitted (--rc / complement N/A; reverse-only is trivial).
    for ds in Genome Transcriptome; do
        fa="${REAL_DATASETS[$ds]:-}"
        [[ -n "$fa" ]] || continue
        region_file="$REGIONS_DIR/${ds}_1kbp_mid.txt"
        [[ -f "$region_file" ]] || continue
        region="$(tr -d '\n' < "$region_file")"
        workload="${ds}/1kbp_mid"
        json="$out_dir/${ds}_1kbp_mid.json"
        qf="$(quote_arg "$fa")"
        qrgn="$(quote_arg "$region")"
        nbytes="$(file_size_bytes "$fa")"
        out_bases="$(region_bases "$region" "$fa")"

        zebrac_clear_commands
        get_add_command perf_rc "$workload" z-fasta-default z-fasta "$json" \
            "$qz get $qf $qrgn > /dev/null" "$nbytes" "$out_bases"
        get_add_command perf_rc "$workload" z-fasta-rc z-fasta "$json" \
            "$qz get $qf $qrgn --rc > /dev/null" "$nbytes" "$out_bases"
        get_add_command perf_rc "$workload" z-fasta-complement-only z-fasta "$json" \
            "$qz get $qf $qrgn --complement-only > /dev/null" "$nbytes" "$out_bases"
        get_add_command perf_rc "$workload" z-fasta-reverse-only z-fasta "$json" \
            "$qz get $qf $qrgn --reverse-only > /dev/null" "$nbytes" "$out_bases"
        bench_has_tool seqtk && get_add_command perf_rc "$workload" seqtk-reference seqtk "$json" \
            "printf '%s\\n' $qrgn | $qt subseq $qf - | $qt seq -r > /dev/null" "$nbytes" "$out_bases"
        bench_has_tool noodles && get_add_command perf_rc "$workload" noodles-rc noodles "$json" \
            "$qn get $qf $qrgn --rc > /dev/null" "$nbytes" "$out_bases"
        bench_has_tool rustbio && get_add_command perf_rc "$workload" rustbio-custom-get-rc rustbio "$json" \
            "$qr get $qf $qrgn --rc > /dev/null" "$nbytes" "$out_bases"
        bench_group "$json"
        echo "  perf_rc ${workload}"
    done
}

run_perf_messy() {
    local out_dir="$RESULTS_DIR/messy_${TIMESTAMP}"
    mkdir -p "$out_dir"
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
            "$ZFASTA" validate --fix -o "$uniform_fasta" "$messy_fasta" > /dev/null 2>&1
            "$ZFASTA" index "$messy_fasta" > /dev/null
            "$ZFASTA" index "$uniform_fasta" > /dev/null
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
        out_bases="$(region_out_bases "$region_full")"
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
    if [[ -n "${GET_RUN_TIMESTAMP:-}" ]] && [[ -f "$METADATA_JSONL" ]]; then
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
    preload_real_indexes

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

    export_manifest_tool_versions
    write_run_manifest "$RUN_MANIFEST" "$TIMESTAMP" "$(basename "$METADATA_JSONL")"
    printf '%s\n' "$TIMESTAMP" > "$RESULTS_DIR/LATEST"
    echo "  manifest: $RUN_MANIFEST"
}

run_report() {
    local py; py="$(report_python)"
    local args=()
    $ALLOW_INCOMPLETE && args+=(--allow-incomplete)
    echo ""
    echo "================================================================"
    echo "  Generating REPORT.md"
    echo "================================================================"
    "$py" "$SCRIPT_DIR/generate_report.py" "${args[@]}"
}

TIMESTAMP="${GET_RUN_TIMESTAMP:-$(date +%Y%m%d_%H%M%S)}"
mkdir -p "$RESULTS_DIR"

$DO_VERIFY && run_verify
$DO_PERF && run_perf
$DO_REPORT && run_report

echo ""
echo "================================================================"
echo "  GET suite complete"
echo "================================================================"
