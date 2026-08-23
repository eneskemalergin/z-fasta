#!/usr/bin/env bash
# Index benchmark master runner: correctness, performance, messy FASTA, and report.
#
# Usage:
#   bash bench/index/run.sh [options]
#
# Default: run every section, then write REPORT.md.
#
# Sampling:
#   --runs 5              measured samples for real and scaling benchmarks
#   --warmup 5            unmeasured warmups before sampling
#   --duration 5000       sampling budget per command, capped by --runs
#
# Correctness does one exit and byte comparison per tool and fixture. Messy performance
# uses 3 samples, 1 warmup, a 1000 ms budget, and --allow-failures.
#
# Common invocations:
#   bash bench/index/run.sh
#   bash bench/index/run.sh --runs 10 --warmup 2
#   bash bench/index/run.sh --skip-report
#   bash bench/index/run.sh --skip-tests --skip-messy
#   bash bench/index/run.sh --skip-benchmarks --skip-messy --skip-report
#   bash bench/index/run.sh --skip-real
#   bash bench/index/run.sh --skip-scaling
#   bash bench/index/run.sh --skip-size
#   bash bench/index/run.sh --scaling-only
#   bash bench/index/run.sh --regenerate-fixtures
#   bash bench/index/run.sh --clean-legacy
#
# Skip flags:
#   --skip-tests --skip-benchmarks --skip-messy --skip-report
#   --skip-verify and --skip-perf are deprecated aliases.
#
# Outputs:
#   results/run_<ts>.json          result-file index and run settings
#   results/metadata_<ts>.jsonl    command and workload metadata
#   results/LATEST                 newest generated run
#   REPORT.md                      benchmark report
#
# Scaling fixtures under bench/shared/cache/scaling/:
#   size_{N}mb.fasta         1 to 1000 MiB, 100 sequences each
#   seqs_budget_{N}.fasta    about 50 MiB total, N in 1k 10k 100k 250k
#   seqs_fixed_{N}.fasta     1024 bp per sequence, N in 100k 250k 500k 1M
#
# Headline and scaling comparisons use z-fasta .fai. Transcriptome and Proteome also use
# all four format and dedup product lanes.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH_ROOT="$(dirname "$SCRIPT_DIR")"
PROJECT_ROOT="$(dirname "$BENCH_ROOT")"
RESULTS_DIR="$SCRIPT_DIR/results"
SCALING_DIR="$BENCH_ROOT/shared/cache/scaling"
DATA_DIR="$BENCH_ROOT/shared/data"
EDGE_DIR="$SCRIPT_DIR/edge_cases"
MESSY_TEST_DIR="$BENCH_ROOT/shared/cache/messy_fixtures"
MESSY_ZEBRAC_DIR="$BENCH_ROOT/shared/cache/messy_perf"

source "$BENCH_ROOT/shared/runner_common.sh"

SIZE_MBS=(1 5 10 25 50 100 250 500 1000)
SEQ_BUDGET_COUNTS=(1000 10000 100000 250000)
SEQ_FIXED_COUNTS=(100000 250000 500000 1000000)
PRODUCT_DATASETS=(Transcriptome Proteome)

RUNS=5
WARMUP=5
ZEBRAC_DURATION_MS="${ZEBRAC_DURATION_MS:-5000}"
DO_TESTS=true
DO_BENCHMARKS=true
DO_MESSY=true
DO_REPORT=true
SKIP_REAL=false
SKIP_SCALING=false
SKIP_SIZE=false
SCALING_ONLY=false
REGENERATE_FIXTURES=false
CLEAN_LEGACY=false
MERGE_BASE=""
INDEX_CORRECTNESS_STATUS="not_run"
INDEX_CORRECTNESS_CHECKS=0
INDEX_CORRECTNESS_REVIEWS=0
RUN_MANIFEST=""
METADATA_JSONL=""
SPOOL_DIR=""
MESSY_DIR=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --runs)              RUNS="$2"; shift 2 ;;
        --warmup)            WARMUP="$2"; shift 2 ;;
        --duration)          ZEBRAC_DURATION_MS="$2"; shift 2 ;;
        --skip-tests|--skip-verify) DO_TESTS=false; shift ;;
        --skip-benchmarks|--skip-perf) DO_BENCHMARKS=false; shift ;;
        --skip-messy)        DO_MESSY=false; shift ;;
        --skip-report)       DO_REPORT=false; shift ;;
        --skip-real)         SKIP_REAL=true; shift ;;
        --skip-scaling)      SKIP_SCALING=true; shift ;;
        --skip-size)         SKIP_SIZE=true; shift ;;
        --scaling-only)
            SCALING_ONLY=true
            DO_TESTS=false
            DO_MESSY=false
            SKIP_REAL=true
            SKIP_SIZE=true
            REGENERATE_FIXTURES=true
            shift
            ;;
        --merge-base)        MERGE_BASE="$2"; shift 2 ;;
        --regenerate-fixtures) REGENERATE_FIXTURES=true; shift ;;
        --clean-legacy)      CLEAN_LEGACY=true; REGENERATE_FIXTURES=true; shift ;;
        -h|--help)
            awk '
                NR == 1 { next }
                /^set -euo pipefail$/ { exit }
                { sub(/^# ?/, ""); print }
            ' "$0"
            exit 0
            ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

if $SCALING_ONLY; then
    DO_BENCHMARKS=true
    if [[ -z "$MERGE_BASE" && -f "$RESULTS_DIR/LATEST" ]]; then
        MERGE_BASE="$(tr -d '[:space:]' < "$RESULTS_DIR/LATEST")"
    fi
    if [[ -z "$MERGE_BASE" ]]; then
        echo "error: --scaling-only needs --merge-base or results/LATEST from a full run" >&2
        exit 1
    fi
    if [[ ! -f "$RESULTS_DIR/run_${MERGE_BASE}.json" ]]; then
        echo "error: base manifest not found: run_${MERGE_BASE}.json" >&2
        exit 1
    fi
fi

# Helpers: scaling setup, manifest, and zebrac benchmark.

ensure_scaling() {
    local args=(--mode all)
    if $CLEAN_LEGACY; then
        args+=(--clean-legacy --force)
    elif $REGENERATE_FIXTURES; then
        args+=(--force)
    fi
    bench_ensure_scaling "${args[@]}"
}

prepare_benchmark_binary() {
    echo "Building ReleaseFast benchmark subject..."
    (cd "$PROJECT_ROOT" && zig build -Doptimize=ReleaseFast)
}

require_real_datasets() {
    local missing=false file
    for file in \
        "$DATA_DIR/REAL_Genome.fa" \
        "$DATA_DIR/REAL_Transcriptome.fa" \
        "$DATA_DIR/REAL_Proteome.fasta"; do
        [[ -f "$file" ]] || missing=true
    done
    if $missing; then
        echo "Fetching missing real datasets..."
        bash "$BENCH_ROOT/shared/download_data.sh"
    fi
    for file in \
        "$DATA_DIR/REAL_Genome.fa" \
        "$DATA_DIR/REAL_Transcriptome.fa" \
        "$DATA_DIR/REAL_Proteome.fasta"; do
        if [[ ! -f "$file" ]]; then
            echo "error: required real dataset is missing: $file" >&2
            return 1
        fi
    done
}

write_run_manifest() {
    local manifest="$1" timestamp="$2" metadata="$3"
    local zig_version zig_target optimize
    zig_version="$(zig version)"
    zig_target="$(zig env | sed -n 's/^[[:space:]]*\.target = "\([^"]*\)",/\1/p')"
    optimize="ReleaseFast"
    python3 - "$manifest" "$timestamp" "$metadata" "$RUNS" "$WARMUP" "$ZEBRAC_DURATION_MS" \
        "$SKIP_REAL" "$SKIP_SCALING" "$SKIP_SIZE" \
        "$zig_version" "$zig_target" "$optimize" \
        "$(IFS=,; echo "${PRODUCT_DATASETS[*]}")" \
        "$INDEX_CORRECTNESS_STATUS" "$INDEX_CORRECTNESS_CHECKS" \
        "$INDEX_CORRECTNESS_REVIEWS" "$MESSY_DIR" <<'PY'
import json, os, sys
from pathlib import Path

manifest, ts, metadata, runs, warmup, duration = sys.argv[1:7]
skip_real, skip_scaling, skip_size = (a == "true" for a in sys.argv[7:10])
zig_version, zig_target, optimize = sys.argv[10:13]
product_workloads = tuple(filter(None, sys.argv[13].split(",")))
correctness_status, correctness_checks, correctness_reviews = sys.argv[14:17]
messy_dir = sys.argv[17]

sections = {}
if not skip_real:
    sections["real"] = f"perf_{ts}"
if not skip_scaling:
    sections["scale_seqs_budget"] = f"scale_seqs_budget_{ts}"
    sections["scale_seqs_fixed"] = f"scale_seqs_fixed_{ts}"
    if not skip_size:
        sections["scale_size"] = f"scale_size_{ts}"
if messy_dir:
    sections["messy"] = messy_dir
if correctness_status == "pass":
    sections["correctness"] = f"tests_{ts}.csv"

tools = {"samtools": os.environ.get("BENCH_VER_SAMTOOLS", "")}
for key in ("seqkit", "fastahack", "pyfaidx", "noodles", "rustbio"):
    val = os.environ.get(f"BENCH_VER_{key.upper()}")
    if val:
        tools[key] = val

out = {
    "schema_version": "index-run.v3",
    "timestamp": ts,
    "runner": "zebrac",
    "mode": "warm",
    "zebrac": os.environ["BENCH_VER_ZEBRAC"],
    "z_fasta": os.environ["BENCH_VER_ZFASTA"],
    "runs": int(runs),
    "warmup": int(warmup),
    "duration_ms": int(duration),
    "skip_real": skip_real,
    "skip_scaling": skip_scaling,
    "skip_size": skip_size,
    "product_workloads": list(product_workloads),
    "metadata": metadata,
    "tools": tools,
    "sections": sections,
    "build": {
        "zig_version": zig_version,
        "target": zig_target,
        "optimize": optimize,
    },
    "correctness": {
        "status": correctness_status,
        "checks": int(correctness_checks),
        "reviews": int(correctness_reviews),
    },
}
Path(manifest).write_text(json.dumps(out, indent=2) + "\n")
PY
}

merge_scaling_manifest() {
    local base_ts="$1" new_ts="$2"
    python3 - "$RESULTS_DIR" "$base_ts" "$new_ts" <<'PY'
import json
from pathlib import Path
import sys

results = Path(sys.argv[1])
base_ts, new_ts = sys.argv[2], sys.argv[3]
base = json.loads((results / f"run_{base_ts}.json").read_text())
new_path = results / f"run_{new_ts}.json"
if not new_path.is_file():
    raise SystemExit(f"error: scaling manifest not found: {new_path}")
new = json.loads(new_path.read_text())
if base.get("schema_version") != "index-run.v3" or new.get("schema_version") != "index-run.v3":
    raise SystemExit("error: scaling merge requires two index-run.v3 manifests")
for key in (
    "build",
    "tools",
    "z_fasta",
    "zebrac",
    "product_workloads",
    "runs",
    "warmup",
    "duration_ms",
):
    if base.get(key) != new.get(key):
        raise SystemExit(f"error: scaling merge identity mismatch: {key}")
sections = dict(base.get("sections") or {})
sections.update(new.get("sections") or {})
sections.pop("scale_seqs", None)
for key, prefix in (("scale_seqs_budget", "scale_seqs_budget"), ("scale_seqs_fixed", "scale_seqs_fixed")):
    candidate = results / f"{prefix}_{new_ts}"
    if candidate.is_dir():
        sections[key] = candidate.name
merged_metadata = results / f"metadata_{new_ts}.jsonl"
parts = []
for name in (base.get("metadata"), new.get("metadata")):
    if not isinstance(name, str):
        raise SystemExit("error: scaling merge manifest has no metadata file")
    path = results / name
    if not path.is_file():
        raise SystemExit(f"error: scaling merge metadata not found: {path}")
    parts.append(path.read_text())
merged_metadata.write_text("".join(parts))
merged = {
    **base,
    **new,
    "timestamp": new_ts,
    "runs": base["runs"],
    "warmup": base["warmup"],
    "skip_real": bool(base.get("skip_real", False)),
    "skip_scaling": False,
    "skip_size": bool(base.get("skip_size", False)),
    "metadata": merged_metadata.name,
    "sections": sections,
    "build": base["build"],
    "correctness": base["correctness"],
}
new_path.write_text(json.dumps(merged, indent=2) + "\n")
print(new_path)
PY
}

bench_add_command() {
    local section="$1" workload="$2" tool="$3" family="$4"
    local json_out="$5" script="$6" input_bytes="$7"
    zebrac_add_command "index" "$section" "$workload" "$tool" "$family" \
        "$input_bytes" "" "$json_out" "$(shell_command "$script")"
}

bench_file() {
    local file="$1" json_out="$2" metadata_jsonl="$3"
    local section="${4:-index}" workload="${5:-$(basename "$json_out" .json)}"
    local mode="${6:-product}"

    local bench_mode="$mode"
    if [[ "$mode" == "messy" ]]; then
        # Messy matrix uses production `.zfi` only; treat like headline for peer set.
        bench_mode="headline"
    fi

    local qf qz qs qk qh qp qn qr qt
    qf="$(quote_arg "$file")"
    qz="$(quote_arg "$ZFASTA")"
    qs="$(quote_arg "$SAMTOOLS")"
    qk="$(quote_arg "$SEQKIT")"
    qh="$(quote_arg "$FASTAHACK")"
    qp="$(quote_arg "$PYFAIDX")"
    qn="$(quote_arg "$NOODLES")"
    qr="$(quote_arg "$RUSTBIO")"
    qt="$(quote_arg "$SPOOL_DIR")"
    local clean="rm -f ${qf}.fai ${qf}.zfi"
    local nbytes
    nbytes="$(file_size_bytes "$file")"

    zebrac_clear_commands
    if [[ "$mode" == "messy" ]]; then
        # Production `.zfi` is the messy contract. `--emit-fai` must refuse non-uniform
        # layouts, so timing emit-fai here falsely marks z-fasta as fail in the matrix.
        bench_add_command "$section" "$workload" "z-fasta" "z-fasta" "$json_out" \
            "$clean; $qz index $qf > /dev/null" "$nbytes"
    elif [[ "$bench_mode" == "headline" ]]; then
        # Peer-comparable lane only (scaling sweeps).
        bench_add_command "$section" "$workload" "z-fasta-fai" "z-fasta" "$json_out" \
            "$clean; TMPDIR=$qt $qz index --emit-fai $qf > /dev/null" "$nbytes"
    else
        # Product matrix: {.fai, .zfi} x {dedup, --no-dedup}.
        bench_add_command "$section" "$workload" "z-fasta-fai" "z-fasta" "$json_out" \
            "$clean; TMPDIR=$qt $qz index --emit-fai $qf > /dev/null" "$nbytes"
        bench_add_command "$section" "$workload" "z-fasta-fai-nodedup" "z-fasta" "$json_out" \
            "$clean; TMPDIR=$qt $qz index --emit-fai --no-dedup $qf > /dev/null" "$nbytes"
        bench_add_command "$section" "$workload" "z-fasta-zfi" "z-fasta" "$json_out" \
            "$clean; $qz index $qf > /dev/null" "$nbytes"
        bench_add_command "$section" "$workload" "z-fasta-zfi-nodedup" "z-fasta" "$json_out" \
            "$clean; $qz index --no-dedup $qf > /dev/null" "$nbytes"
    fi

    bench_add_command "$section" "$workload" "samtools" "samtools" "$json_out" \
        "$clean; $qs faidx $qf" "$nbytes"
    bench_has_tool seqkit && bench_add_command "$section" "$workload" "seqkit" "seqkit" "$json_out" \
        "$clean; $qk faidx $qf > /dev/null 2>&1" "$nbytes"
    bench_has_tool fastahack && bench_add_command "$section" "$workload" "fastahack" "fastahack" "$json_out" \
        "$clean; $qh -i $qf > /dev/null 2>&1" "$nbytes"
    bench_has_tool pyfaidx && bench_add_command "$section" "$workload" "pyfaidx" "pyfaidx" "$json_out" \
        "$clean; $qp $qf --no-output > /dev/null 2>&1" "$nbytes"
    bench_has_tool noodles && bench_add_command "$section" "$workload" "noodles" "noodles" "$json_out" \
        "$clean; $qn index $qf" "$nbytes"
    bench_has_tool rustbio && bench_add_command "$section" "$workload" "rustbio-custom-index" "rustbio" "$json_out" \
        "$clean; $qr index $qf" "$nbytes"

    zebrac_run_current_group "$json_out" "$metadata_jsonl"
    zebrac_clear_commands
}

export_manifest_tool_versions() {
    export_manifest_core_versions
    export_manifest_optional_tool_versions seqkit fastahack pyfaidx noodles rustbio
}

# Correctness: edge cases and lightweight messy variants.

run_tests() {
    bench_require_tool z-fasta
    bench_require_tool samtools
    bench_ensure_messy --fixtures
    mkdir -p "$RESULTS_DIR" "$EDGE_DIR"

    local csv="$RESULTS_DIR/tests_${TIMESTAMP}.csv"
    local fail_log="$RESULTS_DIR/failures_${TIMESTAMP}.log"
    local zf_fai="$RESULTS_DIR/.zf_edge_${TIMESTAMP}_$$.fai"
    local pass=0 fail=0 review=0 total=0
    local file name contract_class match mark zfi_status
    local zf_exit zfi_exit sam_exit seq_exit fh_exit pyfaidx_exit noodles_exit rustbio_exit

    echo ""
    echo "Correctness: edge cases"
    echo "Generating edge case files..."

    > "$EDGE_DIR/empty.fasta"
    echo -n ">" > "$EDGE_DIR/single_byte.fasta"
    echo -n ">seq1" > "$EDGE_DIR/header_no_newline.fasta"
    echo ">seq1" > "$EDGE_DIR/header_only.fasta"
    printf ">%s description\nACGT\n" "$(head -c 10000 < /dev/zero | tr '\0' 'A')" > "$EDGE_DIR/long_header.fasta"
    printf ">longseq\n%s\n" "$(head -c 1000000 < /dev/zero | tr '\0' 'A')" > "$EDGE_DIR/no_wrap_1mb.fasta"
    printf ">special\nACGT*-NXYZ\n" > "$EDGE_DIR/special_chars.fasta"
    printf ">seq1\r\nACGT\r\nGGGG\r\n" > "$EDGE_DIR/crlf.fasta"
    printf ">seq1\nACGT\r\nGGGG\n" > "$EDGE_DIR/mixed_endings.fasta"
    printf ">seq1\nACGT" > "$EDGE_DIR/no_final_newline.fasta"
    printf ">binary\nACGT\x00\x01\x02GGGG\n" > "$EDGE_DIR/binary_data.fasta"
    printf ">seq_émoji_🧬\nACGT\n" > "$EDGE_DIR/unicode_header.fasta"
    printf ">seq1\tmore\ttabs\nACGT\n" > "$EDGE_DIR/tab_in_name.fasta"
    printf ">seq1 description with > symbols > here\nACGT\n" > "$EDGE_DIR/gt_in_description.fasta"
    printf ">seq1\nACGT>GGGG\n" > "$EDGE_DIR/gt_in_sequence.fasta"
    printf " >seq1\nACGT\n" > "$EDGE_DIR/space_before_gt.fasta"
    printf ">all_n\nNNNNNNNNNN\n" > "$EDGE_DIR/all_n.fasta"
    printf ">lower\nacgtacgt\n" > "$EDGE_DIR/lowercase.fasta"
    printf ">sp|P12345|PROT_HUMAN description\nMKWVTFISLL\n" > "$EDGE_DIR/pipe_name.fasta"
    printf ">dup\nAAAA\n>dup\nCCCC\n>dup\nGGGG\n" > "$EDGE_DIR/triple_dup.fasta"
    echo "  $(ls "$EDGE_DIR"/*.fasta 2>/dev/null | wc -l) edge case files"
    echo ""

    echo "test_case,contract_class,zfasta_exit,zfasta_zfi_exit,samtools_exit,seqkit_exit,fastahack_exit,pyfaidx_exit,noodles_exit,rustbio_exit,output_match" > "$csv"
    : > "$fail_log"

    run_index_tools() {
        local file="$1"
        rm -f "${file}.fai" "${file}.zfi" 2>/dev/null || true
        zf_exit=0; "$ZFASTA" index --emit-fai "$file" > "$zf_fai" 2>/dev/null || zf_exit=$?
        rm -f "${file}.fai" 2>/dev/null || true
        sam_exit=0; "$SAMTOOLS" faidx "$file" 2>/dev/null || sam_exit=$?
        seq_exit=127
        if bench_has_tool seqkit; then
            seq_exit=0
            rm -f "${file}.fai" 2>/dev/null || true
            "$SEQKIT" faidx "$file" > /dev/null 2>&1 || seq_exit=$?
        fi
        fh_exit=127
        if bench_has_tool fastahack; then
            fh_exit=0
            rm -f "${file}.fai" 2>/dev/null || true
            "$FASTAHACK" -i "$file" > /dev/null 2>&1 || fh_exit=$?
        fi
        pyfaidx_exit=127
        if bench_has_tool pyfaidx; then
            pyfaidx_exit=0
            rm -f "${file}.fai" 2>/dev/null || true
            "$PYFAIDX" "$file" --no-output > /dev/null 2>&1 || pyfaidx_exit=$?
        fi
        noodles_exit=127
        if bench_has_tool noodles; then
            noodles_exit=0
            rm -f "${file}.fai" 2>/dev/null || true
            "$NOODLES" index "$file" > /dev/null 2>&1 || noodles_exit=$?
        fi
        rustbio_exit=127
        if bench_has_tool rustbio; then
            rustbio_exit=0
            rm -f "${file}.fai" 2>/dev/null || true
            "$RUSTBIO" index "$file" > /dev/null 2>&1 || rustbio_exit=$?
        fi
        zfi_exit=0
        rm -f "${file}.fai" "${file}.zfi" 2>/dev/null || true
        "$ZFASTA" index "$file" > /dev/null 2>&1 || zfi_exit=$?
    }

    score_edge_case() {
        local file="$1" name="$2"
        rm -f "${file}.fai" "${file}.zfi" 2>/dev/null || true
        "$SAMTOOLS" faidx "$file" 2>/dev/null || true
        if [[ $zf_exit -ne 0 && $sam_exit -ne 0 ]]; then
            match="MATCH"
        elif [[ $zf_exit -ne 0 || $sam_exit -ne 0 ]]; then
            match="DIFF"
            echo "FAILURE: $name (exit codes differ: zf=$zf_exit sam=$sam_exit)" >> "$fail_log"
        elif [[ -f "${file}.fai" ]]; then
            if diff -q "$zf_fai" "${file}.fai" &>/dev/null; then
                match="MATCH"
            else
                match="DIFF"
                echo "FAILURE: $name (output differs)" >> "$fail_log"
                diff -u "${file}.fai" "$zf_fai" >> "$fail_log" || true
            fi
        elif [[ ! -s "$zf_fai" ]]; then
            match="MATCH"
        else
            match="DIFF"
            echo "FAILURE: $name (z-fasta produced output, samtools did not)" >> "$fail_log"
        fi
    }

    for file in "$EDGE_DIR"/*.fasta; do
        [[ -f "$file" ]] || continue
        name="$(basename "$file" .fasta)"
        contract_class="fai_parity"
        [[ "$name" == "binary_data" ]] && contract_class="invalid_input"
        total=$((total + 1))
        run_index_tools "$file"
        if [[ "$contract_class" == "invalid_input" ]]; then
            match="REVIEW"
        else
            score_edge_case "$file" "$name"
        fi
        echo "$name,$contract_class,$zf_exit,$zfi_exit,$sam_exit,$seq_exit,$fh_exit,$pyfaidx_exit,$noodles_exit,$rustbio_exit,$match" >> "$csv"
        if [[ "$match" == "MATCH" ]]; then
            pass=$((pass + 1))
            mark="pass"
        elif [[ "$match" == "REVIEW" ]]; then
            review=$((review + 1))
            mark="!"
        else
            fail=$((fail + 1))
            mark="fail"
        fi
        echo "  $mark  $name  (fai=$zf_exit zfi=$zfi_exit sam=$sam_exit) $match"
        rm -f "${file}.fai" "${file}.zfi" "$zf_fai" 2>/dev/null || true
    done

    check_zfi_side_table() {
        python3 - "$1" "$2" <<'PY'
from pathlib import Path
import sys
zfi = Path(sys.argv[1]).read_bytes()
expected = sys.argv[2]
if len(zfi) < 16 or zfi[:4] != b"ZFI\x01":
    raise SystemExit(1)
count = int.from_bytes(zfi[4:8], "little")
side = any(zfi[16 + i * 40 + 10] & 1 for i in range(count))
raise SystemExit(0 if (expected == "side-table") == side else 1)
PY
    }

    run_messy_case() {
        local file="$1" behavior="$2" zfi_expect="$3"
        local name; name=$(basename "$file" .fasta)
        contract_class="fai_parity"
        [[ "$behavior" == "zfasta-only" ]] && contract_class="zfi_messy"
        total=$((total + 1))
        run_index_tools "$file"
        match="DIFF"
        zfi_status="ok"
        # Production `.zfi` is the contract for messy layouts; `--emit-fai` must refuse them.
        if [[ $zfi_exit -ne 0 ]]; then
            zfi_status="index-failed"
        elif ! check_zfi_side_table "${file}.zfi" "$zfi_expect"; then
            zfi_status="bad-zfi"
        fi

        if [[ "$behavior" == "match-samtools" ]]; then
            rm -f "${file}.fai" "${file}.zfi" 2>/dev/null || true
            "$SAMTOOLS" faidx "$file" 2>/dev/null || true
            [[ $zf_exit -eq 0 && $sam_exit -eq 0 && -f "${file}.fai" ]] \
                && diff -q "$zf_fai" "${file}.fai" &>/dev/null \
                && [[ "$zfi_status" == "ok" ]] && match="MATCH"
        elif [[ "$behavior" == "zfasta-only" ]]; then
            # emit-fai nonzero (non-uniform), samtools fails, `.zfi` + side-table ok
            [[ $zf_exit -ne 0 && $sam_exit -ne 0 && "$zfi_status" == "ok" ]] && match="MATCH"
        fi
        echo "$name,$contract_class,$zf_exit,$zfi_exit,$sam_exit,$seq_exit,$fh_exit,$pyfaidx_exit,$noodles_exit,$rustbio_exit,$match" >> "$csv"
        if [[ "$match" == "MATCH" ]]; then
            pass=$((pass + 1))
            mark="pass"
        else
            fail=$((fail + 1))
            mark="fail"
        fi
        echo "  $mark  $name  (fai=$zf_exit zfi=$zfi_exit sam=$sam_exit) $match"
        rm -f "${file}.fai" "${file}.zfi" "$zf_fai" 2>/dev/null || true
    }

    echo ""
    echo "Messy FASTA variants (correctness):"
    run_messy_case "$MESSY_TEST_DIR/uniform.fasta" match-samtools uniform
    run_messy_case "$MESSY_TEST_DIR/mixed_widths.fasta" zfasta-only side-table
    run_messy_case "$MESSY_TEST_DIR/trailing_whitespace.fasta" zfasta-only side-table
    run_messy_case "$MESSY_TEST_DIR/blank_lines.fasta" zfasta-only side-table
    run_messy_case "$MESSY_TEST_DIR/mixed_crlf.fasta" zfasta-only side-table

    echo ""
    echo "  Correctness: $pass / $((total - review)) match ($fail diff, $review review)"
    echo "  CSV: $csv"
    [[ -s "$fail_log" ]] && echo "  Failures: $fail_log" || rm -f "$fail_log"
    if [[ "$fail" -gt 0 ]]; then
        echo "error: index run_tests failed ($fail diff)" >&2
        return 1
    fi
    INDEX_CORRECTNESS_STATUS="pass"
    INDEX_CORRECTNESS_CHECKS="$pass"
    INDEX_CORRECTNESS_REVIEWS="$review"
    return 0
}

# Zebrac performance: real datasets and scaling.

run_benchmarks() {
    bench_require_tool zebrac
    bench_require_tool z-fasta
    bench_require_tool samtools
    mkdir -p "$RESULTS_DIR" "$SCALING_DIR"

    METADATA_JSONL="$RESULTS_DIR/metadata_${TIMESTAMP}.jsonl"
    RUN_MANIFEST="$RESULTS_DIR/run_${TIMESTAMP}.json"
    if $SCALING_ONLY; then
        SPOOL_DIR="$RESULTS_DIR/spool_${MERGE_BASE}"
    else
        SPOOL_DIR="$RESULTS_DIR/spool_${TIMESTAMP}"
    fi
    mkdir -p "$SPOOL_DIR"
    : > "$METADATA_JSONL"
    ZEBRAC_MIN_SAMPLES="$RUNS"
    ZEBRAC_MAX_SAMPLES="$RUNS"
    ZEBRAC_WARMUP="$WARMUP"

    echo ""
    echo "Zebrac performance"
    echo "  Runs: $RUNS | Warmup: $WARMUP | Sampling budget: ${ZEBRAC_DURATION_MS}ms"
    echo "  FAI spool: $SPOOL_DIR"
    echo "  metadata: $METADATA_JSONL"
    echo ""

    if ! $SKIP_REAL; then
        echo "Real datasets"
        require_real_datasets
        local names=(Genome Transcriptome Proteome)
        local files=(
            "$DATA_DIR/REAL_Genome.fa"
            "$DATA_DIR/REAL_Transcriptome.fa"
            "$DATA_DIR/REAL_Proteome.fasta"
        )
        local perf_dir="$RESULTS_DIR/perf_${TIMESTAMP}"
        mkdir -p "$perf_dir"
        local i name file bench_mode product_dataset
        for i in "${!names[@]}"; do
            name="${names[$i]}"
            file="${files[$i]}"
            bench_mode="headline"
            for product_dataset in "${PRODUCT_DATASETS[@]}"; do
                if [[ "$name" == "$product_dataset" ]]; then
                    bench_mode="product"
                    break
                fi
            done
            echo "  $name ($(du -h "$file" | cut -f1))"
            bench_file "$file" "$perf_dir/${name}.json" "$METADATA_JSONL" real "$name" "$bench_mode"
        done
        echo ""
    fi

    if ! $SKIP_SCALING; then
        if $REGENERATE_FIXTURES; then
            shopt -s nullglob
            legacy=("$RESULTS_DIR"/scale_seqs_*)
            ((${#legacy[@]})) && rm -rf "${legacy[@]}"
            shopt -u nullglob
        fi
        ensure_scaling
        # Fail closed if any required family is still missing (stamp/path bugs).
        local missing=false mb c
        for mb in "${SIZE_MBS[@]}"; do [[ -f "$SCALING_DIR/size_${mb}mb.fasta" ]] || missing=true; done
        for c in "${SEQ_BUDGET_COUNTS[@]}"; do [[ -f "$SCALING_DIR/seqs_budget_${c}.fasta" ]] || missing=true; done
        for c in "${SEQ_FIXED_COUNTS[@]}"; do [[ -f "$SCALING_DIR/seqs_fixed_${c}.fasta" ]] || missing=true; done
        if $missing; then
            echo "error: scaling fixtures missing under $SCALING_DIR" >&2
            exit 1
        fi

        if ! $SKIP_SIZE; then
            echo "File-size scaling"
            local size_dir="$RESULTS_DIR/scale_size_${TIMESTAMP}"
            mkdir -p "$size_dir"
            for mb in "${SIZE_MBS[@]}"; do
                echo "  ${mb}MB"
                bench_file "$SCALING_DIR/size_${mb}mb.fasta" "$size_dir/${mb}mb.json" \
                    "$METADATA_JSONL" scale_size "${mb}mb" headline
            done
            echo ""
        fi

        echo "Sequence scaling (constant ~50 MiB)"
        local budget_dir="$RESULTS_DIR/scale_seqs_budget_${TIMESTAMP}"
        mkdir -p "$budget_dir"
        for count in "${SEQ_BUDGET_COUNTS[@]}"; do
            local f="$SCALING_DIR/seqs_budget_${count}.fasta"
            [[ -f "$f" ]] || { echo "error: missing $f" >&2; exit 1; }
            echo "  budget $count ($(du -h "$f" | cut -f1))"
            bench_file "$f" "$budget_dir/${count}.json" "$METADATA_JSONL" \
                scale_seqs_budget "$count" headline
        done
        echo ""

        echo "Sequence scaling (fixed 1024 bp)"
        local fixed_dir="$RESULTS_DIR/scale_seqs_fixed_${TIMESTAMP}"
        mkdir -p "$fixed_dir"
        for count in "${SEQ_FIXED_COUNTS[@]}"; do
            local f="$SCALING_DIR/seqs_fixed_${count}.fasta"
            [[ -f "$f" ]] || { echo "error: missing $f" >&2; exit 1; }
            echo "  fixed $count ($(du -h "$f" | cut -f1))"
            bench_file "$f" "$fixed_dir/${count}.json" "$METADATA_JSONL" \
                scale_seqs_fixed "$count" headline
        done
        echo ""
    fi

}

# Messy FASTA zebrac using proteome-derived cache fixtures.

ensure_messy_perf() {
    bench_ensure_messy --perf
}

run_messy_zebrac() {
    bench_require_tool zebrac
    bench_require_tool z-fasta
    bench_require_tool samtools
    ensure_messy_perf

    local out_dir="$RESULTS_DIR/messy_${TIMESTAMP}"
    local metadata="$out_dir/metadata.jsonl"
    mkdir -p "$out_dir"
    MESSY_DIR="$(basename "$out_dir")"

    local saved_duration="$ZEBRAC_DURATION_MS"
    local saved_min_samples="$ZEBRAC_MIN_SAMPLES"
    local saved_max_samples="$ZEBRAC_MAX_SAMPLES"
    local saved_warmup="$ZEBRAC_WARMUP"
    local saved_allow_failures="$ZEBRAC_ALLOW_FAILURES"
    ZEBRAC_DURATION_MS=1000
    ZEBRAC_MIN_SAMPLES=3
    ZEBRAC_MAX_SAMPLES=3
    ZEBRAC_WARMUP=1
    ZEBRAC_ALLOW_FAILURES=true

    echo ""
    echo ""
    echo "Messy FASTA zebrac: $out_dir"
    for fasta in "$MESSY_ZEBRAC_DIR"/*.fasta; do
        [[ -f "$fasta" ]] || continue
        echo "  $(basename "$fasta" .fasta)"
        bench_file "$fasta" "$out_dir/$(basename "$fasta" .fasta).json" "$metadata" \
            messy "$(basename "$fasta" .fasta)" messy
    done
    # Peers may leave empty/failed .fai beside proteome fixtures; drop sidecars only.
    rm -f "$MESSY_ZEBRAC_DIR"/*.fai "$MESSY_ZEBRAC_DIR"/*.zfi

    ZEBRAC_DURATION_MS="$saved_duration"
    ZEBRAC_MIN_SAMPLES="$saved_min_samples"
    ZEBRAC_MAX_SAMPLES="$saved_max_samples"
    ZEBRAC_WARMUP="$saved_warmup"
    ZEBRAC_ALLOW_FAILURES="$saved_allow_failures"
}

# Report generation.

run_report() {
    local py
    py="$(report_python)"
    local args=("$RESULTS_DIR")
    [[ -n "$RUN_MANIFEST" ]] && args+=(--manifest "$RUN_MANIFEST")
    echo ""
    echo "Generating index report"
    "$py" "$SCRIPT_DIR/generate_report.py" "${args[@]}"
    if [[ -n "$RUN_MANIFEST" ]]; then
        printf '%s\n' "$TIMESTAMP" > "$RESULTS_DIR/LATEST"
    fi
}

# Main.

TIMESTAMP=$(date +%Y%m%d_%H%M%S)

if $DO_BENCHMARKS; then
    prepare_benchmark_binary
fi
if $DO_TESTS; then
    run_tests || exit 1
fi
$DO_BENCHMARKS && run_benchmarks
$DO_MESSY && run_messy_zebrac
if $DO_BENCHMARKS; then
    export_manifest_tool_versions
    write_run_manifest "$RUN_MANIFEST" "$TIMESTAMP" "$(basename "$METADATA_JSONL")"
    if $SCALING_ONLY; then
        merge_scaling_manifest "$MERGE_BASE" "$TIMESTAMP"
    fi
    echo "  manifest: $RUN_MANIFEST"
fi
$DO_REPORT && run_report

echo ""
echo "Index suite complete"
echo ""
