#!/usr/bin/env bash
# Index benchmark master runner (correctness, zebrac performance, messy zebrac, report).
#
# Usage:
#   bash bench/index/run.sh [options]
#
# Default (no options): run every section in order, then write REPORT.md.
#
# ── Defaults ────────────────────────────────────────────────────────
#   --runs 5              zebrac measured samples (real + scaling benchmarks)
#   --warmup 1            zebrac warmup samples before each measured sample
#   --duration 5000       zebrac max ms per sample (override: ZEBRAC_DURATION_MS env)
#
#   Sections enabled: tests, benchmarks, messy zebrac, report (all true).
#   Scaling fixtures: generated on demand under bench/shared/cache/scaling/ if missing.
#
# ── What --runs / --warmup apply to ─────────────────────────────────
#   YES: zebrac performance (real datasets + all three scaling sweeps).
#   NO:  correctness tests (single exit-code pass per tool per file).
#   NO:  messy zebrac (fixed: 3 samples, 1 warmup, 1000 ms, --allow-failures).
#
# ── Common invocations ──────────────────────────────────────────────
#   bash bench/index/run.sh
#   bash bench/index/run.sh --runs 10 --warmup 2
#   bash bench/index/run.sh --skip-report              # data only; report later
#   bash bench/index/run.sh --skip-tests --skip-messy  # zebrac + report
#   bash bench/index/run.sh --skip-verify --skip-perf  # same via aliases
#   bash bench/index/run.sh --skip-benchmarks --skip-messy --skip-report  # correctness only
#   bash bench/index/run.sh --skip-real                # scaling sweeps only
#   bash bench/index/run.sh --skip-scaling             # real datasets only
#   bash bench/index/run.sh --skip-size                # skip file-size sweep
#   bash bench/index/run.sh --scaling-only --merge-base 20260630_231053
#   bash bench/index/run.sh --regenerate-fixtures      # force new scaling FASTA
#   bash bench/index/run.sh --clean-legacy             # drop old seqs_*.fasta + regen
#   bash bench/index/run.sh --allow-incomplete         # report from partial runs
#
# ── Skip flags ──────────────────────────────────────────────────────
# Canonical: --skip-tests --skip-benchmarks --skip-messy --skip-report
# Deprecated aliases (one release cycle): --skip-verify --skip-perf
# --skip-messy skips messy zebrac perf only (never skips messy cases inside run_tests).
#
#   1. Correctness   edge_cases/ + cache/messy_fixtures/ → results/tests_<ts>.csv
#   2. Benchmarks    REAL_* + scaling → perf_*, scale_*, metadata_*, run_*.json
#   3. Messy zebrac  shared/cache/messy_perf/ → results/messy_<ts>/
#   4. Report        generate_report.py → REPORT.md + results/figures/
#
# ── Scaling fixtures (bench/shared/cache/scaling/) ──────────────────
#   size_{N}mb.fasta         1–1000 MB, 100 sequences each
#   seqs_budget_{N}.fasta    ~50 MiB total, N ∈ 1k 10k 100k 250k
#   seqs_fixed_{N}.fasta     1024 bp/seq, N ∈ 100k 250k 500k 1M
#   Materialize: python3 bench/shared/generate_scaling.py [--force]
#   --regenerate-fixtures / --clean-legacy force rebuild (legacy name cleanup too).
#   --scaling-only: skip real+size, regen seq fixtures, re-bench seq sweeps,
#                   merge into --merge-base manifest (defaults to results/LATEST).
#
# ── Outputs ─────────────────────────────────────────────────────────
#   results/LATEST           pointer to newest run_<timestamp>.json
#   results/run_<ts>.json    manifest (sections, tool versions, skip flags)
#   results/metadata_<ts>.jsonl
#
# ── Notes ───────────────────────────────────────────────────────────
#   Partial runs overwrite LATEST; use --allow-incomplete for report drafts.
#   generate_report.py refuses incomplete manifests unless --allow-incomplete.
#   Headline scaling uses z-fasta default only; all three modes on real data.
#
#   -h|--help  print this header

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

# ── Scaling sweep constants ─────────────────────────────────────────
SIZE_MBS=(1 5 10 25 50 100 250 500 1000)
SEQ_BUDGET_COUNTS=(1000 10000 100000 250000)
SEQ_FIXED_COUNTS=(100000 250000 500000 1000000)

# ── Defaults ───────────────────────────────────────────────────────
RUNS=5
WARMUP=1
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
ALLOW_INCOMPLETE=false

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
        --allow-incomplete)  ALLOW_INCOMPLETE=true; shift ;;
        -h|--help)
            sed -n '2,/^set -euo pipefail$/p' "$0" | head -n -1
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

# ══════════════════════════════════════════════════════════════════════
#  Helpers: scaling ensure, manifest, zebrac bench
# ══════════════════════════════════════════════════════════════════════

ensure_scaling() {
    local args=(--mode all)
    if $CLEAN_LEGACY; then
        args+=(--clean-legacy --force)
    elif $REGENERATE_FIXTURES; then
        args+=(--force)
    fi
    bench_ensure_scaling "${args[@]}"
}

write_run_manifest() {
    local manifest="$1" timestamp="$2" metadata="$3"
    python3 - "$manifest" "$timestamp" "$metadata" "$RUNS" "$WARMUP" "$ZEBRAC_DURATION_MS" \
        "$SKIP_REAL" "$SKIP_SCALING" "$SKIP_SIZE" <<'PY'
import json, os, sys
from pathlib import Path

manifest, ts, metadata, runs, warmup, duration = sys.argv[1:7]
skip_real, skip_scaling, skip_size = (a == "true" for a in sys.argv[7:10])

sections = {}
if not skip_real:
    sections["real"] = f"perf_{ts}"
if not skip_scaling:
    sections["scale_seqs_budget"] = f"scale_seqs_budget_{ts}"
    sections["scale_seqs_fixed"] = f"scale_seqs_fixed_{ts}"
    if not skip_size:
        sections["scale_size"] = f"scale_size_{ts}"

tools = {"samtools": os.environ.get("BENCH_VER_SAMTOOLS", "")}
for key in ("seqkit", "fastahack", "pyfaidx", "noodles", "rustbio"):
    val = os.environ.get(f"BENCH_VER_{key.upper()}")
    if val:
        tools[key] = val

out = {
    "schema_version": "index-run.v1",
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
    "metadata": metadata,
    "tools": tools,
    "sections": sections,
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
new = json.loads(new_path.read_text()) if new_path.is_file() else {
    "runs": base.get("runs"),
    "warmup": base.get("warmup"),
    "duration_ms": base.get("duration_ms"),
    "zebrac": base.get("zebrac"),
    "z_fasta": base.get("z_fasta"),
    "tools": base.get("tools", {}),
}
sections = dict(base.get("sections") or {})
sections.update(new.get("sections") or {})
sections.pop("scale_seqs", None)
for key, prefix in (("scale_seqs_budget", "scale_seqs_budget"), ("scale_seqs_fixed", "scale_seqs_fixed")):
    candidate = results / f"{prefix}_{new_ts}"
    if candidate.is_dir():
        sections[key] = candidate.name
merged_metadata = results / f"metadata_{new_ts}.jsonl"
parts = []
for name in (f"metadata_{base_ts}.jsonl", f"metadata_{new_ts}.jsonl"):
    path = results / name
    if path.is_file():
        parts.append(path.read_text())
merged_metadata.write_text("".join(parts))
merged = {
    **base,
    **new,
    "timestamp": new_ts,
    "base_run": base_ts,
    "seq_scaling_revision": "budget+fixed.v1",
    "runs": new.get("runs", base.get("runs")),
    "warmup": new.get("warmup", base.get("warmup")),
    "skip_real": bool(base.get("skip_real", False)),
    "skip_scaling": False,
    "metadata": merged_metadata.name,
    "sections": sections,
}
new_path.write_text(json.dumps(merged, indent=2) + "\n")
(results / "LATEST").write_text(new_ts + "\n")
print(new_path)
PY
}

bench_add_command() {
    local section="$1" workload="$2" tool="$3" family="$4"
    local json_out="$5" script="$6" input_bytes="$7"
    zebrac_add_command "index" "$section" "$workload" "$tool" "$family" \
        "$input_bytes" "" "$json_out" "$(shell_command "$script")"
}

# Write `.zfi` (z-fasta) and `.fai` (samtools) on REAL_* fixtures after zebrac.
# Zebrac lanes delete sidecars between commands; the report size table reads these files.
preserve_real_index_sidecars() {
    bench_require_tool z-fasta
    echo "  Preserving REAL_* index sidecars for report (.zfi + .fai)..."
    local fa
    for fa in \
        "$DATA_DIR/REAL_Genome.fa" \
        "$DATA_DIR/REAL_Transcriptome.fa" \
        "$DATA_DIR/REAL_Proteome.fasta"; do
        [[ -f "$fa" ]] || continue
        "$ZFASTA" index "$fa" > /dev/null
        if bench_has_tool samtools; then
            samtools faidx "$fa" > /dev/null 2>&1 || true
        fi
    done
}

bench_file() {
    local file="$1" json_out="$2" metadata_jsonl="$3"
    local section="${4:-index}" workload="${5:-$(basename "$json_out" .json)}"
    local mode="${6:-all}"

    local zfasta_tool="z-fasta-default" bench_mode="$mode"
    if [[ "$mode" == "messy" ]]; then
        zfasta_tool="z-fasta"
        bench_mode="headline"
    fi

    local qf qz qs qk qh qp qn qr
    qf="$(quote_arg "$file")"
    qz="$(quote_arg "$ZFASTA")"
    qs="$(quote_arg "$SAMTOOLS")"
    qk="$(quote_arg "$SEQKIT")"
    qh="$(quote_arg "$FASTAHACK")"
    qp="$(quote_arg "$PYFAIDX")"
    qn="$(quote_arg "$NOODLES")"
    qr="$(quote_arg "$RUSTBIO")"
    local clean="rm -f ${qf}.fai ${qf}.zfi"
    local clean_zfi="rm -f ${qf}.zfi"
    local nbytes
    nbytes="$(file_size_bytes "$file")"

    zebrac_clear_commands
    bench_add_command "$section" "$workload" "$zfasta_tool" "z-fasta" "$json_out" \
        "$clean; $qz index --emit-fai $qf > /dev/null" "$nbytes"

    if [[ "$bench_mode" != "headline" ]]; then
        bench_add_command "$section" "$workload" "z-fasta-nodedup" "z-fasta" "$json_out" \
            "$clean; $qz index --emit-fai --no-dedup $qf > /dev/null" "$nbytes"
        bench_add_command "$section" "$workload" "z-fasta-lowmem" "z-fasta" "$json_out" \
            "$clean; $qz index --low-mem --emit-fai $qf > /dev/null" "$nbytes"
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

    # Last: production `.zfi` lane. Only remove `.zfi` so competitor `.fai` remains for
    # the report on-disk size table (`bench/index/generate_report.py`).
    if [[ "$bench_mode" != "headline" ]]; then
        bench_add_command "$section" "$workload" "z-fasta-zfi" "z-fasta" "$json_out" \
            "$clean_zfi; $qz index $qf > /dev/null" "$nbytes"
    fi

    zebrac_run_current_group "$json_out" "$metadata_jsonl"
    zebrac_clear_commands
}

export_manifest_tool_versions() {
    export BENCH_VER_ZEBRAC="$(bench_tool_version zebrac)"
    export BENCH_VER_ZFASTA="$(bench_tool_version z-fasta 2>/dev/null || echo unknown)"
    export BENCH_VER_SAMTOOLS="$(bench_tool_version samtools 2>/dev/null || echo unknown)"
    for tool in seqkit fastahack pyfaidx noodles rustbio; do
        if bench_has_tool "$tool"; then
            upper="${tool^^}"
            export "BENCH_VER_${upper}=$(bench_tool_version "$tool")"
        fi
    done
}

# ══════════════════════════════════════════════════════════════════════
#  Correctness: edge cases + lightweight messy variants
# ══════════════════════════════════════════════════════════════════════

run_tests() {
    bench_require_tool z-fasta
    bench_require_tool samtools
    bench_ensure_messy --fixtures
    mkdir -p "$RESULTS_DIR" "$EDGE_DIR"

    local csv="$RESULTS_DIR/tests_${TIMESTAMP}.csv"
    local fail_log="$RESULTS_DIR/failures.log"
    local pass=0 fail=0 total=0

    echo "════════════════════════════════════════════════════════════════"
    echo "  Correctness: edge cases"
    echo "════════════════════════════════════════════════════════════════"
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

    echo "test_case,zfasta_exit,samtools_exit,seqkit_exit,fastahack_exit,noodles_exit,rustbio_exit,output_match" > "$csv"
    : > "$fail_log"

    run_index_tools() {
        local file="$1"
        rm -f "${file}.fai" "${file}.zfi" 2>/dev/null || true
        zf_exit=0; "$ZFASTA" index --emit-fai "$file" > /tmp/zf_edge_test.fai 2>/dev/null || zf_exit=$?
        rm -f "${file}.fai" 2>/dev/null || true
        sam_exit=0; $SAMTOOLS faidx "$file" 2>/dev/null || sam_exit=$?
        seq_exit=0
        if bench_has_tool seqkit; then
            rm -f "${file}.fai" 2>/dev/null || true
            "$SEQKIT" faidx "$file" > /dev/null 2>&1 || seq_exit=$?
        fi
        fh_exit=0
        if bench_has_tool fastahack; then
            rm -f "${file}.fai" 2>/dev/null || true
            "$FASTAHACK" -i "$file" > /dev/null 2>&1 || fh_exit=$?
        fi
        noodles_exit=0
        if bench_has_tool noodles; then
            rm -f "${file}.fai" 2>/dev/null || true
            "$NOODLES" index "$file" > /dev/null 2>&1 || noodles_exit=$?
        fi
        rustbio_exit=0
        if bench_has_tool rustbio; then
            rm -f "${file}.fai" 2>/dev/null || true
            "$RUSTBIO" index "$file" > /dev/null 2>&1 || rustbio_exit=$?
        fi
    }

    score_edge_case() {
        local file="$1" name="$2"
        rm -f "${file}.fai" "${file}.zfi" 2>/dev/null || true
        $SAMTOOLS faidx "$file" 2>/dev/null || true
        if [[ $zf_exit -ne 0 && $sam_exit -ne 0 ]]; then
            match="MATCH"
        elif [[ $zf_exit -ne 0 && $sam_exit -eq 0 ]]; then
            # `--emit-fai` refuses non-representable layouts; default `.zfi` may still succeed.
            rm -f "${file}.zfi" 2>/dev/null || true
            if "$ZFASTA" index "$file" >/dev/null 2>&1; then
                match="MATCH"
            else
                match="DIFF"
                echo "--- FAILURE: $name (emit-fai rejected and .zfi index failed; sam=$sam_exit) ---" >> "$fail_log"
            fi
        elif [[ $zf_exit -ne 0 || $sam_exit -ne 0 ]]; then
            match="DIFF"
            echo "--- FAILURE: $name (exit codes differ: zf=$zf_exit sam=$sam_exit) ---" >> "$fail_log"
        elif [[ -f "${file}.fai" ]]; then
            if diff -q /tmp/zf_edge_test.fai "${file}.fai" &>/dev/null; then
                match="MATCH"
            else
                match="DIFF"
                echo "--- FAILURE: $name (output differs) ---" >> "$fail_log"
                diff -u "${file}.fai" /tmp/zf_edge_test.fai >> "$fail_log" || true
            fi
        elif [[ ! -s /tmp/zf_edge_test.fai ]]; then
            match="MATCH"
        else
            match="DIFF"
            echo "--- FAILURE: $name (z-fasta produced output, samtools did not) ---" >> "$fail_log"
        fi
    }

    for file in "$EDGE_DIR"/*.fasta; do
        [[ -f "$file" ]] || continue
        name=$(basename "$file" .fasta)
        total=$((total + 1))
        run_index_tools "$file"
        score_edge_case "$file" "$name"
        echo "$name,$zf_exit,$sam_exit,$seq_exit,$fh_exit,$noodles_exit,$rustbio_exit,$match" >> "$csv"
        [[ "$match" == "MATCH" ]] && pass=$((pass + 1)) || fail=$((fail + 1))
        echo "  $([[ "$match" == MATCH ]] && echo ✓ || echo ✗)  $name  (zf=$zf_exit sam=$sam_exit) $match"
        rm -f "${file}.fai" "${file}.zfi" /tmp/zf_edge_test.fai 2>/dev/null || true
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
        total=$((total + 1))
        run_index_tools "$file"
        match="DIFF"
        zfi_status="ok"
        # Production `.zfi` is the contract for messy layouts; `--emit-fai` must refuse them.
        rm -f "${file}.fai" "${file}.zfi" 2>/dev/null || true
        "$ZFASTA" index "$file" > /dev/null 2>&1 || zfi_status="index-failed"
        [[ "$zfi_status" == "ok" ]] && check_zfi_side_table "${file}.zfi" "$zfi_expect" || zfi_status="bad-zfi"

        if [[ "$behavior" == "match-samtools" ]]; then
            rm -f "${file}.fai" "${file}.zfi" 2>/dev/null || true
            $SAMTOOLS faidx "$file" 2>/dev/null || true
            [[ $zf_exit -eq 0 && $sam_exit -eq 0 && -f "${file}.fai" ]] \
                && diff -q /tmp/zf_edge_test.fai "${file}.fai" &>/dev/null \
                && [[ "$zfi_status" == "ok" ]] && match="MATCH"
        elif [[ "$behavior" == "zfasta-only" ]]; then
            # emit-fai nonzero (non-uniform), samtools fails, `.zfi` + side-table ok
            [[ $zf_exit -ne 0 && $sam_exit -ne 0 && "$zfi_status" == "ok" ]] && match="MATCH"
        fi
        echo "$name,$zf_exit,$sam_exit,$seq_exit,$fh_exit,$noodles_exit,$rustbio_exit,$match" >> "$csv"
        [[ "$match" == "MATCH" ]] && pass=$((pass + 1)) || fail=$((fail + 1))
        echo "  $([[ "$match" == MATCH ]] && echo ✓ || echo ✗)  $name  (zfi=$zfi_expect) $match"
        rm -f "${file}.fai" "${file}.zfi" /tmp/zf_edge_test.fai 2>/dev/null || true
    }

    echo ""
    echo "Messy FASTA variants (correctness):"
    run_messy_case "$MESSY_TEST_DIR/uniform.fasta" match-samtools uniform
    run_messy_case "$MESSY_TEST_DIR/mixed_widths.fasta" zfasta-only side-table
    run_messy_case "$MESSY_TEST_DIR/trailing_whitespace.fasta" zfasta-only side-table
    run_messy_case "$MESSY_TEST_DIR/blank_lines.fasta" zfasta-only side-table
    run_messy_case "$MESSY_TEST_DIR/mixed_crlf.fasta" zfasta-only side-table

    echo ""
    echo "  Correctness: $pass / $total match ($fail diff)"
    echo "  CSV: $csv"
    [[ -s "$fail_log" ]] && echo "  Failures: $fail_log" || rm -f "$fail_log"
    if [[ "$fail" -gt 0 ]]; then
        echo "error: index run_tests failed ($fail diff)" >&2
        return 1
    fi
    return 0
}

# ══════════════════════════════════════════════════════════════════════
#  Zebrac performance: real datasets + scaling
# ══════════════════════════════════════════════════════════════════════

run_benchmarks() {
    bench_require_tool zebrac
    bench_require_tool z-fasta
    bench_require_tool samtools
    mkdir -p "$RESULTS_DIR" "$SCALING_DIR"

    METADATA_JSONL="$RESULTS_DIR/metadata_${TIMESTAMP}.jsonl"
    RUN_MANIFEST="$RESULTS_DIR/run_${TIMESTAMP}.json"
    : > "$METADATA_JSONL"
    ZEBRAC_MIN_SAMPLES="$RUNS"
    ZEBRAC_MAX_SAMPLES="$RUNS"
    ZEBRAC_WARMUP="$WARMUP"

    echo ""
    echo "════════════════════════════════════════════════════════════════"
    echo "  Zebrac performance"
    echo "════════════════════════════════════════════════════════════════"
    echo "  Runs: $RUNS | Warmup: $WARMUP | Duration: ${ZEBRAC_DURATION_MS}ms"
    echo "  metadata: $METADATA_JSONL"
    echo ""

    if ! $SKIP_REAL; then
        echo "──────────────────────────────────────────────────"
        echo " Real datasets"
        echo "──────────────────────────────────────────────────"
        declare -A DATASETS
        [[ -f "$DATA_DIR/REAL_Genome.fa" ]]        && DATASETS["Genome"]="$DATA_DIR/REAL_Genome.fa"
        [[ -f "$DATA_DIR/REAL_Transcriptome.fa" ]]  && DATASETS["Transcriptome"]="$DATA_DIR/REAL_Transcriptome.fa"
        [[ -f "$DATA_DIR/REAL_Proteome.fasta" ]]    && DATASETS["Proteome"]="$DATA_DIR/REAL_Proteome.fasta"
        if [[ ${#DATASETS[@]} -eq 0 ]]; then
            echo "  Fetching REAL_* datasets..."
            bash "$BENCH_ROOT/shared/download_data.sh"
            [[ -f "$DATA_DIR/REAL_Genome.fa" ]]       && DATASETS["Genome"]="$DATA_DIR/REAL_Genome.fa"
            [[ -f "$DATA_DIR/REAL_Transcriptome.fa" ]] && DATASETS["Transcriptome"]="$DATA_DIR/REAL_Transcriptome.fa"
            [[ -f "$DATA_DIR/REAL_Proteome.fasta" ]]  && DATASETS["Proteome"]="$DATA_DIR/REAL_Proteome.fasta"
        fi
        if [[ ${#DATASETS[@]} -eq 0 ]]; then
            echo "  ERROR: no real datasets available" >&2
        else
            local perf_dir="$RESULTS_DIR/perf_${TIMESTAMP}"
            mkdir -p "$perf_dir"
            for name in "${!DATASETS[@]}"; do
                local file="${DATASETS[$name]}"
                echo "  $name ($(du -h "$file" | cut -f1))"
                bench_file "$file" "$perf_dir/${name}.json" "$METADATA_JSONL" index "$name" all
            done
            preserve_real_index_sidecars
        fi
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
            echo "──────────────────────────────────────────────────"
            echo " File-size scaling"
            echo "──────────────────────────────────────────────────"
            local size_dir="$RESULTS_DIR/scale_size_${TIMESTAMP}"
            mkdir -p "$size_dir"
            for mb in "${SIZE_MBS[@]}"; do
                echo "  ${mb}MB"
                bench_file "$SCALING_DIR/size_${mb}mb.fasta" "$size_dir/${mb}mb.json" \
                    "$METADATA_JSONL" index "${mb}mb" headline
            done
            echo ""
        fi

        echo "──────────────────────────────────────────────────"
        echo " Sequence scaling (bounded ~50 MiB)"
        echo "──────────────────────────────────────────────────"
        local budget_dir="$RESULTS_DIR/scale_seqs_budget_${TIMESTAMP}"
        mkdir -p "$budget_dir"
        for count in "${SEQ_BUDGET_COUNTS[@]}"; do
            local f="$SCALING_DIR/seqs_budget_${count}.fasta"
            [[ -f "$f" ]] || { echo "error: missing $f" >&2; exit 1; }
            echo "  budget $count ($(du -h "$f" | cut -f1))"
            bench_file "$f" "$budget_dir/${count}.json" "$METADATA_JSONL" index "$count" headline
        done
        echo ""

        echo "──────────────────────────────────────────────────"
        echo " Sequence scaling (fixed 1024 bp)"
        echo "──────────────────────────────────────────────────"
        local fixed_dir="$RESULTS_DIR/scale_seqs_fixed_${TIMESTAMP}"
        mkdir -p "$fixed_dir"
        for count in "${SEQ_FIXED_COUNTS[@]}"; do
            local f="$SCALING_DIR/seqs_fixed_${count}.fasta"
            [[ -f "$f" ]] || { echo "error: missing $f" >&2; exit 1; }
            echo "  fixed $count ($(du -h "$f" | cut -f1))"
            bench_file "$f" "$fixed_dir/${count}.json" "$METADATA_JSONL" index "$count" headline
        done
        echo ""
    fi

    export_manifest_tool_versions
    write_run_manifest "$RUN_MANIFEST" "$TIMESTAMP" "$(basename "$METADATA_JSONL")"
    printf '%s\n' "$TIMESTAMP" > "$RESULTS_DIR/LATEST"
    echo "  manifest: $RUN_MANIFEST"

    if $SCALING_ONLY; then
        merge_scaling_manifest "$MERGE_BASE" "$TIMESTAMP"
    fi
}

# ══════════════════════════════════════════════════════════════════════
#  Messy FASTA zebrac (proteome-derived cache fixtures)
# ══════════════════════════════════════════════════════════════════════

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

    ZEBRAC_DURATION_MS=1000
    ZEBRAC_MIN_SAMPLES=3
    ZEBRAC_WARMUP=1
    ZEBRAC_ALLOW_FAILURES=true

    echo ""
    echo "════════════════════════════════════════════════════════════════"
    echo "  Messy FASTA zebrac → $out_dir"
    echo "════════════════════════════════════════════════════════════════"
    for fasta in "$MESSY_ZEBRAC_DIR"/*.fasta; do
        [[ -f "$fasta" ]] || continue
        echo "  → $(basename "$fasta" .fasta)"
        bench_file "$fasta" "$out_dir/$(basename "$fasta" .fasta).json" "$metadata" \
            messy "$(basename "$fasta" .fasta)" messy
    done
    # Peers may leave empty/failed .fai beside proteome fixtures; drop sidecars only.
    rm -f "$MESSY_ZEBRAC_DIR"/*.fai "$MESSY_ZEBRAC_DIR"/*.zfi
}

# ══════════════════════════════════════════════════════════════════════
#  Report
# ══════════════════════════════════════════════════════════════════════

run_report() {
  local py; py="$(report_python)"
  local args=()
  $ALLOW_INCOMPLETE && args+=(--allow-incomplete)
  echo ""
  echo "════════════════════════════════════════════════════════════════"
  echo "  Generating REPORT.md"
  echo "════════════════════════════════════════════════════════════════"
  "$py" "$SCRIPT_DIR/generate_report.py" "${args[@]}"
}

# ══════════════════════════════════════════════════════════════════════
#  Main
# ══════════════════════════════════════════════════════════════════════

TIMESTAMP=$(date +%Y%m%d_%H%M%S)

if $DO_TESTS; then
    run_tests || exit 1
fi
$DO_BENCHMARKS && run_benchmarks
$DO_MESSY && run_messy_zebrac
$DO_REPORT && run_report

echo ""
echo "════════════════════════════════════════════════════════════════"
echo "  Index suite complete"
echo "════════════════════════════════════════════════════════════════"
