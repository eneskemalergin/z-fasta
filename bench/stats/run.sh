#!/usr/bin/env bash
# Stats benchmark runner: correctness (run_tests), then zebrac perf.
#
# Usage:
#   bash bench/stats/run.sh [options]
#
# Defaults: correctness first (89 checks), then perf (--runs 5, --warmup 1, --duration 5000).
# Report generation: python3 bench/stats/generate_report.py
#   (run.sh calls this after perf; --allow-incomplete for drafts)
# Use --skip-tests (alias --skip-verify) to run perf/report only. Correctness failure aborts before perf or report.
# --skip-messy is accepted for suite parity (messy zebrac perf); stats has no messy perf section today.
#
# Canonical flags (index-aligned): --skip-tests --skip-benchmarks --skip-messy --skip-report
# Deprecated aliases (one release cycle): --skip-verify --skip-perf
#
#   bash bench/stats/run.sh
#   bash bench/stats/run.sh --skip-tests --skip-report   # perf only
#   bash bench/stats/run.sh --skip-verify --skip-perf    # same via aliases
#   bash bench/stats/run.sh --skip-full --skip-mode       # scaling only
#   bash bench/stats/run.sh --skip-scale                  # REAL full + mode only
#   bash bench/stats/run.sh --skip-size                   # skip file-size sweep only
#   bash bench/stats/run.sh --skip-seqs                   # skip seq-count sweep only
#   bash bench/stats/run.sh --regenerate-fixtures
#   bash bench/stats/run.sh --allow-incomplete
#   bash bench/stats/run.sh --skip-benchmarks --skip-report   # correctness only
#   STATS_RUN_TIMESTAMP=<ts> bash bench/stats/run.sh --skip-tests ...  # resume same run id
#
# Scaling fixtures (shared cache; gitignored):
#   bench/shared/cache/scaling/size_{N}mb.fasta
#   bench/shared/cache/scaling/seqs_fixed_{N}.fasta
#   Same tree as index (budget fixtures live there too; stats does not bench them).
#   Materialize: python3 bench/shared/generate_scaling.py
#   Indexes (.zfi + .fai) are preloaded once per file; timed stats only loads them.
#
# Not in this runner (deferred):
#   scale_seqs_budget, messy zebrac, layout-twin zebrac
#
# Outputs:
#   results/LATEST           pointer to newest run_<timestamp>.json
#   results/run_<ts>.json    manifest (sections, tool versions, skip flags)
#   results/metadata_<ts>.jsonl
#   results/perf_full_<ts>/  peer + z-fasta full (REAL)
#   results/perf_mode_<ts>/  z-fasta full / full (fai) / indexed (zfi) / indexed (fai)
#   results/scale_size_<ts>/ file-size scaling
#   results/scale_seqs_fixed_<ts>/ sequence-count scaling
#
#   -h|--help  print this header

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH_ROOT="$(dirname "$SCRIPT_DIR")"
PROJECT_ROOT="$(dirname "$BENCH_ROOT")"
RESULTS_DIR="$SCRIPT_DIR/results"
SCALING_DIR="$BENCH_ROOT/shared/cache/scaling"
DATA_DIR="$BENCH_ROOT/shared/data"

source "$BENCH_ROOT/shared/runner_common.sh"

# Scaling sweep constants (match index bench shape)
SIZE_MBS=(1 5 10 25 50 100 250 500 1000)
SEQ_FIXED_COUNTS=(100000 250000 500000 1000000)

# Defaults
RUNS=5
WARMUP=1
ZEBRAC_DURATION_MS="${ZEBRAC_DURATION_MS:-5000}"
DO_TESTS=true
DO_BENCHMARKS=true
DO_FULL=true
DO_MODE=true
DO_SCALE=true
DO_SCALE_SIZE=true
DO_SCALE_SEQS=true
DO_MESSY=true
DO_REPORT=true
REGENERATE_FIXTURES=false
ALLOW_INCOMPLETE=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --runs) RUNS="$2"; shift 2 ;;
        --warmup) WARMUP="$2"; shift 2 ;;
        --duration) ZEBRAC_DURATION_MS="$2"; shift 2 ;;
        --skip-tests|--skip-verify) DO_TESTS=false; shift ;;
        --skip-benchmarks|--skip-perf) DO_BENCHMARKS=false; shift ;;
        --skip-full) DO_FULL=false; shift ;;
        --skip-mode) DO_MODE=false; shift ;;
        --skip-scale) DO_SCALE=false; shift ;;
        --skip-size) DO_SCALE_SIZE=false; shift ;;
        --skip-seqs) DO_SCALE_SEQS=false; shift ;;
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

stats_add_command() {
    local section="$1" workload="$2" tool="$3" family="$4"
    local json_out="$5" script="$6" input_bytes="$7"
    zebrac_add_command "stats" "$section" "$workload" "$tool" "$family" \
        "$input_bytes" "" "$json_out" "$(shell_command "$script")"
}

ensure_scaling_fixtures() {
    local args=()
    if $REGENERATE_FIXTURES; then
        args+=(--force)
    fi
    if ! $DO_SCALE_SIZE && ! $DO_SCALE_SEQS; then
        return 0
    fi
    echo "  Ensuring scaling fixtures under $SCALING_DIR ..."
    if $DO_SCALE_SIZE; then
        bench_ensure_scaling --mode size "${args[@]}"
    fi
    if $DO_SCALE_SEQS; then
        bench_ensure_scaling --mode fixed "${args[@]}"
    fi
}

preload_scaling_indexes() {
    bench_require_tool z-fasta
    echo "  Preloading scaling indexes (.zfi + .fai)..."
    local fa
    if $DO_SCALE_SIZE; then
        local mb
        for mb in "${SIZE_MBS[@]}"; do
            fa="$SCALING_DIR/size_${mb}mb.fasta"
            [[ -f "$fa" ]] || continue
            preload_indexes_for_file "$fa"
        done
    fi
    if $DO_SCALE_SEQS; then
        local count
        for count in "${SEQ_FIXED_COUNTS[@]}"; do
            fa="$SCALING_DIR/seqs_fixed_${count}.fasta"
            [[ -f "$fa" ]] || continue
            preload_indexes_for_file "$fa"
        done
    fi
}

export_manifest_tool_versions() {
    export_manifest_core_versions
    export_manifest_required_tool_versions seqtk
    export_manifest_optional_tool_versions seqkit noodles rustbio
}

write_run_manifest() {
    local manifest="$1" timestamp="$2" metadata="$3"
    python3 - "$manifest" "$timestamp" "$metadata" "$RUNS" "$WARMUP" "$ZEBRAC_DURATION_MS" \
        "$DO_TESTS" \
        "${SECTION_FULL:-}" "${SECTION_MODE:-}" \
        "${SECTION_SCALE_SIZE:-}" "${SECTION_SCALE_SEQS:-}" <<'PY'
import json, os, sys
from pathlib import Path

manifest, ts, metadata, runs, warmup, duration = sys.argv[1:7]
do_verify = sys.argv[7] == "true"
section_full, section_mode, section_size, section_seqs = sys.argv[8:12]

sections = {}
for key, val in (
    ("perf_full", section_full),
    ("perf_mode", section_mode),
    ("scale_size", section_size),
    ("scale_seqs_fixed", section_seqs),
):
    if val:
        sections[key] = val

tools = {}
for name in ("samtools", "seqkit", "seqtk", "noodles", "rustbio"):
    val = os.environ.get(f"BENCH_VER_{name.upper()}")
    if val:
        tools[name] = val

out = {
    "schema_version": "stats-run.v1",
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
    "skip_full": "perf_full" not in sections,
    "skip_mode": "perf_mode" not in sections,
    "skip_scale_size": "scale_size" not in sections,
    "skip_scale_seqs_fixed": "scale_seqs_fixed" not in sections,
    "metadata": metadata,
    "tools": tools,
    "sections": sections,
}

# Resume with --skip-verify must not erase a prior green verify for this timestamp.
manifest_path = Path(manifest)
if not do_verify:
    prior_pass = None
    if manifest_path.is_file():
        try:
            prior = json.loads(manifest_path.read_text())
            if prior.get("verify_pass") is not None and not prior.get("verify_skipped", True):
                prior_pass = str(prior["verify_pass"])
        except Exception:
            prior_pass = None
    if prior_pass is None:
        log = manifest_path.parent / f"verify_{ts}.log"
        if log.is_file():
            text = log.read_text(encoding="utf-8", errors="replace")
            if "ALL PASSED" in text:
                for line in text.splitlines():
                    if line.startswith("Results:"):
                        parts = line.split()
                        if len(parts) >= 2:
                            prior_pass = parts[1]
                        break
    if prior_pass is not None:
        out["verify_skipped"] = False
        out["verify_pass"] = prior_pass

Path(manifest).write_text(json.dumps(out, indent=2) + "\n")
PY
}

# Run one timed tool with live progress (zebrac groups are silent until done).
run_zebrac_tool() {
    local section="$1" workload="$2" tool="$3" family="$4"
    local json_out="$5" script="$6" input_bytes="$7"
    local t0 t1 elapsed
    echo "  >> $section $workload $tool  (runs=$RUNS warmup=$WARMUP)"
    t0="$(date +%s)"
    zebrac_clear_commands
    stats_add_command "$section" "$workload" "$tool" "$family" "$json_out" "$script" "$input_bytes"
    bench_group "$json_out"
    t1="$(date +%s)"
    elapsed=$((t1 - t0))
    echo "  << $section $workload $tool  ${elapsed}s"
}

# Stash .zfi so the loader picks .fai; restore even on failure (shared sidecars).
run_zebrac_tool_fai() {
    local fa="$1" section="$2" workload="$3" tool="$4" family="$5"
    local json_out="$6" script="$7" input_bytes="$8"
    if [[ ! -f "${fa}.fai" ]]; then
        echo "  warn: missing ${fa}.fai; skip $tool for $workload" >&2
        return 0
    fi
    mv -f "${fa}.zfi" "${fa}.zfi.stash" 2>/dev/null || true
    if ! (
        run_zebrac_tool "$section" "$workload" "$tool" "$family" "$json_out" "$script" "$input_bytes"
    ); then
        mv -f "${fa}.zfi.stash" "${fa}.zfi" 2>/dev/null || true
        echo "error: $tool lane failed for $section $workload" >&2
        exit 1
    fi
    mv -f "${fa}.zfi.stash" "${fa}.zfi" 2>/dev/null || true
}

# Timed peers for one FASTA. One zebrac group per tool so progress is visible.
# include_zfasta_modes: also time full(fai), indexed(zfi), indexed(fai).
run_stats_peers() {
    local section="$1" workload="$2" fa="$3" out_dir="$4"
    local include_zfasta_modes="${5:-false}"
    local allow_seqtk="${6:-true}"

    local qf qz qk qn qr qt nbytes json
    qf="$(quote_arg "$fa")"
    qz="$(quote_arg "$ZFASTA")"
    qk="$(quote_arg "$SEQKIT")"
    qn="$(quote_arg "$NOODLES")"
    qr="$(quote_arg "$RUSTBIO")"
    qt="$(quote_arg "$SEQTK")"
    nbytes="$(file_size_bytes "$fa")"

    json="$out_dir/${workload}__z-fasta-full.json"
    run_zebrac_tool "$section" "$workload" z-fasta-full z-fasta "$json" \
        "$qz stats $qf > /dev/null" "$nbytes"

    if [[ "$include_zfasta_modes" == "true" ]]; then
        json="$out_dir/${workload}__z-fasta-full-fai.json"
        run_zebrac_tool_fai "$fa" "$section" "$workload" z-fasta-full-fai z-fasta "$json" \
            "$qz stats $qf > /dev/null" "$nbytes"

        json="$out_dir/${workload}__z-fasta-indexed-zfi.json"
        run_zebrac_tool "$section" "$workload" z-fasta-indexed-zfi z-fasta "$json" \
            "$qz stats --index-only $qf > /dev/null" "$nbytes"

        json="$out_dir/${workload}__z-fasta-indexed-fai.json"
        run_zebrac_tool_fai "$fa" "$section" "$workload" z-fasta-indexed-fai z-fasta "$json" \
            "$qz stats --index-only $qf > /dev/null" "$nbytes"
    fi

    if bench_has_tool noodles; then
        json="$out_dir/${workload}__noodles.json"
        run_zebrac_tool "$section" "$workload" noodles noodles "$json" \
            "$qn stats $qf > /dev/null" "$nbytes"
    fi
    if bench_has_tool rustbio; then
        json="$out_dir/${workload}__rustbio.json"
        run_zebrac_tool "$section" "$workload" rustbio rustbio "$json" \
            "$qr stats $qf > /dev/null" "$nbytes"
    fi
    if bench_has_tool seqkit; then
        json="$out_dir/${workload}__seqkit.json"
        run_zebrac_tool "$section" "$workload" seqkit seqkit "$json" \
            "$qk stats -a -T $qf > /dev/null" "$nbytes"
    fi
    if [[ "$allow_seqtk" == "true" ]] && bench_has_tool seqtk; then
        json="$out_dir/${workload}__seqtk.json"
        run_zebrac_tool "$section" "$workload" seqtk seqtk "$json" \
            "$qt comp $qf > /dev/null" "$nbytes"
    fi
}



# ══════════════════════════════════════════════════════════════════════
#  Correctness: z-fasta stats vs oracles
# ══════════════════════════════════════════════════════════════════════

SKIP_TOOLS=false
SKIP_MESSY=false
SKIP_DEDUP=false
SKIP_EDGE=false
SKIP_LAYOUT=false

# Correctness paths (oracle.py stays external; fixtures under data/verify)
ORACLE="$SCRIPT_DIR/oracle.py"
TEST_DATA="$PROJECT_ROOT/tests/data"
MESSY_DIR="$BENCH_ROOT/shared/cache/messy_fixtures"
# Generated on demand under data/ (gitignored). Do not track FASTA here.
VERIFY_DATA="$SCRIPT_DIR/data/verify"
LAYOUT_TWINS="$VERIFY_DATA/layout_twins"
PYTHON="${PYTHON:-$PROJECT_ROOT/.venv/bin/python}"
[[ -x "$PYTHON" ]] || PYTHON="$(command -v python3)"



FIXTURES=(simple proteome single edge_cases mixed_widths)
# Side-table / whitespace messy: z-fasta .zfi only (peers count raw bytes).
MESSY_FIXTURES=(mixed_widths trailing_whitespace blank_lines mixed_crlf)
# Uniform messy control: samtools .fai works; exercise zfi+fai cross like get.
UNIFORM_MESSY=uniform
# seqkit counts every raw FASTA row; skip edge_cases where index filters empties/dups.
SEQKIT_FIXTURES=(simple proteome single mixed_widths)
# seqtk comp is nucleotide-only; proteome is protein.
SEQTK_FIXTURES=(simple single mixed_widths)
# noodles/rustbio wrappers: clean FASTA comparison peers only (no messy / side-table).
# Richer TSV fields exist so we can compare assembly+composition; they do not gain messy support.
WRAPPER_FIXTURES=(simple proteome single edge_cases mixed_widths)
# Same bases, different wrapping: uniform reference vs messy layouts (z-fasta only).
LAYOUT_TWIN_VARIANTS=(mixed_widths trailing_whitespace blank_lines mixed_crlf)

PASS=0 FAIL=0

pass() { PASS=$((PASS + 1)); printf "  PASS: %s\n" "$1"; }

fail() {
    FAIL=$((FAIL + 1)); printf "  FAIL: %s\n" "$1"
    if [[ -s "${2:-}" ]]; then sed 's/^/    /' "$2"; fi
}

section_hdr() { echo ""; echo "--- [$1] $2 ---"; }

oracle() { "$PYTHON" "$ORACLE" "$@"; }

check_oracle() {
    local label="$1"; shift
    local slug="${label//[^a-zA-Z0-9_]/_}"
    local err="$TMPDIR/oracle_${slug}.err"
    if oracle "$@" 2>"$err"; then
        pass "$label"
    else
        fail "$label" "$err"
    fi
}

expect_fail() {
    local label="$1"; shift
    local err="$TMPDIR/edge_${label//[^a-zA-Z0-9_]/_}.err"
    local out="$TMPDIR/edge_${label//[^a-zA-Z0-9_]/_}.out"
    if "$@" >"$out" 2>"$err"; then
        echo "expected non-zero exit" >"$err"
        fail "$label" "$err"
    else
        pass "$label"
    fi
}

expect_fail_msg() {
    local label="$1" needle="$2"; shift 2
    local err="$TMPDIR/edge_${label//[^a-zA-Z0-9_]/_}.err"
    local out="$TMPDIR/edge_${label//[^a-zA-Z0-9_]/_}.out"
    if "$@" >"$out" 2>"$err"; then
        echo "expected non-zero exit" >"$err"
        fail "$label" "$err"
        return
    fi
    if grep -q "$needle" "$err"; then
        pass "$label"
    else
        echo "stderr missing: $needle" >>"$err"
        fail "$label" "$err"
    fi
}

ensure_index() {
    local f="$1"
    if [[ ! -f "${f}.fai" ]]; then
        if ! "$SAMTOOLS" faidx "$f" 2>"$TMPDIR/faidx.err"; then
            echo "samtools faidx failed for ${f}" >>"$TMPDIR/faidx.err"
            return 1
        fi
    fi
    [[ -f "${f}.zfi" ]] || "$ZFASTA" index "$f" 2>/dev/null
}

prepare_messy() {
    cp "$MESSY_DIR/${1}.fasta" "$2"
    "$ZFASTA" index "$2" >/dev/null 2>&1
}

run_stats() {
    local fasta="$1" mode="$2" out="$3" err="$4"
    local -a args=(stats)
    [[ "$mode" == index-only ]] && args+=(--index-only)
    args+=("$fasta")
    if ! "$ZFASTA" "${args[@]}" >"$out" 2>"$err"; then
        echo "z-fasta stats failed (see stderr)" >"$err"
        return 1
    fi
}

# Same bases, different wrapping. Written under data/verify/layout_twins/ (not tracked).
generate_layout_twins() {
    mkdir -p "$LAYOUT_TWINS"
    "$PYTHON" - "$LAYOUT_TWINS" <<'PY'
from pathlib import Path
import sys

out = Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=True)

(out / "uniform.fasta").write_bytes(
    b">layout_twin uniform wrap\n"
    b"AAAACCCCGGGG\n"
    b"ttttAAAANCCC\n"
    b"G*GGGTTTT\n"
)
(out / "mixed_widths.fasta").write_bytes(
    b">layout_twin mixed widths\n"
    b"AAAACCCCGGGG\n"
    b"ttttAA\n"
    b"AANCCCG*GG\n"
    b"GTTTT\n"
)
(out / "trailing_whitespace.fasta").write_bytes(
    b">layout_twin trailing ws\n"
    b"AAAACCCCGGGG    \n"
    b"ttttAAAANCCC\t\n"
    b"G*GGGTTTT  \n"
)
(out / "blank_lines.fasta").write_bytes(
    b">layout_twin blank lines\n"
    b"AAAACCCCGGGG\n"
    b"\n"
    b"ttttAAAANCCC\n"
    b"\n"
    b"G*GGGTTTT\n"
)
(out / "mixed_crlf.fasta").write_bytes(
    b">layout_twin mixed crlf\r\n"
    b"AAAACCCCGGGG\r\n"
    b"ttttAA\n"
    b"AANCCCG*GGGTTTT\r\n"
)
PY
}

prepare_fixtures() {
    local name fasta
    for name in "${FIXTURES[@]}"; do
        fasta="$TMPDIR/$name.fasta"
        cp "$TEST_DATA/$name.fasta" "$fasta"
        [[ -f "$TEST_DATA/$name.fasta.fai" ]] && cp "$TEST_DATA/$name.fasta.fai" "${fasta}.fai"
        ensure_index "$fasta" || { fail "[oracle] index setup failed for $name" "$TMPDIR/faidx.err"; exit 1; }
        oracle expected "$fasta" >"$TMPDIR/$name.expected.json"
    done
}

verify_modes() {
    local tag="$1" name="$2" fasta="$3" exp="$4" zfi_only="${5:-}"
    local mode out err stash label zfi_txt fai_txt

    section_hdr "$tag" "$name"

    for mode in full index-only; do
        zfi_txt="$TMPDIR/$name.$mode.zfi.txt"
        err="$TMPDIR/$name.$mode.zfi.err"
        if [[ "$zfi_only" == zfi-only ]]; then
            label="[oracle:$name:$mode:zfi] vs BioPython"
            run_stats "$fasta" "$mode" "$zfi_txt" "$err" || { fail "$label" "$err"; continue; }
            check_oracle "$label" check "$mode" zfi "$fasta" "$exp" "$zfi_txt"
            continue
        fi

        # .zfi path vs BioPython
        label="[oracle:$name:$mode:zfi] vs BioPython"
        run_stats "$fasta" "$mode" "$zfi_txt" "$err" || { fail "$label" "$err"; continue; }
        check_oracle "$label" check "$mode" zfi "$fasta" "$exp" "$zfi_txt"

        # .fai selection vs BioPython (stash .zfi so it is absent)
        stash="$TMPDIR/$name.$mode.zfi"
        mv "${fasta}.zfi" "$stash"
        fai_txt="$TMPDIR/$name.$mode.fai.txt"
        err="$TMPDIR/$name.$mode.fai.err"
        label="[oracle:$name:$mode:fai] vs BioPython"
        if ! run_stats "$fasta" "$mode" "$fai_txt" "$err"; then
            fail "$label" "$err"
            mv "$stash" "${fasta}.zfi"
            continue
        fi
        check_oracle "$label" check "$mode" fai "$fasta" "$exp" "$fai_txt"
        mv "$stash" "${fasta}.zfi"

        # Identity: .zfi stats == .fai stats.
        check_oracle "[index:cross] $name $mode .zfi == .fai" same "$mode" "$zfi_txt" "$fai_txt"
    done
}

# One pass per peer tool so failures name the tool (get-style).
# Wrappers: clean fixtures only (WRAPPER_FIXTURES). Never run on messy side-table files.
verify_parity() {
    local name="$1" fasta="$2" exp="$3"
    local idx_txt="$TMPDIR/$name.index-only.zfi.txt"
    local full_txt="$TMPDIR/$name.full.zfi.txt"
    local tool bin out err label
    local ran=0
    local use_wrappers=false
    local sk

    if [[ ! -f "$idx_txt" || ! -f "$full_txt" ]]; then
        err="$TMPDIR/$name.parity_setup.err"
        echo "missing z-fasta stats outputs for $name (run oracle modes first)" >"$err"
        fail "[parity:$name] setup" "$err"
        return
    fi

    for sk in "${WRAPPER_FIXTURES[@]}"; do
        [[ "$sk" == "$name" ]] && use_wrappers=true && break
    done

    if $use_wrappers; then
        for pair in "noodles:$NOODLES" "rustbio:$RUSTBIO"; do
            tool="${pair%%:*}"
            bin="${pair#*:}"
            [[ -x "$bin" ]] || continue
            label="[parity:$tool] $name assembly/composition"
            out="$TMPDIR/$name.$tool.stats"
            err="$TMPDIR/$name.$tool.err"
            if ! "$bin" stats "$fasta" >"$out" 2>"$err"; then
                fail "$label (run)" "$err"
                continue
            fi
            # Stale binary: old wrappers only printed sequences + total_bases.
            if ! grep -q $'^n50\t' "$out" || ! grep -q $'^type\t' "$out"; then
                echo "wrapper output missing n50/type; rebuild tools/${tool}_wrapper (clean-FASTA comparison peer, not messy)" >"$err"
                fail "$label (stale binary)" "$err"
                continue
            fi
            check_oracle "$label" parity "$fasta" "$exp" "$idx_txt" "$full_txt" "$tool" "$out"
            ran=$((ran + 1))
        done
    fi

    for sk in "${SEQKIT_FIXTURES[@]}"; do
        [[ "$sk" == "$name" ]] || continue
        [[ -x "$SEQKIT" ]] || continue
        label="[parity:seqkit] $name assembly stats"
        out="$TMPDIR/$name.seqkit.txt"
        err="$TMPDIR/$name.seqkit.err"
        if ! "$SEQKIT" stats -a -T "$fasta" >"$out" 2>"$err"; then
            fail "$label (run)" "$err"
            break
        fi
        check_oracle "$label" parity "$fasta" "$exp" "$idx_txt" "$full_txt" seqkit "$out"
        ran=$((ran + 1))
        break
    done

    for sk in "${SEQTK_FIXTURES[@]}"; do
        [[ "$sk" == "$name" ]] || continue
        [[ -x "$SEQTK" ]] || continue
        label="[parity:seqtk] $name composition"
        out="$TMPDIR/$name.seqtk.txt"
        err="$TMPDIR/$name.seqtk.err"
        if ! "$SEQTK" comp "$fasta" >"$out" 2>"$err"; then
            fail "$label (run)" "$err"
            break
        fi
        check_oracle "$label" parity "$fasta" "$exp" "$idx_txt" "$full_txt" seqtk "$out"
        ran=$((ran + 1))
        break
    done

    if [[ "$ran" -eq 0 ]]; then
        err="$TMPDIR/$name.parity.err"
        echo "no external tools available for $name" >"$err"
        fail "[parity:$name] external tools" "$err"
    fi
}

verify_file() {
    local name="$1"
    verify_modes oracle "$name" "$TMPDIR/$name.fasta" "$TMPDIR/$name.expected.json"
    if ! $SKIP_TOOLS; then
        verify_parity "$name" "$TMPDIR/$name.fasta" "$TMPDIR/$name.expected.json"
    fi
}

verify_messy() {
    local name fasta exp
    for name in "${MESSY_FIXTURES[@]}"; do
        fasta="$TMPDIR/messy_${name}.fasta"
        [[ -f "$MESSY_DIR/${name}.fasta" ]] || { fail "[extended:messy] missing fixture $name"; continue; }
        prepare_messy "$name" "$fasta" || { fail "[extended:messy] index failed $name"; continue; }
        oracle expected "$fasta" >"$TMPDIR/messy_${name}.expected.json"
        verify_modes "extended:messy" "messy_$name" "$fasta" "$TMPDIR/messy_${name}.expected.json" zfi-only
    done

    # Uniform control: peers can index; exercise .zfi and .fai stats cross-compare.
    name="$UNIFORM_MESSY"
    fasta="$TMPDIR/messy_${name}.fasta"
    [[ -f "$MESSY_DIR/${name}.fasta" ]] || { fail "[extended:messy] missing fixture $name"; return; }
    cp "$MESSY_DIR/${name}.fasta" "$fasta"
    ensure_index "$fasta" || { fail "[extended:messy] index failed $name" "$TMPDIR/faidx.err"; return; }
    oracle expected "$fasta" >"$TMPDIR/messy_${name}.expected.json"
    verify_modes "extended:messy" "messy_$name" "$fasta" "$TMPDIR/messy_${name}.expected.json"
}

verify_dedup_stats() {
    local dedup_fasta="$TMPDIR/dedup_default.fasta"
    local nodedup_fasta="$TMPDIR/dedup_nodedup.fasta"
    local exp_dedup="$TMPDIR/dedup_default.expected.json"
    local exp_nodedup="$TMPDIR/nodedup.expected.json"
    local label="[extended:dedup] --no-dedup indexes more sequences than default"
    local err dedup_seqs nodedup_seqs

    cp "$TEST_DATA/edge_cases.fasta" "$dedup_fasta"
    rm -f "${dedup_fasta}.zfi" "${dedup_fasta}.fai"
    ensure_index "$dedup_fasta" || { fail "[extended:dedup] default index setup" "$TMPDIR/faidx.err"; return; }
    oracle expected "$dedup_fasta" >"$exp_dedup"
    err="$TMPDIR/dedup_default.err"
    run_stats "$dedup_fasta" index-only "$TMPDIR/dedup_default.index-only.txt" "$err" \
        || { fail "[extended:dedup] default index-only stats" "$err"; return; }
    check_oracle "[extended:dedup] default vs BioPython" check index-only zfi \
        "$dedup_fasta" "$exp_dedup" "$TMPDIR/dedup_default.index-only.txt"

    cp "$TEST_DATA/edge_cases.fasta" "$nodedup_fasta"
    rm -f "${nodedup_fasta}.zfi" "${nodedup_fasta}.fai"
    "$ZFASTA" index --no-dedup "$nodedup_fasta" >/dev/null 2>&1 \
        || { fail "[extended:dedup] index --no-dedup"; return; }
    oracle expected "$nodedup_fasta" --no-dedup >"$exp_nodedup"
    verify_modes "extended:dedup" "edge_cases_nodedup" "$nodedup_fasta" "$exp_nodedup" zfi-only

    dedup_seqs="$("$PYTHON" -c "import json; print(json.load(open('$exp_dedup'))['num_seqs'])")"
    nodedup_seqs="$("$PYTHON" -c "import json; print(json.load(open('$exp_nodedup'))['num_seqs'])")"
    if [[ "$nodedup_seqs" -gt "$dedup_seqs" ]]; then
        pass "$label"
    else
        echo "dedup=$dedup_seqs nodedup=$nodedup_seqs" >"$TMPDIR/dedup_cmp.err"
        fail "$label" "$TMPDIR/dedup_cmp.err"
    fi
}

# Same sequence bytes, different on-disk wrapping: stats must match (layout invariant).
verify_layout_twins() {
    section_hdr "extended:layout" "messy == uniform (same bases)"
    local uniform="$TMPDIR/layout_uniform.fasta"
    local variant name fasta err mode
    local uref_full="$TMPDIR/layout_uniform.full.txt"
    local uref_idx="$TMPDIR/layout_uniform.index-only.txt"

    generate_layout_twins
    [[ -f "$LAYOUT_TWINS/uniform.fasta" ]] || {
        fail "[extended:layout] missing fixtures under bench/stats/data/verify/layout_twins"
        return
    }

    cp "$LAYOUT_TWINS/uniform.fasta" "$uniform"
    rm -f "${uniform}.zfi" "${uniform}.fai"
    "$ZFASTA" index "$uniform" >/dev/null 2>&1 \
        || { fail "[extended:layout] index uniform"; return; }
    err="$TMPDIR/layout_uniform.err"
    run_stats "$uniform" full "$uref_full" "$err" \
        || { fail "[extended:layout] uniform full stats" "$err"; return; }
    run_stats "$uniform" index-only "$uref_idx" "$err" \
        || { fail "[extended:layout] uniform index-only stats" "$err"; return; }

    oracle expected "$uniform" >"$TMPDIR/layout_uniform.expected.json"
    check_oracle "[extended:layout] uniform vs BioPython" check full zfi \
        "$uniform" "$TMPDIR/layout_uniform.expected.json" "$uref_full"

    for variant in "${LAYOUT_TWIN_VARIANTS[@]}"; do
        name="layout_$variant"
        fasta="$TMPDIR/${name}.fasta"
        [[ -f "$LAYOUT_TWINS/${variant}.fasta" ]] || {
            fail "[extended:layout] missing $variant fixture"
            continue
        }
        cp "$LAYOUT_TWINS/${variant}.fasta" "$fasta"
        rm -f "${fasta}.zfi" "${fasta}.fai"
        "$ZFASTA" index "$fasta" >/dev/null 2>&1 \
            || { fail "[extended:layout] index $variant"; continue; }

        for mode in full index-only; do
            err="$TMPDIR/${name}.${mode}.err"
            run_stats "$fasta" "$mode" "$TMPDIR/${name}.${mode}.txt" "$err" \
                || { fail "[extended:layout] $variant $mode stats" "$err"; continue; }
            if [[ "$mode" == full ]]; then
                check_oracle "[extended:layout] $variant full == uniform" same full \
                    "$uref_full" "$TMPDIR/${name}.full.txt"
            else
                check_oracle "[extended:layout] $variant index-only == uniform" same index-only \
                    "$uref_idx" "$TMPDIR/${name}.index-only.txt"
            fi
        done

        oracle expected "$fasta" >"$TMPDIR/${name}.expected.json"
        check_oracle "[extended:layout] $variant vs BioPython" check full zfi \
            "$fasta" "$TMPDIR/${name}.expected.json" "$TMPDIR/${name}.full.txt"
    done
}

# Product policy: full stats report source-level extras; index-only never fabricates 0
# on a deduplicated index (n/a). --no-dedup index-only can report retained repeats.
verify_duplicates_policy() {
    section_hdr "extended:dedup" "Duplicates line policy"
    local fasta="$TMPDIR/dup_policy.fasta"
    local out err dups
    err="$TMPDIR/dup_policy.err"

    cp "$TEST_DATA/edge_cases.fasta" "$fasta"
    rm -f "${fasta}.zfi" "${fasta}.fai"
    "$ZFASTA" index "$fasta" >/dev/null 2>&1 \
        || { fail "[extended:dedup] default index for Duplicates policy"; return; }

    out="$TMPDIR/dup_policy.full.txt"
    run_stats "$fasta" full "$out" "$err" \
        || { fail "[extended:dedup] full stats on default index" "$err"; return; }
    dups="$(awk '/^Duplicates:/{print $2; exit}' "$out")"
    if [[ "$dups" == "1" ]]; then
        pass "[extended:dedup] full Duplicates=1 for edge_cases (source extras)"
    else
        echo "expected Duplicates: 1, got ${dups:-missing}" >"$err"
        fail "[extended:dedup] full Duplicates=1 for edge_cases" "$err"
    fi

    out="$TMPDIR/dup_policy.index-only.txt"
    run_stats "$fasta" index-only "$out" "$err" \
        || { fail "[extended:dedup] index-only on default index" "$err"; return; }
    if grep -q '^Duplicates:[[:space:]]*n/a' "$out"; then
        pass "[extended:dedup] index-only Duplicates=n/a on deduped index"
    else
        echo "expected Duplicates: n/a ..., got:" >"$err"
        grep '^Duplicates:' "$out" >>"$err" || true
        fail "[extended:dedup] index-only Duplicates=n/a on deduped index" "$err"
    fi

    rm -f "${fasta}.zfi" "${fasta}.fai"
    "$ZFASTA" index --no-dedup "$fasta" >/dev/null 2>&1 \
        || { fail "[extended:dedup] index --no-dedup for Duplicates policy"; return; }
    out="$TMPDIR/dup_policy.nodedup.index-only.txt"
    run_stats "$fasta" index-only "$out" "$err" \
        || { fail "[extended:dedup] index-only on --no-dedup index" "$err"; return; }
    dups="$(awk '/^Duplicates:/{print $2; exit}' "$out")"
    if [[ "$dups" == "1" ]]; then
        pass "[extended:dedup] index-only Duplicates=1 on --no-dedup index"
    else
        echo "expected Duplicates: 1, got ${dups:-missing}" >"$err"
        fail "[extended:dedup] index-only Duplicates=1 on --no-dedup index" "$err"
    fi

    # --no-dedup with unique names still cannot invent 0 under index-only.
    local uniq="$TMPDIR/dup_policy_uniq.fasta"
    cp "$TEST_DATA/simple.fasta" "$uniq"
    rm -f "${uniq}.zfi" "${uniq}.fai"
    "$ZFASTA" index --no-dedup "$uniq" >/dev/null 2>&1 \
        || { fail "[extended:dedup] index --no-dedup unique names"; return; }
    out="$TMPDIR/dup_policy_uniq.index-only.txt"
    run_stats "$uniq" index-only "$out" "$err" \
        || { fail "[extended:dedup] index-only on --no-dedup unique" "$err"; return; }
    if grep -q '^Duplicates:[[:space:]]*n/a' "$out"; then
        pass "[extended:dedup] index-only Duplicates=n/a on --no-dedup unique names"
    else
        echo "expected Duplicates: n/a ..., got:" >"$err"
        grep '^Duplicates:' "$out" >>"$err" || true
        fail "[extended:dedup] index-only Duplicates=n/a on --no-dedup unique names" "$err"
    fi
}

verify_edge_paths() {
    section_hdr "extended:edge" "CLI and error paths"
    local bare="$TMPDIR/edge_noindex.fasta"
    local empty="$TMPDIR/edge_empty.fasta"
    local junk="$TMPDIR/edge_not_fasta.txt"

    expect_fail_msg "[extended:edge] usage (no file)" "usage:" \
        "$ZFASTA" stats

    expect_fail_msg "[extended:edge] missing file" "file not found" \
        "$ZFASTA" stats "$TMPDIR/does_not_exist.fasta"

    cp "$TEST_DATA/simple.fasta" "$bare"
    rm -f "${bare}.zfi" "${bare}.fai"
    expect_fail_msg "[extended:edge] no index" "no index found" \
        "$ZFASTA" stats "$bare"

    cp "$TEST_DATA/empty.fasta" "$empty"
    rm -f "${empty}.zfi" "${empty}.fai"
    expect_fail_msg "[extended:edge] empty fasta" "file is empty" \
        "$ZFASTA" stats "$empty"

    printf 'this is not fasta\n' >"$junk"
    rm -f "${junk}.zfi" "${junk}.fai"
    expect_fail "[extended:edge] not_fasta (no index)" \
        "$ZFASTA" stats "$junk"

    # Indexed simple: --index-only must omit Composition; full must include it.
    local ok="$TMPDIR/edge_ok.fasta"
    cp "$TEST_DATA/simple.fasta" "$ok"
    ensure_index "$ok" || { fail "[extended:edge] index setup" "$TMPDIR/faidx.err"; return; }
    local out="$TMPDIR/edge_ok.index-only.txt" err="$TMPDIR/edge_ok.err"
    run_stats "$ok" index-only "$out" "$err" || { fail "[extended:edge] index-only run" "$err"; return; }
    if grep -q "Composition:" "$out"; then
        echo "Composition present under --index-only" >"$err"
        fail "[extended:edge] index-only omits Composition" "$err"
    else
        pass "[extended:edge] index-only omits Composition"
    fi
    run_stats "$ok" full "$TMPDIR/edge_ok.full.txt" "$err" || { fail "[extended:edge] full run" "$err"; return; }
    if grep -q "Composition:" "$TMPDIR/edge_ok.full.txt"; then
        pass "[extended:edge] full includes Composition"
    else
        echo "Composition missing" >"$err"
        fail "[extended:edge] full includes Composition" "$err"
    fi
}

run_tests() {
    PASS=0
    FAIL=0
    bench_ensure_messy --fixtures
    mkdir -p "$LAYOUT_TWINS"
    local verify_tmp="$VERIFY_DATA/work"
    rm -rf "$verify_tmp"
    mkdir -p "$verify_tmp"
    TMPDIR="$verify_tmp"


    echo "z-fasta stats verification"
    echo "  z-fasta:  $ZFASTA"
    echo "  skip tools: $SKIP_TOOLS"
    echo "  skip messy: $SKIP_MESSY"
    echo "  skip dedup: $SKIP_DEDUP"
    echo "  skip edge: $SKIP_EDGE"
    echo "  skip layout: $SKIP_LAYOUT"

    cd "$PROJECT_ROOT"
    [[ -x "$ZFASTA" ]] || { echo "Error: run ./zig build first"; rm -rf "$TMPDIR"; return 1; }
    command -v "$SAMTOOLS" &>/dev/null || { echo "Error: samtools not found"; rm -rf "$TMPDIR"; return 1; }
    "$PYTHON" -c "from Bio import SeqIO" 2>/dev/null || { echo "Error: BioPython not available"; rm -rf "$TMPDIR"; return 1; }

    prepare_fixtures
    for name in "${FIXTURES[@]}"; do
        verify_file "$name"
    done

    if ! $SKIP_MESSY; then
        verify_messy
    fi

    if ! $SKIP_LAYOUT; then
        verify_layout_twins
    fi

    if ! $SKIP_DEDUP; then
        verify_dedup_stats
        verify_duplicates_policy
    fi

    if ! $SKIP_EDGE; then
        verify_edge_paths
    fi

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

run_perf_full() {
    local out_dir="$RESULTS_DIR/perf_full_${TIMESTAMP}"
    mkdir -p "$out_dir"
    SECTION_FULL="perf_full_${TIMESTAMP}"

    echo "--------------------------------------------------"
    echo " Full stats (peers + z-fasta full)"
    echo "--------------------------------------------------"
    echo "  Note: each tool is timed alone. On Genome (~3 GB) peers re-parse"
    echo "  the whole file every sample (warmup+runs). Expect long walls."
    echo "  z-fasta full (fai) / indexed lanes live in perf_mode + scaling."

    local ds fa allow_seqtk
    for ds in Genome Transcriptome Proteome; do
        fa="${REAL_DATASETS[$ds]:-}"
        [[ -n "$fa" ]] || continue
        allow_seqtk=true
        [[ "$ds" == "Proteome" ]] && allow_seqtk=false
        echo "  -- dataset $ds ($(du -h "$fa" | cut -f1)) --"
        run_stats_peers perf_full "$ds" "$fa" "$out_dir" false "$allow_seqtk"
        echo "  done perf_full $ds"
    done
}

run_perf_mode() {
    local out_dir="$RESULTS_DIR/perf_mode_${TIMESTAMP}"
    mkdir -p "$out_dir"
    SECTION_MODE="perf_mode_${TIMESTAMP}"

    echo "--------------------------------------------------"
    echo " z-fasta modes (full zfi/fai + indexed zfi/fai)"
    echo "--------------------------------------------------"
    echo "  Policy: .zfi preferred; .fai is compatibility for uniform records."
    echo "  Same stats surface for both (messy side tables need .zfi)."

    local ds fa qf qz nbytes json
    qz="$(quote_arg "$ZFASTA")"
    for ds in Genome Transcriptome Proteome; do
        fa="${REAL_DATASETS[$ds]:-}"
        [[ -n "$fa" ]] || continue
        qf="$(quote_arg "$fa")"
        nbytes="$(file_size_bytes "$fa")"
        echo "  -- dataset $ds --"

        json="$out_dir/${ds}__z-fasta-full.json"
        run_zebrac_tool perf_mode "$ds" z-fasta-full z-fasta "$json" \
            "$qz stats $qf > /dev/null" "$nbytes"

        json="$out_dir/${ds}__z-fasta-full-fai.json"
        run_zebrac_tool_fai "$fa" perf_mode "$ds" z-fasta-full-fai z-fasta "$json" \
            "$qz stats $qf > /dev/null" "$nbytes"

        json="$out_dir/${ds}__z-fasta-indexed-zfi.json"
        run_zebrac_tool perf_mode "$ds" z-fasta-indexed-zfi z-fasta "$json" \
            "$qz stats --index-only $qf > /dev/null" "$nbytes"

        json="$out_dir/${ds}__z-fasta-indexed-fai.json"
        run_zebrac_tool_fai "$fa" perf_mode "$ds" z-fasta-indexed-fai z-fasta "$json" \
            "$qz stats --index-only $qf > /dev/null" "$nbytes"

        echo "  done perf_mode $ds"
    done
}

run_perf_scale_size() {
    local out_dir="$RESULTS_DIR/scale_size_${TIMESTAMP}"
    mkdir -p "$out_dir"
    SECTION_SCALE_SIZE="scale_size_${TIMESTAMP}"

    echo "--------------------------------------------------"
    echo " File-size scaling"
    echo "--------------------------------------------------"

    local mb fa
    for mb in "${SIZE_MBS[@]}"; do
        fa="$SCALING_DIR/size_${mb}mb.fasta"
        [[ -f "$fa" ]] || { echo "error: missing $fa" >&2; exit 1; }
        echo "  -- size ${mb}mb ($(du -h "$fa" | cut -f1)) --"
        run_stats_peers scale_size "${mb}mb" "$fa" "$out_dir" true true
        echo "  done scale_size ${mb}mb"
    done
}

run_perf_scale_seqs() {
    local out_dir="$RESULTS_DIR/scale_seqs_fixed_${TIMESTAMP}"
    mkdir -p "$out_dir"
    SECTION_SCALE_SEQS="scale_seqs_fixed_${TIMESTAMP}"

    echo "--------------------------------------------------"
    echo " Sequence scaling (fixed 1024 bp)"
    echo "--------------------------------------------------"

    local count fa
    for count in "${SEQ_FIXED_COUNTS[@]}"; do
        fa="$SCALING_DIR/seqs_fixed_${count}.fasta"
        [[ -f "$fa" ]] || { echo "error: missing $fa" >&2; exit 1; }
        echo "  -- seqs $count ($(du -h "$fa" | cut -f1)) --"
        run_stats_peers scale_seqs_fixed "$count" "$fa" "$out_dir" true true
        echo "  done scale_seqs_fixed $count"
    done
}

run_perf() {
    bench_require_tool zebrac
    bench_require_tool z-fasta
    METADATA_JSONL="$RESULTS_DIR/metadata_${TIMESTAMP}.jsonl"
    RUN_MANIFEST="$RESULTS_DIR/run_${TIMESTAMP}.json"
    if [[ -n "${STATS_RUN_TIMESTAMP:-}" ]] && [[ -f "$METADATA_JSONL" ]]; then
        # Drop sections we are about to re-run so resume does not duplicate metadata rows.
        python3 - "$METADATA_JSONL" "$DO_FULL" "$DO_MODE" "$DO_SCALE" "$DO_SCALE_SIZE" "$DO_SCALE_SEQS" <<'PY'
import json, sys
from pathlib import Path
path = Path(sys.argv[1])
do_full, do_mode, do_scale, do_size, do_seqs = (a == "true" for a in sys.argv[2:7])
drop = set()
if do_full:
    drop.add("perf_full")
if do_mode:
    drop.add("perf_mode")
if do_scale and do_size:
    drop.add("scale_size")
if do_scale and do_seqs:
    drop.add("scale_seqs_fixed")
rows = [json.loads(l) for l in path.read_text().splitlines() if l.strip()]
kept = [r for r in rows if r.get("section") not in drop]
path.write_text("".join(json.dumps(r, separators=(",", ":")) + "\n" for r in kept))
PY
    else
        : > "$METADATA_JSONL"
    fi
    ZEBRAC_MIN_SAMPLES="$RUNS"
    ZEBRAC_MAX_SAMPLES="$RUNS"
    ZEBRAC_WARMUP="$WARMUP"

    echo ""
    echo "================================================================"
    echo "  Stats zebrac performance"
    echo "================================================================"
    echo "  Runs: $RUNS | Warmup: $WARMUP | Duration: ${ZEBRAC_DURATION_MS}ms"
    echo "  metadata: $METADATA_JSONL"
    echo ""

    ensure_real_data
    preload_real_indexes

    $DO_FULL && run_perf_full
    $DO_MODE && run_perf_mode

    if $DO_SCALE; then
        ensure_scaling_fixtures
        preload_scaling_indexes
        $DO_SCALE_SIZE && run_perf_scale_size
        $DO_SCALE_SEQS && run_perf_scale_seqs
    fi

    : "${SECTION_FULL:=$(existing_section_dir perf_full)}"
    : "${SECTION_MODE:=$(existing_section_dir perf_mode)}"
    : "${SECTION_SCALE_SIZE:=$(existing_section_dir scale_size)}"
    : "${SECTION_SCALE_SEQS:=$(existing_section_dir scale_seqs_fixed)}"

    export_manifest_tool_versions
    write_run_manifest "$RUN_MANIFEST" "$TIMESTAMP" "$(basename "$METADATA_JSONL")"
    printf '%s\n' "$TIMESTAMP" > "$RESULTS_DIR/LATEST"
    echo "  manifest: $RUN_MANIFEST"
}

run_report() {
    if [[ ! -f "$SCRIPT_DIR/generate_report.py" ]]; then
        echo ""
        echo "================================================================"
        echo "  Report skipped (generate_report.py not present yet)"
        echo "================================================================"
    else
        local py; py="$(report_python)"
        local args=()
        $ALLOW_INCOMPLETE && args+=(--allow-incomplete)
        echo ""
        echo "================================================================"
        echo "  Generating REPORT.md"
        echo "================================================================"
        "$py" "$SCRIPT_DIR/generate_report.py" "${args[@]}"
    fi
}

TIMESTAMP="${STATS_RUN_TIMESTAMP:-$(date +%Y%m%d_%H%M%S)}"
mkdir -p "$RESULTS_DIR"

if $DO_TESTS; then
    echo ""
    echo "================================================================"
    echo "  run_tests (correctness)"
    echo "================================================================"
    _verify_log="$RESULTS_DIR/verify_${TIMESTAMP}.log"
    mkdir -p "$RESULTS_DIR"
    # pipefail: tee must not mask run_tests failure
    if run_tests | tee "$_verify_log"; then
        _verify_n="$(awk '/^Results:/{print $2; exit}' "$_verify_log")"
        export BENCH_VERIFY_PASS="${_verify_n:-unknown}"
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
echo "  Stats suite complete"
echo "================================================================"
