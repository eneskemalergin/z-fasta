#!/usr/bin/env bash
# Stats benchmark runner: verify (L2) then zebrac perf (L3).
#
# Usage:
#   bash bench/stats/run.sh [options]
#
# Defaults: verify first (92 checks), then perf (--runs 5, --warmup 1, --duration 5000).
# Report generation: python3 bench/stats/generate_report.py
#   (run.sh calls this after perf; --allow-incomplete for drafts)
#
#   bash bench/stats/run.sh
#   bash bench/stats/run.sh --skip-verify --skip-report   # perf only
#   bash bench/stats/run.sh --skip-full --skip-mode       # scaling only
#   bash bench/stats/run.sh --skip-scale                  # REAL full + mode only
#   bash bench/stats/run.sh --skip-size                   # skip file-size sweep only
#   bash bench/stats/run.sh --skip-seqs                   # skip seq-count sweep only
#   bash bench/stats/run.sh --regenerate-fixtures
#   bash bench/stats/run.sh --allow-incomplete
#   STATS_RUN_TIMESTAMP=<ts> bash bench/stats/run.sh --skip-verify ...  # resume same run id
#
# Scaling fixtures (generated; bench/stats/data/ is gitignored):
#   size_{N}mb.fasta         1-1000 MB, 100 sequences each
#   seqs_fixed_{N}.fasta     1024 bp/seq, N in 100k 250k 500k 1M
#   Created on demand. --regenerate-fixtures overwrites all.
#   Prefer hardlink from bench/index/data/ when present (FASTA + .fai).
#   Indexes (.zfi + .fai) are preloaded once per file; timed stats only loads them.
#
# Not in this runner (deferred / L2-only):
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
SCALING_DIR="$SCRIPT_DIR/data"
DATA_DIR="$BENCH_ROOT/shared/data"

source "$BENCH_ROOT/shared/zebrac_runner.sh"

# Scaling sweep constants (match index bench shape)
SIZE_MBS=(1 5 10 25 50 100 250 500 1000)
SEQ_FIXED_COUNTS=(100000 250000 500000 1000000)

# Defaults
RUNS=5
WARMUP=1
ZEBRAC_DURATION_MS="${ZEBRAC_DURATION_MS:-5000}"
DO_VERIFY=true
DO_PERF=true
DO_FULL=true
DO_MODE=true
DO_SCALE=true
DO_SCALE_SIZE=true
DO_SCALE_SEQS=true
DO_REPORT=true
REGENERATE_FIXTURES=false
ALLOW_INCOMPLETE=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --runs) RUNS="$2"; shift 2 ;;
        --warmup) WARMUP="$2"; shift 2 ;;
        --duration) ZEBRAC_DURATION_MS="$2"; shift 2 ;;
        --skip-verify) DO_VERIFY=false; shift ;;
        --skip-perf) DO_PERF=false; shift ;;
        --skip-full) DO_FULL=false; shift ;;
        --skip-mode) DO_MODE=false; shift ;;
        --skip-scale) DO_SCALE=false; shift ;;
        --skip-size) DO_SCALE_SIZE=false; shift ;;
        --skip-seqs) DO_SCALE_SEQS=false; shift ;;
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

stats_add_command() {
    local section="$1" workload="$2" tool="$3" family="$4"
    local json_out="$5" script="$6" input_bytes="$7"
    zebrac_add_command "stats" "$section" "$workload" "$tool" "$family" \
        "$input_bytes" "" "$json_out" "$(shell_command "$script")"
}

bench_group() {
    local json_out="$1"
    zebrac_run_current_group "$json_out" "$METADATA_JSONL"
    zebrac_clear_commands
}

existing_section_dir() {
    local prefix="$1"
    local dir="$RESULTS_DIR/${prefix}_${TIMESTAMP}"
    if [[ -d "$dir" ]] && compgen -G "${dir}/*.json" > /dev/null; then
        echo "${prefix}_${TIMESTAMP}"
    fi
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

# Rebuild when FASTA is newer, sidecar is missing, z-fasta binary is newer than .zfi,
# or .zfi predates the embedded name-table footer (ZFNM). Same policy as GET.
preload_indexes_for_file() {
    local fa="$1"
    [[ -f "$fa" ]] || return 1
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
}

preload_real_indexes() {
    bench_require_tool z-fasta
    echo "  Preloading REAL_* indexes (.zfi + .fai)..."
    local fa
    for fa in "${REAL_DATASETS[@]}"; do
        preload_indexes_for_file "$fa"
    done
}

generate_scaling_fixtures() {
    local mode="${1:-all}"  # all | size | seq
    mkdir -p "$SCALING_DIR"
    python3 - "$SCALING_DIR" "$mode" "$REGENERATE_FIXTURES" <<'PY'
import sys
from pathlib import Path

data_dir = Path(sys.argv[1])
mode = sys.argv[2]
regenerate = sys.argv[3] == "true"
LINE = "ACGTACGTAC" * 8
MIN_SEQ_LEN = 80
FIXED_SEQ_LEN = 1024
SIZE_MBS = (1, 5, 10, 25, 50, 100, 250, 500, 1000)
FIXED_COUNTS = (100_000, 250_000, 500_000, 1_000_000)
SIZE_SEQ_COUNT = 100


def write_wrapped(handle, seq_len):
    remaining = seq_len
    while remaining > 0:
        chunk = min(80, remaining)
        handle.write(LINE[:chunk] + "\n")
        remaining -= chunk


def write_fixed(count, path):
    with path.open("w", encoding="ascii") as handle:
        for i in range(1, count + 1):
            handle.write(f">seq{i}\n")
            write_wrapped(handle, FIXED_SEQ_LEN)


def write_size(mb, path):
    total = mb * 1024 * 1024
    with path.open("w", encoding="ascii") as handle:
        for i in range(1, SIZE_SEQ_COUNT + 1):
            handle.write(f">seq{i}\n")
            write_wrapped(handle, max(MIN_SEQ_LEN, total // SIZE_SEQ_COUNT - 50))


data_dir.mkdir(parents=True, exist_ok=True)
if mode in ("all", "size"):
    for mb in SIZE_MBS:
        path = data_dir / f"size_{mb}mb.fasta"
        if regenerate or not path.is_file():
            write_size(mb, path)
if mode in ("all", "seq"):
    for count in FIXED_COUNTS:
        path = data_dir / f"seqs_fixed_{count}.fasta"
        if regenerate or not path.is_file():
            write_fixed(count, path)
PY
}

# Prefer hardlinks from bench/index/data/ when present (same generator shape).
# Falls back to generating under bench/stats/data/ so stats does not depend on index.
link_or_skip_from_index() {
    local name="$1"
    local dest="$SCALING_DIR/$name"
    local src="$BENCH_ROOT/index/data/$name"
    [[ -f "$dest" ]] && return 0
    [[ -f "$src" ]] || return 1
    mkdir -p "$SCALING_DIR"
    if ln "$src" "$dest" 2>/dev/null; then
        echo "  linked $name from bench/index/data/"
        # Reuse .fai when available; .zfi is rebuilt by preload (index bench often has no .zfi).
        if [[ ! -f "${dest}.fai" && -f "${src}.fai" ]]; then
            ln "${src}.fai" "${dest}.fai" 2>/dev/null || true
        fi
        return 0
    fi
    return 1
}

ensure_scaling_fixtures() {
    local need_size=false need_seq=false
    mkdir -p "$SCALING_DIR"

    if ! $REGENERATE_FIXTURES; then
        local mb count
        if $DO_SCALE_SIZE; then
            for mb in "${SIZE_MBS[@]}"; do
                link_or_skip_from_index "size_${mb}mb.fasta" || true
            done
        fi
        if $DO_SCALE_SEQS; then
            for count in "${SEQ_FIXED_COUNTS[@]}"; do
                link_or_skip_from_index "seqs_fixed_${count}.fasta" || true
            done
        fi
    fi

    if $DO_SCALE_SIZE; then
        local mb
        for mb in "${SIZE_MBS[@]}"; do
            [[ -f "$SCALING_DIR/size_${mb}mb.fasta" ]] || need_size=true
        done
    fi
    if $DO_SCALE_SEQS; then
        local count
        for count in "${SEQ_FIXED_COUNTS[@]}"; do
            [[ -f "$SCALING_DIR/seqs_fixed_${count}.fasta" ]] || need_seq=true
        done
    fi
    if $REGENERATE_FIXTURES; then
        echo "  Regenerating scaling fixtures under $SCALING_DIR ..."
        generate_scaling_fixtures all
    else
        if $need_size && $need_seq; then
            echo "  Generating missing scaling fixtures under $SCALING_DIR ..."
            generate_scaling_fixtures all
        elif $need_size; then
            echo "  Generating missing size scaling fixtures..."
            generate_scaling_fixtures size
        elif $need_seq; then
            echo "  Generating missing seq scaling fixtures..."
            generate_scaling_fixtures seq
        fi
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
    export BENCH_VER_ZEBRAC="$(bench_tool_version zebrac)"
    export BENCH_VER_ZFASTA="$(bench_tool_version z-fasta 2>/dev/null || echo unknown)"
    export BENCH_VER_SAMTOOLS="$(bench_tool_version samtools 2>/dev/null || echo unknown)"
    local seqtk_ver
    seqtk_ver="$(bench_tool_version seqtk 2>/dev/null | head -1 || true)"
    export BENCH_VER_SEQTK="${seqtk_ver:-unknown}"
    for tool in seqkit noodles rustbio; do
        if bench_has_tool "$tool"; then
            upper="${tool^^}"
            export "BENCH_VER_${upper}=$(bench_tool_version "$tool")"
        fi
    done
}

write_run_manifest() {
    local manifest="$1" timestamp="$2" metadata="$3"
    python3 - "$manifest" "$timestamp" "$metadata" "$RUNS" "$WARMUP" "$ZEBRAC_DURATION_MS" \
        "$DO_VERIFY" \
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

run_verify() {
    echo ""
    echo "================================================================"
    echo "  verify.sh"
    echo "================================================================"
    local log="$RESULTS_DIR/verify_${TIMESTAMP}.log"
    mkdir -p "$RESULTS_DIR"
    if bash "$SCRIPT_DIR/verify.sh" | tee "$log"; then
        local n
        n="$(awk '/^Results:/{print $2; exit}' "$log")"
        export BENCH_VERIFY_PASS="${n:-unknown}"
    else
        export BENCH_VERIFY_PASS=0
        echo "error: verify.sh failed" >&2
        exit 1
    fi
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
    echo "  Policy: .fai is a first-class index for users who skip .zfi."
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

$DO_VERIFY && run_verify
$DO_PERF && run_perf
$DO_REPORT && run_report

echo ""
echo "================================================================"
echo "  Stats suite complete"
echo "================================================================"
