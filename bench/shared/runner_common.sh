#!/usr/bin/env bash
# Shared helpers for bench/*/run.sh (after zebrac_runner / tools).
#
# Source this instead of zebrac_runner.sh directly from suite runners.

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    echo "runner_common.sh is a library; source it from a benchmark script." >&2
    exit 1
fi

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/zebrac_runner.sh"

report_python() {
    if [[ -x "$PROJECT_ROOT/.venv/bin/python" ]]; then
        echo "$PROJECT_ROOT/.venv/bin/python"
    else
        echo python3
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
# or .zfi predates the embedded name-table footer (ZFNM).
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
            "$SAMTOOLS" faidx "$fa" > /dev/null 2>&1 || true
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

# Core versions every suite manifest records.
export_manifest_core_versions() {
    export BENCH_VER_ZEBRAC="$(bench_tool_version zebrac)"
    export BENCH_VER_ZFASTA="$(bench_tool_version z-fasta 2>/dev/null || echo unknown)"
    export BENCH_VER_SAMTOOLS="$(bench_tool_version samtools 2>/dev/null || echo unknown)"
}

# Always export BENCH_VER_<TOOL> (unknown if missing). seqtk uses the first line only.
export_manifest_required_tool_versions() {
    local tool upper ver
    for tool in "$@"; do
        case "$tool" in
            seqtk)
                # Match prior stats behavior: empty stdout becomes "unknown".
                ver="$(bench_tool_version seqtk 2>/dev/null | head -1 || true)"
                export BENCH_VER_SEQTK="${ver:-unknown}"
                ;;
            *)
                upper="${tool^^}"
                export "BENCH_VER_${upper}=$(bench_tool_version "$tool" 2>/dev/null || echo unknown)"
                ;;
        esac
    done
}

# Export BENCH_VER_<TOOL> only when the peer is installed.
export_manifest_optional_tool_versions() {
    local tool upper
    for tool in "$@"; do
        if bench_has_tool "$tool"; then
            upper="${tool^^}"
            export "BENCH_VER_${upper}=$(bench_tool_version "$tool")"
        fi
    done
}
