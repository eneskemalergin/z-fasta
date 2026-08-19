#!/usr/bin/env bash
# Shared benchmark tool discovery helpers.
#
# Source this file from bench/*/run.sh, install_tools.sh, or
# zebrac_runner.sh / download_data.sh. It defines tool paths, availability
# checks, version helpers, and small path utilities.

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    echo "tools.sh is a library; source it from a benchmark script." >&2
    exit 1
fi

BENCH_SHARED_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH_ROOT="$(cd "$BENCH_SHARED_DIR/.." && pwd)"
PROJECT_ROOT="$(cd "$BENCH_ROOT/.." && pwd)"
TOOLS_DIR="$PROJECT_ROOT/tools"
TOOLS_BIN_DIR="$TOOLS_DIR/bin"
TOOLS_LIB_DIR="$TOOLS_DIR/lib"
TOOLS_VENV_DIR="$TOOLS_DIR/venv"
# shellcheck source=../../tools/versions.sh
source "$TOOLS_DIR/versions.sh"

ZFASTA="${ZFASTA:-$PROJECT_ROOT/zig-out/bin/z-fasta}"
ZEBRAC="${ZEBRAC:-$TOOLS_DIR/zebrac}"
SAMTOOLS="${SAMTOOLS:-$TOOLS_BIN_DIR/samtools}"
BEDTOOLS="${BEDTOOLS:-$TOOLS_BIN_DIR/bedtools}"
SEQKIT="${SEQKIT:-$TOOLS_BIN_DIR/seqkit}"
SEQTK="${SEQTK:-$TOOLS_BIN_DIR/seqtk}"
FASTAHACK="${FASTAHACK:-$TOOLS_BIN_DIR/fastahack}"
NOODLES="${NOODLES:-$TOOLS_BIN_DIR/noodles}"
RUSTBIO="${RUSTBIO:-$TOOLS_BIN_DIR/rustbio}"
PYFAIDX="${PYFAIDX:-$TOOLS_BIN_DIR/faidx}"

bench_tool_path() {
    local name="$1"
    case "$name" in
        z-fasta) echo "$ZFASTA" ;;
        zebrac) echo "$ZEBRAC" ;;
        samtools) echo "$SAMTOOLS" ;;
        bedtools) echo "$BEDTOOLS" ;;
        seqkit) echo "$SEQKIT" ;;
        seqtk) echo "$SEQTK" ;;
        fastahack) echo "$FASTAHACK" ;;
        pyfaidx) echo "$PYFAIDX" ;;
        noodles) echo "$NOODLES" ;;
        rustbio) echo "$RUSTBIO" ;;
        *) return 1 ;;
    esac
}

bench_has_tool() {
    local path
    path="$(bench_tool_path "$1")" || return 1
    [[ -n "$path" && -x "$path" ]]
}

bench_require_tool() {
    local name="$1"
    if ! bench_has_tool "$name"; then
        echo "error: required benchmark tool not found or not executable: $name" >&2
        echo "       resolved path: $(bench_tool_path "$name" 2>/dev/null || echo '<unknown>')" >&2
        return 1
    fi
}

bench_tool_version() {
    local name="$1"
    local path
    path="$(bench_tool_path "$name")" || return 1
    case "$name" in
        z-fasta|zebrac|bedtools|pyfaidx)
            [[ -x "$path" ]] && "$path" --version 2>&1 | awk 'NR==1{print; exit}'
            ;;
        samtools)
            [[ -x "$path" ]] && "$path" --version 2>&1 | awk 'NR==1{print; exit}'
            ;;
        seqkit)
            [[ -x "$path" ]] && "$path" version 2>&1 | awk 'NR==1{print; exit}'
            ;;
        seqtk)
            # seqtk prints Version on stderr/stdout and exits nonzero with no args.
            [[ -x "$path" ]] || return 1
            { "$path" 2>&1 || true; } | awk '/^Version:/{print "seqtk " $2; exit}'
            ;;
        fastahack)
            echo "fastahack $FASTAHACK_VERSION (local tool bundle)"
            ;;
        noodles|rustbio)
            [[ -x "$path" ]] && "$path" --version 2>&1 | awk 'NR==1{print; exit}'
            ;;
        *)
            return 1
            ;;
    esac
}

# Portable file size in bytes (GNU and BSD stat).
file_size_bytes() {
    local path="$1"
    stat --printf='%s' "$path" 2>/dev/null || stat -f '%z' "$path" 2>/dev/null || echo 0
}

# Generate messy layouts into bench/shared/cache/ (fixtures and/or perf).
# Usage: bench_ensure_messy [--force] --fixtures|--perf|--fixtures --perf
bench_ensure_messy() {
    local py
    if [[ -x "$TOOLS_VENV_DIR/bin/python" ]]; then
        py="$TOOLS_VENV_DIR/bin/python"
    else
        py=python3
    fi
    "$py" "$BENCH_SHARED_DIR/generate_messy.py" "$@"
}

# Generate synthetic scaling FASTAs into bench/shared/cache/scaling/.
# Usage: bench_ensure_scaling [--force] [--clean-legacy] [--mode all|size|seq]
bench_ensure_scaling() {
    local py
    if [[ -x "$TOOLS_VENV_DIR/bin/python" ]]; then
        py="$TOOLS_VENV_DIR/bin/python"
    else
        py=python3
    fi
    "$py" "$BENCH_SHARED_DIR/generate_scaling.py" "$@"
}
