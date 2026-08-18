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
PEER_TOOLS_DIR="$PROJECT_ROOT/.peer-tools"

ZFASTA="${ZFASTA:-$PROJECT_ROOT/zig-out/bin/z-fasta}"
ZEBRAC="${ZEBRAC:-$TOOLS_DIR/zebrac}"
if [[ -n "${SAMTOOLS:-}" ]]; then
    : # Respect an explicit SAMTOOLS override.
elif [[ -x "$PEER_TOOLS_DIR/bin/samtools" ]]; then
    SAMTOOLS="${SAMTOOLS:-$PEER_TOOLS_DIR/bin/samtools}"
    export LD_LIBRARY_PATH="$PEER_TOOLS_DIR/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
else
    SAMTOOLS="samtools"
fi
if [[ -n "${BEDTOOLS:-}" ]]; then
    : # Respect an explicit BEDTOOLS override.
elif [[ -x "$PEER_TOOLS_DIR/bin/bedtools" ]]; then
    BEDTOOLS="$PEER_TOOLS_DIR/bin/bedtools"
else
    BEDTOOLS="bedtools"
fi
SEQKIT="${SEQKIT:-$TOOLS_DIR/seqkit}"
SEQTK="${SEQTK:-$TOOLS_DIR/seqtk/seqtk}"
FASTAHACK="${FASTAHACK:-$TOOLS_DIR/fastahack-1.0.0/fastahack}"
NOODLES="${NOODLES:-$TOOLS_DIR/noodles_wrapper/target/release/noodles_wrapper}"
RUSTBIO="${RUSTBIO:-$TOOLS_DIR/rustbio_wrapper/target/release/rustbio_wrapper}"

if [[ -x "$PROJECT_ROOT/.venv/bin/faidx" ]]; then
    PYFAIDX="${PYFAIDX:-$PROJECT_ROOT/.venv/bin/faidx}"
elif command -v faidx >/dev/null 2>&1; then
    PYFAIDX="${PYFAIDX:-$(command -v faidx)}"
else
    PYFAIDX="${PYFAIDX:-}"
fi

bench_tool_path() {
    local name="$1"
    case "$name" in
        z-fasta) echo "$ZFASTA" ;;
        zebrac) echo "$ZEBRAC" ;;
        samtools) command -v "$SAMTOOLS" 2>/dev/null || true ;;
        bedtools) command -v "$BEDTOOLS" 2>/dev/null || true ;;
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
            command -v "$SAMTOOLS" >/dev/null 2>&1 && "$SAMTOOLS" --version 2>&1 | awk 'NR==1{print; exit}'
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
            echo "fastahack 1.0.0 (directory pin)"
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
    if [[ -x "$PROJECT_ROOT/.venv/bin/python" ]]; then
        py="$PROJECT_ROOT/.venv/bin/python"
    else
        py=python3
    fi
    "$py" "$BENCH_SHARED_DIR/generate_messy.py" "$@"
}

# Generate synthetic scaling FASTAs into bench/shared/cache/scaling/.
# Usage: bench_ensure_scaling [--force] [--clean-legacy] [--mode all|size|seq]
bench_ensure_scaling() {
    local py
    if [[ -x "$PROJECT_ROOT/.venv/bin/python" ]]; then
        py="$PROJECT_ROOT/.venv/bin/python"
    else
        py=python3
    fi
    "$py" "$BENCH_SHARED_DIR/generate_scaling.py" "$@"
}
