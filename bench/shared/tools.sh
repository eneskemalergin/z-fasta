#!/usr/bin/env bash
# Shared benchmark tool discovery helpers.
#
# Source this file from bench/*/run.sh, verify.sh, install_tools.sh, or
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

ZFASTA="${ZFASTA:-$PROJECT_ROOT/zig-out/bin/z-fasta}"
ZEBRAC="${ZEBRAC:-$TOOLS_DIR/zebrac}"
SAMTOOLS="${SAMTOOLS:-samtools}"
BEDTOOLS="${BEDTOOLS:-bedtools}"
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
        noodles)
            echo "noodles-fasta 0.61 (wrapper)"
            ;;
        rustbio)
            echo "rust-bio 2.2 (custom indexer wrapper)"
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
