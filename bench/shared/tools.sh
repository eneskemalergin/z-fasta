#!/usr/bin/env bash
# Shared benchmark tool discovery helpers.
#
# Source this file from bench/*/run_benchmarks.sh scripts. It defines paths,
# availability checks, and lightweight version helpers for first-party,
# Tier 1, Tier 2, and runner tools.

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

bench_tool_tier() {
    case "$1" in
        z-fasta) echo "first-party" ;;
        zebrac) echo "runner" ;;
        noodles|rustbio) echo "tier2" ;;
        samtools|bedtools|seqkit|seqtk|fastahack|pyfaidx) echo "tier1" ;;
        *) echo "unknown" ;;
    esac
}

bench_tool_label() {
    case "$1" in
        z-fasta) echo "z-fasta" ;;
        zebrac) echo "zebrac" ;;
        noodles) echo "noodles-fasta 0.61" ;;
        rustbio) echo "rust-bio 2.3 (custom indexer)" ;;
        *) echo "$1" ;;
    esac
}

bench_tool_supports_suite() {
    local name="$1"
    local suite="$2"
    case "$name:$suite" in
        z-fasta:index|z-fasta:get|z-fasta:stats) return 0 ;;
        samtools:index|samtools:get) return 0 ;;
        bedtools:get) return 0 ;;
        seqkit:index|seqkit:get|seqkit:stats) return 0 ;;
        seqtk:get|seqtk:stats) return 0 ;;
        fastahack:index|fastahack:get) return 0 ;;
        pyfaidx:index|pyfaidx:get) return 0 ;;
        noodles:index|noodles:get|noodles:stats) return 0 ;;
        rustbio:index|rustbio:get|rustbio:stats) return 0 ;;
        *) return 1 ;;
    esac
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
            [[ -x "$path" ]] && "$path" 2>&1 | awk '/^Version:/{print "seqtk " $2; exit}'
            ;;
        fastahack)
            echo "fastahack 1.0.0 (directory pin)"
            ;;
        noodles)
            echo "noodles-fasta 0.61 (wrapper)"
            ;;
        rustbio)
            echo "rust-bio 2.3 (custom indexer wrapper)"
            ;;
        *)
            return 1
            ;;
    esac
}
