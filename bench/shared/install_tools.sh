#!/usr/bin/env bash
# Verify benchmark tools used by bench/*/run.sh.
#
# Paths and version helpers come from tools.sh (single source of truth).
# This script verifies the local installation created by tools/install.sh.
# It does not install from PATH or mutate a system environment.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=tools.sh
source "$SCRIPT_DIR/tools.sh"
source "$TOOLS_DIR/versions.sh"

VENV="$TOOLS_DIR/venv"
PYFAIDX_PIN="$PYFAIDX_VERSION"
MATPLOTLIB_PIN="$MATPLOTLIB_VERSION"
PANDAS_PIN="$PANDAS_VERSION"
TABULATE_PIN="$TABULATE_VERSION"
SEQTK_PIN="$SEQTK_DISPLAY_VERSION"
SEQKIT_PIN="$SEQKIT_VERSION"
FASTAHACK_PIN="$FASTAHACK_VERSION"
SAMTOOLS_PIN="$SAMTOOLS_VERSION"
BEDTOOLS_PIN="$BEDTOOLS_VERSION"
NOODLES_PIN="$NOODLES_VERSION"
RUSTBIO_PIN="$RUSTBIO_VERSION"

ok() { printf '  [ok] %s\n' "$*"; }
warn() { printf '  [warn] %s\n' "$*"; }
fail() { printf '  [fail] %s\n' "$*"; }

# Tools every full suite path expects. Missing any exits nonzero.
REQUIRED_TOOLS=(
    z-fasta
    zebrac
    samtools
    bedtools
    seqkit
    seqtk
    fastahack
    pyfaidx
    noodles
    rustbio
)

ensure_pyfaidx() {
    local ver=""
    if [[ -x "$VENV/bin/python" ]]; then
        ver="$("$VENV/bin/python" -c "import pyfaidx; print(pyfaidx.__version__)" 2>/dev/null || true)"
    fi
    if [[ "$ver" == "$PYFAIDX_PIN" && -x "$PYFAIDX" ]]; then
        return 0
    fi
    fail "pyfaidx $PYFAIDX_PIN not found in $VENV (run tools/install.sh)"
    return 1
}

ensure_report_python() {
    local versions=""
    if [[ -x "$VENV/bin/python" ]]; then
        versions="$("$VENV/bin/python" -c 'import matplotlib, pandas, tabulate; print(matplotlib.__version__, pandas.__version__, tabulate.__version__)' 2>/dev/null || true)"
    fi
    if [[ "$versions" == "$MATPLOTLIB_PIN $PANDAS_PIN $TABULATE_PIN" ]]; then
        return 0
    fi
    fail "report packages are not $MATPLOTLIB_PIN/$PANDAS_PIN/$TABULATE_PIN in $VENV (run tools/install.sh)"
    return 1
}

check_tool() {
    local name="$1"
    local path ver

    echo "-- $name --"
    if [[ "$name" == "pyfaidx" ]]; then
        ensure_pyfaidx || return 1
    fi

    if ! bench_has_tool "$name"; then
        path="$(bench_tool_path "$name" 2>/dev/null || echo '<unknown>')"
        fail "not found or not executable (resolved: $path)"
        case "$name" in
            z-fasta) echo "     Build: ./zig build -Doptimize=ReleaseFast" ;;
            zebrac) echo "     Expected at: $TOOLS_DIR/zebrac" ;;
            noodles|rustbio|fastahack|seqkit|seqtk|samtools|bedtools)
                echo "     Build: tools/install.sh"
                ;;
        esac
        return 1
    fi

    path="$(bench_tool_path "$name")"
    ver="$(bench_tool_version "$name" 2>/dev/null || true)"
    if [[ -n "$ver" ]]; then
        ok "$ver ($path)"
    else
        ok "$path"
    fi

    case "$name" in
        seqtk)
            local got
            got="$(bench_tool_version seqtk 2>/dev/null | awk '{print $2}' || true)"
            if [[ -n "$got" && "$got" != "$SEQTK_PIN" ]]; then
                fail "seqtk pin is $SEQTK_PIN; got $got"
                return 1
            fi
            ;;
        seqkit)
            local got
            got="$(bench_tool_version seqkit 2>/dev/null || true)"
            if [[ -n "$got" && "$got" != *"v$SEQKIT_PIN"* ]]; then
                fail "seqkit pin is v$SEQKIT_PIN; got $got"
                return 1
            fi
            ;;
        fastahack)
            if [[ "$(bench_tool_version fastahack 2>/dev/null || true)" != *"$FASTAHACK_PIN"* ]]; then
                fail "fastahack pin is $FASTAHACK_PIN"
                return 1
            fi
            ;;
        samtools)
            local got
            got="$(bench_tool_version samtools 2>/dev/null | awk '{print $2}' || true)"
            if [[ -n "$got" && "$got" != "$SAMTOOLS_PIN" ]]; then
                fail "samtools latest pin is $SAMTOOLS_PIN; got $got"
                return 1
            fi
            ;;
        bedtools)
            local got
            got="$(bench_tool_version bedtools 2>/dev/null | sed -n 's/.*v\([0-9][0-9.]*\).*/\1/p' || true)"
            if [[ -n "$got" && "$got" != "$BEDTOOLS_PIN" ]]; then
                fail "bedtools latest pin is $BEDTOOLS_PIN; got $got"
                return 1
            fi
            ;;
        noodles)
            local got
            got="$(bench_tool_version noodles 2>/dev/null | awk '{print $2}' || true)"
            if [[ -n "$got" && "$got" != "$NOODLES_PIN" ]]; then
                fail "noodles-fasta pin is $NOODLES_PIN; got $got"
                return 1
            fi
            ;;
        rustbio)
            local got
            got="$(bench_tool_version rustbio 2>/dev/null | awk '{print $2}' || true)"
            if [[ -n "$got" && "$got" != "$RUSTBIO_PIN" ]]; then
                fail "rust-bio pin is $RUSTBIO_PIN; got $got"
                return 1
            fi
            ;;
    esac
    return 0
}

main() {
    echo ""
    echo "z-fasta benchmark tool verification"
    echo "==================================="
    echo "Paths from: $SCRIPT_DIR/tools.sh"
    echo ""

    local fail_count=0
    local name
    for name in "${REQUIRED_TOOLS[@]}"; do
        if ! check_tool "$name"; then
            fail_count=$((fail_count + 1))
        fi
        echo ""
    done

    echo "-- report Python --"
    if ensure_report_python; then
        ok "matplotlib $MATPLOTLIB_PIN, pandas $PANDAS_PIN, tabulate $TABULATE_PIN ($VENV/bin/python)"
    else
        fail_count=$((fail_count + 1))
    fi
    echo ""

    echo "Summary"
    echo "-------"
    for name in "${REQUIRED_TOOLS[@]}"; do
        if bench_has_tool "$name"; then
            ver="$(bench_tool_version "$name" 2>/dev/null || true)"
            printf '  %-10s %s\n' "$name" "${ver:-present}"
        else
            printf '  %-10s missing\n' "$name"
        fi
    done
    echo ""

    if [[ "$fail_count" -eq 0 ]]; then
        echo "All required tools and report dependencies verified."
    else
        echo "$fail_count tool(s) missing or misconfigured; fix above before running benchmarks."
        exit 1
    fi
}

main "$@"
