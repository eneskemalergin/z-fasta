#!/usr/bin/env bash
# Verify benchmark tools used by bench/*/run.sh.
#
# Paths and version helpers come from tools.sh (single source of truth).
# Optional: install the pinned pyfaidx into .venv when missing.
# Does not check hyperfine (legacy; suites use zebrac).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=tools.sh
source "$SCRIPT_DIR/tools.sh"

VENV="$PROJECT_ROOT/.venv"
PYFAIDX_PIN="0.9.0.4"
SEQTK_PIN="1.5-r133"
SEQKIT_PIN="2.13.0"
FASTAHACK_PIN="1.0.0"
SAMTOOLS_PIN="1.24"
BEDTOOLS_PIN="2.31.1"
NOODLES_PIN="0.66.0"
RUSTBIO_PIN="4.0.1"

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

    # Prefer an already-resolved faidx CLI at the pin (venv or PATH).
    if bench_has_tool pyfaidx; then
        ver="$(bench_tool_version pyfaidx 2>/dev/null || true)"
        if [[ "$ver" == *"$PYFAIDX_PIN"* ]]; then
            return 0
        fi
    fi

    if [[ -x "$VENV/bin/python" ]]; then
        ver="$("$VENV/bin/python" -c "import pyfaidx; print(pyfaidx.__version__)" 2>/dev/null || true)"
        if [[ "$ver" == "$PYFAIDX_PIN" && -x "$VENV/bin/faidx" ]]; then
            PYFAIDX="$VENV/bin/faidx"
            return 0
        fi
        if [[ "$ver" != "$PYFAIDX_PIN" ]]; then
            warn "pyfaidx not at pin $PYFAIDX_PIN in .venv; installing"
            if "$VENV/bin/pip" install --quiet "pyfaidx==$PYFAIDX_PIN"; then
                ver="$("$VENV/bin/python" -c "import pyfaidx; print(pyfaidx.__version__)" 2>/dev/null || true)"
                if [[ "$ver" == "$PYFAIDX_PIN" && -x "$VENV/bin/faidx" ]]; then
                    PYFAIDX="$VENV/bin/faidx"
                    return 0
                fi
            fi
        fi
    fi
    ver="$(python3 -c "import pyfaidx; print(pyfaidx.__version__)" 2>/dev/null || true)"
    if [[ "$ver" == "$PYFAIDX_PIN" ]] && command -v faidx >/dev/null 2>&1; then
        PYFAIDX="$(command -v faidx)"
        return 0
    fi

    fail "pyfaidx $PYFAIDX_PIN not found (pip install pyfaidx==$PYFAIDX_PIN into .venv or PATH)"
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
            noodles) echo "     Build: (cd tools/noodles_wrapper && cargo build --release)" ;;
            rustbio) echo "     Build: (cd tools/rustbio_wrapper && cargo build --release)" ;;
            seqtk) echo "     Build: make -C tools/seqtk" ;;
            fastahack) echo "     Build: make -C tools/fastahack-1.0.0" ;;
            seqkit) echo "     Place binary at tools/seqkit" ;;
            samtools|bedtools) echo "     Install via apt/conda and ensure on PATH" ;;
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
            if [[ "$path" != *"/fastahack-${FASTAHACK_PIN}/"* ]]; then
                fail "fastahack path does not include fastahack-$FASTAHACK_PIN"
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
        echo "All required tools verified."
    else
        echo "$fail_count tool(s) missing or misconfigured; fix above before running benchmarks."
        exit 1
    fi
}

main "$@"
