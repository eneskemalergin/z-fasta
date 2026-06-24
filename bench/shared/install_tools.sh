#!/usr/bin/env bash
# bench/shared/install_tools.sh
#
# Installs / verifies all v0.2.5 benchmark tools at pinned versions and prints
# a canonical VERSION_PINS block for inclusion in REPORT.md.
#
# Pinned versions (last verified 2026-03-30):
#   pyfaidx   0.9.0.3    (PyPI)
#   seqtk     1.5-r133   (pre-built in tools/seqtk/)
#   fastahack 1.0.0      (pre-built in tools/fastahack-1.0.0/)
#   seqkit    2.13.0     (pre-built in tools/seqkit)
#   samtools  1.13       (system)
#   hyperfine 1.12.0     (system)
#
# Usage:
#   bash bench/shared/install_tools.sh
#
# After running, all tools are callable via the helper functions in this script:
#   seqtk_bin, fastahack_bin, seqkit_bin, faidx_bin

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
TOOLS_DIR="$REPO_ROOT/tools"
VENV="$REPO_ROOT/.venv"

# ── Pinned versions ────────────────────────────────────────────────────────────
PYFAIDX_PIN="0.9.0.3"
SEQTK_PIN="1.5-r133"
FASTAHACK_PIN="1.0.0"
SEQKIT_PIN="2.13.0"

RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'; NC='\033[0m'
ok()   { echo -e "${GREEN}  ✓${NC} $*"; }
warn() { echo -e "${YELLOW}  ⚠${NC} $*"; }
fail() { echo -e "${RED}  ✗${NC} $*"; }

# ── Helper: resolve tool paths ─────────────────────────────────────────────────
seqtk_bin()     { echo "$TOOLS_DIR/seqtk/seqtk"; }
fastahack_bin() { echo "$TOOLS_DIR/fastahack-1.0.0/fastahack"; }
seqkit_bin()    { echo "$TOOLS_DIR/seqkit"; }
faidx_bin()     {
    # pyfaidx ships a 'faidx' CLI entry-point; prefer venv, then system python
    if "$VENV/bin/faidx" --version &>/dev/null 2>&1; then
        echo "$VENV/bin/faidx"
    else
        echo "faidx"
    fi
}

# ── Check / install pyfaidx ────────────────────────────────────────────────────
check_pyfaidx() {
    echo "── pyfaidx ──────────────────────────────────────────────────────────────"
    # Check base-env python first (that's where it's currently installed)
    local ver
    ver=$(python -c "import pyfaidx; print(pyfaidx.__version__)" 2>/dev/null || true)
    if [[ "$ver" == "$PYFAIDX_PIN" ]]; then
        ok "pyfaidx $ver (base env python)"
        return 0
    fi

    # Check venv
    ver=$("$VENV/bin/python" -c "import pyfaidx; print(pyfaidx.__version__)" 2>/dev/null || true)
    if [[ "$ver" == "$PYFAIDX_PIN" ]]; then
        ok "pyfaidx $ver (venv)"
        return 0
    fi

    warn "pyfaidx not found at pin $PYFAIDX_PIN; installing into venv"
    "$VENV/bin/pip" install --quiet "pyfaidx==$PYFAIDX_PIN"
    ver=$("$VENV/bin/python" -c "import pyfaidx; print(pyfaidx.__version__)" 2>/dev/null || true)
    if [[ "$ver" == "$PYFAIDX_PIN" ]]; then
        ok "pyfaidx $ver installed"
    else
        fail "pyfaidx install failed (got: '${ver:-none}')"
        return 1
    fi
}

# ── Check seqtk (pre-built) ────────────────────────────────────────────────────
check_seqtk() {
    echo "── seqtk ────────────────────────────────────────────────────────────────"
    local bin
    bin="$(seqtk_bin)"
    if [[ ! -x "$bin" ]]; then
        fail "seqtk binary not found at $bin"
        echo "     Build it: cd tools/seqtk && make"
        return 1
    fi
    local ver
    ver=$("$bin" 2>&1 | awk '/^Version:/{print $2}')
    if [[ "$ver" == "$SEQTK_PIN" ]]; then
        ok "seqtk $ver ($bin)"
    else
        warn "seqtk version mismatch: got '$ver', expected '$SEQTK_PIN'"
        warn "Pin may need updating; tool still usable"
    fi
}

# ── Check fastahack (pre-built) ────────────────────────────────────────────────
check_fastahack() {
    echo "── fastahack ────────────────────────────────────────────────────────────"
    local bin
    bin="$(fastahack_bin)"
    if [[ ! -x "$bin" ]]; then
        fail "fastahack binary not found at $bin"
        echo "     Build it: cd tools/fastahack-1.0.0 && make"
        return 1
    fi
    # fastahack has no --version flag; version is embedded in the directory name
    ok "fastahack $FASTAHACK_PIN ($bin)"
    ok "  Note: no standalone version flag; version derived from directory name"
}

# ── Check seqkit (pre-built) ───────────────────────────────────────────────────
check_seqkit() {
    echo "── seqkit ───────────────────────────────────────────────────────────────"
    local bin
    bin="$(seqkit_bin)"
    if [[ ! -x "$bin" ]]; then
        fail "seqkit binary not found at $bin"
        echo "     Download from https://github.com/shenwei356/seqkit/releases"
        return 1
    fi
    local ver
    ver=$("$bin" version 2>&1 | sed 's/seqkit //')
    if [[ "$ver" == "v$SEQKIT_PIN" ]]; then
        ok "seqkit $ver ($bin)"
    else
        warn "seqkit version mismatch: got '$ver', expected 'v$SEQKIT_PIN'"
        warn "Pin may need updating; tool still usable"
    fi
}

# ── Check system tools ─────────────────────────────────────────────────────────
check_samtools() {
    echo "── samtools ─────────────────────────────────────────────────────────────"
    if ! command -v samtools &>/dev/null; then
        fail "samtools not found in PATH"
        echo "     Install: apt install samtools  OR  conda install -c bioconda samtools"
        return 1
    fi
    local ver
    ver=$(samtools --version 2>&1 | awk 'NR==1{print $2}')
    ok "samtools $ver ($(command -v samtools))"
}

check_hyperfine() {
    echo "── hyperfine ────────────────────────────────────────────────────────────"
    if ! command -v hyperfine &>/dev/null; then
        fail "hyperfine not found in PATH"
        echo "     Install: cargo install hyperfine  OR  apt install hyperfine"
        return 1
    fi
    local ver
    ver=$(hyperfine --version 2>&1 | awk '{print $2}')
    ok "hyperfine $ver ($(command -v hyperfine))"
}

# ── Print version-pins block (for copy-paste into REPORT.md) ──────────────────
print_pins() {
    echo ""
    echo "## VERSION_PINS (copy into REPORT.md)"
    echo "| Tool       | Version         | Source                              |"
    echo "|------------|-----------------|-------------------------------------|"

    local seqtk_ver fastahack_ver seqkit_ver samtools_ver hyperfine_ver pyfaidx_ver
    seqtk_ver=$( { "$(seqtk_bin)" 2>&1; true; } | awk '/^Version:/{print $2}' )
    fastahack_ver="$FASTAHACK_PIN (dir)"
    seqkit_ver=$("$(seqkit_bin)" version 2>&1 | sed 's/seqkit //') || seqkit_ver="unknown"
    samtools_ver=$(samtools --version 2>&1 | awk 'NR==1{print $2}') || samtools_ver="unknown"
    hyperfine_ver=$(hyperfine --version 2>&1 | awk '{print $2}') || hyperfine_ver="unknown"
    pyfaidx_ver=$(python -c "import pyfaidx; print(pyfaidx.__version__)" 2>/dev/null \
                  || "$VENV/bin/python" -c "import pyfaidx; print(pyfaidx.__version__)" 2>/dev/null \
                  || echo "not found")

    printf "| %-10s | %-15s | %-35s |\n" "pyfaidx"   "$pyfaidx_ver"    "PyPI / base-env python"
    printf "| %-10s | %-15s | %-35s |\n" "seqtk"     "$seqtk_ver"      "tools/seqtk/ (make)"
    printf "| %-10s | %-15s | %-35s |\n" "fastahack" "$fastahack_ver"  "tools/fastahack-1.0.0/ (make)"
    printf "| %-10s | %-15s | %-35s |\n" "seqkit"    "$seqkit_ver"     "tools/seqkit (binary)"
    printf "| %-10s | %-15s | %-35s |\n" "samtools"  "$samtools_ver"   "system"
    printf "| %-10s | %-15s | %-35s |\n" "hyperfine" "$hyperfine_ver"  "system"
}

# ── Main ───────────────────────────────────────────────────────────────────────
main() {
    echo ""
    echo "z-fasta v0.2.5 - benchmark tool verification"
    echo "=============================================="
    local fail_count=0

    check_pyfaidx  || (( fail_count++ )) || true
    echo ""
    check_seqtk    || (( fail_count++ )) || true
    echo ""
    check_fastahack || (( fail_count++ )) || true
    echo ""
    check_seqkit   || (( fail_count++ )) || true
    echo ""
    check_samtools || (( fail_count++ )) || true
    echo ""
    check_hyperfine || (( fail_count++ )) || true

    print_pins

    echo ""
    if [[ "$fail_count" -eq 0 ]]; then
        echo -e "${GREEN}All tools verified; ready to run benchmarks.${NC}"
    else
        echo -e "${RED}${fail_count} tool(s) missing or misconfigured; fix above before running benchmarks.${NC}"
        exit 1
    fi
}

main "$@"
