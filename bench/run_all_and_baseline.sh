#!/usr/bin/env bash
# Run index + get + stats, regenerate REPORT.md files, and save a baseline snapshot.
#
# This wraps the per-suite run_benchmarks.sh scripts and generate_report.py.
# Use it when you want one coherent results bundle for save_baseline.py.
#
# Usage:
#   bash bench/run_all_and_baseline.sh [--runs N] [--skip-rc] [--skip-scaling] [--label NAME]
#
# Prerequisites: ./zig build -Doptimize=ReleaseFast, bench/shared/data/, hyperfine

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
RUNS=5
SKIP_RC=false
SKIP_SCALING=false
LABEL=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --runs) RUNS="$2"; shift 2 ;;
        --skip-rc) SKIP_RC=true; shift ;;
        --skip-scaling) SKIP_SCALING=true; shift ;;
        --label) LABEL="$2"; shift 2 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

cd "$PROJECT_ROOT"

echo "=== Building ReleaseFast ==="
./zig build -Doptimize=ReleaseFast

echo ""
echo "=== Index benchmarks ==="
INDEX_ARGS=(--runs "$RUNS")
if $SKIP_SCALING; then
    INDEX_ARGS+=(--skip-scaling)
fi
bash bench/index/run_benchmarks.sh "${INDEX_ARGS[@]}"

echo ""
echo "=== GET benchmarks ==="
GET_ARGS=(--runs "$RUNS")
if $SKIP_SCALING; then
    GET_ARGS+=(--skip-scaling)
fi
if $SKIP_RC; then
    GET_ARGS+=(--skip-rc)
fi
bash bench/get/run_benchmarks.sh "${GET_ARGS[@]}"

echo ""
echo "=== STATS benchmarks ==="
STATS_ARGS=(--runs "$RUNS")
if $SKIP_SCALING; then
    STATS_ARGS+=(--skip-scaling)
fi
bash bench/stats/run_benchmarks.sh "${STATS_ARGS[@]}"

echo ""
echo "=== Generating reports ==="
if [[ -x .venv/bin/python ]]; then
    PY=.venv/bin/python
else
    PY=python3
fi
"$PY" bench/index/generate_report.py
"$PY" bench/get/generate_report.py
"$PY" bench/stats/generate_report.py

echo ""
echo "=== Saving baseline snapshot ==="
BASELINE_ARGS=()
if [[ -n "$LABEL" ]]; then
    BASELINE_ARGS+=(--label "$LABEL")
fi
"$PY" bench/save_baseline.py "${BASELINE_ARGS[@]}"

echo ""
echo "Done. Reports updated; baseline under bench/baselines/ (gitignored)."
