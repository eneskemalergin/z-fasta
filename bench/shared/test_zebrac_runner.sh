#!/usr/bin/env bash
# Smoke test for shared zebrac Bash utilities.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/zebrac_runner.sh"

TMPDIR_LOCAL="$(mktemp -d)"
trap 'rm -rf "$TMPDIR_LOCAL"' EXIT

raw_json="$TMPDIR_LOCAL/raw.json"
metadata_jsonl="$TMPDIR_LOCAL/metadata.jsonl"

ZEBRAC_DURATION_MS=100
ZEBRAC_MIN_SAMPLES=1
ZEBRAC_MAX_SAMPLES=1
ZEBRAC_WARMUP=0

zebrac_run_single \
    "$raw_json" \
    "$metadata_jsonl" \
    "smoke" \
    "runner" \
    "true" \
    "true" \
    "coreutils" \
    "1" \
    "" \
    "/bin/true"

[[ -s "$raw_json" ]] || { echo "raw zebrac JSON was not written" >&2; exit 1; }
[[ -s "$metadata_jsonl" ]] || { echo "metadata JSONL was not written" >&2; exit 1; }

grep -q '"zebrac_version"' "$raw_json" || { echo "raw JSON does not look like zebrac output" >&2; exit 1; }
grep -q '"schema_version":"bench.meta.v1"' "$metadata_jsonl" || { echo "metadata schema missing" >&2; exit 1; }
grep -q '"suite":"smoke"' "$metadata_jsonl" || { echo "metadata suite missing" >&2; exit 1; }
grep -q '"tool":"true"' "$metadata_jsonl" || { echo "metadata tool missing" >&2; exit 1; }
grep -q '"command":"/bin/true"' "$metadata_jsonl" || { echo "metadata command missing" >&2; exit 1; }

echo "shared zebrac runner smoke test passed"
