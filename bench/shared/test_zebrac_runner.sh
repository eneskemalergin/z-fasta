#!/usr/bin/env bash
# Smoke test for shared zebrac Bash utilities.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/runner_common.sh"

# Metadata tests use controlled version probes and never launch a measurement.
(
    bench_tool_version() { return 1; }
    export_manifest_core_versions
    [[ "$BENCH_VER_ZEBRAC" == unknown && "$BENCH_VER_ZFASTA" == unknown && "$BENCH_VER_SAMTOOLS" == unknown ]]
    bench_tool_version() { return 0; }
    export_manifest_core_versions
    [[ "$BENCH_VER_ZEBRAC" == unknown && "$BENCH_VER_ZFASTA" == unknown && "$BENCH_VER_SAMTOOLS" == unknown ]]
    export_manifest_required_tool_versions seqtk seqkit
    [[ "$BENCH_VER_SEQTK" == unknown && "$BENCH_VER_SEQKIT" == unknown ]]
    bench_tool_version() { printf 'partial'; return 1; }
    export_manifest_core_versions
    [[ "$BENCH_VER_ZEBRAC" == unknown && "$BENCH_VER_ZFASTA" == unknown && "$BENCH_VER_SAMTOOLS" == unknown ]]
    bench_tool_version() { printf '%s version\n' "$1"; }
    export_manifest_core_versions
    [[ "$BENCH_VER_ZEBRAC" == 'zebrac version' && "$BENCH_VER_ZFASTA" == 'z-fasta version' && "$BENCH_VER_SAMTOOLS" == 'samtools version' ]]
)
python3 - "$SCRIPT_DIR/zebrac_runner.sh" <<'PYTEST'
import json
import subprocess
import sys

payload = "".join(chr(byte) for byte in range(1, 32)) + 'quote" slash\\ plain'
encoded = subprocess.run(
    ["bash", "-c", 'source "$1"; zebrac_json_string "$2"', "bash", sys.argv[1], payload],
    check=True,
    capture_output=True,
).stdout
assert json.loads(encoded) == payload
PYTEST

if [[ "${1:-}" == --metadata-only ]]; then
    echo "shared runner metadata tests passed"
    exit 0
fi

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
