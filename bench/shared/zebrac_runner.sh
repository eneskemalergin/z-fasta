#!/usr/bin/env bash
# Shared zebrac runner helpers for benchmark suite scripts.
#
# This layer is intentionally Bash. It runs tools/zebrac, preserves raw JSON,
# and writes command/workload metadata as JSONL. Python remains in the existing
# report and plotting layer, where raw zebrac JSON can be joined with metadata.

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    echo "zebrac_runner.sh is a library; source it from a benchmark script." >&2
    exit 1
fi

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/tools.sh"

ZEBRAC_DURATION_MS="${ZEBRAC_DURATION_MS:-5000}"
ZEBRAC_MIN_SAMPLES="${ZEBRAC_MIN_SAMPLES:-5}"
ZEBRAC_MAX_SAMPLES="${ZEBRAC_MAX_SAMPLES:-}"
ZEBRAC_WARMUP="${ZEBRAC_WARMUP:-3}"
ZEBRAC_ALLOW_FAILURES="${ZEBRAC_ALLOW_FAILURES:-false}"

declare -a ZEBRAC_BENCH_COMMANDS=()
declare -a ZEBRAC_BENCH_METADATA=()

zebrac_json_escape() {
    local value="$1"
    value="${value//\\/\\\\}"
    value="${value//\"/\\\"}"
    value="${value//$'\n'/\\n}"
    value="${value//$'\r'/\\r}"
    value="${value//$'\t'/\\t}"
    printf '%s' "$value"
}

zebrac_json_string() {
    printf '"%s"' "$(zebrac_json_escape "$1")"
}

zebrac_json_number_or_null() {
    local value="$1"
    if [[ -z "$value" ]]; then
        printf 'null'
    else
        printf '%s' "$value"
    fi
}

zebrac_clear_commands() {
    ZEBRAC_BENCH_COMMANDS=()
    ZEBRAC_BENCH_METADATA=()
}

zebrac_add_command() {
    if [[ $# -ne 9 ]]; then
        echo "usage: zebrac_add_command <suite> <section> <workload> <tool> <tool_family> <input_bytes> <output_bases> <raw_json> <command>" >&2
        return 1
    fi

    local suite="$1"
    local section="$2"
    local workload="$3"
    local tool="$4"
    local tool_family="$5"
    local input_bytes="$6"
    local output_bases="$7"
    local raw_json="$8"
    local command="$9"

    ZEBRAC_BENCH_COMMANDS+=("$command")
    ZEBRAC_BENCH_METADATA+=("$(
        printf '{'
        printf '"schema_version":"bench.meta.v1",'
        printf '"raw_json":%s,' "$(zebrac_json_string "$raw_json")"
        printf '"suite":%s,' "$(zebrac_json_string "$suite")"
        printf '"section":%s,' "$(zebrac_json_string "$section")"
        printf '"workload":%s,' "$(zebrac_json_string "$workload")"
        printf '"tool":%s,' "$(zebrac_json_string "$tool")"
        printf '"tool_family":%s,' "$(zebrac_json_string "$tool_family")"
        printf '"command":%s,' "$(zebrac_json_string "$command")"
        printf '"input_bytes":%s,' "$(zebrac_json_number_or_null "$input_bytes")"
        printf '"output_bases":%s' "$(zebrac_json_number_or_null "$output_bases")"
        printf '}'
    )")
}

zebrac_write_metadata() {
    local metadata_jsonl="$1"
    mkdir -p "$(dirname "$metadata_jsonl")"
    for row in "${ZEBRAC_BENCH_METADATA[@]}"; do
        printf '%s\n' "$row" >> "$metadata_jsonl"
    done
}

zebrac_run_current_group() {
    if [[ $# -ne 2 ]]; then
        echo "usage: zebrac_run_current_group <raw_json> <metadata_jsonl>" >&2
        return 1
    fi

    local raw_json="$1"
    local metadata_jsonl="$2"

    if [[ "${#ZEBRAC_BENCH_COMMANDS[@]}" -eq 0 ]]; then
        echo "error: no zebrac commands queued" >&2
        return 1
    fi

    bench_require_tool zebrac
    mkdir -p "$(dirname "$raw_json")"

    local args=(
        "$ZEBRAC"
        --quiet
        --duration "$ZEBRAC_DURATION_MS"
        --min-samples "$ZEBRAC_MIN_SAMPLES"
        --warmup "$ZEBRAC_WARMUP"
    )

    if [[ -n "$ZEBRAC_MAX_SAMPLES" ]]; then
        args+=(--max-samples "$ZEBRAC_MAX_SAMPLES")
    fi
    if [[ "$ZEBRAC_ALLOW_FAILURES" == "true" ]]; then
        args+=(--allow-failures)
    fi

    args+=(--json "$raw_json" --)
    args+=("${ZEBRAC_BENCH_COMMANDS[@]}")

    "${args[@]}"
    zebrac_write_metadata "$metadata_jsonl"
}

zebrac_run_single() {
    if [[ $# -ne 10 ]]; then
        echo "usage: zebrac_run_single <raw_json> <metadata_jsonl> <suite> <section> <workload> <tool> <tool_family> <input_bytes> <output_bases> <command>" >&2
        return 1
    fi

    local raw_json="$1"
    local metadata_jsonl="$2"
    shift 2

    zebrac_clear_commands
    zebrac_add_command "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$raw_json" "$8"
    zebrac_run_current_group "$raw_json" "$metadata_jsonl"
    zebrac_clear_commands
}

quote_arg() {
    printf '%q' "$1"
}

shell_command() {
    local script="$1"
    printf 'bash -c %q' "$script"
}
