#!/usr/bin/env bash
# Validate benchmark runner: correctness, then zebrac time and memory, then REPORT.md.
#
# Usage:
#   bash bench/validate/run.sh [options]
#
# Sampling:
#   --runs 5              measured samples per command
#   --warmup 3            unmeasured warmups before sampling
#   --duration 5000       sampling budget per command, capped by --runs
#
# Correctness runs `z-fasta validate --json --summary` on proteome copies with one injected
# defect each and on the shared messy layouts, and compares the exit status and finding kinds
# with values derived from how each input is built. It checks that validate accepts the real
# and messy files that `index` accepts, and compares verdicts with samtools faidx (line layout)
# and SeqKit seq -v (alphabet) where their rules match validate's.
#
# Common invocations:
#   bash bench/validate/run.sh
#   bash bench/validate/run.sh --skip-benchmarks --skip-report   # correctness only
#   bash bench/validate/run.sh --skip-tests --skip-report         # performance only
#
# Skip flags:
#   --skip-tests --skip-benchmarks --skip-messy --skip-report
#   --skip-messy is accepted for suite parity; this suite has no messy performance section.
#   --skip-verify and --skip-perf are aliases.
#
# Outputs:
#   results/run_<ts>.json          run settings, tool versions, and sections
#   results/metadata_<ts>.jsonl    command and workload metadata
#   results/perf_<ts>/             raw zebrac JSON
#   results/tests_<ts>.tsv         correctness checks
#   results/verify_<ts>.log        correctness log
#   results/LATEST                 newest run
#   REPORT.md                      benchmark report
#
# Work files live under bench/validate/data/ (ignored) and are overwritten by each run.
#
#   -h|--help  print this header

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH_ROOT="$(dirname "$SCRIPT_DIR")"
RESULTS_DIR="$SCRIPT_DIR/results"
DATA_DIR="$BENCH_ROOT/shared/data"
WORK_DIR="$SCRIPT_DIR/data/work"
CRAFTED_DIR="$SCRIPT_DIR/data/crafted"
FIXTURE_DIR="$BENCH_ROOT/shared/cache/messy_fixtures"
MESSY_DIR="$BENCH_ROOT/shared/cache/messy_perf"

source "$BENCH_ROOT/shared/runner_common.sh"

RUNS=5
WARMUP=3
ZEBRAC_DURATION_MS="${ZEBRAC_DURATION_MS:-5000}"
DO_TESTS=true
DO_BENCHMARKS=true
DO_REPORT=true

while [[ $# -gt 0 ]]; do
    case $1 in
        --runs) RUNS="$2"; shift 2 ;;
        --warmup) WARMUP="$2"; shift 2 ;;
        --duration) ZEBRAC_DURATION_MS="$2"; shift 2 ;;
        --skip-tests|--skip-verify) DO_TESTS=false; shift ;;
        --skip-benchmarks|--skip-perf) DO_BENCHMARKS=false; shift ;;
        --skip-messy) shift ;;
        --skip-report) DO_REPORT=false; shift ;;
        -h|--help)
            sed -n '2,/^set -euo pipefail$/p' "$0" | head -n -1
            exit 0
            ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

PYTHON="$(report_python)"
LAYOUT_KINDS=",inconsistent_line_widths,trailing_whitespace,mixed_line_endings,"
ERROR_KINDS=",no_sequences,duplicate_name,invalid_character,null_byte,"
FIXTURES=(uniform mixed_widths trailing_whitespace blank_lines mixed_crlf)
VARIANTS=(mixed_widths trailing_whitespace mixed_crlf all_messy)
DATASETS=(Genome Transcriptome Proteome)

PASS=0
FAIL=0
declare -A SUMMARY=()

record() {
    local group="$1" case="$2" check="$3" expected="$4" observed="$5" status
    if [[ "$expected" == "$observed" ]]; then
        PASS=$((PASS + 1))
        status=pass
        printf '  PASS: [%s] %s %s: %s\n' "$group" "$case" "$check" "$observed"
    else
        FAIL=$((FAIL + 1))
        status=fail
        printf '  FAIL: [%s] %s %s: expected %s, got %s\n' "$group" "$case" "$check" "$expected" "$observed"
    fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$group" "$case" "$check" "$expected" "$observed" "$status" >> "$TESTS_TSV"
}

observe() {
    local case="$1" check="$2" first="$3" second="$4"
    printf '  INFO: [difference] %s %s: %s, %s\n' "$case" "$check" "$first" "$second"
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' difference "$case" "$check" "$first" "$second" info >> "$TESTS_TSV"
}

section_hdr() { echo ""; echo "--- [$1] $2 ---"; }

# Prints "<exit>\t<truncated>\t<kinds>" for one validate run; kinds are sorted, or "none".
validate_summary() {
    local file="$1"; shift
    local json status=0
    json="$("$ZFASTA" validate --json --summary "$@" "$file" 2>/dev/null)" || status=$?
    printf '%s' "$json" | "$PYTHON" -c '
import json, sys
text = sys.stdin.read()
data = json.loads(text) if text.strip() else None
if data is None:
    print(sys.argv[1], "unknown", "no-report", sep="\t")
else:
    kinds = sorted(kind for kind, count in data["counts"].items() if count)
    print(sys.argv[1], str(data["truncated"]).lower(), ",".join(kinds) or "none", sep="\t")
' "$status"
}

load_summary() {
    SUMMARY["$1"]="$(validate_summary "$1")"
}

cached_summary() {
    if [[ -z "${SUMMARY[$1]:-}" ]]; then
        echo "error: no validate summary loaded for $1" >&2
        return 1
    fi
    printf '%s\n' "${SUMMARY[$1]}"
}

summary_field() {
    local line="$1" field="$2"
    cut -f "$field" <<<"$line"
}

has_kind() {
    local kinds=",$1," list="$2" kind
    for kind in ${kinds//,/ }; do
        [[ "$list" == *",$kind,"* ]] && return 0
    done
    return 1
}

check_findings() {
    local group="$1" case="$2" file="$3" exit_code="$4" kinds="$5"; shift 5
    local line
    if [[ $# -eq 0 ]]; then line="$(cached_summary "$file")"; else line="$(validate_summary "$file" "$@")"; fi
    record "$group" "$case" findings "exit $exit_code; $kinds" "exit $(summary_field "$line" 1); $(summary_field "$line" 3)"
}

check_capped_findings() {
    local group="$1" case="$2" file="$3" kinds="$4"
    local line truncated exit_code=2
    line="$(cached_summary "$file")"
    truncated="$(summary_field "$line" 2)"
    [[ "$truncated" == true ]] && exit_code=1
    record "$group" "$case" findings "exit $exit_code; $kinds" "exit $(summary_field "$line" 1); $(summary_field "$line" 3)"
}

validate_verdict() {
    if has_kind "$(summary_field "$(cached_summary "$1")" 3)" "$ERROR_KINDS"; then echo reject; else echo accept; fi
}

index_verdict() {
    local file="$1" name="$2"
    ln -sfn "$file" "$WORK_DIR/$name"
    if TMPDIR="$WORK_DIR" "$ZFASTA" index "$WORK_DIR/$name" > /dev/null 2>&1; then echo accept; else echo reject; fi
}

samtools_verdict() {
    if "$SAMTOOLS" faidx --fai-idx "$WORK_DIR/$2.fai" "$1" > /dev/null 2>&1; then echo accept; else echo reject; fi
}

seqkit_verdict() {
    if "$SEQKIT" seq -v -t "$2" "$1" > /dev/null 2>&1; then echo accept; else echo reject; fi
}

layout_verdict() {
    if has_kind "$(summary_field "$(cached_summary "$1")" 3)" "$LAYOUT_KINDS"; then echo reject; else echo accept; fi
}

alphabet_verdict() {
    if has_kind "$(summary_field "$(cached_summary "$1")" 3)" ",invalid_character,"; then echo reject; else echo accept; fi
}

write_crafted() {
    "$PYTHON" - "$DATA_DIR/REAL_Proteome.fasta" "$CRAFTED_DIR" <<'PY'
import sys
from pathlib import Path

source = Path(sys.argv[1]).read_bytes()
out = Path(sys.argv[2])
assert source.startswith(b">") and source.endswith(b"\n")
starts = [0]
at = source.find(b"\n>")
while at != -1:
    starts.append(at + 1)
    at = source.find(b"\n>", at + 1)
records = [source[start:end] for start, end in zip(starts, starts[1:] + [len(source)])]
assert len(records) > 600


def edit_line(record: bytes, line_index: int, edit) -> bytes:
    lines = record.split(b"\n")
    lines[line_index] = edit(lines[line_index])
    return b"\n".join(lines)


def replace(index: int, record: bytes) -> bytes:
    return b"".join(records[:index] + [record] + records[index + 1:])


cases = {
    "duplicate_name": source + records[0],
    "invalid_character": replace(100, edit_line(records[100], 1, lambda line: line[:5] + b"1" + line[6:])),
    "null_byte": replace(200, edit_line(records[200], 0, lambda line: line[:4] + b"\x00" + line[4:])),
    "utf8_bom": b"\xef\xbb\xbf" + source,
    "empty_sequence": b"".join(records[:300] + [b">empty_record no bases\n"] + records[300:]),
    "missing_terminal_newline": source[:-1],
    "interior_space": replace(400, edit_line(records[400], 1, lambda line: line[:10] + b" " + line[10:])),
    "long_header": replace(500, edit_line(records[500], 0, lambda line: line + b" " + b"x" * (2000 - len(line) - 1))),
}
context = b"ACGT" * 30
for name, middle in (("rna_u", b"U"), ("dot_gap", b"."), ("space_inside", b" ")):
    cases[f"difference_{name}"] = b">d\n" + context + middle + context + b"\n"
cases["difference_tab_at_end"] = b">d\n" + context + b"\t\n" + context + b"\n"
out.mkdir(parents=True, exist_ok=True)
for name, data in cases.items():
    (out / f"{name}.fasta").write_bytes(data)
PY
}

run_tests() {
    PASS=0
    FAIL=0
    bench_require_tool z-fasta
    bench_require_tool samtools
    bench_require_tool seqkit
    mkdir -p "$WORK_DIR" "$CRAFTED_DIR"
    printf 'group\tcase\tcheck\texpected\tobserved\tstatus\n' > "$TESTS_TSV"

    echo "z-fasta validate verification"
    echo "  z-fasta:  $ZFASTA"
    echo "  samtools: $SAMTOOLS"
    echo "  seqkit:   $SEQKIT"

    ensure_real_data
    write_crafted

    local c="$CRAFTED_DIR" ds name file
    for name in duplicate_name invalid_character null_byte utf8_bom empty_sequence \
        missing_terminal_newline interior_space long_header difference_rna_u difference_dot_gap \
        difference_space_inside difference_tab_at_end; do
        load_summary "$c/$name.fasta"
    done
    for name in "${FIXTURES[@]}"; do load_summary "$FIXTURE_DIR/$name.fasta"; done
    for name in "${VARIANTS[@]}"; do load_summary "$MESSY_DIR/$name.fasta"; done
    for ds in "${DATASETS[@]}"; do
        file="${REAL_DATASETS[$ds]:-}"
        [[ -n "$file" ]] && load_summary "$file"
    done

    section_hdr crafted "proteome copies with one injected defect each"
    check_findings crafted duplicate_name "$c/duplicate_name.fasta" 1 duplicate_name
    check_findings crafted invalid_character "$c/invalid_character.fasta" 1 invalid_character
    check_findings crafted null_byte "$c/null_byte.fasta" 1 null_byte
    check_findings crafted utf8_bom "$c/utf8_bom.fasta" 2 utf8_bom
    check_findings crafted empty_sequence "$c/empty_sequence.fasta" 2 empty_sequence
    check_findings crafted missing_terminal_newline "$c/missing_terminal_newline.fasta" 2 missing_terminal_newline
    check_findings crafted interior_space "$c/interior_space.fasta" 1 invalid_character
    check_findings crafted long_header "$c/long_header.fasta" 2 long_header
    check_findings crafted strict "$c/empty_sequence.fasta" 1 empty_sequence --strict

    section_hdr messy "shared messy fixtures and proteome variants"
    check_findings fixture uniform "$FIXTURE_DIR/uniform.fasta" 0 none
    check_findings fixture mixed_widths "$FIXTURE_DIR/mixed_widths.fasta" 2 inconsistent_line_widths
    check_findings fixture trailing_whitespace "$FIXTURE_DIR/trailing_whitespace.fasta" 2 trailing_whitespace
    check_findings fixture blank_lines "$FIXTURE_DIR/blank_lines.fasta" 2 inconsistent_line_widths
    check_findings fixture mixed_crlf "$FIXTURE_DIR/mixed_crlf.fasta" 2 inconsistent_line_widths,mixed_line_endings
    check_capped_findings variant mixed_widths "$MESSY_DIR/mixed_widths.fasta" inconsistent_line_widths
    check_capped_findings variant trailing_whitespace "$MESSY_DIR/trailing_whitespace.fasta" trailing_whitespace
    check_capped_findings variant mixed_crlf "$MESSY_DIR/mixed_crlf.fasta" mixed_line_endings
    check_capped_findings variant all_messy "$MESSY_DIR/all_messy.fasta" inconsistent_line_widths,mixed_line_endings,trailing_whitespace

    section_hdr index "validate accepts the files that index accepts"
    for ds in "${DATASETS[@]}"; do
        file="${REAL_DATASETS[$ds]:-}"
        [[ -n "$file" ]] || continue
        record real "$ds" "accepted (index, validate)" "$(index_verdict "$file" "real_$ds.fa")" "$(validate_verdict "$file")"
    done
    for name in "${FIXTURES[@]}"; do
        file="$FIXTURE_DIR/$name.fasta"
        record fixture "$name" "accepted (index, validate)" "$(index_verdict "$file" "fixture_$name.fasta")" "$(validate_verdict "$file")"
    done
    for name in "${VARIANTS[@]}"; do
        file="$MESSY_DIR/$name.fasta"
        record variant "$name" "accepted (index, validate)" "$(index_verdict "$file" "variant_$name.fasta")" "$(validate_verdict "$file")"
    done

    section_hdr samtools "samtools faidx rejects exactly the layouts validate warns about"
    for ds in "${DATASETS[@]}"; do
        file="${REAL_DATASETS[$ds]:-}"
        [[ -n "$file" ]] || continue
        record real "$ds" "layout (samtools, validate)" "$(samtools_verdict "$file" "real_$ds")" "$(layout_verdict "$file")"
    done
    for name in "${FIXTURES[@]}"; do
        file="$FIXTURE_DIR/$name.fasta"
        record fixture "$name" "layout (samtools, validate)" "$(samtools_verdict "$file" "fixture_$name")" "$(layout_verdict "$file")"
    done
    for name in "${VARIANTS[@]}"; do
        file="$MESSY_DIR/$name.fasta"
        record variant "$name" "layout (samtools, validate)" "$(samtools_verdict "$file" "variant_$name")" "$(layout_verdict "$file")"
    done

    section_hdr seqkit "SeqKit seq -v and validate agree where their alphabets match"
    local type
    for ds in "${DATASETS[@]}"; do
        file="${REAL_DATASETS[$ds]:-}"
        [[ -n "$file" ]] || continue
        type=dna
        [[ "$ds" == Proteome ]] && type=protein
        record real "$ds" "alphabet (seqkit, validate)" "$(seqkit_verdict "$file" "$type")" "$(alphabet_verdict "$file")"
    done
    record crafted invalid_character "alphabet (seqkit, validate)" "$(seqkit_verdict "$c/invalid_character.fasta" protein)" "$(alphabet_verdict "$c/invalid_character.fasta")"
    for name in rna_u dot_gap space_inside tab_at_end; do
        file="$c/difference_$name.fasta"
        observe "$name" "alphabet (seqkit, validate)" "$(seqkit_verdict "$file" dna)" "$(alphabet_verdict "$file")"
    done

    echo ""
    echo "Results: $PASS passed, $FAIL failed"
    if [[ "$FAIL" -gt 0 ]]; then
        echo "VERIFICATION FAILED"
        return 1
    fi
    echo "ALL PASSED"
}

run_zebrac_tool() {
    local workload="$1" tool="$2" family="$3" json_out="$4" script="$5" input_bytes="$6"
    echo "  >> perf $workload $tool  (runs=$RUNS warmup=$WARMUP)"
    zebrac_clear_commands
    zebrac_add_command validate perf "$workload" "$tool" "$family" "$input_bytes" "" "$json_out" "$(shell_command "$script")"
    if ! bench_group "$json_out"; then
        echo "error: zebrac failed for $workload $tool" >&2
        echo "  command: $script" >&2
        return 1
    fi
}

write_run_manifest() {
    BENCH_ZFASTA_BUILD="$(bench_build_identity)"
    export BENCH_ZFASTA_BUILD
    python3 - "$RUN_MANIFEST" "$TIMESTAMP" "$(basename "$METADATA_JSONL")" "$RUNS" "$WARMUP" \
        "$ZEBRAC_DURATION_MS" "$DO_TESTS" "perf_${TIMESTAMP}" "tests_${TIMESTAMP}.tsv" <<'PY'
import json, os, sys
from pathlib import Path

manifest, ts, metadata, runs, warmup, duration, did_verify, perf, tests = sys.argv[1:10]
tools = {}
for name in ("samtools", "seqkit"):
    value = os.environ.get(f"BENCH_VER_{name.upper()}")
    if value:
        tools[name] = value
sections = {"perf": perf}
if (Path(manifest).parent / tests).is_file():
    sections["correctness"] = tests
out = {
    "schema_version": "validate-run.v1",
    "timestamp": ts,
    "runner": "zebrac",
    "zebrac": os.environ.get("BENCH_VER_ZEBRAC", ""),
    "z_fasta": os.environ.get("BENCH_VER_ZFASTA", ""),
    "z_fasta_build": json.loads(os.environ["BENCH_ZFASTA_BUILD"]),
    "runs": int(runs),
    "warmup": int(warmup),
    "duration_ms": int(duration),
    "verify_skipped": did_verify != "true",
    "verify_pass": os.environ.get("BENCH_VERIFY_PASS"),
    "metadata": metadata,
    "tools": tools,
    "sections": sections,
}
Path(manifest).write_text(json.dumps(out, indent=2) + "\n")
PY
}

run_perf() {
    bench_require_tool zebrac
    bench_require_tool z-fasta
    bench_require_tool samtools
    bench_require_tool seqkit
    METADATA_JSONL="$RESULTS_DIR/metadata_${TIMESTAMP}.jsonl"
    RUN_MANIFEST="$RESULTS_DIR/run_${TIMESTAMP}.json"
    : > "$METADATA_JSONL"
    ZEBRAC_MIN_SAMPLES="$RUNS"
    ZEBRAC_MAX_SAMPLES="$RUNS"
    ZEBRAC_WARMUP="$WARMUP"
    local out_dir="$RESULTS_DIR/perf_${TIMESTAMP}"
    mkdir -p "$out_dir" "$WORK_DIR"

    echo ""
    echo "================================================================"
    echo "  Validate zebrac performance"
    echo "================================================================"
    echo "  Runs: $RUNS | Warmup: $WARMUP | Duration: ${ZEBRAC_DURATION_MS}ms"
    echo "  metadata: $METADATA_JSONL"

    ensure_real_data
    local ds fa bytes type qf qz qk qs qw
    qz="$(quote_arg "$ZFASTA")"
    qk="$(quote_arg "$SEQKIT")"
    qs="$(quote_arg "$SAMTOOLS")"
    for ds in "${DATASETS[@]}"; do
        fa="${REAL_DATASETS[$ds]:-}"
        [[ -n "$fa" ]] || continue
        bytes="$(file_size_bytes "$fa")"
        type=dna
        [[ "$ds" == Proteome ]] && type=protein
        qf="$(quote_arg "$fa")"
        qw="$(quote_arg "$WORK_DIR/perf_$ds.fai")"
        echo "  -- dataset $ds ($(du -h "$fa" | cut -f1)) --"
        run_zebrac_tool "$ds" z-fasta z-fasta "$out_dir/${ds}__z-fasta.json" "$qz validate --json --summary $qf > /dev/null" "$bytes"
        run_zebrac_tool "$ds" seqkit peer "$out_dir/${ds}__seqkit.json" "$qk seq -j 1 -v -t $type -o /dev/null $qf" "$bytes"
        run_zebrac_tool "$ds" samtools peer "$out_dir/${ds}__samtools.json" "$qs faidx --fai-idx $qw $qf" "$bytes"
    done

    export_manifest_core_versions
    export_manifest_optional_tool_versions seqkit
    write_run_manifest
    printf '%s\n' "$TIMESTAMP" > "$RESULTS_DIR/LATEST"
    echo "  manifest: $RUN_MANIFEST"
}

run_report() {
    echo ""
    echo "================================================================"
    echo "  Generating REPORT.md"
    echo "================================================================"
    "$PYTHON" "$SCRIPT_DIR/generate_report.py"
}

TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
TESTS_TSV="$RESULTS_DIR/tests_${TIMESTAMP}.tsv"
mkdir -p "$RESULTS_DIR"

if $DO_TESTS; then
    echo ""
    echo "================================================================"
    echo "  run_tests (correctness)"
    echo "================================================================"
    bench_ensure_messy --fixtures --perf
    verify_log="$RESULTS_DIR/verify_${TIMESTAMP}.log"
    if run_tests | tee "$verify_log"; then
        BENCH_VERIFY_PASS="$(awk '/^Results:/{print $2; exit}' "$verify_log")"
        export BENCH_VERIFY_PASS
    else
        echo "error: run_tests failed" >&2
        exit 1
    fi
fi
$DO_BENCHMARKS && run_perf
$DO_REPORT && run_report

echo ""
echo "================================================================"
echo "  Validate suite complete"
echo "================================================================"
