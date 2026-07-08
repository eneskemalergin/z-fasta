#!/usr/bin/env bash
# verify.sh - z-fasta stats verification (subject: z-fasta; BioPython + peers are oracles)
#
# Bench layout: one runner (this file), like bench/get/verify.sh. oracle.py is the
# only helper script here; extend it in place, do not add siblings.
#
# Tags: [oracle:biopython] [index:zfi] [index:fai] [index:cross] [index:lowmem]
#       [extended:messy] [extended:dedup] [parity:seqkit] [parity:seqtk]
#       [parity:noodles] [parity:rustbio]
#
# Usage: ./verify.sh [--skip-tools]
set -euo pipefail

SKIP_TOOLS=false

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
ORACLE="$SCRIPT_DIR/oracle.py"
TEST_DATA="$PROJECT_DIR/tests/data"
MESSY_DIR="$PROJECT_DIR/bench/index/messy_variants"
TMPDIR="$SCRIPT_DIR/.verify_tmp"
mkdir -p "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT

ZFASTA="${ZFASTA:-$PROJECT_DIR/zig-out/bin/z-fasta}"
SAMTOOLS="${SAMTOOLS:-samtools}"
SEQKIT="${SEQKIT:-$PROJECT_DIR/tools/seqkit}"
SEQTK="${SEQTK:-$PROJECT_DIR/tools/seqtk/seqtk}"
NOODLES="${NOODLES:-$PROJECT_DIR/tools/noodles_wrapper/target/release/noodles_wrapper}"
RUSTBIO="${RUSTBIO:-$PROJECT_DIR/tools/rustbio_wrapper/target/release/rustbio_wrapper}"
PYTHON="${PYTHON:-$PROJECT_DIR/.venv/bin/python}"
[[ -x "$PYTHON" ]] || PYTHON="$(command -v python3)"

FIXTURES=(simple proteome single edge_cases mixed_widths)
MESSY_FIXTURES=(mixed_line_widths trailing_whitespace blank_lines mixed_crlf_lf)
# seqkit counts every raw FASTA row; skip edge_cases where index filters empties/dups.
SEQKIT_FIXTURES=(simple proteome single mixed_widths)
# seqtk comp is nucleotide-only; proteome is protein.
SEQTK_FIXTURES=(simple single mixed_widths)

PASS=0 FAIL=0

pass() { PASS=$((PASS + 1)); printf "  PASS: %s\n" "$1"; }

fail() {
    FAIL=$((FAIL + 1)); printf "  FAIL: %s\n" "$1"
    if [[ -s "${2:-}" ]]; then sed 's/^/    /' "$2"; fi
}

section_hdr() { echo ""; echo "--- [$1] $2 ---"; }

oracle() { "$PYTHON" "$ORACLE" "$@"; }

check_oracle() {
    local label="$1"; shift
    local slug="${label//[^a-zA-Z0-9_]/_}"
    local err="$TMPDIR/oracle_${slug}.err"
    if oracle "$@" 2>"$err"; then
        pass "$label"
    else
        fail "$label" "$err"
    fi
}

ensure_index() {
    local f="$1"
    if [[ ! -f "${f}.fai" ]]; then
        if ! "$SAMTOOLS" faidx "$f" 2>"$TMPDIR/faidx.err"; then
            echo "samtools faidx failed for ${f}" >>"$TMPDIR/faidx.err"
            return 1
        fi
    fi
    [[ -f "${f}.zfi" ]] || "$ZFASTA" index "$f" 2>/dev/null
}

prepare_messy() {
    cp "$MESSY_DIR/${1}.fasta" "$2"
    "$ZFASTA" index "$2" >/dev/null 2>&1
}

run_stats() {
    local fasta="$1" mode="$2" out="$3" err="$4"
    local -a args=(stats)
    [[ "$mode" == index-only ]] && args+=(--index-only)
    args+=("$fasta")
    if ! "$ZFASTA" "${args[@]}" >"$out" 2>"$err"; then
        echo "z-fasta stats failed (see stderr)" >"$err"
        return 1
    fi
}

prepare_fixtures() {
    local name fasta
    for name in "${FIXTURES[@]}"; do
        fasta="$TMPDIR/$name.fasta"
        cp "$TEST_DATA/$name.fasta" "$fasta"
        [[ -f "$TEST_DATA/$name.fasta.fai" ]] && cp "$TEST_DATA/$name.fasta.fai" "${fasta}.fai"
        ensure_index "$fasta" || { fail "[oracle] index setup failed for $name" "$TMPDIR/faidx.err"; exit 1; }
        oracle expected "$fasta" >"$TMPDIR/$name.expected.json"
    done
}

verify_modes() {
    local tag="$1" name="$2" fasta="$3" exp="$4" zfi_only="${5:-}"
    local mode out err stash label

    section_hdr "$tag" "$name"

    for mode in full index-only; do
        out="$TMPDIR/$name.$mode.zfi.txt"
        err="$TMPDIR/$name.$mode.zfi.err"
        if [[ "$zfi_only" == zfi-only ]]; then
            label="[oracle:$name:$mode:zfi] vs BioPython"
            run_stats "$fasta" "$mode" "$out" "$err" || { fail "$label" "$err"; continue; }
            check_oracle "$label" check "$mode" zfi "$fasta" "$exp" "$out"
            continue
        fi

        label="[oracle:$name:$mode] zfi+fai vs BioPython + cross"
        run_stats "$fasta" "$mode" "$out" "$err" || { fail "$label" "$err"; continue; }

        stash="$TMPDIR/$name.$mode.zfi"
        mv "${fasta}.zfi" "$stash"
        out="$TMPDIR/$name.$mode.fai.txt"
        err="$TMPDIR/$name.$mode.fai.err"
        run_stats "$fasta" "$mode" "$out" "$err" || { fail "$label" "$err"; mv "$stash" "${fasta}.zfi"; continue; }
        mv "$stash" "${fasta}.zfi"

        check_oracle "$label" verify-mode "$mode" "$fasta" "$exp" \
            "$TMPDIR/$name.$mode.zfi.txt" "$TMPDIR/$name.$mode.fai.txt"
    done
}

verify_parity() {
    local name="$1" fasta="$2" exp="$3"
    local -a args=("$fasta" "$exp" "$TMPDIR/$name.index-only.zfi.txt" "$TMPDIR/$name.full.zfi.txt")
    local tool bin pair sk err label tool_count=0

    label="[parity:$name] external tools"
    err="$TMPDIR/$name.parity.err"

    for pair in "noodles:$NOODLES" "rustbio:$RUSTBIO"; do
        tool="${pair%%:*}"
        bin="${pair#*:}"
        [[ -x "$bin" ]] || continue
        if "$bin" stats "$fasta" >"$TMPDIR/$name.$tool.stats" 2>"$TMPDIR/$name.$tool.err"; then
            args+=("$tool" "$TMPDIR/$name.$tool.stats")
            tool_count=$((tool_count + 1))
        else
            fail "$label ($tool run)" "$TMPDIR/$name.$tool.err"
            return
        fi
    done

    for sk in "${SEQKIT_FIXTURES[@]}"; do
        [[ "$sk" == "$name" ]] || continue
        if "$SEQKIT" stats -a "$fasta" >"$TMPDIR/$name.seqkit.txt" 2>"$TMPDIR/$name.seqkit.err"; then
            args+=("seqkit" "$TMPDIR/$name.seqkit.txt")
            tool_count=$((tool_count + 1))
        else
            fail "$label (seqkit run)" "$TMPDIR/$name.seqkit.err"
            return
        fi
    done

    for sk in "${SEQTK_FIXTURES[@]}"; do
        [[ "$sk" == "$name" ]] || continue
        if "$SEQTK" comp "$fasta" >"$TMPDIR/$name.seqtk.txt" 2>"$TMPDIR/$name.seqtk.err"; then
            args+=("seqtk" "$TMPDIR/$name.seqtk.txt")
            tool_count=$((tool_count + 1))
        else
            fail "$label (seqtk run)" "$TMPDIR/$name.seqtk.err"
            return
        fi
    done

    if [[ "$tool_count" -eq 0 ]]; then
        echo "no external tools available for $name" >"$err"
        fail "$label" "$err"
        return
    fi

    if oracle parity "${args[@]}" 2>"$err"; then
        pass "$label"
    else
        fail "$label" "$err"
    fi
}

verify_file() {
    local name="$1"
    verify_modes oracle "$name" "$TMPDIR/$name.fasta" "$TMPDIR/$name.expected.json"
    if ! $SKIP_TOOLS; then
        verify_parity "$name" "$TMPDIR/$name.fasta" "$TMPDIR/$name.expected.json"
    fi
}

verify_messy() {
    # No external-tool parity: wrappers count raw bytes; z-fasta strips whitespace.
    local name fasta exp
    for name in "${MESSY_FIXTURES[@]}"; do
        fasta="$TMPDIR/messy_${name}.fasta"
        [[ -f "$MESSY_DIR/${name}.fasta" ]] || { fail "[extended:messy] missing fixture $name"; continue; }
        prepare_messy "$name" "$fasta" || { fail "[extended:messy] index failed $name"; continue; }
        oracle expected "$fasta" >"$TMPDIR/messy_${name}.expected.json"
        verify_modes "extended:messy" "messy_$name" "$fasta" "$TMPDIR/messy_${name}.expected.json" zfi-only
    done
}

verify_lowmem_stats() {
    local src="$1" name="$2" fasta="$TMPDIR/lowmem_${name}.fasta"
    local err="$TMPDIR/lowmem_${name}.err"

    cp "$src" "$fasta"
    rm -f "${fasta}.zfi" "${fasta}.fai"
    "$ZFASTA" index "$fasta" >/dev/null 2>&1
    run_stats "$fasta" full "$TMPDIR/lowmem_${name}.mmap.full.txt" "$err" \
        || { fail "[index:lowmem] $name full stats (mmap)" "$err"; return; }
    run_stats "$fasta" index-only "$TMPDIR/lowmem_${name}.mmap.idx.txt" "$err" \
        || { fail "[index:lowmem] $name index-only stats (mmap)" "$err"; return; }

    rm -f "${fasta}.zfi"
    "$ZFASTA" index --low-mem "$fasta" >/dev/null 2>&1
    run_stats "$fasta" full "$TMPDIR/lowmem_${name}.low.full.txt" "$err" \
        || { fail "[index:lowmem] $name full stats (--low-mem)" "$err"; return; }
    run_stats "$fasta" index-only "$TMPDIR/lowmem_${name}.low.idx.txt" "$err" \
        || { fail "[index:lowmem] $name index-only stats (--low-mem)" "$err"; return; }

    check_oracle "[index:lowmem] $name full mmap == --low-mem" same full \
        "$TMPDIR/lowmem_${name}.mmap.full.txt" "$TMPDIR/lowmem_${name}.low.full.txt"
    check_oracle "[index:lowmem] $name index-only mmap == --low-mem" same index-only \
        "$TMPDIR/lowmem_${name}.mmap.idx.txt" "$TMPDIR/lowmem_${name}.low.idx.txt"
}

verify_dedup_stats() {
    local dedup_fasta="$TMPDIR/dedup_default.fasta"
    local nodedup_fasta="$TMPDIR/dedup_nodedup.fasta"
    local exp_dedup="$TMPDIR/dedup_default.expected.json"
    local exp_nodedup="$TMPDIR/nodedup.expected.json"
    local label="[extended:dedup] --no-dedup indexes more sequences than default"
    local err dedup_seqs nodedup_seqs

    cp "$TEST_DATA/edge_cases.fasta" "$dedup_fasta"
    rm -f "${dedup_fasta}.zfi" "${dedup_fasta}.fai"
    ensure_index "$dedup_fasta" || { fail "[extended:dedup] default index setup" "$TMPDIR/faidx.err"; return; }
    oracle expected "$dedup_fasta" >"$exp_dedup"
    err="$TMPDIR/dedup_default.err"
    run_stats "$dedup_fasta" index-only "$TMPDIR/dedup_default.index-only.txt" "$err" \
        || { fail "[extended:dedup] default index-only stats" "$err"; return; }
    check_oracle "[extended:dedup] default vs BioPython" check index-only zfi \
        "$dedup_fasta" "$exp_dedup" "$TMPDIR/dedup_default.index-only.txt"

    cp "$TEST_DATA/edge_cases.fasta" "$nodedup_fasta"
    rm -f "${nodedup_fasta}.zfi" "${nodedup_fasta}.fai"
    "$ZFASTA" index --no-dedup "$nodedup_fasta" >/dev/null 2>&1 \
        || { fail "[extended:dedup] index --no-dedup"; return; }
    oracle expected "$nodedup_fasta" --no-dedup >"$exp_nodedup"
    verify_modes "extended:dedup" "edge_cases_nodedup" "$nodedup_fasta" "$exp_nodedup" zfi-only

    dedup_seqs="$("$PYTHON" -c "import json; print(json.load(open('$exp_dedup'))['num_seqs'])")"
    nodedup_seqs="$("$PYTHON" -c "import json; print(json.load(open('$exp_nodedup'))['num_seqs'])")"
    if [[ "$nodedup_seqs" -gt "$dedup_seqs" ]]; then
        pass "$label"
    else
        echo "dedup=$dedup_seqs nodedup=$nodedup_seqs" >"$TMPDIR/dedup_cmp.err"
        fail "$label" "$TMPDIR/dedup_cmp.err"
    fi
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --skip-tools) SKIP_TOOLS=true ;;
        *) echo "Unknown: $1"; exit 1 ;;
    esac
    shift
done

echo "z-fasta stats verification"
echo "  z-fasta:  $ZFASTA"
echo "  skip tools: $SKIP_TOOLS"

cd "$PROJECT_DIR"
[[ -x "$ZFASTA" ]] || { echo "Error: run ./zig build first"; exit 1; }
command -v "$SAMTOOLS" &>/dev/null || { echo "Error: samtools not found"; exit 1; }
"$PYTHON" -c "from Bio import SeqIO" 2>/dev/null || { echo "Error: BioPython not available"; exit 1; }

prepare_fixtures
for name in "${FIXTURES[@]}"; do
    verify_file "$name"
done

verify_messy

section_hdr "index" "low-mem stats parity"
verify_lowmem_stats "$TEST_DATA/mixed_widths.fasta" mixed_widths
verify_lowmem_stats "$MESSY_DIR/mixed_line_widths.fasta" messy_mixed_line_widths

verify_dedup_stats

echo ""
echo "Results: $PASS passed, $FAIL failed"
[[ "$FAIL" -gt 0 ]] && { echo "VERIFICATION FAILED"; exit 1; }
echo "ALL PASSED"
