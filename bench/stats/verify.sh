#!/usr/bin/env bash
# verify.sh - z-fasta stats verification (subject: z-fasta; BioPython + peers are oracles)
#
# Bench layout: one runner (this file), like bench/get/verify.sh. oracle.py is the
# only helper script here; extend it in place, do not add siblings.
#
# Tags: [oracle:biopython] [index:zfi] [index:fai] [index:cross] [index:lowmem]
#       [extended:messy] [extended:dedup] [extended:edge]
#       [parity:seqkit] [parity:seqtk] [parity:noodles] [parity:rustbio]
#
# Usage: ./verify.sh [--skip-tools] [--skip-messy] [--skip-lowmem] [--skip-dedup] [--skip-edge] [--skip-layout]
set -euo pipefail

SKIP_TOOLS=false
SKIP_MESSY=false
SKIP_LOWMEM=false
SKIP_DEDUP=false
SKIP_EDGE=false
SKIP_LAYOUT=false

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
ORACLE="$SCRIPT_DIR/oracle.py"
TEST_DATA="$PROJECT_DIR/tests/data"
MESSY_DIR="$PROJECT_DIR/bench/index/messy_variants"
LAYOUT_TWINS="$SCRIPT_DIR/fixtures/layout_twins"
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
# Side-table / whitespace messy: z-fasta .zfi only (peers count raw bytes).
MESSY_FIXTURES=(mixed_line_widths trailing_whitespace blank_lines mixed_crlf_lf)
# Uniform messy control: samtools .fai works; exercise zfi+fai cross like get.
UNIFORM_MESSY=uniform_control
# seqkit counts every raw FASTA row; skip edge_cases where index filters empties/dups.
SEQKIT_FIXTURES=(simple proteome single mixed_widths)
# seqtk comp is nucleotide-only; proteome is protein.
SEQTK_FIXTURES=(simple single mixed_widths)
# noodles/rustbio wrappers: clean FASTA comparison peers only (no messy / side-table).
# Richer TSV fields exist so we can compare assembly+composition; they do not gain messy support.
WRAPPER_FIXTURES=(simple proteome single edge_cases mixed_widths)
# Same bases, different wrapping: uniform reference vs messy layouts (z-fasta only).
LAYOUT_TWIN_VARIANTS=(mixed_widths trailing_ws blank_lines mixed_crlf)

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

expect_fail() {
    local label="$1"; shift
    local err="$TMPDIR/edge_${label//[^a-zA-Z0-9_]/_}.err"
    local out="$TMPDIR/edge_${label//[^a-zA-Z0-9_]/_}.out"
    if "$@" >"$out" 2>"$err"; then
        echo "expected non-zero exit" >"$err"
        fail "$label" "$err"
    else
        pass "$label"
    fi
}

expect_fail_msg() {
    local label="$1" needle="$2"; shift 2
    local err="$TMPDIR/edge_${label//[^a-zA-Z0-9_]/_}.err"
    local out="$TMPDIR/edge_${label//[^a-zA-Z0-9_]/_}.out"
    if "$@" >"$out" 2>"$err"; then
        echo "expected non-zero exit" >"$err"
        fail "$label" "$err"
        return
    fi
    if grep -q "$needle" "$err"; then
        pass "$label"
    else
        echo "stderr missing: $needle" >>"$err"
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
    local mode out err stash label zfi_txt fai_txt

    section_hdr "$tag" "$name"

    for mode in full index-only; do
        zfi_txt="$TMPDIR/$name.$mode.zfi.txt"
        err="$TMPDIR/$name.$mode.zfi.err"
        if [[ "$zfi_only" == zfi-only ]]; then
            label="[oracle:$name:$mode:zfi] vs BioPython"
            run_stats "$fasta" "$mode" "$zfi_txt" "$err" || { fail "$label" "$err"; continue; }
            check_oracle "$label" check "$mode" zfi "$fasta" "$exp" "$zfi_txt"
            continue
        fi

        # .zfi path vs BioPython
        label="[oracle:$name:$mode:zfi] vs BioPython"
        run_stats "$fasta" "$mode" "$zfi_txt" "$err" || { fail "$label" "$err"; continue; }
        check_oracle "$label" check "$mode" zfi "$fasta" "$exp" "$zfi_txt"

        # .fai fallback vs BioPython (stash .zfi so loader uses .fai)
        stash="$TMPDIR/$name.$mode.zfi"
        mv "${fasta}.zfi" "$stash"
        fai_txt="$TMPDIR/$name.$mode.fai.txt"
        err="$TMPDIR/$name.$mode.fai.err"
        label="[oracle:$name:$mode:fai] vs BioPython"
        if ! run_stats "$fasta" "$mode" "$fai_txt" "$err"; then
            fail "$label" "$err"
            mv "$stash" "${fasta}.zfi"
            continue
        fi
        check_oracle "$label" check "$mode" fai "$fasta" "$exp" "$fai_txt"
        mv "$stash" "${fasta}.zfi"

        # Identity: .zfi stats == .fai stats (same shape as mmap == --low-mem)
        check_oracle "[index:cross] $name $mode .zfi == .fai" same "$mode" "$zfi_txt" "$fai_txt"
    done
}

# One pass per peer tool so failures name the tool (get-style).
# Wrappers: clean fixtures only (WRAPPER_FIXTURES). Never run on messy side-table files.
verify_parity() {
    local name="$1" fasta="$2" exp="$3"
    local idx_txt="$TMPDIR/$name.index-only.zfi.txt"
    local full_txt="$TMPDIR/$name.full.zfi.txt"
    local tool bin out err label
    local ran=0
    local use_wrappers=false
    local sk

    if [[ ! -f "$idx_txt" || ! -f "$full_txt" ]]; then
        err="$TMPDIR/$name.parity_setup.err"
        echo "missing z-fasta stats outputs for $name (run oracle modes first)" >"$err"
        fail "[parity:$name] setup" "$err"
        return
    fi

    for sk in "${WRAPPER_FIXTURES[@]}"; do
        [[ "$sk" == "$name" ]] && use_wrappers=true && break
    done

    if $use_wrappers; then
        for pair in "noodles:$NOODLES" "rustbio:$RUSTBIO"; do
            tool="${pair%%:*}"
            bin="${pair#*:}"
            [[ -x "$bin" ]] || continue
            label="[parity:$tool] $name assembly/composition"
            out="$TMPDIR/$name.$tool.stats"
            err="$TMPDIR/$name.$tool.err"
            if ! "$bin" stats "$fasta" >"$out" 2>"$err"; then
                fail "$label (run)" "$err"
                continue
            fi
            # Stale binary: old wrappers only printed sequences + total_bases.
            if ! grep -q $'^n50\t' "$out" || ! grep -q $'^type\t' "$out"; then
                echo "wrapper output missing n50/type; rebuild tools/${tool}_wrapper (clean-FASTA comparison peer, not messy)" >"$err"
                fail "$label (stale binary)" "$err"
                continue
            fi
            check_oracle "$label" parity "$fasta" "$exp" "$idx_txt" "$full_txt" "$tool" "$out"
            ran=$((ran + 1))
        done
    fi

    for sk in "${SEQKIT_FIXTURES[@]}"; do
        [[ "$sk" == "$name" ]] || continue
        [[ -x "$SEQKIT" ]] || continue
        label="[parity:seqkit] $name assembly stats"
        out="$TMPDIR/$name.seqkit.txt"
        err="$TMPDIR/$name.seqkit.err"
        if ! "$SEQKIT" stats -a -T "$fasta" >"$out" 2>"$err"; then
            fail "$label (run)" "$err"
            break
        fi
        check_oracle "$label" parity "$fasta" "$exp" "$idx_txt" "$full_txt" seqkit "$out"
        ran=$((ran + 1))
        break
    done

    for sk in "${SEQTK_FIXTURES[@]}"; do
        [[ "$sk" == "$name" ]] || continue
        [[ -x "$SEQTK" ]] || continue
        label="[parity:seqtk] $name composition"
        out="$TMPDIR/$name.seqtk.txt"
        err="$TMPDIR/$name.seqtk.err"
        if ! "$SEQTK" comp "$fasta" >"$out" 2>"$err"; then
            fail "$label (run)" "$err"
            break
        fi
        check_oracle "$label" parity "$fasta" "$exp" "$idx_txt" "$full_txt" seqtk "$out"
        ran=$((ran + 1))
        break
    done

    if [[ "$ran" -eq 0 ]]; then
        err="$TMPDIR/$name.parity.err"
        echo "no external tools available for $name" >"$err"
        fail "[parity:$name] external tools" "$err"
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
    local name fasta exp
    for name in "${MESSY_FIXTURES[@]}"; do
        fasta="$TMPDIR/messy_${name}.fasta"
        [[ -f "$MESSY_DIR/${name}.fasta" ]] || { fail "[extended:messy] missing fixture $name"; continue; }
        prepare_messy "$name" "$fasta" || { fail "[extended:messy] index failed $name"; continue; }
        oracle expected "$fasta" >"$TMPDIR/messy_${name}.expected.json"
        verify_modes "extended:messy" "messy_$name" "$fasta" "$TMPDIR/messy_${name}.expected.json" zfi-only
    done

    # Uniform control: peers can index; exercise .zfi and .fai stats cross-compare.
    name="$UNIFORM_MESSY"
    fasta="$TMPDIR/messy_${name}.fasta"
    [[ -f "$MESSY_DIR/${name}.fasta" ]] || { fail "[extended:messy] missing fixture $name"; return; }
    cp "$MESSY_DIR/${name}.fasta" "$fasta"
    ensure_index "$fasta" || { fail "[extended:messy] index failed $name" "$TMPDIR/faidx.err"; return; }
    oracle expected "$fasta" >"$TMPDIR/messy_${name}.expected.json"
    verify_modes "extended:messy" "messy_$name" "$fasta" "$TMPDIR/messy_${name}.expected.json"
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

# Same sequence bytes, different on-disk wrapping: stats must match (layout invariant).
verify_layout_twins() {
    section_hdr "extended:layout" "messy == uniform (same bases)"
    local uniform="$TMPDIR/layout_uniform.fasta"
    local variant name fasta err mode
    local uref_full="$TMPDIR/layout_uniform.full.txt"
    local uref_idx="$TMPDIR/layout_uniform.index-only.txt"

    [[ -f "$LAYOUT_TWINS/uniform.fasta" ]] || {
        fail "[extended:layout] missing fixtures under bench/stats/fixtures/layout_twins"
        return
    }

    cp "$LAYOUT_TWINS/uniform.fasta" "$uniform"
    rm -f "${uniform}.zfi" "${uniform}.fai"
    "$ZFASTA" index "$uniform" >/dev/null 2>&1 \
        || { fail "[extended:layout] index uniform"; return; }
    err="$TMPDIR/layout_uniform.err"
    run_stats "$uniform" full "$uref_full" "$err" \
        || { fail "[extended:layout] uniform full stats" "$err"; return; }
    run_stats "$uniform" index-only "$uref_idx" "$err" \
        || { fail "[extended:layout] uniform index-only stats" "$err"; return; }

    oracle expected "$uniform" >"$TMPDIR/layout_uniform.expected.json"
    check_oracle "[extended:layout] uniform vs BioPython" check full zfi \
        "$uniform" "$TMPDIR/layout_uniform.expected.json" "$uref_full"

    for variant in "${LAYOUT_TWIN_VARIANTS[@]}"; do
        name="layout_$variant"
        fasta="$TMPDIR/${name}.fasta"
        [[ -f "$LAYOUT_TWINS/${variant}.fasta" ]] || {
            fail "[extended:layout] missing $variant fixture"
            continue
        }
        cp "$LAYOUT_TWINS/${variant}.fasta" "$fasta"
        rm -f "${fasta}.zfi" "${fasta}.fai"
        "$ZFASTA" index "$fasta" >/dev/null 2>&1 \
            || { fail "[extended:layout] index $variant"; continue; }

        for mode in full index-only; do
            err="$TMPDIR/${name}.${mode}.err"
            run_stats "$fasta" "$mode" "$TMPDIR/${name}.${mode}.txt" "$err" \
                || { fail "[extended:layout] $variant $mode stats" "$err"; continue; }
            if [[ "$mode" == full ]]; then
                check_oracle "[extended:layout] $variant full == uniform" same full \
                    "$uref_full" "$TMPDIR/${name}.full.txt"
            else
                check_oracle "[extended:layout] $variant index-only == uniform" same index-only \
                    "$uref_idx" "$TMPDIR/${name}.index-only.txt"
            fi
        done

        oracle expected "$fasta" >"$TMPDIR/${name}.expected.json"
        check_oracle "[extended:layout] $variant vs BioPython" check full zfi \
            "$fasta" "$TMPDIR/${name}.expected.json" "$TMPDIR/${name}.full.txt"
    done
}

# Product note: Duplicates line is always 0 (indexed set only). Gate that so it cannot drift silently.
verify_duplicates_policy() {
    section_hdr "extended:dedup" "Duplicates line policy"
    local fasta="$TMPDIR/dup_policy.fasta"
    local out="$TMPDIR/dup_policy.index-only.txt"
    local err="$TMPDIR/dup_policy.err"
    local dups

    cp "$TEST_DATA/edge_cases.fasta" "$fasta"
    rm -f "${fasta}.zfi" "${fasta}.fai"
    "$ZFASTA" index --no-dedup "$fasta" >/dev/null 2>&1 \
        || { fail "[extended:dedup] index --no-dedup for Duplicates policy"; return; }
    run_stats "$fasta" index-only "$out" "$err" \
        || { fail "[extended:dedup] stats on --no-dedup index" "$err"; return; }
    dups="$(awk '/^Duplicates:/{print $2; exit}' "$out")"
    if [[ "$dups" == "0" ]]; then
        pass "[extended:dedup] Duplicates prints 0 on --no-dedup index (indexed set only)"
    else
        echo "expected Duplicates: 0, got ${dups:-missing}" >"$err"
        fail "[extended:dedup] Duplicates prints 0 on --no-dedup index" "$err"
    fi
}

verify_edge_paths() {
    section_hdr "extended:edge" "CLI and error paths"
    local bare="$TMPDIR/edge_noindex.fasta"
    local empty="$TMPDIR/edge_empty.fasta"
    local junk="$TMPDIR/edge_not_fasta.txt"

    expect_fail_msg "[extended:edge] usage (no file)" "usage:" \
        "$ZFASTA" stats

    expect_fail_msg "[extended:edge] missing file" "file not found" \
        "$ZFASTA" stats "$TMPDIR/does_not_exist.fasta"

    cp "$TEST_DATA/simple.fasta" "$bare"
    rm -f "${bare}.zfi" "${bare}.fai"
    expect_fail_msg "[extended:edge] no index" "no index found" \
        "$ZFASTA" stats "$bare"

    cp "$TEST_DATA/empty.fasta" "$empty"
    rm -f "${empty}.zfi" "${empty}.fai"
    expect_fail_msg "[extended:edge] empty fasta" "file is empty" \
        "$ZFASTA" stats "$empty"

    printf 'this is not fasta\n' >"$junk"
    rm -f "${junk}.zfi" "${junk}.fai"
    expect_fail "[extended:edge] not_fasta (no index)" \
        "$ZFASTA" stats "$junk"

    # Indexed simple: --index-only must omit Composition; full must include it.
    local ok="$TMPDIR/edge_ok.fasta"
    cp "$TEST_DATA/simple.fasta" "$ok"
    ensure_index "$ok" || { fail "[extended:edge] index setup" "$TMPDIR/faidx.err"; return; }
    local out="$TMPDIR/edge_ok.index-only.txt" err="$TMPDIR/edge_ok.err"
    run_stats "$ok" index-only "$out" "$err" || { fail "[extended:edge] index-only run" "$err"; return; }
    if grep -q "Composition:" "$out"; then
        echo "Composition present under --index-only" >"$err"
        fail "[extended:edge] index-only omits Composition" "$err"
    else
        pass "[extended:edge] index-only omits Composition"
    fi
    run_stats "$ok" full "$TMPDIR/edge_ok.full.txt" "$err" || { fail "[extended:edge] full run" "$err"; return; }
    if grep -q "Composition:" "$TMPDIR/edge_ok.full.txt"; then
        pass "[extended:edge] full includes Composition"
    else
        echo "Composition missing" >"$err"
        fail "[extended:edge] full includes Composition" "$err"
    fi
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --skip-tools) SKIP_TOOLS=true ;;
        --skip-messy) SKIP_MESSY=true ;;
        --skip-lowmem) SKIP_LOWMEM=true ;;
        --skip-dedup) SKIP_DEDUP=true ;;
        --skip-edge) SKIP_EDGE=true ;;
        --skip-layout) SKIP_LAYOUT=true ;;
        *) echo "Unknown: $1"; exit 1 ;;
    esac
    shift
done

echo "z-fasta stats verification"
echo "  z-fasta:  $ZFASTA"
echo "  skip tools: $SKIP_TOOLS"
echo "  skip messy: $SKIP_MESSY"
echo "  skip lowmem: $SKIP_LOWMEM"
echo "  skip dedup: $SKIP_DEDUP"
echo "  skip edge: $SKIP_EDGE"
echo "  skip layout: $SKIP_LAYOUT"

cd "$PROJECT_DIR"
[[ -x "$ZFASTA" ]] || { echo "Error: run ./zig build first"; exit 1; }
command -v "$SAMTOOLS" &>/dev/null || { echo "Error: samtools not found"; exit 1; }
"$PYTHON" -c "from Bio import SeqIO" 2>/dev/null || { echo "Error: BioPython not available"; exit 1; }

prepare_fixtures
for name in "${FIXTURES[@]}"; do
    verify_file "$name"
done

if ! $SKIP_MESSY; then
    verify_messy
fi

if ! $SKIP_LAYOUT; then
    verify_layout_twins
fi

if ! $SKIP_LOWMEM; then
    section_hdr "index" "low-mem stats parity"
    verify_lowmem_stats "$TEST_DATA/simple.fasta" simple
    verify_lowmem_stats "$TEST_DATA/mixed_widths.fasta" mixed_widths
    verify_lowmem_stats "$MESSY_DIR/mixed_line_widths.fasta" messy_mixed_line_widths
fi

if ! $SKIP_DEDUP; then
    verify_dedup_stats
    verify_duplicates_policy
fi

if ! $SKIP_EDGE; then
    verify_edge_paths
fi

echo ""
echo "Results: $PASS passed, $FAIL failed"
[[ "$FAIL" -gt 0 ]] && { echo "VERIFICATION FAILED"; exit 1; }
echo "ALL PASSED"
