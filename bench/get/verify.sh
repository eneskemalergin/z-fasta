#!/usr/bin/env bash
# verify.sh - z-fasta get verification (subject: z-fasta only; others are oracles)
#
# Single-file runner (like bench/index/run.sh). Messy/BED/RC expected output uses an
# embedded Python oracle below; extend that or use awk/samtools, do not add sibling scripts.
#
# Tags: [parity:samtools] [parity:bedtools] [parity:seq] [extended:messy]
#       [extended:header] [extended:summary] [extended:dedup] [index:zfi] [index:fai]
#       [index:cross] [index:lowmem]
#
# Section 0 asserts .zfi == .fai byte-identical output for all uniform-index GET modes
# (single, multi incl. 15-region scan / 16+ hash boundary, BED incl. large/stranded,
# names, RC/complement/reverse). Messy non-uniform stays .zfi-only.
#
# Usage: ./verify.sh [--skip-index] [--skip-get] [--skip-multi] [--skip-bed] [--skip-rc] [--skip-edge]
set -euo pipefail

# --- config ---
SKIP_INDEX=false SKIP_GET=false SKIP_MULTI=false SKIP_BED=false SKIP_RC=false SKIP_EDGE=false

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
MESSY_DIR="$PROJECT_DIR/bench/index/messy_variants"
FIXTURE_CACHE="$SCRIPT_DIR/.fixture_cache"
TMPDIR="$SCRIPT_DIR/.verify_tmp"
mkdir -p "$TMPDIR" "$FIXTURE_CACHE"
trap 'rm -rf "$TMPDIR"' EXIT

MESSY_POS_CASES=(
    "mixed_line_widths:1-1,1-32,3-24,10-20"
    "trailing_whitespace:1-8,7-20,24-24"
    "blank_lines:1-4,9-12,17-20"
    "mixed_crlf_lf:1-8,3-18,7-14"
)

ZFASTA="${ZFASTA:-$PROJECT_DIR/zig-out/bin/z-fasta}"
SAMTOOLS="${SAMTOOLS:-samtools}"
BEDTOOLS="${BEDTOOLS:-bedtools}"
SEQTK="${SEQTK:-$PROJECT_DIR/tools/seqtk/seqtk}"
NOODLES="${NOODLES:-$PROJECT_DIR/tools/noodles_wrapper/target/release/noodles_wrapper}"
RUSTBIO="${RUSTBIO:-$PROJECT_DIR/tools/rustbio_wrapper/target/release/rustbio_wrapper}"

PASS=0 FAIL=0

# --- primitives ---
pass() { PASS=$((PASS + 1)); printf "  PASS: %s\n" "$1"; }

fail() {
    FAIL=$((FAIL + 1)); printf "  FAIL: %s\n" "$1"
    local ef="${2:-$TMPDIR/expected.tmp}" gf="${3:-$TMPDIR/got.tmp}"
    echo "    expected:"; head -3 "$ef" 2>/dev/null | sed 's/^/      /'
    echo "    got:";      head -3 "$gf" 2>/dev/null | sed 's/^/      /'
}

diff_oracle() {
    diff -q "$1" "$2" >/dev/null 2>&1 && pass "$3" || fail "$3" "${4:-$1}" "$2"
}

section_hdr() { echo ""; echo "--- [$1] $2 ---"; }

ensure_index() {
    local f="$1"
    [[ -f "${f}.fai" ]] || "$SAMTOOLS" faidx "$f" 2>/dev/null
    [[ -f "${f}.zfi" ]] || "$ZFASTA" index "$f" 2>/dev/null
}

# Messy-FASTA oracles (region/rc/rev/comp/bed). argv: MODE FASTA ... OUT
oracle() {
    python3 - "$@" <<'PY'
import sys
from pathlib import Path

RC = str.maketrans(
    "ACGTURYSWKMBDHVNacgturyswkmbdhvn",
    "TGCAAYRSWMKVHDBNtgcaayrswmkvhdbn",
)


def load_seqs(fasta: str) -> dict[str, str]:
    seqs: dict[str, list[str]] = {}
    cur: str | None = None
    for line in Path(fasta).read_text().splitlines():
        if line.startswith(">"):
            cur = line[1:].split()[0]
            seqs[cur] = []
        elif cur is not None:
            seqs[cur].append("".join(c for c in line if not c.isspace()))
    return {name: "".join(parts) for name, parts in seqs.items()}


def write_region(path: str, name: str, start: int, end: int, seq: str) -> None:
    Path(path).write_text(f">{name}:{start}-{end}\n{seq}\n")


def cmd_region(argv: list[str], transform) -> None:
    fasta, name, start, end, out = argv[2], argv[3], int(argv[4]), int(argv[5]), argv[6]
    seq = load_seqs(fasta)[name][start - 1 : end]
    if transform is not None:
        seq = transform(seq)
    write_region(out, name, start, end, seq)


def cmd_bed(argv: list[str], stranded_rc: bool = False) -> None:
    fasta, bed, out = argv[2], argv[3], argv[4]
    seqs = load_seqs(fasta)
    parts: list[str] = []
    for line in Path(bed).read_text().splitlines():
        if not line or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
            continue
        fields = line.split("\t")
        chrom, s0, e0 = fields[0], int(fields[1]), int(fields[2])
        seq = seqs[chrom][s0:e0]
        if stranded_rc and (fields[5] if len(fields) >= 6 else ".") != "-":
            seq = seq.translate(RC)[::-1]
        parts.append(f">{chrom}:{s0 + 1}-{e0}\n{seq}\n")
    Path(out).write_text("".join(parts))


COMMANDS = {
    "region": lambda a: cmd_region(a, None),
    "rc": lambda a: cmd_region(a, lambda s: s.translate(RC)[::-1]),
    "rev": lambda a: cmd_region(a, lambda s: s[::-1]),
    "comp": lambda a: cmd_region(a, lambda s: s.translate(RC)),
    "bed": lambda a: cmd_bed(a),
    "bed-stranded-rc": lambda a: cmd_bed(a, stranded_rc=True),
}

if len(sys.argv) < 2 or sys.argv[1] not in COMMANDS:
    print("usage: oracle <region|rc|rev|comp|bed|bed-stranded-rc> ...", file=sys.stderr)
    sys.exit(2)
COMMANDS[sys.argv[1]](sys.argv)
PY
}

_fixture_paths() {
    local prefix="$1" src="$2" desc="$3"
    local base="${src##*/}"; base="${base%.fasta}"
    local safe="${desc// /_}"; safe="${safe//\//_}"; safe="${safe//|/_}"
    FIXTURE_FASTA="$TMPDIR/${prefix}_${base}_${safe}.fasta"
    FIXTURE_STASH="$TMPDIR/${prefix}_${base}_${safe}.zfi"
}

_zfasta_get() {
    local fasta="$1" stdin_file="${2:-}"; shift 2
    if [[ -n "$stdin_file" ]]; then
        cat "$stdin_file" | "$ZFASTA" get "$fasta" "$@"
    else
        "$ZFASTA" get "$fasta" "$@"
    fi
}

bed_to_regions() {
    awk -F'\t' '!/^[[:space:]]*$|^#|^track|^browser/{sub(/\r$/,""); printf "%s:%d-%d\n",$1,$2+1,$3}' "$1"
}

names_to_regions() { awk '!/^[[:space:]]*$|^#/{sub(/\r$/,""); print}' "$1"; }

# Split name:start-end on the last ':' (handles sp|P12345|PROT_HUMAN:1-10).
# Sets PARSE_NAME, PARSE_START, PARSE_END; returns 0 when coords present.
parse_region_spec() {
    local region="$1"
    PARSE_NAME="" PARSE_START="" PARSE_END=""
    [[ "$region" == *:* ]] || { PARSE_NAME="$region"; return 1; }
    local coords="${region##*:}"
    PARSE_NAME="${region%:"$coords"}"
    [[ "$coords" == *-* ]] || return 1
    PARSE_START="${coords%%-*}"
    PARSE_END="${coords#*-}"
    [[ -n "$PARSE_START" && -n "$PARSE_END" ]] || return 1
    return 0
}

fasta_concat_sequence() {
    awk '/^>/{next} {gsub(/[[:space:]]/, ""); printf "%s", $0} END{print ""}' "$1" > "$2"
}

normalize_bed_output() {
    local raw_in="$1" bed_in="$2" honor_strand="$3" fasta_out="$4"
    awk -v bed_file="$bed_in" -v hs="$honor_strand" '
        function flush(    h, s2) {
            if (!idx) return
            if (!(idx in c)) { printf("record %d exceeds BED rows\n", idx) > "/dev/stderr"; exit 1 }
            h = c[idx] ":" s1[idx] "-" e[idx]
            if (hs == "1" && str[idx] == "-") h = h ":rc"
            print ">" h; s2 = seq
            while (length(s2) > 60) { print substr(s2, 1, 60); s2 = substr(s2, 61) }
            print s2; seq = ""
        }
        BEGIN {
            FS = "\t"
            while ((getline < bed_file) > 0) {
                sub(/\r$/, "", $0)
                if ($0 == "" || $0 ~ /^#/ || $0 ~ /^track/ || $0 ~ /^browser/) continue
                n++; c[n] = $1; s1[n] = $2 + 1; e[n] = $3; str[n] = NF >= 6 ? $6 : "."
            }
            close(bed_file); idx = 0; seq = ""
        }
        /^>/ { flush(); idx++; next }
        { gsub(/[[:space:]]/, "", $0); seq = seq $0 }
        END { flush(); if (idx != n) { printf("mismatch: fasta=%d bed=%d\n", idx, n) > "/dev/stderr"; exit 1 } }
    ' "$raw_in" > "$fasta_out"
}

expect_fail() {
    local label="$1"; shift
    "$@" >/dev/null 2>/dev/null && fail "$label (expected error)" || pass "$label"
}

verify_oracle_transform() {
    local mode="$1" label="$2" fasta="$3" name="$4" start="$5" end="$6"
    shift 6
    oracle "$mode" "$fasta" "$name" "$start" "$end" "$TMPDIR/expected.tmp"
    local zf_args=(get "$fasta" "${name}:${start}-${end}")
    zf_args+=("$@")
    "$ZFASTA" "${zf_args[@]}" > "$TMPDIR/got.tmp" 2>/dev/null \
        || { fail "$label (z-fasta err)"; return; }
    diff_oracle "$TMPDIR/expected.tmp" "$TMPDIR/got.tmp" "$label"
}

verify_exact_get() {
    local label="$1" expected="$2" fasta="$3"; shift 3
    "$ZFASTA" get "$fasta" "$@" > "$TMPDIR/got.tmp" 2>/dev/null \
        && diff_oracle "$expected" "$TMPDIR/got.tmp" "$label" \
        || fail "$label"
}

# --- parity: samtools ---
parity_samtools() {
    local label="$1" fasta="$2"; shift 2
    "$SAMTOOLS" faidx "$fasta" "$@" > "$TMPDIR/expected.tmp" 2>/dev/null \
        || { fail "$label (samtools err)"; return; }
    "$ZFASTA" get "$fasta" "$@" > "$TMPDIR/got.tmp" 2>/dev/null \
        || { fail "$label (z-fasta err)"; return; }
    diff_oracle "$TMPDIR/expected.tmp" "$TMPDIR/got.tmp" "$label"
}

parity_samtools_region() { parity_samtools "[parity:samtools] $3" "$1" "$2"; }

parity_samtools_regions() {
    local fasta="$1" desc="$2"; shift 2
    parity_samtools "[parity:samtools] $desc" "$fasta" "$@"
}

parity_samtools_rc() {
    local label="[parity:samtools] $3" fasta="$1" region="$2"
    "$SAMTOOLS" faidx -i --mark-strand no "$fasta" "$region" > "$TMPDIR/expected.tmp" 2>/dev/null \
        || { fail "$label (samtools err)"; return; }
    "$ZFASTA" get "$fasta" "$region" --rc > "$TMPDIR/got.tmp" 2>/dev/null \
        || { fail "$label (z-fasta err)"; return; }
    diff_oracle "$TMPDIR/expected.tmp" "$TMPDIR/got.tmp" "$label"
}

verify_open_ended_region() {
    local fasta="$1" name="$2" start="$3" len="$4" desc="$5"
    local region="${name}:${start}-"
    local seq_label="[parity:seq] $desc bases vs samtools"
    local hdr_label="[extended:header] $desc open suffix (samtools) vs clamped end (z-fasta)"

    "$SAMTOOLS" faidx "$fasta" "$region" > "$TMPDIR/expected.tmp" 2>/dev/null \
        || { fail "$seq_label (samtools err)"; return; }
    "$ZFASTA" get "$fasta" "$region" > "$TMPDIR/got.tmp" 2>/dev/null \
        || { fail "$seq_label (z-fasta err)"; return; }

    fasta_concat_sequence "$TMPDIR/expected.tmp" "$TMPDIR/exp.seq"
    fasta_concat_sequence "$TMPDIR/got.tmp" "$TMPDIR/got.seq"
    diff_oracle "$TMPDIR/exp.seq" "$TMPDIR/got.seq" "$seq_label"

    local st_hdr zf_hdr
    st_hdr=$(head -1 "$TMPDIR/expected.tmp")
    zf_hdr=$(head -1 "$TMPDIR/got.tmp")
    if [[ "$st_hdr" == ">${name}:${start}-" && "$zf_hdr" == ">${name}:${start}-${len}" ]]; then
        pass "$hdr_label"
    else
        fail "$hdr_label" "$TMPDIR/expected.tmp" "$TMPDIR/got.tmp"
    fi
}

# verify_index_cross: byte-identical get with .zfi present vs .fai-only fallback.
# Optional third arg: existing file → read GET input from stdin (e.g. BED via pipe).
verify_index_cross() {
    local desc="$1" src="$2"
    shift 2
    local stdin_file=""
    [[ $# -gt 0 && -f "$1" ]] && { stdin_file="$1"; shift; }

    _fixture_paths xidx "$src" "$desc"
    local fasta="$FIXTURE_FASTA" stash="$FIXTURE_STASH"
    local tag="[index:cross] $desc"

    cp "$src" "$fasta"
    ensure_index "$fasta"
    [[ -f "${fasta}.zfi" && -f "${fasta}.fai" ]] || { fail "[index:cross] $desc sidecars missing"; return; }

    if ! _zfasta_get "$fasta" "$stdin_file" "$@" > "$TMPDIR/zfi.out" 2>/dev/null; then
        fail "[index:zfi] $desc get failed"
        return
    fi
    mv "${fasta}.zfi" "$stash"
    if ! _zfasta_get "$fasta" "$stdin_file" "$@" > "$TMPDIR/fai.out" 2>/dev/null; then
        mv "$stash" "${fasta}.zfi"
        fail "[index:fai] $desc get via .fai failed"
        return
    fi
    diff_oracle "$TMPDIR/zfi.out" "$TMPDIR/fai.out" "$tag"
    mv "$stash" "${fasta}.zfi"
}

verify_index() {
    local mode="$1" src="$2" target="$3" desc="$4"
    local base="${src##*/}"; base="${base%.fasta}"
    _fixture_paths "idx_${mode}" "$src" "$desc"
    local fasta="$FIXTURE_FASTA" stash="$FIXTURE_STASH"
    local tag="idx $mode $base $desc"
    local zf_get=(get "$fasta") st_get=(faidx)

    cp "$src" "$fasta"
    ensure_index "$fasta"
    [[ -f "${fasta}.zfi" && -f "${fasta}.fai" ]] || { fail "[index:zfi] $tag sidecars missing"; return; }

    if [[ "$mode" == bed ]]; then
        bed_to_regions "$target" > "$TMPDIR/idx_regions.txt"
        zf_get+=(--bed "$target")
        st_get+=(-r "$TMPDIR/idx_regions.txt" "$fasta")
    else
        zf_get+=("$target")
        st_get+=("$fasta" "$target")
    fi

    "$ZFASTA" "${zf_get[@]}" > "$TMPDIR/zfi.out" 2>/dev/null \
        || { fail "[index:zfi] $tag get failed"; return; }
    pass "[index:zfi] $tag get with .zfi present"

    "$SAMTOOLS" "${st_get[@]}" > "$TMPDIR/st.out" 2>/dev/null \
        || { fail "[parity:samtools] $tag (samtools err)"; return; }
    diff_oracle "$TMPDIR/st.out" "$TMPDIR/zfi.out" "[parity:samtools] $tag via .zfi"

    mv "${fasta}.zfi" "$stash"
    "$ZFASTA" "${zf_get[@]}" > "$TMPDIR/fai.out" 2>/dev/null || {
        mv "$stash" "${fasta}.zfi"; fail "[index:fai] $tag get via .fai failed"; return
    }
    pass "[index:fai] $tag get with .zfi absent (.fai fallback)"
    diff_oracle "$TMPDIR/st.out" "$TMPDIR/fai.out" "[parity:samtools] $tag via .fai fallback"
    diff_oracle "$TMPDIR/zfi.out" "$TMPDIR/fai.out" "[index:cross] $tag .zfi == .fai fallback"
    mv "$stash" "${fasta}.zfi"
}

verify_low_mem() {
    local src="$1" target="$2" desc="$3" mode="${4:-positional}"
    _fixture_paths "lowmem_${mode}" "$src" "$desc"
    local fasta="$FIXTURE_FASTA"
    local tag="low-mem $mode ${src##*/} $desc"
    local get_cmd=(get "$fasta")

    [[ "$mode" == bed ]] && get_cmd+=(--bed "$target") || get_cmd+=("$target")

    cp "$src" "$fasta"
    rm -f "${fasta}.zfi" "${fasta}.fai"
    "$ZFASTA" index "$fasta" >/dev/null 2>&1 \
        || { fail "[index:lowmem] $tag mmap index failed"; return; }
    "$ZFASTA" "${get_cmd[@]}" > "$TMPDIR/mmap.out" 2>/dev/null \
        || { fail "[index:lowmem] $tag mmap get failed"; return; }

    rm -f "${fasta}.zfi"
    "$ZFASTA" index --low-mem "$fasta" >/dev/null 2>&1 \
        || { fail "[index:lowmem] $tag streaming index failed"; return; }
    "$ZFASTA" "${get_cmd[@]}" > "$TMPDIR/low.out" 2>/dev/null \
        || { fail "[index:lowmem] $tag streaming get failed"; return; }
    diff_oracle "$TMPDIR/mmap.out" "$TMPDIR/low.out" "[index:lowmem] $tag mmap == --low-mem"
}

verify_messy_zfi_required() {
    local name="$1" start="$2" end="$3"
    local fasta="$TMPDIR/messy_nozfi_${name}.fasta"
    local stash="$TMPDIR/messy_nozfi_${name}.zfi"
    [[ -f "$MESSY_DIR/${name}.fasta" ]] || { fail "[extended:messy] missing fixture $name"; return; }
    prepare_messy_fixture "$name" "$fasta" || { fail "[extended:messy] index failed $name"; return; }
    mv "${fasta}.zfi" "$stash"
    rm -f "${fasta}.fai"
    expect_fail "[extended:messy] $name:${start}-${end} without .zfi" \
        "$ZFASTA" get "$fasta" "${name}:${start}-${end}"
    mv "$stash" "${fasta}.zfi"
}

# --- extended: messy (Python oracle) ---
prepare_messy_fixture() {
    cp "$MESSY_DIR/${1}.fasta" "$2"
    "$ZFASTA" index "$2" >/dev/null 2>&1
}

verify_messy_region() {
    local fasta="$1" name="$2" start="$3" end="$4" desc="$5"
    local mode=region label="[extended:messy] $desc" extra=()
    [[ "${6:-}" == "--rc" ]] && { mode=rc; label="$label --rc"; extra=(--rc); }
    verify_oracle_transform "$mode" "$label" "$fasta" "$name" "$start" "$end" "${extra[@]}"
}

verify_expected_bed() {
    local label="[extended:messy] $1" fasta="$2" bed="$3" expected="$4"
    shift 4
    "$ZFASTA" get "$fasta" --bed "$bed" "$@" > "$TMPDIR/got.tmp" 2>/dev/null \
        || { fail "$label (z-fasta err)"; return; }
    diff_oracle "$expected" "$TMPDIR/got.tmp" "$label"
}

# --- parity: bedtools ---
parity_bedtools() {
    local label="[parity:bedtools] $1" bed="$2" honor_strand="$3" fasta="$4"
    local chunk="${5:-4096}" stdin="${6:-0}" strand_flag="${7:---strand-aware}"
    local bt=(bedtools getfasta -fi "$fasta" -bed "$bed")
    [[ "$honor_strand" == "1" ]] && bt+=(-s)
    local zf=(get "$fasta" --bed "$bed" --chunk-size "$chunk")
    [[ "$honor_strand" == "1" ]] && zf+=("$strand_flag")
    if [[ "$stdin" == "1" ]]; then
        zf[2]="--bed"
        zf[3]="-"
        "${bt[@]}" > "$TMPDIR/bt.raw" 2>/dev/null || { fail "$label (bedtools err)"; return; }
        cat "$bed" | "$ZFASTA" "${zf[@]}" > "$TMPDIR/zf.raw" 2>/dev/null || { fail "$label (z-fasta err)"; return; }
    else
        "${bt[@]}" > "$TMPDIR/bt.raw" 2>/dev/null || { fail "$label (bedtools err)"; return; }
        "$ZFASTA" "${zf[@]}" > "$TMPDIR/zf.raw" 2>/dev/null || { fail "$label (z-fasta err)"; return; }
    fi
    normalize_bed_output "$TMPDIR/bt.raw" "$bed" "$honor_strand" "$TMPDIR/bt.norm"
    normalize_bed_output "$TMPDIR/zf.raw" "$bed" "$honor_strand" "$TMPDIR/zf.norm"
    diff_oracle "$TMPDIR/bt.norm" "$TMPDIR/zf.norm" "$label"
}

parity_bed_samtools() {
    local label="[parity:samtools] $1" bed="$2" fasta="$3"
    local chunk="${4:-4096}" stdin="${5:-0}"
    local zf=(get "$fasta" --chunk-size "$chunk")
    bed_to_regions "$bed" > "$TMPDIR/st_regions.txt"
    if [[ "$stdin" == "1" ]]; then
        zf+=(--bed -)
        cat "$bed" | "$ZFASTA" "${zf[@]}" > "$TMPDIR/got.tmp" 2>/dev/null || { fail "$label (z-fasta err)"; return; }
    else
        zf+=(--bed "$bed")
        "$ZFASTA" "${zf[@]}" > "$TMPDIR/got.tmp" 2>/dev/null || { fail "$label (z-fasta err)"; return; }
    fi
    "$SAMTOOLS" faidx -r "$TMPDIR/st_regions.txt" "$fasta" > "$TMPDIR/expected.tmp" 2>/dev/null \
        || { fail "$label (samtools err)"; return; }
    diff_oracle "$TMPDIR/expected.tmp" "$TMPDIR/got.tmp" "$label"
}

parity_bed_rc_composed() {
    local label="[parity:bedtools] $1" fasta="$2" bed="$3"
    local chunk="${4:-4096}" stdin="${5:-0}"
    local zf=(get "$fasta" --strand-aware --rc --chunk-size "$chunk")
    if ! "$BEDTOOLS" getfasta -fi "$fasta" -bed "$bed" -s 2>/dev/null \
        | "$SEQTK" seq -r > "$TMPDIR/bt.raw"; then
        fail "$label (bedtools/seqtk err)"
        return
    fi
    if [[ "$stdin" == "1" ]]; then
        cat "$bed" | "$ZFASTA" "${zf[@]}" --bed - > "$TMPDIR/zf.raw" 2>/dev/null || { fail "$label (z-fasta err)"; return; }
    else
        "$ZFASTA" "${zf[@]}" --bed "$bed" > "$TMPDIR/zf.raw" 2>/dev/null || { fail "$label (z-fasta err)"; return; }
    fi
    normalize_bed_output "$TMPDIR/bt.raw" "$bed" 0 "$TMPDIR/bt.norm"
    normalize_bed_output "$TMPDIR/zf.raw" "$bed" 0 "$TMPDIR/zf.norm"
    diff_oracle "$TMPDIR/bt.norm" "$TMPDIR/zf.norm" "$label"
}

# --- optional wrapper cross-checks ---
wrapper_vs() {
    local desc="$1" cmp="$2" bin="$3"; shift 3
    [[ -x "$bin" ]] || return 0
    if "$bin" "$@" > "$TMPDIR/wrapper.tmp" 2>/dev/null; then
        diff_oracle "$cmp" "$TMPDIR/wrapper.tmp" "$desc"
    else
        fail "$desc (wrapper err)"
    fi
}

wrapper_region() {
    local desc="$1" fasta="$2" region="$3"
    wrapper_vs "$desc noodles" "$TMPDIR/got.tmp" "$NOODLES" get "$fasta" "$region"
    wrapper_vs "$desc rustbio" "$TMPDIR/got.tmp" "$RUSTBIO" get "$fasta" "$region"
}

wrapper_regions() {
    local desc="$1" fasta="$2"; shift 2
    "$SAMTOOLS" faidx "$fasta" "$@" > "$TMPDIR/expected.tmp" 2>/dev/null || return 0
    wrapper_vs "$desc noodles" "$TMPDIR/expected.tmp" "$NOODLES" get "$fasta" "$@"
    wrapper_vs "$desc rustbio" "$TMPDIR/expected.tmp" "$RUSTBIO" get "$fasta" "$@"
}

wrapper_rc() {
    local desc="$1" fasta="$2" region="$3"
    local st
    st=$("$SAMTOOLS" faidx -i --mark-strand no "$fasta" "$region" 2>/dev/null) || return 0
    printf "%s\n" "$st" > "$TMPDIR/expected.tmp"
    wrapper_vs "$desc noodles" "$TMPDIR/expected.tmp" "$NOODLES" get "$fasta" "$region" --rc
    wrapper_vs "$desc rustbio" "$TMPDIR/expected.tmp" "$RUSTBIO" get "$fasta" "$region" --rc
}

wrapper_bed() {
    local desc="$1" bed="$2" fasta="${3:-$PROJECT_DIR/tests/data/simple.fasta}"
    bed_to_regions "$bed" > "$TMPDIR/st_regions.txt"
    "$SAMTOOLS" faidx -r "$TMPDIR/st_regions.txt" "$fasta" > "$TMPDIR/expected.tmp" 2>/dev/null || return 0
    for pair in "noodles:$NOODLES" "rustbio:$RUSTBIO"; do
        IFS=: read -r name bin <<< "$pair"
        [[ -x "$bin" ]] || continue
        if "$bin" get "$fasta" --bed "$bed" > "$TMPDIR/${name}.raw" 2>/dev/null; then
            normalize_bed_output "$TMPDIR/expected.tmp" "$bed" 0 "$TMPDIR/st.norm"
            normalize_bed_output "$TMPDIR/${name}.raw" "$bed" 0 "$TMPDIR/${name}.norm"
            diff_oracle "$TMPDIR/st.norm" "$TMPDIR/${name}.norm" "$desc $name"
        else
            fail "$desc ($name err)"
        fi
    done
}

# --- fixtures ---
gen_names_file() { printf '# names file\nseq2\n\nseq1\nseq2\n' > "$1"; }

gen_bed_file() {
    local out="$1" count="$2"
    : > "$out"
    printf '# synthetic BED for %s\ntrack name=verify_bed_%s\n' "$count" "$count" >> "$out"
    cat >> "$out" <<'BED'
seq1	0	24	full_seq1	0	+
seq1	0	5	dup_a	0	+
seq1	0	5	dup_b	0	+
seq1	2	8	overlap_minus	0	-
seq2	0	12	full_seq2	0	-
seq2	0	4	short_plus	0	+
BED
    for ((i=6; i<count; i++)); do
        local idx=$i chrom len span start end strand
        ((idx%2==0)) && chrom=seq1 len=24 || chrom=seq2 len=12
        ((len<=4)) && span=1 || span=$((idx%4+1))
        start=$((idx%(len-span+1))); end=$((start+span))
        ((idx%3==0)) && strand=- || strand=+
        printf '%s\t%d\t%d\t%s_%d\t0\t%s\n' "$chrom" "$start" "$end" "$chrom" "$idx" "$strand" >> "$out"
    done
}

gen_iupac_fixture() {
    local out="$FIXTURE_CACHE/iupac.fasta"
    [[ -f "${out}.zfi" ]] && { printf '%s\n' "$out"; return; }
    printf '>iupac_all\nACGTURYSWKMBDHVNacgturyswkmbdhvnu\n' > "$out"
    "$ZFASTA" index "$out" >/dev/null 2>&1
    printf '%s\n' "$out"
}

gen_chrom_fixture() {
    local total="${1:-16384}"
    local out="$FIXTURE_CACHE/chrom_${total}.fasta"
    [[ -f "${out}.zfi" ]] && { printf '%s\n' "$out"; return; }
    {
        echo ">chrSynthetic"
        awk -v total="$total" 'BEGIN{pattern="ACGTNRYWSKMBDHVacgtnrywskmbdhv"; s=""
            while(length(s)<total){s=s pattern} s=substr(s,1,total)
            while(length(s)>0){print substr(s,1,71); s=substr(s,72)}}'
    } > "$out"
    "$ZFASTA" index "$out" >/dev/null 2>&1
    printf '%s\n' "$out"
}

gen_bed_rc_file() {
    local out="$1" stranded_out="$2"
    : > "$out"; : > "$stranded_out"
    for ((i=0; i<1200; i++)); do
        local s=$((1000+i*400)) e strand
        e=$((s+120))
        printf '%s\t%d\t%d\n' "chrSynthetic" "$s" "$e" >> "$out"
        ((i%2==0)) && strand=+ || strand=-
        printf '%s\t%d\t%d\tregion_%05d\t0\t%s\n' "chrSynthetic" "$s" "$e" "$i" "$strand" >> "$stranded_out"
    done
}

REG20=(); for i in {1..10}; do REG20+=("seq1:${i}-${i}" "seq2:${i}-${i}"); done
REG20_REV=(); for i in {10..1}; do REG20_REV+=("seq2:${i}-${i}" "seq1:${i}-${i}"); done

messy_bed_for() {
    case "$1" in
        mixed_line_widths) cat <<'BED'
mixed_line_widths	0	12	line1	0	+
mixed_line_widths	9	20	cross	0	+
BED
            ;;
        trailing_whitespace) cat <<'BED'
trailing_whitespace	0	8	start	0	+
trailing_whitespace	16	24	end	0	+
BED
            ;;
        blank_lines) cat <<'BED'
blank_lines	0	4	start	0	+
blank_lines	8	12	mid	0	+
BED
            ;;
        mixed_crlf_lf) cat <<'BED'
mixed_crlf_lf	0	8	start	0	+
mixed_crlf_lf	6	14	cross	0	+
BED
            ;;
    esac
}

# --- sections ---
section0_index() {
    section_hdr 0 "Index path coverage (.zfi vs .fai fallback)"
    echo "  Uniform fixtures only; messy_variants stay .zfi-only ([extended:messy])"
    local simple="tests/data/simple.fasta" proteome="tests/data/proteome.fasta"
    local edge="tests/data/edge_cases.fasta" mixed="tests/data/mixed_widths.fasta"
    local bed_small="$TMPDIR/idx_bed_small.bed" bed_medium="$TMPDIR/idx_bed_medium.bed"
    local bed_large="$TMPDIR/idx_bed_large.bed"
    local mx_bed="$TMPDIR/idx_mx.bed" names="$TMPDIR/idx_names.txt"
    local bed_rc="$TMPDIR/idx_bed_rc.bed"
    local -a p16=() reg15=() reg1024=()
    local prot1="sp|P12345|PROT_HUMAN" prot2="sp|Q98765|ANOT_MOUSE"

    gen_bed_file "$bed_small" 10
    gen_bed_file "$bed_medium" 100
    gen_bed_file "$bed_large" 1000
    gen_names_file "$names"
    for i in {1..8}; do
        p16+=("${prot1}:${i}-$((i + 2))" "${prot2}:${i}-$((i + 2))")
    done
    reg15=("${REG20[@]:0:15}")
    for _ in {1..64}; do reg1024+=("seq1:1-1"); done

    cat > "$mx_bed" <<'BED'
mixed1	54	75	mixed1_span	0	+
mixed2	74	95	mixed2_span	0	-
mixed3	57	72	mixed3_span	0	+
BED
    cat > "$bed_rc" <<'BED'
seq1	0	5	plus	0	+
seq1	0	5	minus	0	-
seq1	2	8	overlap_minus	0	-
seq2	0	4	short_plus	0	+
BED

    # Spot checks with samtools oracle on both index paths.
    verify_index positional "$simple" "seq1:1-10" "sub-region"
    verify_index positional "$simple" "seq1" "full sequence"
    verify_index positional tests/data/single.fasta "single_sequence:1-4" "single record"
    verify_index positional "$proteome" "sp|P12345|PROT_HUMAN:1-10" "pipe name"
    verify_index positional "$mixed" "mixed1:55-75" "mixed-width record"
    verify_index bed "$simple" "$bed_small" "small BED"

    echo "  [index:cross] single-region"
    verify_index_cross "sub-region" "$simple" "seq1:1-10"
    verify_index_cross "full sequence" "$simple" "seq1"
    verify_index_cross "open-ended region" "$simple" "seq1:10-"
    verify_index_cross "pipe name sub-region" "$proteome" "sp|P12345|PROT_HUMAN:1-10"
    verify_index_cross "mixed-width sub-region" "$mixed" "mixed1:55-75"
    verify_index_cross "duplicate name (last wins)" "$edge" "dupname"
    verify_index_cross "lowercase sub-region" "$edge" "lowercase:1-6"

    echo "  [index:cross] multi-region"
    verify_index_cross "two sub-regions" "$simple" "seq1:1-10" "seq1:13-24"
    verify_index_cross "full + sub-region" "$simple" "seq1" "seq2:3-10"
    verify_index_cross "15 regions (scan path max)" "$simple" "${reg15[@]}"
    verify_index_cross "16 regions (hash map path)" "$proteome" "${p16[@]}"
    verify_index_cross "20 regions (sort path)" "$simple" "${REG20[@]}"
    verify_index_cross "64 regions (hash map path)" "$simple" "${reg1024[@]}"

    echo "  [index:cross] BED and names"
    verify_index_cross "small BED default chunk" "$simple" --bed "$bed_small" --chunk-size 3
    verify_index_cross "medium BED chunk 97" "$simple" --bed "$bed_medium" --chunk-size 97
    verify_index_cross "large BED chunk 257" "$simple" --bed "$bed_large" --chunk-size 257
    verify_index_cross "BED chunk-size -1" "$simple" --bed "$bed_medium" --chunk-size -1
    verify_index_cross "BED chunk-size 1" "$simple" --bed "$bed_small" --chunk-size 1
    verify_index_cross "BED via stdin" "$simple" "$bed_small" --bed - --chunk-size 3
    verify_index_cross "BED stranded" "$simple" --bed "$bed_small" --chunk-size 3 --strand-aware
    verify_index_cross "BED --honor-strand alias" "$simple" --bed "$bed_small" --chunk-size 4096 --honor-strand
    verify_index_cross "large BED stranded chunk 257" "$simple" --bed "$bed_large" --chunk-size 257 --strand-aware
    verify_index_cross "mixed-width BED" "$mixed" --bed "$mx_bed"
    verify_index_cross "names file" "$simple" --names "$names"

    echo "  [index:cross] RC / complement / reverse"
    verify_index_cross "--rc" "$simple" "seq1:1-5" --rc
    verify_index_cross "--complement-only" "$simple" "seq1:1-5" --complement-only
    verify_index_cross "--reverse-only" "$simple" "seq1:1-5" --reverse-only
    verify_index_cross "multi --rc" "$simple" "seq1:1-5" "seq1:10-15" "seq1:20-24" --rc
    verify_index_cross "names --rc" "$simple" --names "$names" --rc
    verify_index_cross "BED --strand-aware --rc" "$simple" --bed "$bed_rc" --strand-aware --rc
    verify_index_cross "--rc --annotate-rc" "$simple" "seq1:1-5" --rc --annotate-rc
    verify_index_cross "--complement-only --annotate-rc" "$simple" "seq1:1-5" --complement-only --annotate-rc
    verify_index_cross "IUPAC --rc" "$(gen_iupac_fixture)" "iupac_all:1-33" --rc
    verify_index_cross "IUPAC --complement-only" "$(gen_iupac_fixture)" "iupac_all:1-33" --complement-only

    verify_low_mem "$simple" "seq1:1-10" "positional"
    verify_low_mem "$mixed" "mixed1:55-75" "mixed-width"
    verify_low_mem "$MESSY_DIR/mixed_line_widths.fasta" "mixed_line_widths:3-24" "messy positional"
    verify_low_mem "$simple" "$bed_small" "small BED" bed
    verify_messy_zfi_required mixed_line_widths 1 8
}

section1() {
    section_hdr 1 "Single-region extraction"

    test_seq() {
        local fasta="$1" name="$2" len="$3"
        local label="${fasta##*/}:${name}"
        parity_samtools_region "$fasta" "$name" "$label full"
        wrapper_region "$label full wrapper" "$fasta" "$name"
        if ((len>=10)); then
            parity_samtools_region "$fasta" "${name}:1-10" "$label :1-10"
            wrapper_region "$label :1-10 wrapper" "$fasta" "${name}:1-10"
            local ls=$((len-9))
            parity_samtools_region "$fasta" "${name}:${ls}-${len}" "$label :${ls}-${len}"
            wrapper_region "$label :${ls}-${len} wrapper" "$fasta" "${name}:${ls}-${len}"
        else
            parity_samtools_region "$fasta" "${name}:1-${len}" "$label :1-${len}"
        fi
        if ((len>=101)); then
            local ms=$((len/2-50)) me=$((len/2+50))
            ((ms<1)) && ms=1; ((me>len)) && me=$len
            parity_samtools_region "$fasta" "${name}:${ms}-${me}" "$label mid-span"
            wrapper_region "$label mid-span wrapper" "$fasta" "${name}:${ms}-${me}"
        fi
        parity_samtools_region "$fasta" "${name}:1-1" "$label :1-1"
        parity_samtools_region "$fasta" "${name}:${len}-${len}" "$label :${len}-${len}"
        parity_samtools_region "$fasta" "${name}:1-${len}" "$label :1-${len}"
        parity_samtools_region "$fasta" "${name}:1-$((len+100))" "$label clamp"
    }

    test_file() {
        local fasta="$1"
        [[ -f "${fasta}.fai" ]] || return
        echo "  ${fasta##*/}"
        while IFS=$'\t' read -r name length _; do
            ((length>0)) && test_seq "$fasta" "$name" "$length"
        done < "${fasta}.fai"
    }

    for f in tests/data/simple.fasta tests/data/proteome.fasta tests/data/single.fasta tests/data/edge_cases.fasta; do
        [[ -f "$f" ]] && { ensure_index "$f"; test_file "$f"; }
    done
    cp tests/data/mixed_widths.fasta "$TMPDIR/mixed_widths.fasta"
    ensure_index "$TMPDIR/mixed_widths.fasta"
    test_file "$TMPDIR/mixed_widths.fasta"

    verify_open_ended_region tests/data/simple.fasta seq1 10 24 "simple open-ended region"

    echo "  messy_variants [extended:messy]"
    for case in "${MESSY_POS_CASES[@]}"; do
        IFS=: read -r name regions <<< "$case"
        local tgt="$TMPDIR/messy_${name}.fasta"
        [[ -f "$MESSY_DIR/${name}.fasta" ]] || { fail "[extended:messy] missing fixture $name"; continue; }
        prepare_messy_fixture "$name" "$tgt" || { fail "[extended:messy] index failed $name"; continue; }
        IFS=, read -ra spans <<< "$regions"
        for span in "${spans[@]}"; do
            IFS=- read -r start end <<< "$span"
            verify_messy_region "$tgt" "$name" "$start" "$end" "$name:${start}-${end}"
        done
    done
}

section2() {
    section_hdr 2 "Multi-region extraction"
    local simple="tests/data/simple.fasta" proteome="tests/data/proteome.fasta" edge="tests/data/edge_cases.fasta"
    ensure_index "$simple"; ensure_index "$proteome"; ensure_index "$edge"
    local LONG="AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"

    echo "  simple.fasta"
    parity_samtools_regions "$simple" "same seq, two sub-regions" "seq1:1-10" "seq1:13-24"
    parity_samtools_regions "$simple" "seqs file order" "seq1:1-12" "seq2:1-6"
    parity_samtools_regions "$simple" "seqs reversed order" "seq2:1-6" "seq1:1-12"
    parity_samtools_regions "$simple" "overlapping" "seq1:1-15" "seq1:10-24"
    parity_samtools_regions "$simple" "duplicate" "seq1:1-12" "seq1:1-12"
    parity_samtools_regions "$simple" "full + sub" "seq1" "seq2:3-10"
    parity_samtools_regions "$simple" "both full" "seq1" "seq2"
    parity_samtools_regions "$simple" "triple duplicate" "seq1:1-5" "seq1:1-5" "seq1:1-5"
    parity_samtools_regions "$simple" "full + 2 sub" "seq1" "seq2:1-6" "seq1:10-20"
    parity_samtools_regions "$simple" "all reversed" "seq2" "seq1"
    parity_samtools_regions "$simple" "20 regions sort path" "${REG20[@]}"
    parity_samtools_regions "$simple" "20 regions reversed" "${REG20_REV[@]}"
    parity_samtools_regions "$simple" "15 regions (scan path max)" "${REG20[@]:0:15}"
    wrapper_regions "multi same seq" "$simple" "seq1:1-10" "seq1:13-24"
    wrapper_regions "multi file order" "$simple" "seq1:1-12" "seq2:1-6"
    wrapper_regions "multi 20 regions" "$simple" "${REG20[@]}"

    echo "  edge_cases.fasta"
    parity_samtools_regions "$edge" "long + normal" "single_line" "$LONG"
    parity_samtools_regions "$edge" "long reversed" "$LONG" "single_line"
    parity_samtools_regions "$edge" "lowercase + nonstandard" "lowercase" "nonstandard"
    parity_samtools_regions "$edge" "overlapping same seq" "single_line:1-5" "single_line:3-8"

    echo "  proteome.fasta"
    parity_samtools_regions "$proteome" "pipe names" "sp|P12345|PROT_HUMAN" "sp|Q98765|ANOT_MOUSE"
    parity_samtools_regions "$proteome" "pipe reversed" "sp|Q98765|ANOT_MOUSE" "sp|P12345|PROT_HUMAN"
    parity_samtools_regions "$proteome" "pipe sub + full" "sp|P12345|PROT_HUMAN:1-10" "sp|Q98765|ANOT_MOUSE"
    parity_samtools_regions "$proteome" "pipe duplicate" "sp|P12345|PROT_HUMAN" "sp|P12345|PROT_HUMAN"
    local prot1="sp|P12345|PROT_HUMAN" prot2="sp|Q98765|ANOT_MOUSE" P16=()
    for i in {1..8}; do P16+=("${prot1}:${i}-$((i+2))" "${prot2}:${i}-$((i+2))"); done
    parity_samtools_regions "$proteome" "16 proteome sub-regions" "${P16[@]}"

    ALL_EDGE=(); while IFS=$'\t' read -r name _; do ALL_EDGE+=("$name"); done < "${edge}.fai"
    (( ${#ALL_EDGE[@]} > 0 )) && {
        parity_samtools_regions "$edge" "all edge seqs" "${ALL_EDGE[@]}"
        REV_EDGE=(); for ((i=${#ALL_EDGE[@]}-1; i>=0; i--)); do REV_EDGE+=("${ALL_EDGE[$i]}"); done
        parity_samtools_regions "$edge" "all edge reversed" "${REV_EDGE[@]}"
    }

    REG1024=()
    for _ in {1..1024}; do REG1024+=("seq1:1-1"); done
    parity_samtools_regions "$simple" "1024 positional regions" "${REG1024[@]}"
    REG1025=("${REG1024[@]}" "seq1:1-1")
    expect_fail "1025 positional regions rejected" "$ZFASTA" get "$simple" "${REG1025[@]}"
}

section3() {
    section_hdr 3 "BED and names extraction"
    local simple="$PROJECT_DIR/tests/data/simple.fasta"
    ensure_index "$simple"
    local BED_SMALL="$TMPDIR/small.bed" BED_MEDIUM="$TMPDIR/medium.bed"
    gen_bed_file "$BED_SMALL" 10; gen_bed_file "$BED_MEDIUM" 100; gen_bed_file "$TMPDIR/large.bed" 1000

    run_size_suite() {
        local label="$1" count="$2" chunk="$3"
        local bed="$TMPDIR/${label}.bed"
        gen_bed_file "$bed" "$count"
        echo "  suite $label ($count rows, chunk=$chunk)"
        parity_bedtools "$label default" "$bed" 0 "$simple" "$chunk"
        parity_bedtools "$label stranded" "$bed" 1 "$simple" "$chunk"
        parity_bed_samtools "$label vs samtools" "$bed" "$simple" "$chunk"
        parity_bedtools "$label via stdin" "$bed" 0 "$simple" "$chunk" 1
        parity_bed_samtools "$label stdin vs samtools" "$bed" "$simple" "$chunk" 1
        parity_bedtools "$label stranded stdin" "$bed" 1 "$simple" "$chunk" 1
    }

    run_size_suite small 10 3
    run_size_suite medium 100 97
    run_size_suite large 1000 257

    echo "  chunk-size 1"
    parity_bed_samtools "small chunk-size 1" "$BED_SMALL" "$simple" 1
    parity_bedtools "small chunk-size 1" "$BED_SMALL" 0 "$simple" 1

    echo "  chunk-size -1"
    parity_bed_samtools "medium all-in-one" "$BED_MEDIUM" "$simple" -1
    parity_bedtools "medium all-in-one" "$BED_MEDIUM" 0 "$simple" -1

    if command -v "$BEDTOOLS" &>/dev/null; then
        parity_bedtools "--honor-strand alias" "$BED_SMALL" 1 "$simple" 4096 0 --honor-strand
    fi

    "$ZFASTA" get "$simple" seq1:1-10 seq2:1-6 --summary > "$TMPDIR/got.tmp" 2>"$TMPDIR/summary.err" \
        && grep -Eq '^summary: regions=[0-9]+ total_bases=[0-9]+ elapsed_s=' "$TMPDIR/summary.err" \
        && pass "[extended:summary] --summary stderr" || fail "[extended:summary] --summary stderr"

    local NAMES="$TMPDIR/names.txt"
    gen_names_file "$NAMES"
    names_to_regions "$NAMES" > "$TMPDIR/st_names.txt"
    "$ZFASTA" get "$simple" --names "$NAMES" > "$TMPDIR/got.tmp" 2>/dev/null \
        && "$SAMTOOLS" faidx -r "$TMPDIR/st_names.txt" "$simple" > "$TMPDIR/expected.tmp" 2>/dev/null \
        && diff_oracle "$TMPDIR/expected.tmp" "$TMPDIR/got.tmp" "[parity:samtools] names file" \
        || fail "[parity:samtools] names file"

    local mx="$TMPDIR/mx.fasta" mx_bed="$TMPDIR/mx.bed"
    cp "$PROJECT_DIR/tests/data/mixed_widths.fasta" "$mx"
    ensure_index "$mx"
    cat > "$mx_bed" <<'BED'
mixed1	54	75	mixed1_span	0	+
mixed2	74	95	mixed2_span	0	-
mixed3	57	72	mixed3_span	0	+
BED
    command -v "$BEDTOOLS" &>/dev/null && {
        parity_bedtools "mixed-width default" "$mx_bed" 0 "$mx"
        parity_bedtools "mixed-width stranded" "$mx_bed" 1 "$mx"
    }
    parity_bed_samtools "mixed-width vs samtools" "$mx_bed" "$mx"

    for messy_name in mixed_line_widths trailing_whitespace blank_lines mixed_crlf_lf; do
        local mf="$TMPDIR/messy_bed_${messy_name}.fasta" bed="$TMPDIR/messy_bed_${messy_name}.bed" exp="$TMPDIR/messy_bed_${messy_name}.fa"
        prepare_messy_fixture "$messy_name" "$mf" || { fail "[extended:messy] bed index $messy_name"; continue; }
        messy_bed_for "$messy_name" > "$bed"
        oracle bed "$mf" "$bed" "$exp"
        verify_expected_bed "messy BED $messy_name" "$mf" "$bed" "$exp"
    done

    local nu="$TMPDIR/nu_rc_bed.fasta" nu_bed="$TMPDIR/nu_rc_bed.bed" nu_exp="$TMPDIR/nu_rc_bed_expected.fa"
    prepare_messy_fixture mixed_line_widths "$nu" || fail "[extended:messy] stranded RC index"
    cat > "$nu_bed" <<'BED'
mixed_line_widths	9	20	plus	0	+
mixed_line_widths	14	28	minus	0	-
BED
    oracle bed-stranded-rc "$nu" "$nu_bed" "$nu_exp"
    verify_expected_bed "messy BED --strand-aware --rc" "$nu" "$nu_bed" "$nu_exp" --strand-aware --rc

    wrapper_bed "small BED" "$BED_SMALL"
    wrapper_bed "medium BED" "$BED_MEDIUM"
    wrapper_bed "mixed-width BED" "$mx_bed" "$mx"

    printf 'seq1:1-5\n' > "$TMPDIR/names_reg.txt"
    expect_fail "names file region syntax rejected" "$ZFASTA" get "$simple" --names "$TMPDIR/names_reg.txt"
    printf 'seq1\t5\t5\tzero_len\t0\t+\n' > "$TMPDIR/bad.bed"
    expect_fail "0-length BED rejected" "$ZFASTA" get "$simple" --bed "$TMPDIR/bad.bed"
    printf 'seq1\t5\n' > "$TMPDIR/bad.bed"
    expect_fail "short BED rejected" "$ZFASTA" get "$simple" --bed "$TMPDIR/bad.bed"
    printf 'seq1\t0\t5\tname\t0\t?\n' > "$TMPDIR/bad.bed"
    expect_fail "invalid BED strand rejected" "$ZFASTA" get "$simple" --bed "$TMPDIR/bad.bed" --strand-aware
}

section4() {
    section_hdr 4 "Reverse complement"
    local simple="$PROJECT_DIR/tests/data/simple.fasta" single="$PROJECT_DIR/tests/data/single.fasta"
    local edge="$PROJECT_DIR/tests/data/edge_cases.fasta"
    local iupac chrom
    iupac=$(gen_iupac_fixture); chrom=$(gen_chrom_fixture 500000)
    ensure_index "$iupac"; ensure_index "$chrom"

    for r in seq1:1-5 seq1:10-15 seq2:1-12; do
        parity_samtools_rc "$simple" "$r" "simple $r"
        wrapper_rc "simple $r wrapper" "$simple" "$r"
    done
    parity_samtools_rc "$single" "single_sequence:1-4" "single span"
    parity_samtools_rc "$edge" "lowercase:1-12" "lowercase preservation"
    parity_samtools_rc "$iupac" "iupac_all:1-33" "IUPAC full"
    parity_samtools_rc "$chrom" "chrSynthetic" "chrom full"

    local mx="$TMPDIR/mx_rc.fasta"
    cp "$PROJECT_DIR/tests/data/mixed_widths.fasta" "$mx"
    ensure_index "$mx"
    parity_samtools_rc "$mx" "mixed1:55-75" "mixed-width region"

    for spec in "mixed_line_widths:3:24" "trailing_whitespace:7:20" "blank_lines:9:12" "mixed_crlf_lf:7:14"; do
        IFS=: read -r n s e <<< "$spec"
        local mf="$TMPDIR/messy_rc_${n}.fasta"
        prepare_messy_fixture "$n" "$mf" || { fail "[extended:messy] RC index $n"; continue; }
        verify_messy_region "$mf" "$n" "$s" "$e" "messy RC $n:${s}-${e}" --rc
    done

    verify_oracle_transform rev "[parity:seq] simple --reverse-only" "$simple" seq1 1 5 --reverse-only
    verify_oracle_transform comp "[parity:seq] simple --complement-only" "$simple" seq1 1 5 --complement-only
    printf '>seq1:1-5 (reverse)\nATGCA\n' > "$TMPDIR/rev_annot.fa"
    verify_exact_get "[extended:header] annotated --reverse-only" "$TMPDIR/rev_annot.fa" "$simple" seq1:1-5 --reverse-only --annotate-rc
    printf '>seq1:1-5 (complement)\nTGCAT\n' > "$TMPDIR/comp_annot.fa"
    verify_exact_get "[extended:header] annotated --complement-only" "$TMPDIR/comp_annot.fa" "$simple" seq1:1-5 --complement-only --annotate-rc

    printf '>iupac_all:1-33 (complement)\nTGCAAYRSWMKVHDBNtgcaayrswmkvhdbna\n' > "$TMPDIR/iupac_comp.fa"
    "$ZFASTA" get "$iupac" iupac_all:1-33 --complement-only --annotate-rc > "$TMPDIR/got.tmp" 2>/dev/null \
        && diff_oracle "$TMPDIR/iupac_comp.fa" "$TMPDIR/got.tmp" "[extended:header] IUPAC --complement-only annotated" \
        || fail "[extended:header] IUPAC --complement-only annotated"

    : > "$TMPDIR/expected.tmp"
    for r in seq1:1-5 seq1:10-15 seq1:20-24; do "$ZFASTA" get "$simple" "$r" --rc >> "$TMPDIR/expected.tmp" 2>/dev/null; done
    "$ZFASTA" get "$simple" seq1:1-5 seq1:10-15 seq1:20-24 --rc > "$TMPDIR/got.tmp" 2>/dev/null
    diff_oracle "$TMPDIR/expected.tmp" "$TMPDIR/got.tmp" "[parity:samtools] multi RC self-consistency"

    local nf="$TMPDIR/names_rc.txt"
    printf 'seq2\nseq1\n' > "$nf"
    names_to_regions "$nf" > "$TMPDIR/st_rc_names.txt"
    "$SAMTOOLS" faidx -i --mark-strand no -r "$TMPDIR/st_rc_names.txt" "$simple" > "$TMPDIR/expected.tmp" 2>/dev/null \
        && "$ZFASTA" get "$simple" --names "$nf" --rc > "$TMPDIR/got.tmp" 2>/dev/null \
        && diff_oracle "$TMPDIR/expected.tmp" "$TMPDIR/got.tmp" "[parity:samtools] names --rc" \
        || fail "[parity:samtools] names --rc"

    local bed_comp="$TMPDIR/bed_comp.bed" bed_exp="$TMPDIR/bed_comp_expected.fa"
    cat > "$bed_comp" <<'BED'
seq1	0	5	plus	0	+
seq1	0	5	minus	0	-
seq1	2	8	overlap_minus	0	-
seq2	0	4	short_plus	0	+
BED
    printf '>seq1:1-5\nTACGT\n>seq1:1-5\nACGTA\n>seq1:3-8\nGTACGT\n>seq2:1-4\nCCCC\n' > "$bed_exp"
    "$ZFASTA" get "$simple" --bed "$bed_comp" --strand-aware --rc > "$TMPDIR/got.tmp" 2>/dev/null \
        && diff_oracle "$bed_exp" "$TMPDIR/got.tmp" "[parity:seq] BED --strand-aware --rc exact" \
        || fail "[parity:seq] BED --strand-aware --rc exact"

    printf '>seq1:1-5 (reverse complement)\nTACGT\n' > "$TMPDIR/annot_rc.fa"
    "$ZFASTA" get "$simple" seq1:1-5 --rc --annotate-rc > "$TMPDIR/got.tmp" 2>/dev/null \
        && diff_oracle "$TMPDIR/annot_rc.fa" "$TMPDIR/got.tmp" "[extended:header] annotated --rc header" \
        || fail "[extended:header] annotated --rc header"

    if command -v "$BEDTOOLS" &>/dev/null && [[ -x "$SEQTK" ]]; then
        parity_bed_rc_composed "BED stranded RC" "$simple" "$bed_comp" 2 0
        parity_bed_rc_composed "BED stranded RC stdin" "$simple" "$bed_comp" 2 1
        gen_bed_rc_file "$TMPDIR/bed_rc.bed" "$TMPDIR/bed_str.bed"
        parity_bed_rc_composed "BED stranded RC chrom" "$chrom" "$TMPDIR/bed_str.bed" 4096 0
    fi

    oracle rc "$iupac" iupac_all 1 33 "$TMPDIR/iupac_rc_py.fa"
    "$ZFASTA" get "$iupac" iupac_all:1-33 --rc > "$TMPDIR/got.tmp" 2>/dev/null \
        && diff_oracle "$TMPDIR/iupac_rc_py.fa" "$TMPDIR/got.tmp" "[parity:seq] IUPAC --rc Python oracle" \
        || fail "[parity:seq] IUPAC --rc Python oracle"

    "$ZFASTA" get tests/data/proteome.fasta 'sp|P12345|PROT_HUMAN:1-10' --rc > "$TMPDIR/got.tmp" 2> "$TMPDIR/got.err" \
        && fail "protein RC rejected (should fail)" \
        || grep -q 'reverse complement is not defined' "$TMPDIR/got.err" \
        && pass "protein RC rejected" || fail "protein RC rejected (wrong err)"
}

section5() {
    section_hdr 5 "Edge cases and error paths"
    local E="$PROJECT_DIR/tests/data"

    verify_edge_case() {
        local desc="$1" fasta="$2" region="$3" exp_zf="$4" exp_st="$5"
        local exp_noo="${6:--}" exp_rb="${7:--}"
        local tag="[parity:samtools]"
        (( exp_st == 1 && exp_zf == 0 )) && tag="[extended:messy]"

        "$ZFASTA" index "$fasta" 2>/dev/null || true
        [[ -f "${fasta}.fai" ]] || "$SAMTOOLS" faidx "$fasta" 2>/dev/null || true

        local zf_exit=0 st_exit=0 zf_out st_out
        zf_out=$("$ZFASTA" get "$fasta" "$region" 2>/dev/null) || zf_exit=$?
        st_out=$("$SAMTOOLS" faidx "$fasta" "$region" 2>/dev/null) || st_exit=$?

        if (( (exp_zf == 0 && zf_exit == 0) || (exp_zf == 1 && zf_exit != 0) )); then pass "$tag $desc z-fasta exit ok"
        else fail "$tag $desc z-fasta (expected exit $exp_zf, got $zf_exit)"; fi
        if (( (exp_st == 0 && st_exit == 0) || (exp_st == 1 && st_exit != 0) )); then pass "$tag $desc samtools exit ok"
        else fail "$tag $desc samtools (expected exit $exp_st, got $st_exit)"; fi

        if (( zf_exit == 0 && st_exit == 0 )); then
            printf "%s" "$zf_out" > "$TMPDIR/zf_edge.tmp"
            printf "%s" "$st_out" > "$TMPDIR/st_edge.tmp"
            diff_oracle "$TMPDIR/st_edge.tmp" "$TMPDIR/zf_edge.tmp" "$tag $desc byte parity"
        fi
        if (( zf_exit == 0 && st_exit != 0 )); then
            if parse_region_spec "$region"; then
                oracle region "$fasta" "$PARSE_NAME" "$PARSE_START" "$PARSE_END" "$TMPDIR/expected.tmp"
                printf "%s\n" "$zf_out" > "$TMPDIR/zf_edge.tmp"
                diff_oracle "$TMPDIR/expected.tmp" "$TMPDIR/zf_edge.tmp" "$tag $desc Python oracle"
            else
                fail "$tag $desc Python oracle skipped (region has no coords)"
            fi
        fi
        if [[ -x "$NOODLES" && "$exp_noo" != "-" ]]; then
            local noo_exit=0; noo_out=$("$NOODLES" get "$fasta" "$region" 2>/dev/null) || noo_exit=$?
            (( (exp_noo == 0 && noo_exit == 0) || (exp_noo == 1 && noo_exit != 0) )) \
                && pass "$tag $desc noodles exit ok" || fail "$tag $desc noodles (expected $exp_noo, got $noo_exit)"
        fi
        if [[ -x "$RUSTBIO" && "$exp_rb" != "-" ]]; then
            local rb_exit=0; rb_out=$("$RUSTBIO" get "$fasta" "$region" 2>/dev/null) || rb_exit=$?
            (( (exp_rb == 0 && rb_exit == 0) || (exp_rb == 1 && rb_exit != 0) )) \
                && pass "$tag $desc rustbio exit ok" || fail "$tag $desc rustbio (expected $exp_rb, got $rb_exit)"
        fi
    }

    verify_edge_case "empty.fasta" "$E/empty.fasta" "nonexistent" 1 1 1 1
    verify_edge_case "not_fasta.txt" "$E/not_fasta.txt" "simple" 1 1 1 1
    verify_edge_case "edge_cases:single_line" "$E/edge_cases.fasta" "single_line" 0 0 0 0
    verify_edge_case "edge_cases:empty_seq" "$E/edge_cases.fasta" "empty_seq" 1 1 1 1
    verify_edge_case "edge_cases:dupname" "$E/edge_cases.fasta" "dupname" 0 0 0 0
    verify_edge_case "edge_cases:lowercase" "$E/edge_cases.fasta" "lowercase:1-6" 0 0 0 0
    verify_edge_case "edge_cases:nonstandard" "$E/edge_cases.fasta" "nonstandard" 0 0 0 0

    for case in "${MESSY_POS_CASES[@]}"; do
        IFS=: read -r name regions <<< "$case"
        IFS=, read -r first _ <<< "$regions"
        IFS=- read -r start end <<< "$first"
        local tgt="$TMPDIR/e_${name}.fasta"
        [[ -f "$MESSY_DIR/${name}.fasta" ]] || continue
        cp "$MESSY_DIR/${name}.fasta" "$tgt"
        verify_edge_case "messy:$name:$start-$end" "$tgt" "${name}:${start}-${end}" 0 1 1 1
    done

    local dup="$TMPDIR/dup_edge.fasta"
    cp "$E/edge_cases.fasta" "$dup"
    "$ZFASTA" index --no-dedup "$dup" >/dev/null 2>&1 || fail "[extended:dedup] index"
    printf '>dupname\nCCCCCCCC\n' > "$TMPDIR/dup_expected.fa"
    "$ZFASTA" get "$dup" dupname > "$TMPDIR/got.tmp" 2>/dev/null \
        && diff_oracle "$TMPDIR/dup_expected.fa" "$TMPDIR/got.tmp" "[extended:dedup] dupname last record" \
        || fail "[extended:dedup] dupname last record"

    expect_fail "get on nonexistent seq" "$ZFASTA" get "$E/simple.fasta" NOSUCHSEQ
    expect_fail "get on invalid region" "$ZFASTA" get "$E/simple.fasta" seq1:invalid
    expect_fail "get with zero args" "$ZFASTA" get "$E/simple.fasta"
    expect_fail "mutually exclusive --rc and --complement-only" "$ZFASTA" get "$E/simple.fasta" seq1:1-5 --rc --complement-only
    expect_fail "mutually exclusive --rc and --reverse-only" "$ZFASTA" get "$E/simple.fasta" seq1:1-5 --rc --reverse-only
    expect_fail "mutually exclusive --complement-only and --reverse-only" "$ZFASTA" get "$E/simple.fasta" seq1:1-5 --complement-only --reverse-only
    expect_fail "--annotate-rc without transform" "$ZFASTA" get "$E/simple.fasta" seq1:1-5 --annotate-rc
}

# --- main ---
while [[ $# -gt 0 ]]; do
    case $1 in
        --skip-index) SKIP_INDEX=true ;;
        --skip-get) SKIP_GET=true ;;
        --skip-multi) SKIP_MULTI=true ;;
        --skip-bed) SKIP_BED=true ;;
        --skip-rc) SKIP_RC=true ;;
        --skip-edge) SKIP_EDGE=true ;;
        *) echo "Unknown: $1"; exit 1 ;;
    esac
    shift
done

echo "z-fasta get verification"
echo "  z-fasta:  $ZFASTA"
echo "  samtools: $(command -v "$SAMTOOLS" 2>/dev/null || echo missing)"
echo "  bedtools: $(command -v "$BEDTOOLS" 2>/dev/null || echo missing)"
echo "  skip: index=$SKIP_INDEX get=$SKIP_GET multi=$SKIP_MULTI bed=$SKIP_BED rc=$SKIP_RC edge=$SKIP_EDGE"

cd "$PROJECT_DIR"
[[ -x "$ZFASTA" ]] || { echo "Error: z-fasta not found"; exit 1; }
command -v "$SAMTOOLS" &>/dev/null || { echo "Error: samtools not found"; exit 1; }

! $SKIP_INDEX && section0_index
! $SKIP_GET   && section1
! $SKIP_MULTI && section2
! $SKIP_BED   && section3
! $SKIP_RC    && section4
! $SKIP_EDGE  && section5

echo ""
echo "Results: $PASS passed, $FAIL failed"
[[ "$FAIL" -gt 0 ]] && { echo "VERIFICATION FAILED"; exit 1; }
echo "ALL PASSED"
