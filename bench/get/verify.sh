#!/usr/bin/env bash
# verify.sh - Unified verification for z-fasta get
#
# Consolidates verify_get, verify_multi_get, verify_bed, verify_rc.
# Tests z-fasta get against samtools faidx, bedtools getfasta, and
# expected values. Validates Tier 2 Rust wrappers (noodles, rustbio)
# where feature sets overlap.
#
# Usage: ./verify.sh [--skip-get] [--skip-multi] [--skip-bed] [--skip-rc] [--skip-edge]
set -euo pipefail

# ══════════════════════════════════════════════════════════════════════
#  Config
# ══════════════════════════════════════════════════════════════════════

SKIP_GET=false; SKIP_MULTI=false; SKIP_BED=false; SKIP_RC=false; SKIP_EDGE=false

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
MESSY_DIR="$PROJECT_DIR/bench/index/messy_variants"
TMPDIR="$SCRIPT_DIR/.verify_tmp"
mkdir -p "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT

ZFASTA="${ZFASTA:-$PROJECT_DIR/zig-out/bin/z-fasta}"
SAMTOOLS="${SAMTOOLS:-samtools}"
BEDTOOLS="${BEDTOOLS:-bedtools}"
SEQTK="${SEQTK:-$PROJECT_DIR/tools/seqtk/seqtk}"
NOODLES="${NOODLES:-$PROJECT_DIR/tools/noodles_wrapper/target/release/noodles_wrapper}"
RUSTBIO="${RUSTBIO:-$PROJECT_DIR/tools/rustbio_wrapper/target/release/rustbio_wrapper}"

PASS=0; FAIL=0

# ══════════════════════════════════════════════════════════════════════
#  Helpers
# ══════════════════════════════════════════════════════════════════════

pass() { PASS=$((PASS + 1)); printf "  PASS: %s\n" "$1"; }

fail() {
    FAIL=$((FAIL + 1)); printf "  FAIL: %s\n" "$1"
    local ef="${2:-$TMPDIR/expected.tmp}" gf="${3:-$TMPDIR/got.tmp}"
    echo "    expected:"; head -3 "$ef" 2>/dev/null | sed 's/^/      /'
    echo "    got:";      head -3 "$gf" 2>/dev/null | sed 's/^/      /'
}

t2_ok() { [[ -x "$NOODLES" || -x "$RUSTBIO" ]]; }

ensure_index() {
    local f="$1"
    [[ -f "${f}.fai" ]] || "$SAMTOOLS" faidx "$f" 2>/dev/null
    [[ -f "${f}.zfi" ]] || "$ZFASTA" index "$f" 2>/dev/null
}

bed_to_regions() {
    awk -F'\t' '!/^[[:space:]]*$|^#|^track|^browser/{sub(/\r$/,""); printf "%s:%d-%d\n",$1,$2+1,$3}' "$1"
}

names_to_regions() {
    awk '!/^[[:space:]]*$|^#/{sub(/\r$/,""); print}' "$1"
}

normalize_bed_output() {
    local fasta_in="$1" bed_in="$2" honor_strand="$3" fasta_out="$4"
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
    ' "$fasta_in" > "$fasta_out"
}

write_expected_region() {
    local fasta="$1" name="$2" start="$3" end="$4" out="$5"
    python3 - "$fasta" "$name" "$start" "$end" "$out" <<'PY'
from pathlib import Path; import sys
f,n,s,e,o=sys.argv[1],sys.argv[2],int(sys.argv[3]),int(sys.argv[4]),sys.argv[5]
seqs={}; cur=None
for l in Path(f).read_text().splitlines():
    if l.startswith('>'): cur=l[1:].split()[0]; seqs[cur]=[]
    elif cur: seqs[cur].append(''.join(c for c in l if not c.isspace()))
seq=''.join(seqs[n])[s-1:e]; Path(o).write_text(f'>{n}:{s}-{e}\n{seq}\n')
PY
}

write_expected_rc_region() {
    local fasta="$1" name="$2" start="$3" end="$4" out="$5"
    python3 - "$fasta" "$name" "$start" "$end" "$out" <<'PY'
from pathlib import Path; import sys
f,n,s,e,o=sys.argv[1],sys.argv[2],int(sys.argv[3]),int(sys.argv[4]),sys.argv[5]
t=str.maketrans("ACGTURYSWKMBDHVNacgturyswkmbdhvn","TGCAAYRSWMKVHDBNtgcaayrswmkvhdbn")
seqs={}; cur=None
for l in Path(f).read_text().splitlines():
    if l.startswith('>'): cur=l[1:].split()[0]; seqs[cur]=[]
    elif cur: seqs[cur].append(''.join(c for c in l if not c.isspace()))
seq=''.join(seqs[n])[s-1:e].translate(t)[::-1]; Path(o).write_text(f'>{n}:{s}-{e}\n{seq}\n')
PY
}

write_expected_bed_regions() {
    local fasta="$1" bed="$2" out="$3"
    python3 - "$fasta" "$bed" "$out" <<'PY'
from pathlib import Path; import sys
fasta,bed,out=sys.argv[1],sys.argv[2],sys.argv[3]
seqs={}; cur=None
for l in Path(fasta).read_text().splitlines():
    if l.startswith('>'): cur=l[1:].split()[0]; seqs[cur]=[]
    elif cur: seqs[cur].append(''.join(c for c in l if not c.isspace()))
parts=[]
for l in Path(bed).read_text().splitlines():
    if not l or l.startswith('#') or l.startswith('track') or l.startswith('browser'): continue
    f=l.split('\t'); chrom=f[0]; s0=int(f[1]); e0=int(f[2])
    seq=''.join(seqs[chrom])[s0:e0]
    parts.append(f'>{chrom}:{s0+1}-{e0}\n{seq}\n')
Path(out).write_text(''.join(parts))
PY
}

# ══════════════════════════════════════════════════════════════════════
#  Tier 2 wrapper helpers
# ══════════════════════════════════════════════════════════════════════

tier2_check() {
    local fasta="$1" region="$2" desc="$3"
    [[ -x "$NOODLES" ]] && {
        "$NOODLES" get "$fasta" "$region" > "$TMPDIR/noodles.tmp" 2>/dev/null && {
            diff -q "$TMPDIR/got.tmp" "$TMPDIR/noodles.tmp" >/dev/null 2>&1 && pass "$desc noodles" \
                || fail "$desc noodles" "$TMPDIR/got.tmp" "$TMPDIR/noodles.tmp"
        } || fail "$desc (noodles err)"
    }
    [[ -x "$RUSTBIO" ]] && {
        "$RUSTBIO" get "$fasta" "$region" > "$TMPDIR/rustbio.tmp" 2>/dev/null && {
            diff -q "$TMPDIR/got.tmp" "$TMPDIR/rustbio.tmp" >/dev/null 2>&1 && pass "$desc rustbio" \
                || fail "$desc rustbio" "$TMPDIR/got.tmp" "$TMPDIR/rustbio.tmp"
        } || fail "$desc (rustbio err)"
    }
}

tier2_multi() {
    local desc="$1" fasta="$2"; shift 2
    "$SAMTOOLS" faidx "$fasta" "$@" > "$TMPDIR/expected.tmp" 2>/dev/null || return
    [[ -x "$NOODLES" ]] && {
        "$NOODLES" get "$fasta" "$@" > "$TMPDIR/got.tmp" 2>/dev/null && {
            diff -q "$TMPDIR/expected.tmp" "$TMPDIR/got.tmp" >/dev/null 2>&1 && pass "$desc noodles" \
                || fail "$desc noodles" "$TMPDIR/expected.tmp" "$TMPDIR/got.tmp"
        } || fail "$desc (noodles err)"
    }
    [[ -x "$RUSTBIO" ]] && {
        "$RUSTBIO" get "$fasta" "$@" > "$TMPDIR/got.tmp" 2>/dev/null && {
            diff -q "$TMPDIR/expected.tmp" "$TMPDIR/got.tmp" >/dev/null 2>&1 && pass "$desc rustbio" \
                || fail "$desc rustbio" "$TMPDIR/expected.tmp" "$TMPDIR/got.tmp"
        } || fail "$desc (rustbio err)"
    }
}

tier2_bed() {
    local desc="$1" bed="$2" fasta="${3:-$PROJECT_DIR/tests/data/simple.fasta}"
    local regions="$TMPDIR/st_regions.txt"
    bed_to_regions "$bed" > "$regions"
    "$SAMTOOLS" faidx -r "$regions" "$fasta" > "$TMPDIR/expected.tmp" 2>/dev/null || return
    [[ -x "$NOODLES" ]] && {
        "$NOODLES" get "$fasta" --bed "$bed" > "$TMPDIR/noodles.raw" 2>/dev/null && {
            normalize_bed_output "$TMPDIR/noodles.raw" "$bed" 0 "$TMPDIR/noodles.norm"
            normalize_bed_output "$TMPDIR/expected.tmp" "$bed" 0 "$TMPDIR/st.norm"
            diff -q "$TMPDIR/st.norm" "$TMPDIR/noodles.norm" >/dev/null 2>&1 && pass "$desc noodles" \
                || fail "$desc noodles" "$TMPDIR/st.norm" "$TMPDIR/noodles.norm"
        } || fail "$desc (noodles err)"
    }
    [[ -x "$RUSTBIO" ]] && {
        "$RUSTBIO" get "$fasta" --bed "$bed" > "$TMPDIR/rustbio.raw" 2>/dev/null && {
            normalize_bed_output "$TMPDIR/rustbio.raw" "$bed" 0 "$TMPDIR/rustbio.norm"
            normalize_bed_output "$TMPDIR/expected.tmp" "$bed" 0 "$TMPDIR/st.norm"
            diff -q "$TMPDIR/st.norm" "$TMPDIR/rustbio.norm" >/dev/null 2>&1 && pass "$desc rustbio" \
                || fail "$desc rustbio" "$TMPDIR/st.norm" "$TMPDIR/rustbio.norm"
        } || fail "$desc (rustbio err)"
    }
}

tier2_names() {
    local desc="$1" names_file="$2" fasta="${3:-$PROJECT_DIR/tests/data/simple.fasta}"
    local regions="$TMPDIR/st_names.txt"
    names_to_regions "$names_file" > "$regions"
    "$SAMTOOLS" faidx -r "$regions" "$fasta" > "$TMPDIR/expected.tmp" 2>/dev/null || return
    [[ -x "$NOODLES" ]] && {
        "$NOODLES" get "$fasta" --names "$names_file" > "$TMPDIR/noodles.tmp" 2>/dev/null && {
            diff -q "$TMPDIR/expected.tmp" "$TMPDIR/noodles.tmp" >/dev/null 2>&1 && pass "$desc noodles" \
                || fail "$desc noodles" "$TMPDIR/expected.tmp" "$TMPDIR/noodles.tmp"
        } || fail "$desc (noodles err)"
    }
    [[ -x "$RUSTBIO" ]] && {
        "$RUSTBIO" get "$fasta" --names "$names_file" > "$TMPDIR/rustbio.tmp" 2>/dev/null && {
            diff -q "$TMPDIR/expected.tmp" "$TMPDIR/rustbio.tmp" >/dev/null 2>&1 && pass "$desc rustbio" \
                || fail "$desc rustbio" "$TMPDIR/expected.tmp" "$TMPDIR/rustbio.tmp"
        } || fail "$desc (rustbio err)"
    }
}

tier2_rc() {
    local fasta="$1" region="$2" desc="$3"
    local st
    st=$("$SAMTOOLS" faidx -i --mark-strand no "$fasta" "$region" 2>/dev/null) || return
    printf "%s\n" "$st" > "$TMPDIR/expected.tmp"
    [[ -x "$NOODLES" ]] && {
        r=$("$NOODLES" get "$fasta" "$region" --rc 2>/dev/null) && {
            printf "%s\n" "$r" > "$TMPDIR/got.tmp"
            diff -q "$TMPDIR/expected.tmp" "$TMPDIR/got.tmp" >/dev/null 2>&1 && pass "$desc noodles" \
                || fail "$desc noodles" "$TMPDIR/expected.tmp" "$TMPDIR/got.tmp"
        } || fail "$desc (noodles err)"
    }
    [[ -x "$RUSTBIO" ]] && {
        r=$("$RUSTBIO" get "$fasta" "$region" --rc 2>/dev/null) && {
            printf "%s\n" "$r" > "$TMPDIR/got.tmp"
            diff -q "$TMPDIR/expected.tmp" "$TMPDIR/got.tmp" >/dev/null 2>&1 && pass "$desc rustbio" \
                || fail "$desc rustbio" "$TMPDIR/expected.tmp" "$TMPDIR/got.tmp"
        } || fail "$desc (rustbio err)"
    }
}

# ══════════════════════════════════════════════════════════════════════
#  Fixture generators
# ══════════════════════════════════════════════════════════════════════

gen_names_file() {
    local out="$1"
    printf '# names file\nseq2\n\nseq1\nseq2\n' > "$out"
}

gen_bed_file() {
    local out="$1" count="$2"
    : > "$out"
    printf '# synthetic BED for %s\n' "$count" >> "$out"
    printf 'track name=verify_bed_%s\n' "$count" >> "$out"
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
    local out="$1"
    [[ -f "${out}.zfi" ]] && return
    printf '>iupac_all\nACGTURYSWKMBDHVNacgturyswkmbdhvnu\n' > "$out"
    "$ZFASTA" index "$out" >/dev/null 2>&1
}

gen_chrom_fixture() {
    local out="$1" total="${2:-16384}"
    [[ -f "${out}.zfi" ]] && return
    {
        echo ">chrSynthetic"
        awk -v total="$total" 'BEGIN{pattern="ACGTNRYWSKMBDHVacgtnrywskmbdhv"; s=""
            while(length(s)<total){s=s pattern} s=substr(s,1,total)
            while(length(s)>0){print substr(s,1,71); s=substr(s,72)}}'
    } > "$out"
    "$ZFASTA" index "$out" >/dev/null 2>&1
}

gen_bed_rc_file() {
    local out="$1" stranded_out="$2"
    : > "$out"; : > "$stranded_out"
    for ((i=0; i<1200; i++)); do
        s=$((1000+i*400)); e=$((s+120))
        printf '%s\t%d\t%d\n' "chrSynthetic" "$s" "$e" >> "$out"
        if ((i%2==0)); then strand='+'; else strand='-'; fi
        printf '%s\t%d\t%d\tregion_%05d\t0\t%s\n' "chrSynthetic" "$s" "$e" "$i" "$strand" >> "$stranded_out"
    done
}

# ══════════════════════════════════════════════════════════════════════
#  Global fixture data
# ══════════════════════════════════════════════════════════════════════

REG20=(); for i in {1..10}; do REG20+=("seq1:${i}-${i}" "seq2:${i}-${i}"); done
REG20_REV=(); for i in {10..1}; do REG20_REV+=("seq2:${i}-${i}" "seq1:${i}-${i}"); done

# ══════════════════════════════════════════════════════════════════════
#  [1] Single-Region Extraction
# ══════════════════════════════════════════════════════════════════════

section1() {
    echo ""; echo "═══════════════════════════════════════════════════"
    echo " [1] Single-Region Extraction"
    echo "═══════════════════════════════════════════════════"

    verify() {
        local fasta="$1" region="$2" desc="$3"
        "$SAMTOOLS" faidx "$fasta" "$region" > "$TMPDIR/expected.tmp" 2>/dev/null || { fail "$desc (samtools err)"; return; }
        "$ZFASTA" get "$fasta" "$region" > "$TMPDIR/got.tmp" 2>/dev/null      || { fail "$desc (z-fasta err)"; return; }
        diff -q "$TMPDIR/expected.tmp" "$TMPDIR/got.tmp" >/dev/null 2>&1 && pass "$desc" || fail "$desc"
    }

    verify_expected() {
        local fasta="$1" region="$2" desc="$3" expected="$4"
        printf "%s" "$expected" > "$TMPDIR/expected.tmp"
        "$ZFASTA" get "$fasta" "$region" > "$TMPDIR/got.tmp" 2>/dev/null || { fail "$desc (z-fasta err)"; return; }
        diff -q "$TMPDIR/expected.tmp" "$TMPDIR/got.tmp" >/dev/null 2>&1 && pass "$desc" || fail "$desc"
    }

    test_seq() {
        local fasta="$1" name="$2" len="$3"
        local label="${fasta##*/}:${name}"
        verify "$fasta" "$name" "$label full"
        tier2_check "$fasta" "$name" "$label full Tier 2"
        if ((len>=10)); then
            verify "$fasta" "${name}:1-10" "$label :1-10"
            tier2_check "$fasta" "${name}:1-10" "$label :1-10 Tier 2"
            local ls=$((len-9))
            verify "$fasta" "${name}:${ls}-${len}" "$label :${ls}-${len}"
            tier2_check "$fasta" "${name}:${ls}-${len}" "$label :${ls}-${len} Tier 2"
        else
            verify "$fasta" "${name}:1-${len}" "$label :1-${len}"
        fi
        if ((len>=101)); then
            local ms=$((len/2-50)) me=$((len/2+50))
            ((ms<1)) && ms=1; ((me>len)) && me=$len
            verify "$fasta" "${name}:${ms}-${me}" "$label mid-span"
            tier2_check "$fasta" "${name}:${ms}-${me}" "$label mid-span Tier 2"
        fi
        verify "$fasta" "${name}:1-1" "$label :1-1"
        verify "$fasta" "${name}:${len}-${len}" "$label :${len}-${len}"
        verify "$fasta" "${name}:1-${len}" "$label :1-${len}"
        verify "$fasta" "${name}:1-$((len+100))" "$label :1-$((len+100)) (clamp)"
    }

    test_file() {
        local fasta="$1"
        [[ -f "${fasta}.fai" ]] || return
        echo "  ${fasta##*/}"; echo ""
        while IFS=$'\t' read -r name length _; do ((length>0)) && test_seq "$fasta" "$name" "$length"; done < "${fasta}.fai"
    }

    for f in tests/data/simple.fasta tests/data/proteome.fasta tests/data/single.fasta tests/data/edge_cases.fasta; do
        [[ -f "$f" ]] && { ensure_index "$f"; test_file "$f"; }
    done

    cp tests/data/mixed_widths.fasta "$TMPDIR/mixed_widths.fasta"
    ensure_index "$TMPDIR/mixed_widths.fasta"
    test_file "$TMPDIR/mixed_widths.fasta"

    for fixture in mixed_line_widths:3:24 trailing_whitespace:7:20 blank_lines:5:16 mixed_crlf_lf:3:18; do
        IFS=: read -r name start end <<< "$fixture"
        local src="$MESSY_DIR/${name}.fasta"; tgt="$TMPDIR/${name}.fasta"
        [[ -f "$src" ]] || continue
        cp "$src" "$tgt"
        "$ZFASTA" index "$tgt" >/dev/null 2>&1 || continue
        write_expected_region "$tgt" "$name" "$start" "$end" "$TMPDIR/expected.tmp"
        "$ZFASTA" get "$tgt" "${name}:${start}-${end}" > "$TMPDIR/got.tmp" 2>/dev/null || continue
        diff -q "$TMPDIR/expected.tmp" "$TMPDIR/got.tmp" >/dev/null 2>&1 && pass "$name ${start}-${end}" || fail "$name ${start}-${end}"
    done
}

# ══════════════════════════════════════════════════════════════════════
#  [2] Multi-Region Extraction
# ══════════════════════════════════════════════════════════════════════

section2() {
    echo ""; echo "═══════════════════════════════════════════════════"
    echo " [2] Multi-Region Extraction"
    echo "═══════════════════════════════════════════════════"

    verify_multi() {
        local desc="$1" fasta="$2"; shift 2
        "$SAMTOOLS" faidx "$fasta" "$@" > "$TMPDIR/expected.tmp" 2>/dev/null || { fail "$desc (samtools err)"; return; }
        "$ZFASTA" get "$fasta" "$@"        > "$TMPDIR/got.tmp"      2>/dev/null || { fail "$desc (z-fasta err)"; return; }
        diff -q "$TMPDIR/expected.tmp" "$TMPDIR/got.tmp" >/dev/null 2>&1 && pass "$desc" || fail "$desc"
    }

    local simple="tests/data/simple.fasta" proteome="tests/data/proteome.fasta" edge="tests/data/edge_cases.fasta"
    local LONG="AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"

    echo "  simple.fasta"
    verify_multi "same seq, two sub-regions" "$simple" "seq1:1-10" "seq1:13-24"
    verify_multi "seqs file order" "$simple" "seq1:1-12" "seq2:1-6"
    verify_multi "seqs reversed order" "$simple" "seq2:1-6" "seq1:1-12"
    verify_multi "overlapping" "$simple" "seq1:1-15" "seq1:10-24"
    verify_multi "duplicate" "$simple" "seq1:1-12" "seq1:1-12"
    verify_multi "full + sub" "$simple" "seq1" "seq2:3-10"
    verify_multi "both full" "$simple" "seq1" "seq2"
    verify_multi "triple duplicate" "$simple" "seq1:1-5" "seq1:1-5" "seq1:1-5"
    verify_multi "full + 2 sub" "$simple" "seq1" "seq2:1-6" "seq1:10-20"
    verify_multi "all reversed" "$simple" "seq2" "seq1"
    verify_multi "20 regions sort path" "$simple" "${REG20[@]}"
    verify_multi "20 regions reversed" "$simple" "${REG20_REV[@]}"

    tier2_multi "multi same seq" "$simple" "seq1:1-10" "seq1:13-24"
    tier2_multi "multi file order" "$simple" "seq1:1-12" "seq2:1-6"
    tier2_multi "multi reversed" "$simple" "seq2:1-6" "seq1:1-12"
    tier2_multi "multi full+sub" "$simple" "seq1" "seq2:3-10"
    tier2_multi "multi both full" "$simple" "seq1" "seq2"
    tier2_multi "multi 20 regions" "$simple" "${REG20[@]}"

    echo "  edge_cases.fasta"
    verify_multi "long + normal" "$edge" "single_line" "$LONG"
    verify_multi "long reversed" "$edge" "$LONG" "single_line"
    verify_multi "lowercase + nonstandard" "$edge" "lowercase" "nonstandard"
    verify_multi "overlapping same seq" "$edge" "single_line:1-5" "single_line:3-8"

    echo "  proteome.fasta"
    verify_multi "pipe names" "$proteome" "sp|P12345|PROT_HUMAN" "sp|Q98765|ANOT_MOUSE"
    verify_multi "pipe reversed" "$proteome" "sp|Q98765|ANOT_MOUSE" "sp|P12345|PROT_HUMAN"
    verify_multi "pipe sub + full" "$proteome" "sp|P12345|PROT_HUMAN:1-10" "sp|Q98765|ANOT_MOUSE"
    verify_multi "pipe duplicate" "$proteome" "sp|P12345|PROT_HUMAN" "sp|P12345|PROT_HUMAN"

    local prot1="sp|P12345|PROT_HUMAN" prot2="sp|Q98765|ANOT_MOUSE" P16=()
    for i in {1..8}; do P16+=("${prot1}:${i}-$((i+2))" "${prot2}:${i}-$((i+2))"); done
    verify_multi "16 proteome sub-regions (sort path)" "$proteome" "${P16[@]}"

    ALL_EDGE=(); while IFS=$'\t' read -r name _; do ALL_EDGE+=("$name"); done < "${edge}.fai"
    (( ${#ALL_EDGE[@]} > 0 )) && {
        verify_multi "all edge seqs" "$edge" "${ALL_EDGE[@]}"
        REV_EDGE=(); for ((i=${#ALL_EDGE[@]}-1; i>=0; i--)); do REV_EDGE+=("${ALL_EDGE[$i]}"); done
        verify_multi "all edge reversed" "$edge" "${REV_EDGE[@]}"
    }
}

# ══════════════════════════════════════════════════════════════════════
#  [3] BED & Names Extraction
# ══════════════════════════════════════════════════════════════════════

section3() {
    echo ""; echo "═══════════════════════════════════════════════════"
    echo " [3] BED & Names Extraction"
    echo "═══════════════════════════════════════════════════"

    local simple_fasta="$PROJECT_DIR/tests/data/simple.fasta"
    ensure_index "$simple_fasta"

    local BED_SMALL="$TMPDIR/small.bed" BED_MEDIUM="$TMPDIR/medium.bed" BED_LARGE="$TMPDIR/large.bed"
    gen_bed_file "$BED_SMALL" 10; gen_bed_file "$BED_MEDIUM" 100; gen_bed_file "$BED_LARGE" 1000

    verify_case() {
        local desc="$1" bed="$2" hs="$3" fasta="${4:-$simple_fasta}" chunk="${5:-4096}" stdin="${6:-0}"
        local bt_cmd="bedtools getfasta -fi $fasta -bed $bed"; [[ "$hs" == "1" ]] && bt_cmd+=" -s"
        local zf_args=(get "$fasta" --bed "$bed" --chunk-size "$chunk"); [[ "$hs" == "1" ]] && zf_args+=(--strand-aware)
        if [[ "$stdin" == "1" ]]; then
            zf_args=(${zf_args[@]/"--bed $bed"/"--bed -"})
            $bt_cmd 2>/dev/null > "$TMPDIR/bt.raw" || { fail "$desc (bedtools err)"; return; }
            cat "$bed" | "$ZFASTA" "${zf_args[@]}" > "$TMPDIR/zf.raw" 2>/dev/null || { fail "$desc (z-fasta err)"; return; }
        else
            $bt_cmd 2>/dev/null > "$TMPDIR/bt.raw" || { fail "$desc (bedtools err)"; return; }
            "$ZFASTA" "${zf_args[@]}" > "$TMPDIR/zf.raw" 2>/dev/null || { fail "$desc (z-fasta err)"; return; }
        fi
        normalize_bed_output "$TMPDIR/bt.raw" "$bed" "$hs" "$TMPDIR/bt.norm" || return
        normalize_bed_output "$TMPDIR/zf.raw" "$bed" "$hs" "$TMPDIR/zf.norm" || return
        diff -q "$TMPDIR/bt.norm" "$TMPDIR/zf.norm" >/dev/null 2>&1 && pass "$desc" || fail "$desc"
    }

    verify_samtools() {
        local desc="$1" bed="$2" fasta="${3:-$simple_fasta}" chunk="${4:-4096}" stdin="${5:-0}"
        local zf_args=(get "$fasta" --chunk-size "$chunk")
        local regions="$TMPDIR/st_regions.txt"
        bed_to_regions "$bed" > "$regions"
        if [[ "$stdin" == "1" ]]; then
            zf_args+=(--bed -)
            cat "$bed" | "$ZFASTA" "${zf_args[@]}" > "$TMPDIR/got.tmp" 2>/dev/null || { fail "$desc (z-fasta err)"; return; }
        else
            zf_args+=(--bed "$bed")
            "$ZFASTA" "${zf_args[@]}" > "$TMPDIR/got.tmp" 2>/dev/null || { fail "$desc (z-fasta err)"; return; }
        fi
        "$SAMTOOLS" faidx -r "$regions" "$fasta" > "$TMPDIR/expected.tmp" 2>/dev/null || { fail "$desc (samtools err)"; return; }
        diff -q "$TMPDIR/expected.tmp" "$TMPDIR/got.tmp" >/dev/null 2>&1 && pass "$desc" || fail "$desc"
    }

    verify_names() {
        local desc="$1" names_file="$2"
        local regions="$TMPDIR/st_names.txt"
        names_to_regions "$names_file" > "$regions"
        "$ZFASTA" get "$simple_fasta" --names "$names_file" > "$TMPDIR/got.tmp" 2>/dev/null || { fail "$desc (z-fasta err)"; return; }
        "$SAMTOOLS" faidx -r "$regions" "$simple_fasta" > "$TMPDIR/expected.tmp" 2>/dev/null || { fail "$desc (samtools err)"; return; }
        diff -q "$TMPDIR/expected.tmp" "$TMPDIR/got.tmp" >/dev/null 2>&1 && pass "$desc" || fail "$desc"
    }

    verify_expected_bed_case() {
        local desc="$1" fasta="$2" bed="$3" expected_file="$4"
        "$ZFASTA" get "$fasta" --bed "$bed" > "$TMPDIR/got.tmp" 2>/dev/null || { fail "$desc (z-fasta err)"; return; }
        diff -q "$expected_file" "$TMPDIR/got.tmp" >/dev/null 2>&1 && pass "$desc" || fail "$desc"
    }

    run_size_suite() {
        local label="$1" count="$2" chunk="$3"
        local bed="$TMPDIR/${label}.bed"
        gen_bed_file "$bed" "$count"
        echo ""
        echo "=== ${label} (${count} BED rows, chunk=${chunk}) ==="
        verify_case "${label} default" "$bed" 0 "$simple_fasta" "$chunk"
        verify_case "${label} stranded" "$bed" 1 "$simple_fasta" "$chunk"
        verify_samtools "${label} vs samtools" "$bed" "$simple_fasta" "$chunk"
        verify_case "${label} via stdin" "$bed" 0 "$simple_fasta" "$chunk" 1
        verify_samtools "${label} via stdin vs samtools" "$bed" "$simple_fasta" "$chunk" 1
        verify_case "${label} stranded via stdin" "$bed" 1 "$simple_fasta" "$chunk" 1
    }

    local NAMES="$TMPDIR/names.txt"
    gen_names_file "$NAMES"

    echo "  simple.fasta sized suites"
    run_size_suite "small" 10 3
    run_size_suite "medium" 100 97
    run_size_suite "large" 1000 257

    echo ""; echo "=== names-file extraction ==="
    verify_names "names file vs samtools" "$NAMES"

    echo ""; echo "=== mixed-width FASTA BED ==="
    local mx_fasta="$TMPDIR/mx.fasta" mx_bed="$TMPDIR/mx.bed"
    cp "$PROJECT_DIR/tests/data/mixed_widths.fasta" "$mx_fasta"
    ensure_index "$mx_fasta"
    cat > "$mx_bed" <<'BED'
mixed1	54	75	mixed1_span	0	+
mixed2	74	95	mixed2_span	0	-
mixed3	57	72	mixed3_span	0	+
BED
    if command -v "$BEDTOOLS" &>/dev/null; then
        verify_case "mixed-width default" "$mx_bed" 0 "$mx_fasta"
        verify_case "mixed-width stranded" "$mx_bed" 1 "$mx_fasta"
    fi
    verify_samtools "mixed-width vs samtools" "$mx_bed" "$mx_fasta"
    tier2_bed "mixed-width vs samtools" "$mx_bed" "$mx_fasta"

    echo ""; echo "=== non-uniform FASTA BED ==="
    local nu_fasta="$TMPDIR/nu.fasta" nu_bed="$TMPDIR/nu.bed" nu_expected="$TMPDIR/nu_expected.fa"
    cp "$MESSY_DIR/mixed_line_widths.fasta" "$nu_fasta"
    "$ZFASTA" index "$nu_fasta" >/dev/null 2>&1
    cat > "$nu_bed" <<'BED'
mixed_line_widths	2	24	span1	0	+
mixed_line_widths	12	28	span2	0	+
BED
    write_expected_bed_regions "$nu_fasta" "$nu_bed" "$nu_expected"
    verify_expected_bed_case "non-uniform BED expected" "$nu_fasta" "$nu_bed" "$nu_expected"

    echo ""; echo "=== Tier 2 Rust wrappers (BED / names) ==="
    tier2_bed "small BED vs samtools" "$BED_SMALL"
    tier2_bed "medium BED vs samtools" "$BED_MEDIUM"
    tier2_bed "mixed-width BED vs samtools" "$mx_bed" "$mx_fasta"
    tier2_names "names vs samtools" "$NAMES"

    # Names file with region syntax should error (--names is name-only, not regions)
    echo ""; echo "=== names-file region syntax ==="
    local names_reg="$TMPDIR/names_reg.txt"
    printf 'seq1:1-5\n' > "$names_reg"
    "$ZFASTA" get "$simple_fasta" --names "$names_reg" > "$TMPDIR/got.tmp" 2>/dev/null && \
        fail "names file region syntax should error" || pass "names file region syntax rejected"

    # BED edge states: error paths
    echo ""; echo "=== BED error paths ==="
    local bad_bed="$TMPDIR/bad.bed"
    printf 'seq1\t5\t5\tzero_len\t0\t+\n' > "$bad_bed"
    "$ZFASTA" get "$simple_fasta" --bed "$bad_bed" > "$TMPDIR/got.tmp" 2>/dev/null && \
        fail "0-length BED interval" || pass "0-length BED rejected"

    printf 'seq1\t5\n' > "$bad_bed"
    "$ZFASTA" get "$simple_fasta" --bed "$bad_bed" > "$TMPDIR/got.tmp" 2>/dev/null && \
        fail "short BED (<3 cols)" || pass "short BED rejected"

    printf 'seq1\t0\t5\tname\t0\t?\n' > "$bad_bed"
    "$ZFASTA" get "$simple_fasta" --bed "$bad_bed" --strand-aware > "$TMPDIR/got.tmp" 2>/dev/null && \
        fail "invalid BED strand" || pass "invalid BED strand rejected"
}

# ══════════════════════════════════════════════════════════════════════
#  [4] Reverse Complement
# ══════════════════════════════════════════════════════════════════════

section4() {
    echo ""; echo "═══════════════════════════════════════════════════"
    echo " [4] Reverse Complement"
    echo "═══════════════════════════════════════════════════"

    local simple="$PROJECT_DIR/tests/data/simple.fasta"
    local single="$PROJECT_DIR/tests/data/single.fasta"
    local edge="$PROJECT_DIR/tests/data/edge_cases.fasta"
    local iupac="$TMPDIR/iupac.fasta" chrom="$TMPDIR/chrom.fasta"
    gen_iupac_fixture "$iupac"; ensure_index "$iupac"
    gen_chrom_fixture "$chrom" 500000; ensure_index "$chrom"

    # --rc single region vs samtools
    verify_rc() {
        local fasta="$1" region="$2" desc="$3"
        "$SAMTOOLS" faidx -i --mark-strand no "$fasta" "$region" > "$TMPDIR/expected.tmp" 2>/dev/null || { fail "$desc (samtools err)"; return; }
        "$ZFASTA" get "$fasta" "$region" --rc > "$TMPDIR/got.tmp" 2>/dev/null || { fail "$desc (z-fasta err)"; return; }
        diff -q "$TMPDIR/expected.tmp" "$TMPDIR/got.tmp" >/dev/null 2>&1 && pass "$desc" || fail "$desc"
    }

    # exact expected output
    verify_exact_output() {
        local desc="$1" expected_file="$2"; shift 2
        "$ZFASTA" "$@" > "$TMPDIR/got.tmp" 2>"$TMPDIR/got.err" || { fail "$desc (z-fasta err)"; return; }
        diff -q "$expected_file" "$TMPDIR/got.tmp" >/dev/null 2>&1 && pass "$desc" || fail "$desc" "$expected_file" "$TMPDIR/got.tmp"
    }

    # multi-region --rc self-consistency (batch vs per-region)
    verify_rc_chain() {
        local desc="$1" fasta="$2"; shift 2
        : > "$TMPDIR/expected.tmp"
        for r in "$@"; do "$ZFASTA" get "$fasta" "$r" --rc >> "$TMPDIR/expected.tmp" 2>/dev/null; done
        "$ZFASTA" get "$fasta" "$@" --rc > "$TMPDIR/got.tmp" 2>/dev/null
        diff -q "$TMPDIR/expected.tmp" "$TMPDIR/got.tmp" >/dev/null 2>&1 && pass "$desc" || fail "$desc"
    }

    # names-file --rc vs samtools
    verify_rc_names() {
        local desc="$1" fasta="$2" nf="$3"
        local regions="$TMPDIR/st_rc_names.txt"
        names_to_regions "$nf" > "$regions"
        "$SAMTOOLS" faidx -i --mark-strand no -r "$regions" "$fasta" > "$TMPDIR/expected.tmp" 2>/dev/null || { fail "$desc (samtools err)"; return; }
        "$ZFASTA" get "$fasta" --names "$nf" --rc > "$TMPDIR/got.tmp" 2>/dev/null || { fail "$desc (z-fasta err)"; return; }
        diff -q "$TMPDIR/expected.tmp" "$TMPDIR/got.tmp" >/dev/null 2>&1 && pass "$desc" || fail "$desc"
    }

    # BED + strand-aware + RC vs bedtools+seqtk
    verify_bed_rc_composed() {
        local desc="$1" fasta="$2" bed="$3" chunk="$4" stdin="$5"
        local zf_args=(get "$fasta" --strand-aware --rc --chunk-size "$chunk")
        local bt=$("$BEDTOOLS" getfasta -fi "$fasta" -bed "$bed" -s 2>/dev/null) || { fail "$desc (bedtools err)"; return; }
        local mapped=$(printf "%s\n" "$bt" | "$SEQTK" seq -r 2>/dev/null) || { fail "$desc (seqtk err)"; return; }
        printf "%s\n" "$mapped" > "$TMPDIR/bt.raw"
        if [[ "$stdin" == "1" ]]; then
            cat "$bed" | "$ZFASTA" "${zf_args[@]}" --bed - > "$TMPDIR/zf.raw" 2>/dev/null || { fail "$desc (z-fasta err)"; return; }
        else
            "$ZFASTA" "${zf_args[@]}" --bed "$bed" > "$TMPDIR/zf.raw" 2>/dev/null || { fail "$desc (z-fasta err)"; return; }
        fi
        normalize_bed_output "$TMPDIR/bt.raw" "$bed" 0 "$TMPDIR/bt.norm"
        normalize_bed_output "$TMPDIR/zf.raw" "$bed" 0 "$TMPDIR/zf.norm"
        diff -q "$TMPDIR/bt.norm" "$TMPDIR/zf.norm" >/dev/null 2>&1 && pass "$desc" || fail "$desc"
    }

    # protein RC rejection
    verify_protein_error() {
        local desc="$1"
        "$ZFASTA" get tests/data/proteome.fasta 'sp|P12345|PROT_HUMAN:1-10' --rc > "$TMPDIR/got.tmp" 2> "$TMPDIR/got.err" \
            && { fail "$desc (should fail)"; return; }
        grep -q 'reverse complement is not defined' "$TMPDIR/got.err" && pass "$desc" || fail "$desc (wrong err)"
    }

    # ── simple.fasta RC ──
    echo "  simple.fasta"
    for r in seq1:1-5 seq1:10-15 seq2:1-12; do
        verify_rc "$simple" "$r" "simple $r"
        tier2_rc "$simple" "$r" "simple $r"
    done

    # ── single.fasta RC ──
    echo "  single.fasta"
    verify_rc "$single" "single_sequence:1-4" "single span"
    tier2_rc "$single" "single_sequence:1-4" "single span"

    # ── edge_cases.fasta lowercase RC ──
    echo "  edge_cases.fasta"
    verify_rc "$edge" "lowercase:1-12" "lowercase preservation"

    # ── IUPAC / chrom fixtures RC ──
    echo "  IUPAC/chrom fixtures"
    verify_rc "$iupac" "iupac_all:1-33" "IUPAC full"
    verify_rc "$chrom" "chrSynthetic" "chrom full"
    tier2_rc "$iupac" "iupac_all:1-33" "IUPAC full"
    tier2_rc "$chrom" "chrSynthetic" "chrom full"

    # ── mixed-width RC ──
    echo "  mixed-width RC"
    local mx="$TMPDIR/mx_rc.fasta"
    cp "$PROJECT_DIR/tests/data/mixed_widths.fasta" "$mx"
    ensure_index "$mx"
    verify_rc "$mx" "mixed1:55-75" "mixed-width region"

    # ── non-uniform RC exact output ──
    echo "  non-uniform RC"
    local nu="$TMPDIR/nu_rc.fasta" nu_exp="$TMPDIR/nu_rc_expected.fa"
    cp "$MESSY_DIR/mixed_line_widths.fasta" "$nu"
    "$ZFASTA" index "$nu" >/dev/null 2>&1
    write_expected_rc_region "$nu" "mixed_line_widths" 3 24 "$nu_exp"
    verify_exact_output "non-uniform --rc exact" "$nu_exp" get "$nu" mixed_line_widths:3-24 --rc

    # ── IUPAC --complement-only exact ──
    local iupac_comp_exp="$TMPDIR/iupac_comp.fa"
    printf '>iupac_all:1-33 (complement)\nTGCAAYRSWMKVHDBNtgcaayrswmkvhdbna\n' > "$iupac_comp_exp"
    verify_exact_output "IUPAC --complement-only exact" "$iupac_comp_exp" \
        get "$iupac" iupac_all:1-33 --complement-only --annotate-rc

    # ── multi-region RC self-consistency ──
    echo "  multi-region RC"
    verify_rc_chain "multi RC same seq" "$simple" seq1:1-5 seq1:10-15 seq1:20-24
    verify_rc_chain "multi RC mixed" "$simple" seq2:1-4 seq1:1-5 seq2:5-12 seq1:10-15
    verify_rc_chain "multi RC 20-region" "$simple" "${REG20_REV[@]}"

    # ── names-file RC ──
    echo "  names-file RC"
    local nf="$TMPDIR/names_rc.txt"
    printf 'seq2\nseq1\n' > "$nf"
    verify_rc_names "names --rc vs samtools" "$simple" "$nf"

    # ── --annotate-rc / --reverse-only / --complement-only exact ──
    local a1="$TMPDIR/annot_rc.fa" a2="$TMPDIR/reverse.fa" a3="$TMPDIR/complement.fa"
    printf '>seq1:1-5 (reverse complement)\nTACGT\n' > "$a1"
    printf '>seq1:1-5 (reverse)\nATGCA\n' > "$a2"
    printf '>seq1:1-5 (complement)\nTGCAT\n' > "$a3"
    verify_exact_output "annotated --rc header" "$a1" get "$simple" seq1:1-5 --rc --annotate-rc
    verify_exact_output "annotated --reverse-only" "$a2" get "$simple" seq1:1-5 --reverse-only --annotate-rc
    verify_exact_output "annotated --complement-only" "$a3" get "$simple" seq1:1-5 --complement-only --annotate-rc

    # ── BED --strand-aware --rc composition exact ──
    echo "  BED stranded RC exact"
    local bed_comp="$TMPDIR/bed_comp.bed" bed_comp_exp="$TMPDIR/bed_comp_expected.fa"
    cat > "$bed_comp" <<'BED'
seq1	0	5	plus	0	+
seq1	0	5	minus	0	-
seq1	2	8	overlap_minus	0	-
seq2	0	4	short_plus	0	+
BED
    printf '>seq1:1-5\nTACGT\n>seq1:1-5\nACGTA\n>seq1:3-8\nGTACGT\n>seq2:1-4\nCCCC\n' > "$bed_comp_exp"
    verify_exact_output "BED --strand-aware --rc" "$bed_comp_exp" \
        get "$simple" --bed "$bed_comp" --strand-aware --rc

    # ── BED stranded RC vs bedtools+seqtk ──
    echo "  BED stranded RC vs bedtools+seqtk"
    if command -v "$BEDTOOLS" &>/dev/null && [[ -x "$SEQTK" ]]; then
        verify_bed_rc_composed "BED stranded RC" "$simple" "$bed_comp" 2 0
        verify_bed_rc_composed "BED stranded RC via stdin" "$simple" "$bed_comp" 2 1
    fi

    # ── BED stranded RC on chrom fixture (large) ──
    echo "  BED stranded RC chrom fixture"
    local bed_rc_file="$TMPDIR/bed_rc.bed" bed_str="$TMPDIR/bed_str.bed"
    gen_bed_rc_file "$bed_rc_file" "$bed_str"
    if command -v "$BEDTOOLS" &>/dev/null && [[ -x "$SEQTK" ]]; then
        verify_bed_rc_composed "BED stranded RC chrom" "$chrom" "$bed_str" 4096 0
    fi

    # ── Tier 2 wrappers: BED + RC (wrappers' --honor-strand is broken on BED) ──
    echo "  Tier 2 BED + RC"
    tier2_bed_rc() {
        local desc="$1" fasta="$2" bed="$3"
        local regions="$TMPDIR/st_regions.txt"
        bed_to_regions "$bed" > "$regions"
        "$SAMTOOLS" faidx -i --mark-strand no -r "$regions" "$fasta" > "$TMPDIR/expected.tmp" 2>/dev/null || return
        normalize_bed_output "$TMPDIR/expected.tmp" "$bed" 0 "$TMPDIR/st.norm"
        [[ -x "$NOODLES" ]] && {
            r=$("$NOODLES" get "$fasta" --bed "$bed" --rc 2>/dev/null) && {
                printf "%s\n" "$r" > "$TMPDIR/noodles.raw"
                normalize_bed_output "$TMPDIR/noodles.raw" "$bed" 0 "$TMPDIR/noodles.norm"
                diff -q "$TMPDIR/st.norm" "$TMPDIR/noodles.norm" >/dev/null 2>&1 && pass "$desc noodles" \
                    || fail "$desc noodles" "$TMPDIR/st.norm" "$TMPDIR/noodles.norm"
            } || pass "$desc (noodles --rc not supported)"
        }
        [[ -x "$RUSTBIO" ]] && {
            r=$("$RUSTBIO" get "$fasta" --bed "$bed" --rc 2>/dev/null) && {
                printf "%s\n" "$r" > "$TMPDIR/rustbio.raw"
                normalize_bed_output "$TMPDIR/rustbio.raw" "$bed" 0 "$TMPDIR/rustbio.norm"
                diff -q "$TMPDIR/st.norm" "$TMPDIR/rustbio.norm" >/dev/null 2>&1 && pass "$desc rustbio" \
                    || fail "$desc rustbio" "$TMPDIR/st.norm" "$TMPDIR/rustbio.norm"
            } || pass "$desc (rustbio --rc not supported)"
        }
    }
    tier2_bed_rc "BED --rc" "$simple" "$bed_comp"

    # ── IUPAC complement cross-validation vs Python ──
    echo "  IUPAC complement cross-validation"
    cross_check_iupac() {
        local region="$1" desc="$2"
        write_expected_rc_region "$iupac" "iupac_all" 1 33 "$TMPDIR/iupac_rc_py.fa"
        "$ZFASTA" get "$iupac" "$region" --rc > "$TMPDIR/got.tmp" 2>/dev/null
        diff -q "$TMPDIR/iupac_rc_py.fa" "$TMPDIR/got.tmp" >/dev/null 2>&1 && pass "$desc z-fasta" \
            || fail "$desc z-fasta" "$TMPDIR/iupac_rc_py.fa" "$TMPDIR/got.tmp"
        [[ -x "$NOODLES" ]] && {
            r=$("$NOODLES" get "$iupac" "$region" --rc 2>/dev/null) && {
                printf "%s\n" "$r" > "$TMPDIR/noodles.tmp"
                diff -q "$TMPDIR/iupac_rc_py.fa" "$TMPDIR/noodles.tmp" >/dev/null 2>&1 && pass "$desc noodles" \
                    || fail "$desc noodles" "$TMPDIR/iupac_rc_py.fa" "$TMPDIR/noodles.tmp"
            }
        }
        [[ -x "$RUSTBIO" ]] && {
            r=$("$RUSTBIO" get "$iupac" "$region" --rc 2>/dev/null) && {
                printf "%s\n" "$r" > "$TMPDIR/rustbio.tmp"
                diff -q "$TMPDIR/iupac_rc_py.fa" "$TMPDIR/rustbio.tmp" >/dev/null 2>&1 && pass "$desc rustbio" \
                    || fail "$desc rustbio" "$TMPDIR/iupac_rc_py.fa" "$TMPDIR/rustbio.tmp"
            }
        }
    }
    cross_check_iupac "iupac_all:1-33" "IUPAC revcomp"

    # ── protein RC rejection ──
    echo "  protein RC rejection"
    verify_protein_error "protein RC rejected"
}

# ══════════════════════════════════════════════════════════════════════
#  [5] Edge Case & Fail-State Matrix
# ══════════════════════════════════════════════════════════════════════

section5() {
    echo ""; echo "═══════════════════════════════════════════════════"
    echo " [5] Edge Case & Fail-State Matrix"
    echo "═══════════════════════════════════════════════════"
    echo "  Tests which tools accept/reject each edge-case FASTA"

    # verify_edge_case <desc> <fasta> <region> <exp_zf> <exp_st> [exp_noo] [exp_rb]
    # exp_*: 0 = expected to pass (exit 0), 1 = expected to fail (exit != 0)
    # When both zf and st pass, output is diffed.
    # When only zf passes (z-fasta-only), output is verified against Python.
    verify_edge_case() {
        local desc="$1" fasta="$2" region="$3" exp_zf="$4" exp_st="$5"
        local exp_noo="${6:--}" exp_rb="${7:--}"

        "$ZFASTA" index "$fasta" 2>/dev/null || true  # ensure .zfi exists if possible
        [[ -f "${fasta}.fai" ]] || "$SAMTOOLS" faidx "$fasta" 2>/dev/null || true

        local zf_exit=0 st_exit=0 zf_out st_out
        zf_out=$("$ZFASTA" get "$fasta" "$region" 2>/dev/null) || zf_exit=$?
        st_out=$("$SAMTOOLS" faidx "$fasta" "$region" 2>/dev/null) || st_exit=$?

        local ok=true
        # z-fasta
        if (( exp_zf == 0 && zf_exit == 0 )) || (( exp_zf == 1 && zf_exit != 0 )); then
            pass "$desc z-fasta"
        else
            fail "$desc z-fasta (expected exit $exp_zf, got $zf_exit)"
            ok=false
        fi
        # samtools
        if (( exp_st == 0 && st_exit == 0 )) || (( exp_st == 1 && st_exit != 0 )); then
            pass "$desc samtools"
        else
            fail "$desc samtools (expected exit $exp_st, got $st_exit)"
            ok=false
        fi
        # diff when both pass
        if (( zf_exit == 0 && st_exit == 0 )); then
            printf "%s" "$zf_out" > "$TMPDIR/zf_edge.tmp"
            printf "%s" "$st_out" > "$TMPDIR/st_edge.tmp"
            diff -q "$TMPDIR/st_edge.tmp" "$TMPDIR/zf_edge.tmp" >/dev/null 2>&1 && pass "$desc output match" \
                || fail "$desc output match" "$TMPDIR/st_edge.tmp" "$TMPDIR/zf_edge.tmp"
        fi
        # z-fasta-only: verify against Python expected
        if (( zf_exit == 0 && st_exit != 0 )); then
            # region is like "name:start-end" - parse it
            local rname="${region%%:*}" rstart="${region#*:}"; rstart="${rstart%%-*}"
            local rend="${region#*-}"
            write_expected_region "$fasta" "$rname" "$rstart" "$rend" "$TMPDIR/expected.tmp"
            printf "%s\n" "$zf_out" > "$TMPDIR/zf_edge.tmp"
            diff -q "$TMPDIR/expected.tmp" "$TMPDIR/zf_edge.tmp" >/dev/null 2>&1 && pass "$desc z-fasta-only verified" \
                || fail "$desc z-fasta-only" "$TMPDIR/expected.tmp" "$TMPDIR/zf_edge.tmp"
        fi
        # noodles (if available)
        if [[ -x "$NOODLES" && "$exp_noo" != "-" ]]; then
            local noo_exit=0
            noo_out=$("$NOODLES" get "$fasta" "$region" 2>/dev/null) || noo_exit=$?
            if (( exp_noo == 0 && noo_exit == 0 )) || (( exp_noo == 1 && noo_exit != 0 )); then
                pass "$desc noodles"
            else
                fail "$desc noodles (expected exit $exp_noo, got $noo_exit)"
            fi
        fi
        # rustbio (if available)
        if [[ -x "$RUSTBIO" && "$exp_rb" != "-" ]]; then
            local rb_exit=0
            rb_out=$("$RUSTBIO" get "$fasta" "$region" 2>/dev/null) || rb_exit=$?
            if (( exp_rb == 0 && rb_exit == 0 )) || (( exp_rb == 1 && rb_exit != 0 )); then
                pass "$desc rustbio"
            else
                fail "$desc rustbio (expected exit $exp_rb, got $rb_exit)"
            fi
        fi
    }

    verify_error_case() {
        local desc="$1"; shift
        local got_exit=0
        "$@" >/dev/null 2>/dev/null || got_exit=$?
        if (( got_exit != 0 )); then
            pass "$desc"
        else
            fail "$desc (expected error, exit $got_exit)"
        fi
    }

    echo ""
    echo "  --- Extraction from edge-case files ---"
    E=$PROJECT_DIR/tests/data
    verify_edge_case "empty.fasta" "$E/empty.fasta" "nonexistent" 1 1 1 1
    verify_edge_case "not_fasta.txt" "$E/not_fasta.txt" "simple" 1 1 1 1
    verify_edge_case "edge_cases:single_line" "$E/edge_cases.fasta" "single_line" 0 0 0 0
    verify_edge_case "edge_cases:empty_seq" "$E/edge_cases.fasta" "empty_seq" 1 1 1 1
    verify_edge_case "edge_cases:dupname" "$E/edge_cases.fasta" "dupname" 0 0 0 0
    verify_edge_case "edge_cases:lowercase" "$E/edge_cases.fasta" "lowercase:1-6" 0 0 0 0
    verify_edge_case "edge_cases:nonstandard" "$E/edge_cases.fasta" "nonstandard" 0 0 0 0

    echo ""
    echo "  --- Messy FASTA (z-fasta-only files) ---"
    for fixture in mixed_line_widths:3:24 trailing_whitespace:7:20 blank_lines:5:16 mixed_crlf_lf:3:18; do
        IFS=: read -r name start end <<< "$fixture"
        local src="$MESSY_DIR/${name}.fasta" tgt="$TMPDIR/e_${name}.fasta"
        [[ -f "$src" ]] || continue
        cp "$src" "$tgt"
        verify_edge_case "messy:$name:$start-$end" "$tgt" "${name}:${start}-${end}" 0 1 1 1
    done

    echo ""
    echo "  --- Error paths ---"
    verify_error_case "get on nonexistent seq" "$ZFASTA" get "$E/simple.fasta" NOSUCHSEQ
    verify_error_case "get on invalid region" "$ZFASTA" get "$E/simple.fasta" seq1:invalid
    verify_error_case "get with zero args" "$ZFASTA" get "$E/simple.fasta"
}

# ══════════════════════════════════════════════════════════════════════
#  Main
# ══════════════════════════════════════════════════════════════════════

while [[ $# -gt 0 ]]; do
    case $1 in --skip-get)   SKIP_GET=true   ;; --skip-multi) SKIP_MULTI=true ;;
              --skip-bed)    SKIP_BED=true    ;; --skip-rc)    SKIP_RC=true    ;;
              --skip-edge)   SKIP_EDGE=true   ;;
              *) echo "Unknown: $1"; exit 1 ;;
    esac; shift
done

echo "z-fasta get unified verification"
echo "═══════════════════════════════════"
echo "  z-fasta: $ZFASTA"
echo "  samtools: $(command -v "$SAMTOOLS" 2>/dev/null || echo not found)"
echo "  bedtools: $(command -v "$BEDTOOLS" 2>/dev/null || echo not found)"
echo "  seqtk:    $([ -x "$SEQTK" ] && echo found || echo not found)"
echo "  noodles:  $([ -x "$NOODLES" ] && echo found || echo not found)"
echo "  rustbio:  $([ -x "$RUSTBIO" ] && echo found || echo not found)"
echo "  skip: get=$SKIP_GET multi=$SKIP_MULTI bed=$SKIP_BED rc=$SKIP_RC edge=$SKIP_EDGE"
echo ""

cd "$PROJECT_DIR"

[[ -x "$ZFASTA" ]] || { echo "Error: z-fasta not found at $ZFASTA"; exit 1; }
command -v "$SAMTOOLS" &>/dev/null || { echo "Error: samtools not found"; exit 1; }

! $SKIP_GET   && section1
! $SKIP_MULTI && section2
! $SKIP_BED   && section3
! $SKIP_RC    && section4
! $SKIP_EDGE  && section5

echo ""
echo "═══════════════════════════════════════════════════════════════"
echo "  Results: $PASS passed, $FAIL failed"
echo "═══════════════════════════════════════════════════════════════"
[[ "$FAIL" -gt 0 ]] && { echo "VERIFICATION FAILED"; exit 1; }
echo "ALL PASSED"
