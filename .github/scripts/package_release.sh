#!/usr/bin/env bash
set -euo pipefail

version=${1:?usage: package_release.sh VERSION SYSTEM}
system=${2:?usage: package_release.sh VERSION SYSTEM}
exe=z-fasta

bin="zig-out/bin/${exe}"
[[ -f "$bin" ]] || {
    echo "missing release binary: $bin" >&2
    exit 1
}

got=$("$bin" --version)
[[ "$got" == "z-fasta ${version}" ]] || {
    echo "binary --version: got '$got'" >&2
    exit 1
}

size=$(wc -c < "$bin" | tr -d ' ')
if [[ "$size" -gt 1048576 ]]; then
    echo "release binary too large (${size} bytes); expected stripped ReleaseFast (<=1048576 bytes)" >&2
    exit 1
fi
if [[ "$system" == linux_* ]] && command -v file >/dev/null 2>&1; then
    info=$(file -b "$bin" || true)
    if echo "$info" | grep -qiE 'with debug_info|not stripped'; then
        echo "linux release binary still has debug info: $info" >&2
        exit 1
    fi
fi

work_dir=.release-work
dist="${work_dir}/dist"
smoke="${work_dir}/smoke"
trap 'rm -rf "$work_dir"' EXIT
rm -rf "$work_dir" release-out
mkdir -p "$dist" "$smoke" release-out
cp "$bin" "${dist}/${exe}"
cp LICENSE "${dist}/LICENSE"

archive="z-fasta_${version}_${system}.tar.gz"
chmod +x "${dist}/${exe}"
tar -C "$dist" -czf "release-out/${archive}" "$exe" LICENSE

tar -xf "release-out/${archive}" -C "$smoke"
(
    cd "$smoke"
    export PATH="/usr/bin:/bin:/usr/local/bin"
    command -v zig >/dev/null 2>&1 && {
        echo "zig still on PATH" >&2
        exit 1
    }

    printf '>seq1\nACGT\n>seq2\nGGGG\n' > test.fasta
    "./${exe}" --help >/dev/null
    "./${exe}" --version
    "./${exe}" validate test.fasta >/dev/null
    "./${exe}" index test.fasta
    "./${exe}" index --emit-fai test.fasta > test.fai
    printf 'seq1\t4\t6\t4\t5\nseq2\t4\t17\t4\t5\n' > expected.fai
    cmp expected.fai test.fai
    "./${exe}" get test.fasta seq1:1-4 > get.out
    printf '>seq1:1-4\nACGT\n' > expected.get
    cmp expected.get get.out
    "./${exe}" stats test.fasta > stats.out
    grep -Fx '  indexed_records: 2' stats.out
    grep -Fx '  total_symbols: 8' stats.out
    grep -Fx '  n50: 4' stats.out
    grep -Fx '  gc: 75.00%' stats.out
    status=0
    "./${exe}" stats --error-format=json > failure.out 2> failure.json || status=$?
    [[ "$status" -eq 1 ]]
    test ! -s failure.out
    grep -F '"kind":"failure"' failure.json >/dev/null
    grep -F '"code":"usage.missing_argument"' failure.json >/dev/null
)
