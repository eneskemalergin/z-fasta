#!/usr/bin/env bash
set -euo pipefail

version=${1:?usage: package_release.sh VERSION SYSTEM EXE}
system=${2:?usage: package_release.sh VERSION SYSTEM EXE}
exe=${3:?usage: package_release.sh VERSION SYSTEM EXE}

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
    echo "release binary too large (${size} bytes); expected stripped ReleaseFast (<1 MiB)" >&2
    exit 1
fi
if [[ "$system" == linux_* && "$size" -gt 750000 ]]; then
    echo "linux release binary too large (${size} bytes); expected <=750 KiB stripped RF" >&2
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

if [[ "$exe" == *.exe ]]; then
    archive="z-fasta_${version}_${system}.zip"
    tar -C "$dist" -a -cf "release-out/${archive}" "$exe" LICENSE
else
    archive="z-fasta_${version}_${system}.tar.gz"
    chmod +x "${dist}/${exe}"
    tar -C "$dist" -czf "release-out/${archive}" "$exe" LICENSE
fi

tar -xf "release-out/${archive}" -C "$smoke"
(
    cd "$smoke"
    if [[ "${RUNNER_OS:-}" == "Windows" ]]; then
        export PATH="/usr/bin:/bin:${SYSTEMROOT}/System32"
    else
        export PATH="/usr/bin:/bin:/usr/local/bin"
    fi
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
    test -s test.fai
    "./${exe}" get test.fasta seq1:1-4 > get.out
    test -s get.out
    "./${exe}" stats test.fasta > stats.out
    test -s stats.out
)
