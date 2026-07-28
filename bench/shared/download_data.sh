#!/usr/bin/env bash
#
# Download real datasets for benchmarking z-fasta (~4.1 GB uncompressed).
#
# Reads bench/shared/datasets.manifest (URL, expected size, sha256).
# Skip-if-present: existing file with matching byte size is left alone.
# Fresh downloads always check size and sha256. Pass --verify to hash-check
# files that were skipped by size.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=tools.sh
source "$SCRIPT_DIR/tools.sh"

DATA_DIR="$SCRIPT_DIR/data"
MANIFEST="$SCRIPT_DIR/datasets.manifest"
VERIFY_EXISTING=false

usage() {
    cat <<'EOF'
Usage: bash bench/shared/download_data.sh [--verify]

  --verify   Also sha256-check files that already match the expected size.
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --verify) VERIFY_EXISTING=true; shift ;;
        -h|--help) usage; exit 0 ;;
        *)
            echo "error: unknown option: $1" >&2
            usage >&2
            exit 1
            ;;
    esac
done

if [[ ! -f "$MANIFEST" ]]; then
    echo "error: manifest missing: $MANIFEST" >&2
    exit 1
fi

mkdir -p "$DATA_DIR"
cd "$DATA_DIR"

file_sha256() {
    local path="$1"
    if command -v sha256sum >/dev/null 2>&1; then
        sha256sum -- "$path" | awk '{print $1}'
    elif command -v shasum >/dev/null 2>&1; then
        shasum -a 256 -- "$path" | awk '{print $1}'
    else
        echo "error: need sha256sum or shasum to verify datasets" >&2
        exit 1
    fi
}

verify_file() {
    local file="$1" expect_bytes="$2" expect_sha="$3"
    local got_bytes got_sha
    got_bytes="$(file_size_bytes "$file")"
    if [[ "$got_bytes" != "$expect_bytes" ]]; then
        echo "error: $file size $got_bytes (expected $expect_bytes)" >&2
        return 1
    fi
    got_sha="$(file_sha256 "$file")"
    if [[ "$got_sha" != "$expect_sha" ]]; then
        echo "error: $file sha256 $got_sha (expected $expect_sha)" >&2
        return 1
    fi
    return 0
}

download_one() {
    local id="$1" file="$2" expect_bytes="$3" expect_sha="$4" url="$5"
    local gz="${file}.gz"
    local got_bytes

    echo "[$id] $file"
    if [[ -f "$file" ]]; then
        got_bytes="$(file_size_bytes "$file")"
        if [[ "$got_bytes" != "$expect_bytes" ]]; then
            echo "error: $file exists with size $got_bytes (expected $expect_bytes)." >&2
            echo "       Delete it and re-run download_data.sh to fetch the pinned revision." >&2
            exit 1
        fi
        if $VERIFY_EXISTING; then
            echo "  Verifying sha256..."
            verify_file "$file" "$expect_bytes" "$expect_sha" || exit 1
            echo "  ok (size + sha256)"
        else
            echo "  Skipping: already present ($got_bytes bytes)"
        fi
        return 0
    fi

    echo "  Downloading..."
    curl -fL --retry 3 --retry-delay 2 -o "$gz" "$url"
    gunzip -f "$gz"
    if [[ ! -f "$file" ]]; then
        echo "error: gunzip did not produce $file" >&2
        exit 1
    fi
    echo "  Verifying size + sha256..."
    if ! verify_file "$file" "$expect_bytes" "$expect_sha"; then
        rm -f -- "$file" "$gz"
        exit 1
    fi
    # Prior sidecars may belong to a different revision; rebuild below.
    rm -f -- "${file}.zfi" "${file}.fai"
    echo "  Downloaded: $file ($(du -h "$file" | cut -f1))"
}

echo "=== Downloading Real Datasets ==="
echo "Target:   $DATA_DIR"
echo "Manifest: $MANIFEST"
echo

entry_count=0
declare -a MANIFEST_FILES=()
while IFS=$'\t' read -r id filename bytes sha256 url || [[ -n "${id:-}" ]]; do
    [[ -z "${id:-}" || "$id" == \#* ]] && continue
    if [[ -z "${filename:-}" || -z "${bytes:-}" || -z "${sha256:-}" || -z "${url:-}" ]]; then
        echo "error: incomplete manifest row for id=${id:-?}" >&2
        exit 1
    fi
    download_one "$id" "$filename" "$bytes" "$sha256" "$url"
    MANIFEST_FILES+=("$filename")
    entry_count=$((entry_count + 1))
    echo
done <"$MANIFEST"

if [[ "$entry_count" -eq 0 ]]; then
    echo "error: no dataset rows in $MANIFEST" >&2
    exit 1
fi

echo "=== Download Complete ==="
echo
du -h REAL_* 2>/dev/null | grep -E '\.(fa|fasta)$' || echo "No REAL_* FASTA files found"

if [[ -x "$ZFASTA" ]]; then
    echo
    echo "=== Reference indexes (.zfi + .fai) ==="
    for fa in "${MANIFEST_FILES[@]}"; do
        [[ -f "$fa" ]] || continue
        if [[ ! -f "${fa}.zfi" ]]; then
            echo "  indexing $fa -> .zfi"
            "$ZFASTA" index "$fa"
        fi
        if [[ ! -f "${fa}.fai" ]] && bench_has_tool samtools; then
            echo "  indexing $fa -> .fai (samtools)"
            "$SAMTOOLS" faidx "$fa"
        fi
    done
fi
