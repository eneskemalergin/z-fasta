#!/usr/bin/env bash
# Print CHANGELOG.md body for "## [X.Y.Z]". Fails if missing, empty, or Unreleased.
# Usage: changelog_section.sh <changelog> <version>
set -euo pipefail

changelog=${1:?}
version=${2:?}

[[ "$version" =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]] || {
  echo "version must be X.Y.Z: $version" >&2
  exit 1
}
[[ -f "$changelog" ]] || {
  echo "missing $changelog" >&2
  exit 1
}

header=$(grep -E "^## \[${version}\]" "$changelog" || true)
[[ -n "$header" ]] || {
  echo "missing ## [${version}] in $changelog" >&2
  exit 1
}
[[ "$header" != *[Uu]nreleased* ]] || {
  echo "still Unreleased: $header" >&2
  exit 1
}

body=$(awk -v ver="$version" '
  /^## \[/ {
    if (on) exit
    if (index($0, "[" ver "]")) { on = 1; next }
  }
  on { print }
' "$changelog")
body=$(printf '%s\n' "$body" | sed '/./,$!d')

[[ -n "$body" ]] || {
  echo "empty section for ${version}" >&2
  exit 1
}
printf '%s\n' "$body"
