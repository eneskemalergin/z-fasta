#!/usr/bin/env bash
# Builds the repository's social preview: the README logo over the light tagline on white, 1280x640, 8-bit PNG.
#   tools/social_preview.sh [OUT]    (default assets/social-preview.png; needs ImageMagick 7's magick)
# GitHub's API cannot set the preview: upload the file under Settings > General > Social preview.
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
out=${1:-$root/assets/social-preview.png}
work=$(mktemp -d /var/tmp/z-fasta-preview.XXXXXX)
magick -background none -density 384 "$root/assets/logo-readme.svg" -resize x330 "$work/logo.png"
magick -background none -density 192 "$root/assets/tagline-light.svg" -resize 1080x "$work/tagline.png"
magick -size 1280x640 xc:white "$work/logo.png" -gravity north -geometry +0+95 -composite \
    "$work/tagline.png" -gravity north -geometry +0+470 -composite -strip -depth 8 "$out"
rm "$work/logo.png" "$work/tagline.png"
rmdir "$work"
