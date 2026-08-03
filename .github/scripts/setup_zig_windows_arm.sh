#!/usr/bin/env bash
set -euo pipefail

archive=zig-x86_64-windows-0.16.0.zip

curl -fsSL -o "$archive" "https://ziglang.org/download/0.16.0/${archive}"
unzip -q "$archive" -d zig-sdk
echo "$GITHUB_WORKSPACE/zig-sdk/zig-x86_64-windows-0.16.0" >> "$GITHUB_PATH"
