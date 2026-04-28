#!/bin/bash
# z-fasta startup-floor benchmark harness.
# Builds tiny Zig 0.16 probes in /tmp and compares them with z-fasta CLI startup.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH_ROOT="$(dirname "$SCRIPT_DIR")"
PROJECT_ROOT="$(dirname "$BENCH_ROOT")"
RESULTS_DIR="$SCRIPT_DIR/results"
ZIG="$PROJECT_ROOT/zig-0.16.0/zig"
ZFASTA="$PROJECT_ROOT/zig-out/bin/z-fasta"

RUNS=30
WARMUP=5

while [[ $# -gt 0 ]]; do
    case "$1" in
        --runs) RUNS="$2"; shift 2 ;;
        --warmup) WARMUP="$2"; shift 2 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

command -v hyperfine &>/dev/null || { echo "Error: hyperfine not found"; exit 1; }
[[ -x "$ZIG" ]] || { echo "Error: Zig 0.16 binary not found at $ZIG"; exit 1; }

if [[ ! -x "$ZFASTA" ]]; then
    echo "z-fasta binary not found; building ReleaseFast..."
    (cd "$PROJECT_ROOT" && "$ZIG" build -Doptimize=ReleaseFast) >/dev/null
fi

mkdir -p "$RESULTS_DIR"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
JSON_OUT="$RESULTS_DIR/startup_${TIMESTAMP}.json"

TMP_DIR="${TMPDIR:-/tmp}/z-fasta-startup-probes"
rm -rf "$TMP_DIR"
mkdir -p "$TMP_DIR"

cat > "$TMP_DIR/empty_main.zig" <<'ZIG'
pub fn main() void {}
ZIG

cat > "$TMP_DIR/process_init_empty.zig" <<'ZIG'
const std = @import("std");

pub fn main(init: std.process.Init) void {
    _ = init;
}
ZIG

cat > "$TMP_DIR/process_args.zig" <<'ZIG'
const std = @import("std");

pub fn main(init: std.process.Init) void {
    var args = std.process.Args.Iterator.init(init.minimal.args);
    _ = args.skip();
    while (args.next()) |_| {}
}
ZIG

cat > "$TMP_DIR/std_io_write.zig" <<'ZIG'
const std = @import("std");

pub fn main(init: std.process.Init) void {
    std.Io.File.writeStreamingAll(.stdout(), init.io, "probe\n") catch {};
}
ZIG

cat > "$TMP_DIR/posix_write.zig" <<'ZIG'
const std = @import("std");

pub fn main() void {
    _ = std.os.linux.write(1, "probe\n".ptr, "probe\n".len);
}
ZIG

build_probe() {
    local name="$1"
    "$ZIG" build-exe "$TMP_DIR/${name}.zig" -O ReleaseFast -femit-bin="$TMP_DIR/$name" >/dev/null
}

build_probe empty_main
build_probe process_init_empty
build_probe process_args
build_probe std_io_write
build_probe posix_write

echo "Startup benchmark"
echo "  Runs: $RUNS | Warmup: $WARMUP"
echo "  Output: $JSON_OUT"
echo ""

hyperfine \
    --warmup "$WARMUP" \
    --runs "$RUNS" \
    --export-json "$JSON_OUT" \
    -n true '/bin/true' \
    -n zig-empty "$TMP_DIR/empty_main" \
    -n zig-process-init "$TMP_DIR/process_init_empty" \
    -n zig-process-args "$TMP_DIR/process_args --version get stats index" \
    -n zig-std-io-write "$TMP_DIR/std_io_write >/dev/null" \
    -n zig-posix-write "$TMP_DIR/posix_write >/dev/null" \
    -n z-fasta-version "$ZFASTA --version >/dev/null" \
    -n z-fasta-help "$ZFASTA --help >/dev/null"

echo ""
echo "Wrote $JSON_OUT"
