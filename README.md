# z-fasta ⚡

[![CI](https://github.com/eneskemalergin/z-fasta/actions/workflows/ci.yml/badge.svg)](https://github.com/eneskemalergin/z-fasta/actions/workflows/ci.yml)

A ruthless, zero-allocation, SIMD-accelerated FASTA indexer written in Zig. `z-fasta` is designed to be a drop-in replacement for `samtools faidx`, providing the exact same byte-for-byte output at up to **17x the speed**.

## Why z-fasta?

Modern bioinformatics workflows are bottlenecked by legacy text parsers. `z-fasta` bypasses standard I/O overhead by memory-mapping (`mmap`) the entire FASTA file and using explicit SIMD instructions to scan for sequence headers at the theoretical limit of your NVMe drive.

- **Drop-in replacement:** Emits standard `.fai` format, byte-identical to `samtools faidx`.
- **Single binary:** No dependencies, no `conda` environments, no `glibc` version errors.
- **Safe:** Built with Zig's `ArenaAllocator` to guarantee zero memory leaks.

## Installation

```bash
# Download Zig 0.14.0 (if needed)
curl -L https://ziglang.org/download/0.14.0/zig-linux-x86_64-0.14.0.tar.xz | tar xJ

# Build
zig build -Doptimize=ReleaseFast

# The executable is now at ./zig-out/bin/z-fasta
```

## Usage

```bash
z-fasta index [options] <file.fasta>

Options:
  --emit-fai    Output FAI format to stdout (default: create .zfi binary file)
  --no-dedup    Disable duplicate name filtering (maximizes speed)
  --low-mem     Use chunked reader instead of mmap (limits RAM to 4 MB)
  --help        Show help message
  --version     Print version
```

### Examples

```bash
# Create .zfi binary index (default, compact binary format)
z-fasta index genome.fa

# Output .fai to stdout (samtools-compatible)
z-fasta index --emit-fai genome.fa > genome.fai

# Ultra-low memory mode for massive files on small machines
z-fasta index --low-mem large_metagenome.fa > large_metagenome.fai

# Maximum speed mode (allows duplicate headers)
z-fasta index --emit-fai --no-dedup transcriptome.fa
```

## Performance & Correctness

Tested on an AMD Ryzen 9 3950X using real biological datasets (warm cache).

| Dataset | Size | z-fasta (no-dedup) | samtools | Speedup |
| --- | --- | --- | --- | --- |
| Human Genome | 3.0 GB | 0.57s | 9.15s | **15.9x** |
| Transcriptome | 972 MB | 0.10s | 1.79s | **17.5x** |
| Proteome | 66 MB | 0.005s | 0.05s | **9.4x** |

> See [bench/REPORT.md](bench/REPORT.md) for full results including scaling curves and memory analysis.

### Memory / Execution Modes

| Mode | Speed | Heap Memory | Architecture |
| --- | --- | --- | --- |
| `--no-dedup` | **Fastest** | **< 1 MB** | `mmap` + SIMD. No hash map tracking. |
| `default` | Fast | ~45 MB | `mmap` + SIMD. Tracks duplicate headers. |
| `--low-mem` | Slower | **4.0 MB** | `read()` + 4 MB buffer. Bypasses `mmap`. |

> *Note: `time` utilities report VmRSS equal to the file size for mmap modes because the OS maps the file to virtual memory. Actual private heap allocation is minimal.*

### Correctness

`z-fasta` has been rigorously tested against `samtools faidx` across 20 distinct edge cases (zero-byte files, missing trailing newlines, mixed `\r\n` endings, unicode headers, binary garbage, etc).

**Result: 20 / 20 edge cases match `samtools` behavior (exit codes and outputs).**

## Benchmarking

```bash
# Download real test data (~4 GB)
./bench/download_data.sh

# Run benchmarks (requires hyperfine, samtools)
./bench/run_benchmarks.sh

# Run edge case correctness tests
./bench/run_tests.sh

# Generate report from results
python3 bench/generate_report.py
```

## Output Formats

| Format | Flag | Description |
| --- | --- | --- |
| `.zfi` | *(default)* | Compact binary index. Fast to read/write programmatically. |
| `.fai` | `--emit-fai` | Tab-separated text, identical to `samtools faidx` output. |

## Development

```bash
# Build (debug)
zig build

# Run tests
zig build test

# Build optimized binary
zig build -Doptimize=ReleaseFast
```

## Roadmap

This tool is the foundation for a larger high-performance genomics/proteomics suite:

- `z-fasta get` — O(1) sub-sequence extraction.
- `z-fasta digest` — High-speed in-silico Trypsin digestion for mass spectrometry.

## License

MIT

---

> *Aligned life in bytes,*
> 
> *FASTA sings through mirrored streams*
> 
> *humans bloom as code.*
