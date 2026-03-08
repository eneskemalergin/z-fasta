# z-fasta ⚡

[![CI](https://github.com/eneskemalergin/z-fasta/actions/workflows/ci.yml/badge.svg)](https://github.com/eneskemalergin/z-fasta/actions/workflows/ci.yml)

A ruthless, zero-allocation, SIMD-accelerated FASTA toolkit written in Zig. `z-fasta` indexes, extracts, and summarizes FASTA files — a drop-in replacement for `samtools faidx` at up to **17x the speed**.

## Why z-fasta?

Modern bioinformatics workflows are bottlenecked by legacy text parsers. `z-fasta` bypasses standard I/O overhead by memory-mapping (`mmap`) the entire FASTA file and using explicit SIMD instructions to scan for sequence headers at the theoretical limit of your NVMe drive.

- **Drop-in replacement:** `z-fasta get` output is byte-identical to `samtools faidx`.
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

### Index

```bash
z-fasta index [options] <file.fasta>

Options:
  --emit-fai    Output FAI format to stdout (default: create .zfi binary file)
  --no-dedup    Disable duplicate name filtering (maximizes speed)
  --low-mem     Use chunked reader instead of mmap (limits RAM to 4 MB)
  --help        Show help message
  --version     Print version
```

### Get (sequence extraction)

```bash
z-fasta get <file.fasta> <region>
```

Extract sequences or sub-regions from an indexed FASTA file. Output is **byte-identical** to `samtools faidx`.

Requires an index — either `.zfi` (preferred) or `.fai`. If `.zfi` is not found, falls back to `.fai` automatically.

**Region formats:**

| Format | Description |
| --- | --- |
| `NAME` | Full sequence |
| `NAME:START-END` | 1-based, inclusive sub-region |
| `NAME:START-` | From START to end of sequence |

Handles Ensembl-style names containing colons (e.g., `chromosome:GRCh38:1:1:248956422:1`).

### Stats

```bash
z-fasta stats [options] <file.fasta>

Options:
  --index-only  Compute stats from index only (no FASTA scan, <500 μs)
```

Compute assembly/proteome statistics. Automatically detects nucleotide vs. protein sequences.

**Tier 1 (index-only):** sequence count, total bases, min/max/mean/median lengths, N50, L50, N90, L90, AU, duplicate count.

**Tier 2 (default):** full composition scan — nucleotide frequencies, GC content (N excluded), GC skew, soft-masked fraction. For proteins: top 3 amino acids with full names.

### Examples

```bash
# Create .zfi binary index (default, compact binary format)
z-fasta index genome.fa

# Output .fai to stdout (samtools-compatible)
z-fasta index --emit-fai genome.fa > genome.fai

# Extract a full sequence
z-fasta get genome.fa chr1

# Extract a sub-region (1-based, inclusive)
z-fasta get genome.fa chr1:1000000-2000000

# Assembly stats (full composition scan)
z-fasta stats genome.fa

# Quick stats from index only (sub-millisecond)
z-fasta stats --index-only genome.fa
```

## Performance & Correctness

Tested on an AMD Ryzen 9 3950X using real biological datasets (warm cache).

| Dataset | Size | z-fasta (no-dedup) | samtools | Speedup |
| --- | --- | --- | --- | --- |
| Human Genome | 3.0 GB | 0.57s | 9.15s | **15.9x** |
| Transcriptome | 972 MB | 0.10s | 1.79s | **17.5x** |
| Proteome | 66 MB | 0.005s | 0.05s | **9.4x** |

> See [bench/index/REPORT.md](bench/index/REPORT.md) for full results including scaling curves and memory analysis.

### Memory / Execution Modes

| Mode | Speed | Heap Memory | Architecture |
| --- | --- | --- | --- |
| `--no-dedup` | **Fastest** | **< 1 MB** | `mmap` + SIMD. No hash map tracking. |
| `default` | Fast | ~45 MB | `mmap` + SIMD. Tracks duplicate headers. |
| `--low-mem` | Slower | **4.0 MB** | `read()` + 4 MB buffer. Bypasses `mmap`. |

> *Note: `time` utilities report VmRSS equal to the file size for mmap modes because the OS maps the file to virtual memory. Actual private heap allocation is minimal.*

### Correctness

`z-fasta` has been rigorously tested:

- **Index:** 20/20 edge cases match `samtools` behavior (exit codes and outputs).
- **Get:** 90/90 samtools diff tests pass — byte-identical extraction across 5 test files, covering full sequences, sub-regions, single bases, line-boundary spans, and clamped ranges.
- **Stats:** 107/107 BioPython verification tests pass — exact agreement on all Tier 1 and Tier 2 values across nucleotide and protein files.
- **Unit tests:** 80/80 Zig unit tests (19 index + 29 get + 32 stats).

## Benchmarking

```bash
# Download real test data (~4 GB)
bash bench/shared/download_data.sh

# Run indexer benchmarks (requires hyperfine, samtools)
bash bench/index/run_benchmarks.sh

# Run indexer edge case correctness tests
bash bench/index/run_tests.sh

# Verify get command against samtools (byte-identical diff)
bash bench/get/verify_get.sh

# Verify stats against BioPython
.venv/bin/python bench/stats/verify_stats.py

# Generate indexer benchmark report
python3 bench/index/generate_report.py
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

# Run all tests (index + get + stats)
zig build test --summary all

# Build optimized binary
zig build -Doptimize=ReleaseFast
```

## Roadmap

- [x] `z-fasta index` — SIMD-accelerated FASTA indexing (v0.1)
- [x] `z-fasta get` — O(1) byte-offset sequence extraction (v0.2)
- [x] `z-fasta stats` — Assembly/proteome statistics (v0.2)
- [ ] Multi-region queries, reverse complement, BED input (v0.3)
- [ ] `z-fasta digest` — In-silico trypsin digestion for mass spectrometry

## License

MIT

---

> *Aligned life in bytes,*
> 
> *FASTA sings through mirrored streams*
> 
> *humans bloom as code.*
