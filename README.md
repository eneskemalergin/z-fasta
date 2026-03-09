# z-fasta ⚡

[![CI](https://github.com/eneskemalergin/z-fasta/actions/workflows/ci.yml/badge.svg)](https://github.com/eneskemalergin/z-fasta/actions/workflows/ci.yml)

A ruthless, arena-allocated FASTA toolkit written in Zig. `z-fasta` indexes, extracts, and summarizes FASTA files — SIMD-accelerated indexing up to **17× faster** than `samtools faidx`, sub-millisecond region extraction, and instantaneous assembly statistics from a compact binary index.

## Why z-fasta?

Modern bioinformatics workflows are bottlenecked by legacy text parsers. `z-fasta` bypasses standard I/O overhead by memory-mapping (`mmap`) the entire FASTA file, using explicit SIMD instructions in the indexer to scan for sequence headers at the theoretical limit of your NVMe drive.

- **Drop-in replacement:** Both `z-fasta index --emit-fai` and `z-fasta get` produce output byte-identical to `samtools faidx`. Falls back from `.zfi` to `.fai` with mtime + file-size staleness validation.
- **Single binary:** No dependencies, no `conda` environments, no `glibc` version errors.
- **Arena-allocated:** Uses Zig's `ArenaAllocator` — zero memory leaks, minimal heap overhead in all modes.

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

| Format           | Description                   |
| ---------------- | ----------------------------- |
| `NAME`           | Full sequence                 |
| `NAME:START-END` | 1-based, inclusive sub-region |
| `NAME:START-`    | From START to end of sequence |

Handles Ensembl-style names containing colons (e.g., `chromosome:GRCh38:1:1:248956422:1`).

### Stats

```bash
z-fasta stats [options] <file.fasta>

Options:
  --index-only  Compute stats from index only (no FASTA scan, < 1 ms)
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

All timings on AMD Ryzen 9 3950X, warm cache.

### Index — SIMD-Accelerated Indexing

| Dataset       | Size   | z-fasta (no-dedup) | samtools | Speedup   |
| ------------- | ------ | ------------------ | -------- | --------- |
| Human Genome  | 3.0 GB | 0.57s              | 9.15s    | **15.9×** |
| Transcriptome | 972 MB | 0.10s              | 1.79s    | **17.5×** |
| Proteome      | 66 MB  | 0.005s             | 0.05s    | **9.4×**  |

| Mode         | Heap Memory | Notes                                                                          |
| ------------ | ----------- | ------------------------------------------------------------------------------ |
| `--no-dedup` | **< 1 MB**  | Fastest — `mmap` + SIMD, no deduplication hash map.                            |
| `default`    | ~45 MB      | `mmap` + SIMD, deduplicates sequence names.                                    |
| `--low-mem`  | **4 MB**    | `read()` + fixed 4 MB buffer — no `mmap`, for memory-constrained environments. |

> *`mmap` modes show VmRSS ≈ file size (OS-mapped pages); actual private heap is < 1 MB or ~45 MB as above.*  
> See [bench/index/REPORT.md](bench/index/REPORT.md) for full scaling curves and memory analysis.

### Get — O(1) Region Extraction

| Dataset                | Region          | z-fasta     | samtools | Speedup      |
| ---------------------- | --------------- | ----------- | -------- | ------------ |
| Any (warm cache)       | 100 bp – 10 kbp | **~0.6 ms** | ~1.5 ms  | **2.3–2.5×** |
| Proteome (14 MB)       | 1 kbp region    | 4.0 ms      | 11.3 ms  | **2.9×**     |
| Transcriptome (459 MB) | 1 kbp region    | 128 ms      | 284 ms   | **2.2×**     |

> Region extraction is O(1) regardless of file size — a single `pread` to the precomputed byte offset. Note: fastahack is faster than z-fasta for large (≥50 MB) single full-sequence extraction due to a simpler write path; z-fasta leads on multi-sequence real datasets.  
> See [bench/get/REPORT.md](bench/get/REPORT.md) for full results.

### Stats — Assembly/Proteome Statistics

| Mode       | Dataset              | z-fasta    | seqkit  | Speedup    |
| ---------- | -------------------- | ---------- | ------- | ---------- |
| Index-only | Genome (3.0 GB)      | **0.7 ms** | 2.1 s   | **~3000×** |
| Index-only | Proteome (14 MB)     | **4.9 ms** | 24.6 ms | **5×**     |
| Full scan  | 1 GB single-seq file | 0.89 s     | 0.41 s  | 0.46×      |
| Full scan  | Proteome (14 MB)     | **15 ms**  | 25 ms   | **1.6×**   |

> Index-only time is constant (< 1 ms) regardless of file size — reads only the binary `.zfi` header. Full-scan throughput is ~1.1 GB/s. seqkit is faster on large single-sequence synthetic files; z-fasta leads on multi-sequence files (proteomes, transcriptomes) and computes richer statistics (N50, GC composition, skew, amino acid breakdown) that seqkit does not provide.  
> See [bench/stats/REPORT.md](bench/stats/REPORT.md) for full results.

### Correctness

- **Index:** 20/20 edge cases match `samtools faidx` (exit codes and output).
- **Get:** 90/90 byte-identical diff tests pass across 5 test files — full sequences, sub-regions, single bases, line-boundary spans, clamped ranges.
- **Stats:** 107/107 BioPython verification tests pass — exact agreement on all Tier 1 and Tier 2 values across nucleotide and protein files.
- **Unit tests:** 63/63 Zig unit tests (19 index · 12 get · 32 stats).

## Benchmarking

```bash
# Download real test data (~4 GB, one-time)
bash bench/shared/download_data.sh

# ── Index ─────────────────────────────────────────────────────────
bash bench/index/run_benchmarks.sh       # timing + memory
bash bench/index/run_tests.sh            # 20 edge-case correctness tests
python3 bench/index/generate_report.py   # → bench/index/REPORT.md

# ── Get ───────────────────────────────────────────────────────────
bash bench/get/run_benchmarks.sh         # latency, scaling, real datasets
bash bench/get/verify_get.sh             # 90 byte-identical diff tests vs samtools
python3 bench/get/generate_report.py     # → bench/get/REPORT.md

# ── Stats ─────────────────────────────────────────────────────────
bash bench/stats/run_benchmarks.sh       # full/index-only, scaling, throughput
.venv/bin/python bench/stats/verify_stats.py  # 107 BioPython verification tests
python3 bench/stats/generate_report.py   # → bench/stats/REPORT.md
```

Add `--skip-real` to the `get` / `stats` scripts to skip real dataset runs (~3 GB downloads required otherwise). See [bench/README.md](bench/README.md) for prerequisites and full instructions.

## Output Formats

| Format | Flag         | Description                                                |
| ------ | ------------ | ---------------------------------------------------------- |
| `.zfi` | *(default)*  | Compact binary index. Fast to read/write programmatically. |
| `.fai` | `--emit-fai` | Tab-separated text, identical to `samtools faidx` output.  |

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

**Delivered**

- [x] `z-fasta index` — SIMD-accelerated FASTA indexing (v0.1)
- [x] `z-fasta get` — O(1) byte-offset sequence extraction (v0.2)
- [x] `z-fasta stats` — Assembly/proteome statistics with index-only mode (v0.2)
- [x] Unified benchmark suite with per-module reports and figures (v0.2.2)

**Near-term**

- [ ] Expanded tool comparison across all subcommands in benchmark reports (v0.2.3)
- [ ] Multi-region queries: multiple `NAME:START-END` args or BED file input (v0.3)
- [ ] Reverse complement output flag for `z-fasta get` (v0.3)

**Long-term / Exploratory**

- [ ] `z-fasta validate` — FASTA format validator with detailed error reporting (v0.3+)
- [ ] `z-fasta digest` — In-silico trypsin digestion for mass spectrometry (v0.3+)
- [ ] Zig version upgrade to 0.15+ for async I/O and improved SIMD support (v0.4+)
- [ ] Parallel mmap scanning for multi-threaded indexing on NVMe arrays
- [ ] Native BGZF / gzip streaming read support

## License

MIT — see [LICENSE](LICENSE)

---

> *Aligned life in bytes,*
>
>*FASTA sings through mirrored streams.*  
>
>*Humans bloom as code.*
