<!-- markdownlint-disable MD033 MD036 MD041 -->
<div align="center">
  <h1>z-fasta ⚡</h1>
  <p>
    Fast, modular FASTA/FASTQ toolkit built in Zig.<br/>
    SIMD-accelerated indexing, O(1) region extraction, and instant assembly stats.<br/>
    a drop-in replacement for <code>samtools faidx</code>, <code>seqkit</code>, and <code>fastahack</code>.
  </p>
  <p>Current version: <strong>v0.2.3</strong></p>
  <br/>
  <a href="https://github.com/eneskemalergin/z-fasta/actions/workflows/ci.yml"><img src="https://img.shields.io/badge/CI-passing-22c55e?style=for-the-badge" alt="CI" /></a>
  <a href="https://ziglang.org/download/0.14.0/"><img src="https://img.shields.io/badge/Zig-0.14.0-F7A41D?style=for-the-badge&logo=zig&logoColor=white" alt="Zig 0.14.0" /></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/license-MIT-6366f1?style=for-the-badge" alt="License: MIT" /></a>
  <a href="bench/"><img src="https://img.shields.io/badge/samtools-20×_faster-0ea5e9?style=for-the-badge" alt="20× faster than samtools" /></a>
</div>
<!-- markdownlint-enable MD041 -->

---

Quick links: [Installation](#installation) · [Usage](#usage) · [Performance & Correctness](#performance--correctness) · [Benchmarking](#benchmarking) · [Roadmap](#roadmap)

## Why z-fasta?

Modern bioinformatics workflows are bottlenecked by legacy text parsers. `z-fasta` bypasses standard I/O overhead by memory-mapping (`mmap`) the entire FASTA file, using explicit SIMD instructions in the indexer to scan for sequence headers at the theoretical limit of your NVMe drive.

- **Drop-in replacement:** Both `z-fasta index --emit-fai` and `z-fasta get` produce output byte-identical to `samtools faidx`. Falls back from `.zfi` to `.fai` with mtime + file-size staleness validation.
- **Single binary:** No dependencies, no `conda` environments, no `glibc` version errors.
- **Arena-allocated:** Uses Zig's `ArenaAllocator`. Zero memory leaks, minimal heap overhead in all modes.

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
z-fasta get <file.fasta> <region> [region ...]
```

Extract one or more sequences or sub-regions from an indexed FASTA file. Output is **byte-identical** to `samtools faidx`. Multiple regions are accepted in a single call; the index loads once and results stream in CLI order.

Requires an index: either `.zfi` (preferred) or `.fai`. If `.zfi` is not found, falls back to `.fai` automatically.

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

**Tier 2 (default):** full composition scan: nucleotide frequencies, GC content (N excluded), GC skew, soft-masked fraction. For proteins: top 3 amino acids with full names.

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

# Extract multiple regions in one call (index loads once)
z-fasta get genome.fa chr1:1000-2000 chr2:5000-6000 chrX:100-200

# Assembly stats (full composition scan)
z-fasta stats genome.fa

# Quick stats from index only (sub-millisecond)
z-fasta stats --index-only genome.fa
```

## Performance & Correctness

All timings on AMD Ryzen 9 3950X, warm cache.

### Index: SIMD-Accelerated Indexing

| Dataset       | Size   | z-fasta (no-dedup) | samtools | fastahack | pyfaidx | Speedup vs samtools |
| ------------- | ------ | ------------------ | -------- | --------- | ------- | ------------------- |
| Human Genome  | 3.0 GB | 0.45s              | 9.02s    | 21.79s    | 27.04s  | **20.2×**           |
| Transcriptome | 972 MB | 0.11s              | 1.79s    | 5.60s     | 6.36s   | **16.5×**           |
| Proteome      | 66 MB  | 0.006s             | 0.054s   | 0.253s    | 0.360s  | **8.6×**            |

| Mode         | Heap Memory | Notes                                                                         |
| ------------ | ----------- | ----------------------------------------------------------------------------- |
| `--no-dedup` | **< 1 MB**  | Fastest. `mmap` + SIMD, no deduplication hash map.                            |
| `default`    | ~45 MB      | `mmap` + SIMD, deduplicates sequence names.                                   |
| `--low-mem`  | **4 MB**    | `read()` + fixed 4 MB buffer, no `mmap`. For memory-constrained environments. |

> _`mmap` modes show VmRSS ≈ file size (OS-mapped pages); actual private heap is < 1 MB or ~45 MB as above._
See [bench/index/REPORT.md](bench/index/REPORT.md) for full scaling curves and memory analysis.

### Get: O(1) Region Extraction

| Dataset                | Region          | z-fasta     | samtools | seqtk  | pyfaidx | Speedup vs samtools |
| ---------------------- | --------------- | ----------- | -------- | ------ | ------- | ------------------- |
| Any (warm cache)       | 100 bp – 10 kbp | **~0.6 ms** | ~1.5 ms  | ~4 ms  | ~60 ms  | **2.2–2.9×**        |
| Proteome (14 MB)       | 1 kbp region    | 3.9 ms      | 11.1 ms  | 7.1 ms | 116 ms  | **2.8×**            |
| Transcriptome (459 MB) | 1 kbp region    | 126 ms      | 275 ms   | 220 ms | 1098 ms | **2.2×**            |

> Region extraction is O(1) regardless of file size. The index resolves a direct byte offset into the FASTA, then `z-fasta` streams bases from the mapped file while skipping line breaks. Note: fastahack is faster than z-fasta for large (>=50 MB) single full-sequence extraction due to a simpler write path; z-fasta leads on multi-sequence real datasets.

**Multi-region (v0.2.4):** `z-fasta get` accepts multiple regions per call, loading the index once and streaming all results in CLI order.

| Regions | z-fasta  | samtools  | seqtk    | Speedup vs samtools |
| ------- | -------- | --------- | -------- | ------------------- |
| 1       | 137 ms   | 280 ms    | 218 ms   | **2.0×**            |
| 10      | 140 ms   | 269 ms    | 221 ms   | **1.9×**            |
| 50      | 135 ms   | 272 ms    | 219 ms   | **2.0×**            |
| 100     | 138 ms   | 277 ms    | 218 ms   | **2.0×**            |

> Benchmarked on REAL_Transcriptome.fa (972 MB, 254,070 sequences). Latency is dominated by mmap startup cost and stays constant across region counts — each additional region costs ~1 µs of index lookup. seqtk performs a full-file scan per call regardless of region count and is listed for reference only.

See [bench/get/REPORT.md](bench/get/REPORT.md) for full results.

### Stats: Assembly/Proteome Statistics

| Mode       | Dataset              | z-fasta     | seqkit -a | seqtk comp | Speedup vs seqkit -a |
| ---------- | -------------------- | ----------- | --------- | ---------- | -------------------- |
| Index-only | Genome (3.0 GB)      | **0.8 ms**  | 17.4 s    | N/A        | **~21,000x**         |
| Index-only | Proteome (14 MB)     | **4.9 ms**  | 57.5 ms   | N/A        | **~12x**             |
| Full scan  | 1 GB single-seq file | 0.89 s      | 0.41 s    | N/A        | 0.46x                |
| Full scan  | Proteome (14 MB)     | **14.5 ms** | 57.5 ms   | 91.9 ms    | **~4×**              |

> Index-only time is constant (< 1 ms) regardless of file size. It reads only the binary `.zfi` header. Full-scan throughput is ~1.1 GB/s. seqkit -a is faster on large single-sequence synthetic files; z-fasta leads on multi-sequence files (proteomes, transcriptomes) and computes richer statistics (N50, GC composition, skew, amino acid breakdown) that seqkit does not provide.
See [bench/stats/REPORT.md](bench/stats/REPORT.md) for full results.

### Correctness

- **Index:** 20/20 edge cases match `samtools faidx` (exit codes and output).
- **Get:** 90/90 single-region and 22/22 multi-region byte-identical diff tests pass vs samtools across 5+ test files: full sequences, sub-regions, single bases, line-boundary spans, clamped ranges, duplicate regions, reversed CLI order, sort-path (≥16 regions).
- **Stats:** 107/107 BioPython verification tests pass: exact agreement on all Tier 1 and Tier 2 values across nucleotide and protein files.
- **Unit tests:** 85/85 Zig unit tests (23 index · 30 get · 32 stats).
- **Messy FASTA:** z-fasta is the only tool tested that correctly indexes mixed-width and trailing-whitespace FASTA files. samtools, fastahack, and pyfaidx all reject them. See [bench/index/REPORT.md](bench/index/REPORT.md) for the full compatibility matrix.

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
| `.zfi` | _(default)_  | Compact binary index. Fast to read/write programmatically. |
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

- [x] `z-fasta index`: SIMD-accelerated FASTA indexing (v0.1)
- [x] `z-fasta get`: O(1) byte-offset sequence extraction (v0.2)
- [x] `z-fasta stats`: Assembly/proteome statistics with index-only mode (v0.2)
- [x] Unified benchmark suite with per-module reports and figures (v0.2.2)
- [x] Expanded tool comparison: pyfaidx, seqtk added across all benchmark modules; messy FASTA compatibility matrix (v0.2.3)
- [x] Multi-region `get`: single call with N regions, index loads once, results stream in CLI order; ~2× faster than samtools across 1–100 regions (v0.2.4)

**Near-term**

- [ ] v0.2.5: BED file input
    - [ ] `--bed regions.bed` flag for batch extraction from BED files
    - [ ] BED coordinates are 0-based half-open; z-fasta converts to 1-based inclusive internally
    - [ ] Mix `--bed` with positional `NAME:START-END` args in one call
    - [ ] Enables direct comparison with `bedtools getfasta`
- [ ] v0.2.6: Reverse complement
    - [ ] `--rc` flag for `z-fasta get` to output the reverse complement of any extracted region
    - [ ] Comptime 256-element complement table, zero-cost lookup baked into the binary
    - [ ] Works with single regions, multi-region, and `--bed` batch calls
- [ ] v0.3.0: Validate + Tier 2 benchmarks + release polish
    - [ ] `z-fasta validate`: single-pass FASTA format checker with line-numbered error/warning output
    - [ ] Checks: duplicate names, inconsistent line widths, invalid characters, empty sequences, missing terminal newline
    - [ ] `--strict` flag treats warnings as errors
    - [ ] Tier 2 benchmark suite: noodles, rust-bio, Fusta, htslib, bedtools comparisons
    - [ ] Fix GET on messy FASTA (mixed-width and trailing-whitespace files indexed but not retrievable)

**Long-term / Exploratory**

- [ ] `z-fasta digest`: In-silico trypsin digestion for mass spectrometry (v0.4+)
- [ ] Zig version upgrade to 0.15+ for async I/O and improved SIMD support (v0.4+)
- [ ] Parallel mmap scanning for multi-threaded indexing on NVMe arrays
- [ ] Native BGZF / gzip streaming read support

## License

MIT. See [LICENSE](LICENSE)

---

<p align="center"><em>Aligned life in bytes,<br>
FASTA sings through mirrored streams.<br>
Humans bloom as code.</em></p>
