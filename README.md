<!-- markdownlint-disable MD033 MD036 MD041 -->
<div align="center">
  <h1>z-fasta ⚡</h1>
  <p>
    Fast, modular FASTA toolkit built in Zig.<br/>
    SIMD-accelerated indexing, O(1) region extraction, and instant assembly stats.<br/>
    samtools-compatible FASTA indexing and extraction, benchmarked against <code>seqkit</code>, <code>fastahack</code>, and <code>pyfaidx</code>.
  </p>
  <p>Current version: <strong>v0.2.6</strong></p>
  <br/>
  <a href="https://github.com/eneskemalergin/z-fasta/actions/workflows/ci.yml"><img src="https://img.shields.io/badge/CI-passing-22c55e?style=for-the-badge" alt="CI" /></a>
  <a href="https://ziglang.org/download/0.16.0/"><img src="https://img.shields.io/badge/Zig-0.16.0-F7A41D?style=for-the-badge&logo=zig&logoColor=white" alt="Zig 0.16.0" /></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/license-MIT-6366f1?style=for-the-badge" alt="License: MIT" /></a>
  <a href="bench/"><img src="https://img.shields.io/badge/indexing-22%C3%97_faster_than_samtools-0ea5e9?style=for-the-badge" alt="22x faster than samtools for indexing" /></a>
</div>
<!-- markdownlint-enable MD041 -->

---

Quick links: [Supported Today](#supported-today) · [Installation](#installation) · [Usage](#usage) · [Performance & Correctness](#performance--correctness) · [Benchmarking](#benchmarking) · [Roadmap](#roadmap)

## Supported Today

`z-fasta` is focused on uncompressed FASTA workflows: building indexes, extracting one or many indexed regions, and computing assembly/proteome statistics. It supports compact `.zfi` indexes and samtools-compatible `.fai` output. FASTQ is outside the current scope; BED input, reverse complement output, and compressed FASTA/BGZF streams remain roadmap items rather than current features.

## Why z-fasta?

Modern bioinformatics workflows are often bottlenecked by legacy text parsers. `z-fasta` keeps the hot paths close to the data: memory-mapped FASTA input for the default indexer, explicit SIMD header scanning, compact binary indexes, and startup-conscious CLI dispatch for tiny commands.

- **samtools-compatible output:** Both `z-fasta index --emit-fai` and `z-fasta get` produce output byte-identical to `samtools faidx` for the verified cases. Lookup falls back from `.zfi` to `.fai` with mtime + file-size staleness validation.
- **Single binary:** No dependencies, no `conda` environments, no `glibc` version errors.
- **Arena-scoped allocations:** Uses Zig's `ArenaAllocator` for short-lived command state, keeping heap overhead low and cleanup simple.

## Installation

```bash
# Download Zig 0.16.0 if you are not using the vendored toolchain
curl -L https://ziglang.org/download/0.16.0/zig-linux-x86_64-0.16.0.tar.xz | tar xJ

# Build with the vendored Zig 0.16.0 toolchain
./zig-0.16.0/zig build -Doptimize=ReleaseFast

# The executable is now at ./zig-out/bin/z-fasta
./zig-out/bin/z-fasta --help
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
  --index-only  Compute stats from index only (no FASTA scan; startup-dominated)
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

# Quick stats from index only (does not scan FASTA sequence bytes)
z-fasta stats --index-only genome.fa
```

## Performance & Correctness

All timings on AMD Ryzen 9 3950X, warm cache.

### Index: SIMD-Accelerated Indexing

| Dataset       | Size   | z-fasta (no-dedup) | samtools | fastahack | pyfaidx | Speedup vs samtools |
| ------------- | ------ | ------------------ | -------- | --------- | ------- | ------------------- |
| Human Genome  | 3.0 GB | 0.39s              | 9.03s    | 21.73s    | 27.48s  | **22.9×**           |
| Transcriptome | 972 MB | 0.093s             | 1.79s    | 5.72s     | 6.50s   | **19.3×**           |
| Proteome      | 66 MB  | 0.0056s            | 0.055s   | 0.275s    | 0.368s  | **10.0×**           |

| Mode         | Genome timing | Memory behavior                                                                  |
| ------------ | ------------- | -------------------------------------------------------------------------------- |
| `--no-dedup` | **0.39s**     | Fastest on repeated-name-free inputs. mmap-backed; MaxRSS reflects mapped pages. |
| `default`    | 0.40s         | Deduplicates names while staying in the same mmap-backed performance class.      |
| `--low-mem`  | 2.46s         | Streaming path; measured at 4.5 MB MaxRSS on the genome benchmark.               |

> _`mmap` modes show RSS close to the mapped FASTA size because `/usr/bin/time -v` counts mapped pages, not just private heap._
See [bench/index/REPORT.md](bench/index/REPORT.md) for full scaling curves and memory analysis.

### Get: O(1) Region Extraction

| Dataset                | Region          | z-fasta        | samtools   | seqtk     | pyfaidx | Speedup vs samtools |
| ---------------------- | --------------- | -------------- | ---------- | --------- | ------- | ------------------- |
| Any (warm cache)       | 100 bp – 10 kbp | **0.7–0.9 ms** | 1.5–1.6 ms | 4–34 ms   | ~60 ms  | **1.8–2.1×**        |
| Proteome (14 MB)       | 1 kbp region    | 1.3 ms         | 10.9 ms    | 7.2 ms    | 119 ms  | **8.4×**            |
| Transcriptome (972 MB) | 1 kbp region    | 25.3 ms        | 278.7 ms   | 220.3 ms  | 1103 ms | **11.0×**           |

> Small-region extraction is O(1), but on this host the end-to-end CLI path is startup-dominated below roughly 10 kbp. The Zig 0.16 minimal startup path keeps those calls under 1 ms on synthetic warm-cache fixtures. For very large full-sequence extraction, fastahack can still win on raw write-path overhead; z-fasta stays ahead of samtools across the real-dataset GET cases.

**Multi-region (v0.2.4):** `z-fasta get` accepts multiple regions per call, loading the index once and streaming all results in CLI order.

| Regions | z-fasta | samtools | seqtk  | Speedup vs samtools |
| ------- | ------- | -------- | ------ | ------------------- |
| 1       | 25.6 ms | 289 ms   | 221 ms | **11.3×**           |
| 10      | 33.8 ms | 283 ms   | 226 ms | **8.4×**            |
| 50      | 66.7 ms | 292 ms   | 225 ms | **4.4×**            |
| 100     | 66.7 ms | 279 ms   | 222 ms | **4.2×**            |

> Benchmarked on REAL_Transcriptome.fa (972 MB, 254,070 sequences). Latency is dominated by index resolution and output setup rather than region byte count. seqtk performs a full-file scan per call regardless of region count and is listed for reference only.

See [bench/get/REPORT.md](bench/get/REPORT.md) for full results.

### Stats: Assembly/Proteome Statistics

| Mode       | Dataset              | z-fasta     | seqkit -a | seqtk comp | Speedup vs seqkit -a |
| ---------- | -------------------- | ----------- | --------- | ---------- | -------------------- |
| Index-only | Genome (3.0 GB)      | **0.9 ms**  | 17.45 s   | N/A        | **~19,000×**         |
| Index-only | Proteome (14 MB)     | **2.9 ms**  | 57.8 ms   | N/A        | **~20×**             |
| Full scan  | 1 GB single-seq file | **0.78 s**  | 5.62 s    | 2.65 s     | **~7×**              |
| Full scan  | Proteome (14 MB)     | **11.8 ms** | 57.8 ms   | 93.0 ms    | **~4.9×**            |

> Index-only time is effectively constant with file size and is best described as startup-dominated. It reads `.zfi` index data and computes length-derived metrics without scanning FASTA sequence bytes. Full-scan throughput on synthetic files is ~1.3 GB/s, and the latest benchmark report has z-fasta ahead of seqkit on the real genome/proteome/transcriptome stats cases while still computing richer statistics.
See [bench/stats/REPORT.md](bench/stats/REPORT.md) for full results.

### Correctness

- **Index:** 20/20 edge cases match `samtools faidx` (exit codes and output).
- **Get:** 90/90 single-region and 22/22 multi-region byte-identical diff tests pass vs samtools across 5+ test files: full sequences, sub-regions, single bases, line-boundary spans, clamped ranges, duplicate regions, reversed CLI order, sort-path (≥16 regions).
- **Stats:** 107/107 BioPython verification tests pass: exact agreement on all Tier 1 and Tier 2 values across nucleotide and protein files.
- **Unit tests:** 86/86 Zig unit tests (24 index · 30 get · 32 stats).
- **Messy FASTA:** z-fasta is the only tool tested that correctly indexes mixed-width and trailing-whitespace FASTA files. samtools, fastahack, and pyfaidx all reject them. See [bench/index/REPORT.md](bench/index/REPORT.md) for the full compatibility matrix.

## Benchmarking

```bash
# Download real test data (~4 GB, one-time)
bash bench/shared/download_data.sh

# ── Index ─────────────────────────────────────────────────────────
./zig-0.16.0/zig build -Doptimize=ReleaseFast
bash bench/index/run_benchmarks.sh       # timing + memory
bash bench/index/run_tests.sh            # 20 edge-case correctness tests
.venv/bin/python bench/index/generate_report.py   # → bench/index/REPORT.md

# ── Get ───────────────────────────────────────────────────────────
bash bench/get/run_benchmarks.sh         # latency, scaling, real datasets
bash bench/get/verify_get.sh             # 90 byte-identical diff tests vs samtools
.venv/bin/python bench/get/generate_report.py     # → bench/get/REPORT.md

# ── Stats ─────────────────────────────────────────────────────────
bash bench/stats/run_benchmarks.sh       # full/index-only, scaling, throughput
.venv/bin/python bench/stats/verify_stats.py  # 107 BioPython verification tests
.venv/bin/python bench/stats/generate_report.py   # → bench/stats/REPORT.md
```

Full local refresh, in the same order used before publishing benchmark updates:

```bash
./zig-0.16.0/zig build -Doptimize=ReleaseFast && bash bench/index/run_tests.sh && bash bench/get/verify_get.sh && bash bench/get/verify_multi_get.sh && .venv/bin/python bench/stats/verify_stats.py && bash bench/index/run_benchmarks.sh --runs 5 && .venv/bin/python bench/index/generate_report.py && bash bench/get/run_benchmarks.sh --runs 5 && .venv/bin/python bench/get/generate_report.py && bash bench/stats/run_benchmarks.sh --runs 5 && .venv/bin/python bench/stats/generate_report.py && bash bench/perf-recovery/run_startup.sh
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
./zig-0.16.0/zig build

# Run all tests (index + get + stats)
./zig-0.16.0/zig build test --summary all

# Build optimized binary
./zig-0.16.0/zig build -Doptimize=ReleaseFast
```

## Roadmap

**Delivered**

- [x] `z-fasta index`: SIMD-accelerated FASTA indexing (v0.1)
- [x] `z-fasta get`: O(1) byte-offset sequence extraction (v0.2)
- [x] `z-fasta stats`: Assembly/proteome statistics with index-only mode (v0.2)
- [x] Unified benchmark suite with per-module reports and figures (v0.2.2)
- [x] Expanded tool comparison: pyfaidx, seqtk added across all benchmark modules; messy FASTA compatibility matrix (v0.2.3)
- [x] Multi-region `get`: single call with N regions, index loads once, results stream in CLI order; ~2× faster than samtools across 1–100 regions (v0.2.4)
- [x] Zig 0.16.0 migration plus benchmark/report refresh for v0.2.5
- [x] v0.2.6 performance recovery: lower startup overhead, faster index loading, buffered GET emission, fixed-width stats/index fast paths, and refreshed benchmark reports

**Near-term**

- [ ] v0.2.7: BED file input
    - [ ] `--bed regions.bed` flag for batch extraction from BED files
    - [ ] BED coordinates are 0-based half-open; z-fasta converts to 1-based inclusive internally
    - [ ] Mix `--bed` with positional `NAME:START-END` args in one call
    - [ ] Enables direct comparison with `bedtools getfasta`
- [ ] v0.2.8: Reverse complement
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
- [ ] Parallel mmap scanning for multi-threaded indexing on NVMe arrays
- [ ] Native BGZF / gzip streaming read support

## License

MIT. See [LICENSE](LICENSE)

---

<p align="center"><em>Aligned life in bytes,<br>
FASTA sings through mirrored streams.<br>
Humans bloom as code.</em></p>
