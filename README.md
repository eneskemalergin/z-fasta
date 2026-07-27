<!-- markdownlint-disable MD033 MD036 MD041 -->
<div align="center">
  <h1>z-fasta ⚡</h1>
  <p>
    Fast, modular FASTA toolkit built in Zig.<br/>
    SIMD-accelerated indexing, O(1) region extraction, and instant assembly stats.<br/>
    samtools-compatible FASTA indexing and extraction, benchmarked against <code>seqkit</code>, <code>fastahack</code>, and <code>pyfaidx</code>.
  </p>
  <p>Current release: <strong>v0.2.9</strong></p>
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

`z-fasta` is focused on uncompressed FASTA workflows: building indexes, extracting one or many indexed regions, and computing assembly/proteome statistics. It supports compact `.zfi` indexes and samtools-compatible `.fai` output. `get` accepts positional regions, BED files, BED from stdin, names files, strand-aware extraction, and explicit orientation transforms through `--rc`, `--reverse-only`, `--complement-only`, and `--annotate-rc`. FASTQ and compressed FASTA/BGZF streams remain outside the current scope.

## Why z-fasta?

Modern bioinformatics workflows are often bottlenecked by legacy text parsers. `z-fasta` keeps the hot paths close to the data: memory-mapped FASTA input for the default indexer, explicit SIMD header scanning, compact binary indexes, and startup-conscious CLI dispatch for tiny commands.

- **samtools-compatible output:** Both `z-fasta index --emit-fai` (uniform FASTA only) and `z-fasta get` produce output byte-identical to `samtools faidx` for the verified cases. `--emit-fai` rejects non-uniform layouts and directs you to default `.zfi` indexing. Lookup prefers `.zfi` (size + embedded source mtime) and falls back to `.fai` (mtime-only identity).
- **Single binary:** No dependencies, no `conda` environments, no `glibc` version errors.
- **Arena-scoped allocations:** Uses Zig's `ArenaAllocator` for short-lived command state, keeping heap overhead low and cleanup simple.

## Installation

```bash
# Download Zig 0.16.0 if you are not using the vendored toolchain
curl -L https://ziglang.org/download/0.16.0/zig-linux-x86_64-0.16.0.tar.xz | tar xJ

# Build with the repo-local Zig wrapper (uses ./zig-0.16.0/zig)
./zig build -Doptimize=ReleaseFast

# The executable is now at ./zig-out/bin/z-fasta
./zig-out/bin/z-fasta --help
```

## Usage

### Index

```bash
z-fasta index [options] <file.fasta>

Options:
  --emit-fai    Output FAI to stdout when every record has fixed line geometry;
                otherwise fails and directs you to default .zfi indexing
  --no-dedup    Keep duplicate sequence names in the index (default: first wins at
                index time). get resolves duplicate names to the last record.
  --low-mem     Stream input with bounded RAM; same outputs as default index
  --help        Show help message
  --version     Print version
```

### Get (sequence extraction)

```bash
z-fasta get <file.fasta> [--bed file.bed|-] [--names file.txt]
            [--strand-aware] [--summary]
            [--rc|--complement-only|--reverse-only] [--annotate-rc]
            [--chunk-size N|-1] <region> [region ...]
```

Extract one or more sequences or sub-regions from an indexed FASTA file. Output is **byte-identical** to `samtools faidx` for the positional-region path, and the BED / names batch flows are verified against `bedtools getfasta` and `samtools faidx -r`. Multiple regions are accepted in a single call; the index loads once and results stream in CLI order. BED rows and names-file entries are appended in source order ahead of later positional arguments.

Requires an index: either `.zfi` (preferred) or `.fai`. If `.zfi` is not found, falls back to `.fai` automatically.

**Region formats:** `NAME` (full sequence); `NAME:START-END` (1-based inclusive sub-region); `NAME:START-` (from START through end of sequence). Ensembl-style names with colons work (for example `chromosome:GRCh38:1:1:248956422:1`).

**GET flags:**

- **`--bed file.bed`** / **`--bed -`**: BED regions from a file or stdin. Coordinates are 0-based half-open; z-fasta converts to 1-based inclusive internally.
- **`--names file.txt`**: one full-sequence name per line for long batch lists.
- **`--strand-aware`** (alias **`--honor-strand`**): read BED column 6; `-` applies reverse-complement before any global orientation flag.
- **`--rc`**: reverse-complement output. Verified against `samtools faidx -i --mark-strand no`. Mutually exclusive with **`--complement-only`** and **`--reverse-only`**.
- **`--complement-only`**: complement without reverse. Nucleotide records only.
- **`--reverse-only`**: reverse without complement.
- **`--annotate-rc`**: append a transform suffix to headers (for example `(reverse complement)`). Default output stays samtools-style.
- **`--summary`**: print region count, total bases, elapsed time, and regions/sec to stderr.
- **`--chunk-size N`**: BED batch size (default `4096`). Use `1` only for debugging. **`--chunk-size -1`**: load all BED rows in one batch when memory allows (names/BED inputs over 512 MiB are rejected with `-1`).

Positional CLI regions are capped at **1024** per invocation (`error: too many regions`). BED and **`--names`** lists are not subject to that cap.

Complement-based transforms error on protein FASTA input so **`--rc`** and **`--complement-only`** stay nucleotide-only.

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

# Reverse-complement a region
z-fasta get genome.fa chr1:1000-2000 --rc

# Reverse without complementing
z-fasta get genome.fa chr1:1000-2000 --reverse-only

# Complement without reversing
z-fasta get genome.fa chr1:1000-2000 --complement-only

# Add explicit transform text to the FASTA header
z-fasta get genome.fa chr1:1000-2000 --rc --annotate-rc

# Extract regions from BED
z-fasta get genome.fa --bed regions.bed

# Read BED from stdin
awk '$5 > 100' raw.bed | z-fasta get genome.fa --bed -

# Extract whole sequences from a names file
z-fasta get genome.fa --names ids.txt

# Respect BED strand and print a stderr summary
z-fasta get genome.fa --bed regions.bed --strand-aware --summary

# Compose BED strand handling with a global reverse-complement flip
z-fasta get genome.fa --bed regions.bed --honor-strand --rc

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
| Human Genome  | 3.0 GB | 0.39s              | 9.03s    | 21.73s    | 27.48s  | **22.9x**           |
| Transcriptome | 972 MB | 0.093s             | 1.79s    | 5.72s     | 6.50s   | **19.3x**           |
| Proteome      | 66 MB  | 0.0056s            | 0.055s   | 0.275s    | 0.368s  | **10.0x**           |

**Index modes** on Genome (warm cache; [bench/index/REPORT.md](bench/index/REPORT.md) run `20260706_134943`): **default** mmap ~0.40s, RSS ≈ mapped FASTA size; **`--low-mem`** stream ~1.61s FAI / ~3.4 MB RSS, same `.zfi` bytes as default; **`--no-dedup`** ~0.38s; **`--emit-fai`** writes FAI to stdout only (no on-disk `.zfi`).

Both default and **`--low-mem`** write the same `.zfi` unless **`--emit-fai`** is set. The gap is how the FASTA is read during the build, not the on-disk format.

> _`mmap` modes show RSS close to the mapped FASTA size because `/usr/bin/time -v` counts mapped pages, not just private heap._
> See [bench/index/REPORT.md](bench/index/REPORT.md) for full scaling curves and memory analysis.

### Get: O(1) Region Extraction

| Dataset                | Region          | z-fasta        | samtools   | seqtk    | pyfaidx | Speedup vs samtools |
| ---------------------- | --------------- | -------------- | ---------- | -------- | ------- | ------------------- |
| Any (warm cache)       | 100 bp - 10 kbp | **0.7-0.9 ms** | 1.5-1.6 ms | 4-34 ms  | ~60 ms  | **1.8-2.1x**        |
| Proteome (14 MB)       | 1 kbp region    | 1.3 ms         | 10.9 ms    | 7.2 ms   | 119 ms  | **8.4x**            |
| Transcriptome (972 MB) | 1 kbp region    | 25.3 ms        | 278.7 ms   | 220.3 ms | 1103 ms | **11.0x**           |

> Small-region extraction is O(1), but on this host the end-to-end CLI path is startup-dominated below roughly 10 kbp. The historical checked-in benchmark report for v0.2.6 was generated under a faster local benchmark environment than the current reruns; direct side-by-side rebuilds of v0.2.6, v0.2.7, and current `main` on the same machine do not reproduce a material no-flag `get` regression. For very large full-sequence extraction, fastahack can still win on raw write-path overhead; z-fasta stays ahead of samtools across the real-dataset GET cases.

Orientation note: **`--rc`** uses the same mmap-backed extraction path and applies reverse traversal plus complement lookup during emission instead of materializing a second copy of the region.

**Multi-region (v0.2.4):** `z-fasta get` accepts multiple regions per call, loading the index once and streaming all results in CLI order.

| Regions | z-fasta | samtools | seqtk  | Speedup vs samtools |
| ------- | ------- | -------- | ------ | ------------------- |
| 1       | 25.6 ms | 289 ms   | 221 ms | **11.3x**           |
| 10      | 33.8 ms | 283 ms   | 226 ms | **8.4x**            |
| 50      | 66.7 ms | 292 ms   | 225 ms | **4.4x**            |
| 100     | 66.7 ms | 279 ms   | 222 ms | **4.2x**            |

> Benchmarked on REAL_Transcriptome.fa (972 MB, 254,070 sequences). Latency is dominated by index resolution and output setup rather than region byte count. seqtk performs a full-file scan per call regardless of region count and is listed for reference only.
>
> GET and stats benchmark reports (`bench/get/REPORT.md`, `bench/stats/REPORT.md`) are from pre-cleanup runs; perf harnesses are pending rebuild for v0.3.0. Index numbers above are from [bench/index/REPORT.md](bench/index/REPORT.md) run `20260706_134943`.

### Stats: Assembly/Proteome Statistics

| Mode       | Dataset              | z-fasta     | seqkit -a | seqtk comp | Speedup vs seqkit -a |
| ---------- | -------------------- | ----------- | --------- | ---------- | -------------------- |
| Index-only | Genome (3.0 GB)      | **0.9 ms**  | 17.45 s   | N/A        | **~19,000x**         |
| Index-only | Proteome (14 MB)     | **2.9 ms**  | 57.8 ms   | N/A        | **~20x**             |
| Full scan  | 1 GB single-seq file | **0.78 s**  | 5.62 s    | 2.65 s     | **~7x**              |
| Full scan  | Proteome (14 MB)     | **11.8 ms** | 57.8 ms   | 93.0 ms    | **~4.9x**            |

> Historical stats benchmarks; harness and `bench/stats/REPORT.md` pending rebuild. Index-only mode reads `.zfi` without scanning sequence bytes.

### Correctness

- **Index:** edge-case and messy-variant correctness via `bench/index/run.sh` (`edge_cases/` generated on the fly; `messy_variants/` checked in).
- **Get:** verification in `bench/get/verify.sh` (single/multi-region, BED, names, RC, messy FASTA) against samtools, bedtools, and seqtk where applicable.
- **Stats:** verification harness pending rebuild (removed in v0.3.0 bench cleanup).
- **Unit tests:** `./zig build test` (index, get, stats, complement, BED parser, validator).
- **Messy FASTA:** z-fasta indexes mixed-width and trailing-whitespace FASTA files that samtools, fastahack, and pyfaidx reject. See [bench/index/REPORT.md](bench/index/REPORT.md) for the compatibility matrix.

## Benchmarking

```bash
# Download real test data (~4 GB, one-time)
bash bench/shared/download_data.sh

# Index (full suite: correctness, zebrac perf, messy zebrac, report)
./zig build -Doptimize=ReleaseFast
bash bench/index/run.sh
bash bench/index/run.sh --skip-report         # benchmarks only; report separately

# Get verification (bench perf report pending rebuild)
bash bench/get/verify.sh

# Regenerate index report after a benchmark run
.venv/bin/python bench/index/generate_report.py   # -> bench/index/REPORT.md
```

See [bench/index/README.md](bench/index/README.md) for index runner flags. GET/stats zebrac suites and baseline snapshots (`bench/save_baseline.py`) were removed in the July 2026 bench cleanup and are planned to return in v0.3.0.

## Output Formats

- **`.zfi`** (default): compact binary index for fast programmatic read/write.
- **`.fai`** (`--emit-fai`): tab-separated text for uniform, FAI-representable records only; byte-identical to `samtools faidx` for those cases. Non-uniform layouts require default `.zfi` indexing.

## Development

```bash
# Build (debug)
./zig build

# Run all tests (index + get + stats)
./zig build test --summary all

# Build optimized binary
./zig build -Doptimize=ReleaseFast
```

## Roadmap

**Delivered**

- [x] `z-fasta index`: SIMD-accelerated FASTA indexing (v0.1)
- [x] `z-fasta get`: O(1) byte-offset sequence extraction (v0.2)
- [x] `z-fasta stats`: Assembly/proteome statistics with index-only mode (v0.2)
- [x] Unified benchmark suite with per-module reports and figures (v0.2.2)
- [x] Expanded tool comparison: pyfaidx, seqtk added across all benchmark modules; messy FASTA compatibility matrix (v0.2.3)
- [x] Multi-region `get`: single call with N regions, index loads once, results stream in CLI order; ~2x faster than samtools across 1-100 regions (v0.2.4)
- [x] Zig 0.16.0 migration plus benchmark/report refresh for v0.2.5
- [x] v0.2.6 performance recovery: lower startup overhead, faster index loading, buffered GET emission, fixed-width stats/index fast paths, and refreshed benchmark reports
- [x] v0.2.7 BED batch extraction: `--bed`, `--bed -`, `--names`, `--strand-aware`, bounded chunked processing, and verification/benchmark coverage
- [x] v0.2.8 reverse/complement extraction: `--rc`, `--reverse-only`, `--complement-only`, `--annotate-rc`, RC verification, and integrated RC benchmark/report coverage
- [x] v0.2.9 overall quality improvements around memory safety, performance optimization, and code cleanup.

**Near-term**

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
