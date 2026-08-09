<!-- markdownlint-disable MD033 MD036 MD041 -->
<div align="center">
  <img src="assets/logo-readme.svg" alt="z-fasta" width="220">
  <p>
    Fast, modular FASTA toolkit built in Zig.<br/>
    SIMD-accelerated indexing, O(1) region extraction, validation, and assembly stats.<br/>
    samtools-compatible indexing and extraction, benchmarked against <code>seqkit</code>, <code>fastahack</code>, <code>pyfaidx</code>, and other peers.
  </p>
  <p>Current release: <strong>v0.3.1</strong></p>
  <br/>
  <a href="https://github.com/eneskemalergin/z-fasta/actions/workflows/ci.yml"><img src="https://img.shields.io/badge/CI-passing-22c55e?style=for-the-badge" alt="CI" /></a>
  <a href="https://ziglang.org/download/0.16.0/"><img src="https://img.shields.io/badge/Zig-0.16.0-F7A41D?style=for-the-badge&logo=zig&logoColor=white" alt="Zig 0.16.0" /></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/license-MIT-6366f1?style=for-the-badge" alt="License: MIT" /></a>
</div>

---

Quick links: [Supported Today](#supported-today) | [Installation](#installation) | [Usage](#usage) | [Performance & Correctness](#performance--correctness) | [Benchmarking](#benchmarking) | [Roadmap](#roadmap)

## Supported Today

`z-fasta` targets uncompressed FASTA: `index`, `get`, `stats`, and `validate`. Default indexes are compact `.zfi` files (embedded names, optional side tables, source identity). `--emit-fai` writes samtools-compatible `.fai` only when every record has fixed line geometry.

`get` accepts positional regions, BED files, BED from stdin, names files, strand-aware extraction, and orientation transforms (`--rc`, `--reverse-only`, `--complement-only`, `--annotate-rc`). FASTQ and compressed FASTA/BGZF remain out of scope.

CI builds and smokes the binary on Linux, macOS, and Windows. GET reads requested spans through the same portable file path on every platform. Commands that still map data use the portable mapping layer; POSIX memory advice remains optional.

## Why z-fasta?

Modern bioinformatics workflows are often bottlenecked by legacy text parsers. `z-fasta` keeps the hot paths close to the data: bounded buffered indexing, SIMD scanning, compact binary indexes, and a startup-conscious CLI for tiny commands.

- **samtools-compatible output:** `z-fasta index --emit-fai` (uniform FASTA only) and positional `z-fasta get` match `samtools faidx` on the verified cases. `--emit-fai` rejects non-uniform layouts and points you at default `.zfi` indexing.
- **Single binary:** no runtime dependencies, no `conda` environments, no `glibc` version traps.
- **Arena-scoped allocations:** short-lived command state uses Zig arenas so heap overhead stays low and cleanup stays simple.

## Installation

```bash
# Download Zig 0.16.0 if you are not using the vendored toolchain
curl -L https://ziglang.org/download/0.16.0/zig-linux-x86_64-0.16.0.tar.xz | tar xJ

# Build with the repo-local Zig wrapper (uses ./zig-0.16.0/zig)
./zig build -Doptimize=ReleaseFast

# The executable is now at ./zig-out/bin/z-fasta
./zig-out/bin/z-fasta --help
./zig-out/bin/z-fasta --version
```

## Usage

### Index

```bash
z-fasta index [options] <file.fasta>

Options:
  --emit-fai    Output FAI to stdout when every record has fixed line geometry;
                otherwise fails and directs you to default .zfi indexing
  --no-dedup    Keep duplicate sequence names in the index (default: first wins at
                index time). get resolves duplicate names to the first record.
  --help        Show help message
  --version     Print version
```

`index` reads the FASTA through one bounded buffer and writes `{file}.zfi`. Sequence names longer than 65535 bytes are rejected (`HeaderTooLong`).

### Get (sequence extraction)

```bash
z-fasta get <file.fasta> [options] <region> [region ...]
z-fasta get <file.fasta> [options] --names file.txt|-
z-fasta get <file.fasta> [options] --bed file.bed|- [--strand-aware]
```

Extract one or more sequences or sub-regions from an indexed FASTA. Choose exactly one request source: positional regions, `--names`, or `--bed`. Positional output is **byte-identical** to `samtools faidx` on the verified path. BED and names flows are checked against `bedtools getfasta` and `samtools faidx -r`. Every source preserves its request order.

Requires an index. A present `.zfi` is authoritative, including when it is invalid. `.fai` is used only when `.zfi` is absent. Messy (non-uniform) records need a side-table `.zfi`; a plain `.fai` cannot describe them.

**Region formats:** `NAME` (full sequence); `NAME:START-END` (1-based inclusive); `NAME:START-` (from START through end of sequence). Ensembl-style names with colons work (for example `chromosome:GRCh38:1:1:248956422:1`).

**GET flags:**

- **`--bed file.bed`** / **`--bed -`**: BED regions from a file or stdin. Coordinates are 0-based half-open; z-fasta converts to 1-based inclusive internally.
- **`--names file.txt`** / **`--names -`**: one literal full-sequence name per line from a file or stdin. Blank and comment lines are skipped.
- **`--strand-aware`** (alias **`--honor-strand`**): BED only. Read column 6; `-` applies reverse-complement before any global orientation flag.
- **`--rc`**: reverse-complement output. Verified against `samtools faidx -i --mark-strand no`. Mutually exclusive with **`--complement-only`** and **`--reverse-only`**.
- **`--complement-only`**: complement without reverse. Nucleotide records only.
- **`--reverse-only`**: reverse without complement.
- **`--annotate-rc`**: append the final composed orientation to transformed headers. It requires a global transform or strand-aware BED input. Identity output has no suffix.
- **`--summary`**: print region count, total bases, elapsed time, and regions/sec to stderr after successful output. Timing includes index loading and request acquisition.

Positional CLI regions are capped at **1024** per invocation (`error: too many regions`). Names and BED input have no row-count or whole-input size cap. Both stream through reusable request storage instead of retaining the complete input.

Complement-based transforms error on protein FASTA so **`--rc`** and **`--complement-only`** stay nucleotide-only. That guard samples up to the first **256** bases of each requested record.

### Validate

```bash
z-fasta validate [options] <file.fasta>

Options:
  --strict                 Treat warnings as errors
  --json                   Emit JSON Lines instead of text
  --summary                With --json, emit one summary object
  --fix -o <file.fasta>    Write a fixed FASTA to a separate output path
  --fix-format-only        With --fix, proceed despite alphabet errors
  --schema NAME            Header schema: uniprot or refseq
  --custom-alphabet CHARS  Override nucleotide/protein alphabet checks
  --max-header-len N       Warn on headers longer than N bytes (default: 1024)
```

Checks structure, alphabets, and headers (duplicate names, invalid characters, null bytes, UTF-8 BOM, inconsistent line widths, trailing whitespace, empty sequences, missing terminal newlines, mixed line endings, long headers, schema violations). The event list stops at 10000 issues with a deterministic error. `--fix` streams rewritten FASTA to `-o` when the fix is safe.

### Stats

```bash
z-fasta stats <file.fasta>
```

Reports indexed record count, total symbols, shortest and longest records and names, mean, quartiles, median, range, N50/L50, N90/L90, auN, and complete composition. Nucleotide output distinguishes T, U, mixed T/U, ambiguity, invalid symbols, GC, GC skew, and lowercase input. Protein output reports every standard and ambiguous residue, stop and invalid symbols, and lowercase input. Percentages use total symbols as their denominator.

### Sequence type classification

All commands share `stats.detectType` (IUPAC nucleotide alphabet; nucleotide if those letters are strictly more than 90% of counted bases). Sampling differs on purpose:

- `stats`: full-file composition scan
- `get --rc` / `--complement-only`: up to 256 bases per record
- `validate`: up to the first 100000 sequence bases; `--json --summary` reports `sequence_type`, `type_bases_sampled`, and `type_sample_cap`

### Examples

```bash
# Create .zfi binary index (default)
z-fasta index genome.fa

# Output .fai to stdout (uniform FASTA only)
z-fasta index --emit-fai genome.fa > genome.fai

# Keep duplicate names in the index (get resolves to the first record)
z-fasta index --no-dedup genome.fa

# Extract a full sequence / sub-region / multi-region
z-fasta get genome.fa chr1
z-fasta get genome.fa chr1:1000000-2000000
z-fasta get genome.fa chr1:1000-2000 chr2:5000-6000

# Orientation transforms
z-fasta get genome.fa chr1:1000-2000 --rc
z-fasta get genome.fa chr1:1000-2000 --reverse-only
z-fasta get genome.fa chr1:1000-2000 --complement-only
z-fasta get genome.fa chr1:1000-2000 --rc --annotate-rc

# BED / names / strand
z-fasta get genome.fa --bed regions.bed
awk '$5 > 100' raw.bed | z-fasta get genome.fa --bed -
z-fasta get genome.fa --names ids.txt
z-fasta get genome.fa --bed regions.bed --strand-aware --summary

# Validate and optionally rewrite a clean copy
z-fasta validate genome.fa
z-fasta validate --json --summary genome.fa
z-fasta validate --fix -o genome.clean.fa genome.fa

# Assembly stats
z-fasta stats genome.fa
```

## Formats and behavior

### Index formats

- **`.zfi` (default, authoritative when present):** binary index with embedded names, optional per-line side tables for messy layout, and source identity (size plus embedded FASTA mtime).
- **`.fai` (compatibility):** text index for uniform, FAI-representable records only. `--emit-fai` refuses variable widths and other non-representable layouts instead of writing a misleading file. Stale `.fai` checks are mtime-only (weaker than `.zfi`).

### Messy FASTA

Variable line widths, trailing spaces or tabs, blank lines, mixed CRLF/LF, and a final line without a terminal newline are indexed in `.zfi` via side tables. `get` uses those tables for correct extraction. Uniform records keep O(1) byte-offset math. Peers that require classic FAI geometry often reject these files; that is expected. `validate --fix` can rewrite a clean uniform copy when you need `.fai` or stricter tooling.

### Duplicate names

Duplicate _names_ are not the same as identical _sequence contents_. Default `index` keeps the first name; `index --no-dedup` keeps all. `get` resolves a repeated name to the **first** exact match. Duplicate lookup is ambiguous, so use unique identifiers when records must be addressable individually. `validate` reports `duplicate_name`. `stats` prints source-level extras (`sum(k-1)`).

## Support and limits

**Formats:** `.zfi` is the default. It handles uniform records and messy layout (variable widths, trailing whitespace, blank lines, mixed CRLF/LF, missing final newline) via side tables. `.fai` is compatibility only for uniform, FAI-representable records; `--emit-fai` refuses messy layout. A present `.zfi` blocks `.fai` selection when it is stale or invalid. Compressed FASTA, BGZF, and FASTQ are out of scope.

**Names:** index hard-rejects sequence names longer than 65535 bytes. Empty identifiers are valid but difficult to address and should be avoided. An explicit empty positional argument retrieves one; blank names-file lines stay skipped and BED requires a non-empty chromosome. `validate --max-header-len` warns above N bytes (default 1024) and does not raise that index limit.

**Get / validate caps:** at most 1024 positional regions per `get` call. Names and BED request names may be at most 65535 bytes. Streaming request storage retains up to 4 MiB of unique name bytes plus one final name per batch; names batches hold at most 65536 requests and BED batches at most 4096. `validate` stops after 10000 retained events.

**Memory:** index reads sequence payload through a bounded buffer. GET reads requested FASTA spans through one descriptor and retains only index metadata, active requests, side tables, and fixed buffers. Stats and validate still map FASTA data, so full scans can approach the mapped file size.

**Platforms:** Linux, macOS, and Windows share the same commands and GET file-access path. Memory mapping and advice used by other commands stay behind the portable platform layer. CI builds and smokes six native platform lanes: x86_64 and arm64 on all three operating systems. Tagged releases publish the same six archives; v0.3.0 was the first complete six-archive release.

## Performance & Correctness

All timings below are on AMD Ryzen 9 3950X with warm cache. See the linked benchmark reports for the full methods and results.

### Index: SIMD-Accelerated Indexing

| Dataset       | Size on disk | z-fasta (.fai) | noodles  | rust-bio | samtools | vs noodles | vs samtools |
| ------------- | ------------ | -------------- | -------- | -------- | -------- | ---------- | ----------- |
| Human Genome  | ~2.9 GiB     | 0.3714 s       | 1.3337 s | 3.6716 s | 9.0976 s | **3.6x**   | **24.5x**   |
| Transcriptome | ~459 MiB     | 0.2124 s       | 0.3658 s | 0.6995 s | 1.8315 s | **1.72x**  | **8.6x**    |
| Proteome      | ~13 MiB      | 0.0105 s       | 0.0174 s | 0.0255 s | 0.0598 s | **1.67x**  | **5.7x**    |

See the [detailed index benchmark report](bench/index/REPORT.md) for tool definitions, methodology, memory, deduplication, and scaling results.

### Get: O(1) Region Extraction

| Dataset                  | Region       | z-fasta (.zfi) | z-fasta (.fai) | noodles | rust-bio | samtools | vs noodles | vs samtools |
| ------------------------ | ------------ | -------------- | -------------- | ------- | -------- | -------- | ---------- | ----------- |
| Genome (~2.9 GiB)        | 1 kbp region | 2.0 ms         | 2.2 ms         | 2.6 ms  | 2.7 ms   | 3.3 ms   | **1.28x**  | **1.62x**   |
| Transcriptome (~459 MiB) | 1 kbp region | 4.8 ms         | 27.1 ms        | 88.1 ms | 530.5 ms | 292.0 ms | **18.5x**  | **61.2x**   |
| Proteome (~13 MiB)       | 1 kbp region | 2.4 ms         | 3.6 ms         | 6.8 ms  | 20.0 ms  | 12.7 ms  | **2.9x**   | **5.4x**    |

> Small-region extraction is O(1), but on this host the end-to-end CLI path is startup-dominated below roughly 10 kbp. For very large full-sequence extraction, fastahack can still win on raw write-path overhead; z-fasta stays ahead of samtools across the real-dataset GET cases.

**`--rc`** applies reverse traversal and complement lookup during output without materializing a second copy of the region.

**Multi-region:** one call loads the index once and streams results in request order. Transcriptome, 1 kbp per region:

| Regions | z-fasta (.zfi) | z-fasta (.fai) | noodles  | samtools | vs noodles | vs samtools |
| ------- | -------------- | -------------- | -------- | -------- | ---------- | ----------- |
| 1       | 4.9 ms         | 26.9 ms        | 86.5 ms  | 283.1 ms | **17.7x**  | **57.9x**   |
| 10      | 24.5 ms        | 30.6 ms        | 89.4 ms  | 283.9 ms | **3.6x**   | **11.6x**   |
| 100     | 24.2 ms        | 32.2 ms        | 124.5 ms | 290.9 ms | **5.1x**   | **12.0x**   |
| 1,000   | 30.0 ms        | 37.1 ms        | 407.7 ms | 295.8 ms | **13.6x**  | **9.9x**    |

> Benchmarked on REAL_Transcriptome.fa. Latency is dominated by index resolution and output setup rather than region byte count.

**BED:** one invocation streams 1 kbp intervals from a BED file. Transcriptome:

| BED rows | z-fasta (.zfi) | z-fasta (.fai) | noodles  | bedtools | vs noodles | vs bedtools |
| -------- | -------------- | -------------- | -------- | -------- | ---------- | ----------- |
| 10       | 26.0 ms        | 57.7 ms        | 91.8 ms  | 585.6 ms | **3.5x**   | **22.5x**   |
| 100      | 25.5 ms        | 55.9 ms        | 130.6 ms | 582.4 ms | **5.1x**   | **22.9x**   |
| 1,000    | 27.4 ms        | 59.9 ms        | 469.4 ms | 575.9 ms | **17.1x**  | **21.0x**   |
| 10,000   | 43.6 ms        | 75.4 ms        | 3.700 s  | 670.3 ms | **84.8x**  | **15.4x**   |

See the [detailed GET benchmark report](bench/get/REPORT.md) for tool definitions, methodology, memory, page faults, BED, reverse-complement, and messy-layout results.

### Stats: Assembly/Proteome Statistics

| Dataset                  | z-fasta (.zfi) | z-fasta (.fai) | noodles | SeqKit `stats -a` | vs noodles | vs SeqKit |
| ------------------------ | -------------- | -------------- | ------- | ----------------- | ---------- | --------- |
| Genome (~2.9 GiB)        | **2.752 s**    | **2.757 s**    | 6.119 s | 17.777 s          | **2.22x**  | **6.46x** |
| Transcriptome (~459 MiB) | **0.382 s**    | **0.416 s**    | 1.107 s | 2.403 s           | **2.90x**  | **6.29x** |
| Proteome (~13 MiB)       | **12.3 ms**    | **15.5 ms**    | 37.0 ms | 59.6 ms           | **3.00x**  | **4.83x** |

> Noodles is the project benchmark wrapper that computes the same agreed complete field set. SeqKit reports a smaller assembly and composition field set, so its column is a familiar ecosystem latency reference rather than an equivalent-work comparison. Rust-bio and Seqtk remain in the detailed report but are omitted here to keep the summary readable.

See the [detailed stats benchmark report](bench/stats/REPORT.md) for field coverage, seven-sample variation, RSS, page faults, correctness, and `.zfi`/`.fai` scaling results.

### Correctness

- **Index:** `bench/index/run.sh` edge cases **24/24 matching cases plus one `binary_data` review**; messy layouts from `python3 bench/shared/generate_messy.py` into `bench/shared/cache/messy_fixtures/` (correctness) and `messy_perf/` (proteome perf).
- **Get:** `bench/get/run.sh` correctness **418/418** (positional, multi-region, BED, names, RC, messy) against samtools, bedtools, and seqtk where applicable.
- **Stats:** `bench/stats/run.sh` checks an independent exact oracle, `.zfi`/`.fai` parity, layout twins, messy fixtures, deduplication, CLI errors, and complete noodles/rust-bio peer parity.
- **Unit tests:** `./zig build test` (index, get, stats, complement, BED parser, validator).
- **Messy FASTA:** z-fasta indexes and extracts mixed-width and trailing-whitespace FASTA that samtools-style FAI tools reject. Details: [bench/index/REPORT.md](bench/index/REPORT.md).

## Benchmarking

```bash
# Download real test data (~4 GB, one-time; checksummed via datasets.manifest)
bash bench/shared/download_data.sh

./zig build -Doptimize=ReleaseFast

# Index: correctness, zebrac perf, messy zebrac, report
bash bench/index/run.sh
# Correctness only:
bash bench/index/run.sh --skip-benchmarks --skip-messy --skip-report

# GET: correctness, then optional perf + report
bash bench/get/run.sh
# Correctness only:
bash bench/get/run.sh --skip-benchmarks --skip-report

# Stats: correctness, then optional perf + report
bash bench/stats/run.sh
# Correctness only:
bash bench/stats/run.sh --skip-benchmarks --skip-report
```

Details: [bench/README.md](bench/README.md). Shared helpers live in `bench/shared/`.

## Development

```bash
./zig build
./zig build test --summary all
./zig build test -Doptimize=ReleaseFast --summary all
./zig build -Doptimize=ReleaseFast
./zig-out/bin/z-fasta --version
```

## Roadmap

**Implemented in current development**

- [x] `index`, `get`, `stats`, and `validate`
- [x] Preferred `.zfi` with side tables; `.fai` only for representable uniform records
- [x] Messy FASTA GET via side tables; `validate --fix` for safe rewrites
- [x] Bounded index reader shared by `.zfi` and `.fai` output
- [x] Portable Linux / macOS / Windows mapping and CI smoke
- [x] Zebrac benchmark suites for index, GET, and stats with verify gates
- [x] Stripped ReleaseFast artifacts, ship-mode tests, and six-platform CI smoke

**Deferred**

- [ ] Compressed or BGZF FASTA
- [ ] Parallel indexing
- [ ] `z-fasta digest` and other domain-specific analyses
- [ ] Broader paper-style baseline collection

## License

MIT. See [LICENSE](LICENSE)

---

<p align="center"><em>Aligned life in bytes,<br>
FASTA sings through mirrored streams.<br>
Humans bloom as code.</em></p>
