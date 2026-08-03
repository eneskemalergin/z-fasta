<!-- markdownlint-disable MD033 MD036 MD041 -->
<div align="center">
  <h1>z-fasta</h1>
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
<!-- markdownlint-enable MD041 -->

---

Quick links: [Supported Today](#supported-today) | [Installation](#installation) | [Usage](#usage) | [Performance & Correctness](#performance--correctness) | [Benchmarking](#benchmarking) | [Roadmap](#roadmap)

## Supported Today

`z-fasta` targets uncompressed FASTA: `index`, `get`, `stats`, and `validate`. Default indexes are compact `.zfi` files (embedded names, optional side tables, source identity). `--emit-fai` writes samtools-compatible `.fai` only when every record has fixed line geometry.

`get` accepts positional regions, BED files, BED from stdin, names files, strand-aware extraction, and orientation transforms (`--rc`, `--reverse-only`, `--complement-only`, `--annotate-rc`). FASTQ and compressed FASTA/BGZF remain out of scope.

CI builds and smokes the binary on Linux, macOS, and Windows. Mapping-backed commands go through a portable `MemoryMap` layer; POSIX memory advice is an optimization and a no-op on Windows.

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
                index time). get resolves duplicate names to the last record.
  --help        Show help message
  --version     Print version
```

`index` reads the FASTA through one bounded buffer and writes `{file}.zfi`. Sequence names longer than 65535 bytes are rejected (`HeaderTooLong`).

### Get (sequence extraction)

```bash
z-fasta get <file.fasta> [--bed file.bed|-] [--names file.txt]
            [--strand-aware] [--summary]
            [--rc|--complement-only|--reverse-only] [--annotate-rc]
            [--chunk-size N|-1] <region> [region ...]
```

Extract one or more sequences or sub-regions from an indexed FASTA. Positional output is **byte-identical** to `samtools faidx` on the verified path. BED and names flows are checked against `bedtools getfasta` and `samtools faidx -r`. Multiple regions load the index once and stream in CLI order. BED rows and names-file entries are appended in source order ahead of later positional arguments.

Requires an index: `.zfi` preferred, then `.fai`. Messy (non-uniform) records need a side-table `.zfi`; a plain `.fai` cannot describe them.

**Region formats:** `NAME` (full sequence); `NAME:START-END` (1-based inclusive); `NAME:START-` (from START through end of sequence). Ensembl-style names with colons work (for example `chromosome:GRCh38:1:1:248956422:1`).

**GET flags:**

- **`--bed file.bed`** / **`--bed -`**: BED regions from a file or stdin. Coordinates are 0-based half-open; z-fasta converts to 1-based inclusive internally.
- **`--names file.txt`**: one full-sequence name per line for long batch lists.
- **`--strand-aware`** (alias **`--honor-strand`**): read BED column 6; `-` applies reverse-complement before any global orientation flag.
- **`--rc`**: reverse-complement output. Verified against `samtools faidx -i --mark-strand no`. Mutually exclusive with **`--complement-only`** and **`--reverse-only`**.
- **`--complement-only`**: complement without reverse. Nucleotide records only.
- **`--reverse-only`**: reverse without complement.
- **`--annotate-rc`**: append a transform suffix to headers (for example `(reverse complement)`). Default output stays samtools-style.
- **`--summary`**: print region count, total bases, elapsed time, and regions/sec to stderr.
- **`--chunk-size N`**: BED batch size (default `4096`). Use `1` only for debugging. **`--chunk-size -1`**: load all BED rows in one batch when memory allows (BED inputs over 512 MiB are rejected with `-1`).
- **`--names`**: always loads the whole names file into memory, capped at 512 MiB. `--chunk-size` does not stream `--names`.

Positional CLI regions are capped at **1024** per invocation (`error: too many regions`). BED row count and **`--names`** line count are not subject to that cap (size limits above still apply).

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
z-fasta stats [options] <file.fasta>

Options:
  --index-only  Compute stats from index only (no FASTA scan; startup-dominated)
```

**Index-derived:** sequence count, total bases, min/max/mean/median lengths, N50, L50, N90, L90, AU, shortest/longest names, and duplicate reporting.

**Full scan (default):** composition on top of the index metrics: nucleotide frequencies, GC content (N excluded), GC skew, soft-masked fraction; for proteins, top 3 amino acids with full names.

### Sequence type classification

All commands share `stats.detectType` (IUPAC nucleotide alphabet; nucleotide if those letters are strictly more than 90% of counted bases). Sampling differs on purpose:

- `stats` (default): full-file composition scan
- `stats --index-only`: no type (composition skipped)
- `get --rc` / `--complement-only`: up to 256 bases per record
- `validate`: up to the first 100000 sequence bases; `--json --summary` reports `sequence_type`, `type_bases_sampled`, and `type_sample_cap`

### Examples

```bash
# Create .zfi binary index (default)
z-fasta index genome.fa

# Output .fai to stdout (uniform FASTA only)
z-fasta index --emit-fai genome.fa > genome.fai

# Keep duplicate names in the index (get still resolves to the last record)
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
z-fasta stats --index-only genome.fa
```

## Formats and behavior

### Index formats

- **`.zfi` (default, preferred):** binary index with embedded names, optional per-line side tables for messy layout, and source identity (size plus embedded FASTA mtime). Load prefers `.zfi` over `.fai`.
- **`.fai` (compatibility):** text index for uniform, FAI-representable records only. `--emit-fai` refuses variable widths and other non-representable layouts instead of writing a misleading file. Stale `.fai` checks are mtime-only (weaker than `.zfi`).

### Messy FASTA

Variable line widths, trailing spaces or tabs, blank lines, mixed CRLF/LF, and a final line without a terminal newline are indexed in `.zfi` via side tables. `get` uses those tables for correct extraction. Uniform records keep O(1) byte-offset math. Peers that require classic FAI geometry often reject these files; that is expected. `validate --fix` can rewrite a clean uniform copy when you need `.fai` or stricter tooling.

### Duplicate names

Duplicate _names_ are not the same as identical _sequence contents_. Default `index` keeps the first name; `index --no-dedup` keeps all. `get` resolves a repeated name to the **last** matching record in the loaded index. `validate` reports `duplicate_name`. Full `stats` prints source-level extras (`sum(k-1)`); `stats --index-only` prints `n/a (run without --index-only)` on a deduplicated index and never fabricates `0` (with an `--no-dedup` index it reports repeats kept in the index).

## Support and limits

**Formats:** `.zfi` is the default. It handles uniform records and messy layout (variable widths, trailing whitespace, blank lines, mixed CRLF/LF, missing final newline) via side tables. `.fai` is compatibility only for uniform, FAI-representable records; `--emit-fai` refuses messy layout. Load prefers `.zfi`. Compressed FASTA, BGZF, and FASTQ are out of scope.

**Names:** index hard-rejects sequence names longer than 65535 bytes. `validate --max-header-len` warns above N bytes (default 1024) and does not raise that index limit.

**Get / validate caps:** at most 1024 positional regions per `get` call. `--names` loads the whole file (max 512 MiB; `--chunk-size` does not stream it). BED defaults to 4096-row batches; `--chunk-size -1` loads the whole BED up to 512 MiB. `validate` stops after 10000 retained events.

**Memory:** index reads sequence payload through a bounded buffer. Get, stats, and validate retain mapped FASTA paths, so their RSS can approach the mapped file size. `stats --index-only` skips the sequence scan. Dense BED releases mapped pages behind the cursor; sparse gets on large FASTAs prefer positional reads.

**Platforms:** Linux, macOS, and Windows share the same commands through portable mapping. Memory-advice hints are POSIX-only (unused on Windows). CI builds and smokes six native platform lanes: x86_64 and arm64 on all three operating systems. Tagged releases publish the same six archives; v0.3.0 was the first complete six-archive release.

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

| Dataset                  | Region       | z-fasta | noodles | rust-bio | samtools | Speedup vs samtools |
| ------------------------ | ------------ | ------- | ------- | -------- | -------- | ------------------- |
| Genome (~2.9 GiB)        | 1 kbp region | 2.1 ms  | 2.5 ms  | 2.5 ms   | 3.2 ms   | **1.5x**            |
| Proteome (~13 MiB)       | 1 kbp region | 2.3 ms  | 6.7 ms  | 19.8 ms  | 12.7 ms  | **5.5x**            |
| Transcriptome (~459 MiB) | 1 kbp region | 4.4 ms  | 87.6 ms | 544.1 ms | 289.3 ms | **65.8x**           |

> Small-region extraction is O(1), but on this host the end-to-end CLI path is startup-dominated below roughly 10 kbp. For very large full-sequence extraction, fastahack can still win on raw write-path overhead; z-fasta stays ahead of samtools across the real-dataset GET cases.

**`--rc`** uses the same mmap-backed extraction path and applies reverse traversal plus complement lookup during emission instead of materializing a second copy of the region.

**Multi-region:** one call loads the index once and streams results in CLI order ([bench/get/REPORT.md](bench/get/REPORT.md) run `20260801_093130`, Transcriptome, 1 kbp per region):

| Regions | z-fasta | samtools | noodles  | Speedup vs samtools |
| ------- | ------- | -------- | -------- | ------------------- |
| 1       | 5.5 ms  | 294.2 ms | 88.1 ms  | **53.1x**           |
| 10      | 28.3 ms | 294.3 ms | 91.6 ms  | **10.4x**           |
| 100     | 28.2 ms | 298.3 ms | 129.0 ms | **10.6x**           |
| 1,000   | 37.2 ms | 297.0 ms | 462.0 ms | **8.0x**            |

> Benchmarked on REAL_Transcriptome.fa. Latency is dominated by index resolution and output setup rather than region byte count.

### Stats: Assembly/Proteome Statistics

| Mode       | Dataset                  | z-fasta     | seqkit -a | Speedup vs seqkit -a |
| ---------- | ------------------------ | ----------- | --------- | -------------------- |
| Index-only | Genome (~2.9 GiB)        | **2.11 ms** | 17.418 s  | **~8,250x**          |
| Index-only | Proteome (~13 MiB)       | **6.08 ms** | 63.3 ms   | **~10.4x**           |
| Full scan  | Genome (~2.9 GiB)        | **5.149 s** | 17.418 s  | **3.38x**            |
| Full scan  | Transcriptome (~459 MiB) | **0.870 s** | 2.390 s   | **2.75x**            |
| Full scan  | Proteome (~13 MiB)       | **25.7 ms** | 63.3 ms   | **2.46x**            |

> See [bench/stats/REPORT.md](bench/stats/REPORT.md). Index-only mode reads the sidecar without scanning sequence bytes.

### Correctness

- **Index:** `bench/index/run.sh` edge cases **24/24 matching cases plus one `binary_data` review**; messy layouts from `python3 bench/shared/generate_messy.py` into `bench/shared/cache/messy_fixtures/` (correctness) and `messy_perf/` (proteome perf).
- **Get:** `bench/get/run.sh` correctness **405/405** (positional, multi-region, BED, names, RC, messy) against samtools, bedtools, and seqtk where applicable.
- **Stats:** `bench/stats/run.sh` correctness **89/89** (BioPython oracle, index formats, layout twins, messy fixtures, duplicates policy, peer parity).
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
