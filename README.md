<!-- markdownlint-disable MD033 MD036 MD041 -->
<div align="center">
  <h1>z-fasta ⚡</h1>
  <p>
    Fast, modular FASTA toolkit built in Zig.<br/>
    SIMD-accelerated indexing, O(1) region extraction, validation, and assembly stats.<br/>
    samtools-compatible indexing and extraction, benchmarked against <code>seqkit</code>, <code>fastahack</code>, <code>pyfaidx</code>, and other peers.
  </p>
  <p>Current release: <strong>v0.3.0</strong></p>
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

`z-fasta` targets uncompressed FASTA: `index`, `get`, `stats`, and `validate`. Default indexes are compact `.zfi` files (embedded names, optional side tables, source identity). `--emit-fai` writes samtools-compatible `.fai` only when every record has fixed line geometry.

`get` accepts positional regions, BED files, BED from stdin, names files, strand-aware extraction, and orientation transforms (`--rc`, `--reverse-only`, `--complement-only`, `--annotate-rc`). FASTQ and compressed FASTA/BGZF remain out of scope.

CI builds and smokes the binary on Linux, macOS, and Windows. File views go through a portable `MemoryMap` layer; POSIX memory advice is an optimization and a no-op on Windows.

## Why z-fasta?

Modern bioinformatics workflows are often bottlenecked by legacy text parsers. `z-fasta` keeps the hot paths close to the data: memory-mapped FASTA for the default indexer, SIMD header scanning, compact binary indexes, and a startup-conscious CLI for tiny commands.

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
  --low-mem     Stream input with bounded RAM; same on-disk outputs as default index
  --help        Show help message
  --version     Print version
```

Default `index` and `--low-mem` both write `{file}.zfi` with matching bytes on the supported fixtures. `--low-mem` streams the FASTA instead of mapping it; use it when RSS must stay near a few megabytes. Sequence names longer than 65535 bytes are rejected (`HeaderTooLong`).

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

# Bounded-RAM index (same .zfi bytes as default on supported fixtures)
z-fasta index --low-mem genome.fa

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

Duplicate *names* are not the same as identical *sequence contents*. Default `index` keeps the first name; `index --no-dedup` keeps all. `get` resolves a repeated name to the **last** matching record in the loaded index. `validate` reports `duplicate_name`. Full `stats` prints source-level extras (`sum(k-1)`); `stats --index-only` prints `n/a (run without --index-only)` on a deduplicated index and never fabricates `0` (with an `--no-dedup` index it reports repeats kept in the index).

## Support and limits

**Formats:** `.zfi` is the default. It handles uniform records and messy layout (variable widths, trailing whitespace, blank lines, mixed CRLF/LF, missing final newline) via side tables. `.fai` is compatibility only for uniform, FAI-representable records; `--emit-fai` refuses messy layout. Load prefers `.zfi`. Compressed FASTA, BGZF, and FASTQ are out of scope.

**Names:** index hard-rejects sequence names longer than 65535 bytes. `validate --max-header-len` warns above N bytes (default 1024) and does not raise that index limit.

**Get / validate caps:** at most 1024 positional regions per `get` call. `--names` loads the whole file (max 512 MiB; `--chunk-size` does not stream it). BED defaults to 4096-row batches; `--chunk-size -1` loads the whole BED up to 512 MiB. `validate` stops after 10000 retained events.

**Memory:** default index/get/stats/validate use mapped file views (RSS often near mapped FASTA size). `index --low-mem` streams with bounded RAM and the same `.zfi` bytes on supported inputs. `stats --index-only` skips the sequence scan. Dense BED releases mapped pages behind the cursor; sparse gets on large FASTAs prefer positional reads.

**Platforms:** Linux, macOS, and Windows share the same commands through portable mapping. Memory-advice hints are POSIX-only (unused on Windows). CI smokes all three. Tagged multi-triple release archives are not published yet.

## Performance & Correctness

All timings below are on AMD Ryzen 9 3950X, warm cache, from the checked-in suite reports (subject binary was labeled 0.2.9 in those runs; behavior matches current gates). Regenerate reports after the v0.3.0 tag if you need subject strings to match exactly.

### Index: SIMD-Accelerated Indexing

| Dataset       | Size (on disk) | z-fasta (no-dedup) | samtools | fastahack | pyfaidx | Speedup vs samtools |
| ------------- | -------------- | ------------------ | -------- | --------- | ------- | ------------------- |
| Human Genome  | ~2.9 GiB       | 0.39s              | 9.03s    | 21.73s    | 27.48s  | **22.9x**           |
| Transcriptome | ~459 MiB       | 0.093s             | 1.79s    | 5.72s     | 6.50s   | **19.3x**           |
| Proteome      | ~13 MiB        | 0.0056s            | 0.055s   | 0.275s    | 0.368s  | **10.0x**           |

**Index modes** on Genome (warm cache; [bench/index/REPORT.md](bench/index/REPORT.md) run `20260706_134943`): **default** mmap ~0.40s, RSS near mapped FASTA size; **`--low-mem`** stream ~1.61s with ~3.4 MB RSS and the same `.zfi` bytes as default; **`--no-dedup`** ~0.38s; **`--emit-fai`** writes FAI to stdout only (no on-disk `.zfi`).

> Mapped modes show RSS close to the FASTA size because `/usr/bin/time -v` counts mapped pages, not only private heap. Full curves: [bench/index/REPORT.md](bench/index/REPORT.md).

### Get: O(1) Region Extraction

| Dataset                  | Region          | z-fasta        | samtools   | seqtk    | pyfaidx | Speedup vs samtools |
| ------------------------ | --------------- | -------------- | ---------- | -------- | ------- | ------------------- |
| Any (warm cache)         | 100 bp - 10 kbp | **0.7-0.9 ms** | 1.5-1.6 ms | 4-34 ms  | ~60 ms  | **1.8-2.1x**        |
| Proteome (~13 MiB)       | 1 kbp region    | 1.3 ms         | 10.9 ms    | 7.2 ms   | 119 ms  | **8.4x**            |
| Transcriptome (~459 MiB) | 1 kbp region    | 25.3 ms        | 278.7 ms   | 220.3 ms | 1103 ms | **11.0x**           |

> Small-region extraction is O(1), but on this host the end-to-end CLI path is startup-dominated below roughly 10 kbp. For very large full-sequence extraction, fastahack can still win on raw write-path overhead; z-fasta stays ahead of samtools across the real-dataset GET cases.

**`--rc`** uses the same mmap-backed extraction path and applies reverse traversal plus complement lookup during emission instead of materializing a second copy of the region.

**Multi-region:** one call loads the index once and streams results in CLI order ([bench/get/REPORT.md](bench/get/REPORT.md) run `20260709_094913`):

| Regions | z-fasta | samtools | seqtk  | Speedup vs samtools |
| ------- | ------- | -------- | ------ | ------------------- |
| 1       | 25.6 ms | 289 ms   | 221 ms | **11.3x**           |
| 10      | 33.8 ms | 283 ms   | 226 ms | **8.4x**            |
| 50      | 66.7 ms | 292 ms   | 225 ms | **4.4x**            |
| 100     | 66.7 ms | 279 ms   | 222 ms | **4.2x**            |

> Benchmarked on REAL_Transcriptome.fa. Latency is dominated by index resolution and output setup rather than region byte count. seqtk performs a full-file scan per call and is listed for reference only.

### Stats: Assembly/Proteome Statistics

| Mode       | Dataset              | z-fasta     | seqkit -a | seqtk comp | Speedup vs seqkit -a |
| ---------- | -------------------- | ----------- | --------- | ---------- | -------------------- |
| Index-only | Genome (~2.9 GiB)    | **0.9 ms**  | 17.45 s   | N/A        | **~19,000x**         |
| Index-only | Proteome (~13 MiB)   | **2.9 ms**  | 57.8 ms   | N/A        | **~20x**             |
| Full scan  | 1 GB single-seq file | **0.78 s**  | 5.62 s    | 2.65 s     | **~7x**              |
| Full scan  | Proteome (~13 MiB)   | **11.8 ms** | 57.8 ms   | 93.0 ms    | **~4.9x**            |

> See [bench/stats/REPORT.md](bench/stats/REPORT.md). Index-only mode reads the sidecar without scanning sequence bytes.

### Correctness

- **Index:** `bench/index/run.sh` edge cases **25/25**; messy layouts from `python3 bench/shared/generate_messy.py` into `bench/shared/cache/messy_fixtures/` (correctness) and `messy_perf/` (proteome perf).
- **Get:** `bench/get/run.sh` correctness **409/409** (positional, multi-region, BED, names, RC, messy, low-mem parity) against samtools, bedtools, and seqtk where applicable.
- **Stats:** `bench/stats/run.sh` correctness **95/95** (BioPython oracle, index formats, layout twins, messy fixtures, duplicates policy, peer parity).
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
./zig build -Doptimize=ReleaseFast
./zig-out/bin/z-fasta --version
```

## Roadmap

**Delivered through v0.3.0**

- [x] `index`, `get`, `stats`, and `validate`
- [x] Preferred `.zfi` with side tables; `.fai` only for representable uniform records
- [x] Messy FASTA GET via side tables; `validate --fix` for safe rewrites
- [x] `index --low-mem` streaming path with mmap parity
- [x] Portable Linux / macOS / Windows mapping and CI smoke
- [x] Zebrac benchmark suites for index, GET, and stats with verify gates

**Deferred (not v0.3.0)**

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
