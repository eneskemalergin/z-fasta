<!-- markdownlint-disable MD024 -->
# Changelog

All notable changes to z-fasta will be documented in this file.

## [0.3.0] - Unreleased

### Added

- `src/validator.zig`: new `validate` subcommand for FASTA structure, alphabet, and header validation with human-readable text, JSON Lines, and JSON summary output modes
- `z-fasta validate`: checks for duplicate names, invalid characters, null bytes, UTF-8 BOM, inconsistent line widths, trailing whitespace, empty sequences, missing terminal newlines, mixed line endings, long headers, and schema violations
- `z-fasta validate --strict`: promotes all warnings to errors (exit code 1)
- `z-fasta validate --json`: streaming JSON Lines output with `schema_version: "v1"` on every event
- `z-fasta validate --json --summary`: single JSON object with counts and first examples per kind
- `z-fasta validate --fix -o <file>`: writes a fixed FASTA with format-level issues resolved (BOM stripped, lines rewrapped to modal width, line endings normalized, trailing whitespace removed, terminal newline appended)
- `z-fasta validate --fix-format-only`: proceeds with format fixes even when character-level errors exist
- `z-fasta validate --schema uniprot` and `--schema refseq`: header schema validation for UniProt (`sp|`, `tr|`, `db|`) and RefSeq (`NC_`, `NM_`, `NR_`, etc.) formats
- `z-fasta validate --custom-alphabet "ACGTN..."`: overrides the built-in nucleotide or protein alphabet
- `z-fasta validate --max-header-len N`: warns when headers exceed N bytes (default: 1024)
- `src/index_format.zig`: `.zfi` v2 format with per-record `is_uniform_width` flag and side-table offset for non-uniform FASTA records; backward-compatible with v0.2.x (uniform records are byte-identical)
- `src/indexer.zig`: `scanZfiIndex()` detects non-uniform-width records and writes side-tables with per-line `(base_start, byte_offset, line_bytes, line_bases)` triples
- `src/getter.zig`: binary search over side-table for non-uniform records; O(log L) per base lookup where L is the number of lines in the sequence
- `bench/index/run.sh`: unified index suite runner combining correctness tests, zebrac performance benchmarks, messy FASTA zebrac, and report generation
- `bench/shared/zebrac_runner.sh`: shared zebrac harness with configurable runs, warmup, duration, and metadata JSONL output
- `bench/shared/tools.sh`: source of truth for local tool paths and Tier 1/Tier 2 labels
- `tests/test_index.zig`: side-table validation tests for v0.2 to v0.3 index compatibility

### Changed

- `src/index_format.zig`: `.zfi` magic bumped from `ZFI\x01` to `ZFI\x02`; v0.3 reads both v1 and v2 formats
- `src/indexer.zig`: `scanFastaRecords` now passes sequence data and uniform-width flag to the emit callback for side-table construction
- `src/main.zig`: index writing uses `scanZfiIndex()` with in-memory record and side-table arrays instead of streaming write with a dummy header
- `src/stats.zig`: composition scan handles non-uniform records via side-table; whitespace check uses `byte > ' '` instead of explicit `\n`/`\r` comparison
- `build.zig`: added `test_validator` target for `src/validator.zig` unit tests
- `bench/save_baseline.py` and `bench/compare_baseline.py`: replace the old hyperfine-only baseline path with normalized `bench.baseline.v2` snapshots that ingest zebrac JSON, hyperfine JSON, metadata JSONL, memory CSVs, and throughput CSVs
- `bench/README.md`: rewritten with new layout, runner documentation, and verification workflow
- `bench/get/verify_get.sh` and `bench/get/verify_bed.sh`: added messy FASTA verification cases
- `bench/stats/verify_stats.py`: added messy FASTA stats verification
- `README.md`: updated benchmark commands to use `bench/index/run.sh` and the new verification workflow

### Fixed

- `src/getter.zig`: `get` on messy (non-uniform-width) FASTAs now retrieves correctly using the v0.2 index format extension; previously emitted garbage or failed silently on mixed-width files
- `src/index_format.zig`: v0.2.x binaries emit a clear upgrade error when encountering v0.3 `.zfi` files with the new magic

### Removed

- `bench/wrappers/`: retire the obsolete standalone Tier 2 hyperfine lane now that Rust wrappers are integrated into the main index/get/stats suites
- `bench/index/run_benchmarks.sh` and `bench/index/run_tests.sh`: replaced by `bench/index/run.sh`
- `bench/index/results/figures/scaling_seqs.png` and `speedup.png`: removed legacy figure files

## [0.2.9] - 2026-06-24

Various memory safety, optimizations and re-building benchmarking framework to work with less dependencies. Also some code cleanup, standardization.

### Changed

- `src/main.zig`: write `.zfi` indexes via a temp file and rename on success
- `src/getter.zig`: cap multi-region sort buffers at 64 MiB per region and 256 MiB total; reject BED/names inputs over 512 MiB when `--chunk-size -1`
- `src/getter.zig`: use `MADV.SEQUENTIAL` on multi-region sort-path reads; keep `MADV.RANDOM` for BED, reverse reads, and batches under 16 regions
- `src/index_format.zig`: reject `.zfi` records whose `seq_offset` is past EOF
- `src/indexer.zig`: restore fast `seq_len` counting on uniform records while keeping v0.2.8 validation on messy records
- `src/stats.zig`: SIMD `countCompositionSlice` on the fixed-width full-scan path
- `.gitignore`: ignore Python bytecode, caches, and `.venv/`
- `src/main.zig`: route CLI errors through `index_format.printErrorAndExit`; flush stdout before `--emit-fai` error exit; VERSION bumped to 0.2.9; help text for `--low-mem`, `--no-dedup`, region cap, and `--chunk-size 1`
- `build.zig.zon`: package version bumped to 0.2.9
- `src/indexer.zig`: flush stdout before `--low-mem` error exit
- `src/getter.zig`: clearer BED `StreamTooLong` error (line number and 4096-byte buffer limit); `detectRecordType` documents base-not-byte sampling
- `src/stats.zig`: reuse one 64-byte `formatComma` buffer (byte-identical output)
- `README.md`: index/get flag docs aligned with CLI help

### Added

- `bench/save_baseline.py`, `bench/compare_baseline.py`, `bench/run_all_and_baseline.sh`: local baseline snapshots under `bench/baselines/` for regression diffing
- `src/getter.zig`: unit tests for sort-buffer caps and sequential `madvise` gating
- `src/stats.zig`: unit test for `countCompositionSlice`
- `tests/test_index.zig`: `seq_len` edge-case tests and `seq_offset` EOF rejection test

### Removed

- Tracked `__pycache__/*.pyc` files
- `src/getter.zig`: dead `resolveRegionsByRecordScan`

### Validation

- `./zig build test --summary all` (108/108 tests passed)
- `./zig build -Doptimize=ReleaseFast`
- `bash bench/index/run_tests.sh` (20/20 passed)
- `bash bench/get/verify_get.sh` (90/90 passed)
- `bash bench/get/verify_multi_get.sh` (22/22 passed)
- `bash bench/get/verify_bed.sh` (16/16 passed)
- `bash bench/get/verify_rc.sh` (19/19 passed)
- `.venv/bin/python bench/stats/verify_stats.py` (107/107 passed)

## [0.2.8] - 2026-05-15

### Added

- `z-fasta get --rc`: reverse-complement extraction for positional regions, names-file batches, and BED-driven extraction
- `z-fasta get --complement-only`: complement extraction without reversing
- `z-fasta get --reverse-only`: reverse extraction without complementing
- `z-fasta get --annotate-rc`: optional transform annotation for FASTA headers
- `bench/get/verify_rc.sh`: RC / reverse / complement verification against `samtools`, `bedtools`, and `seqtk`
- `bench/get/RC_STRATEGY.md`: checked-in note for the shipped reverse strategy and the current measurement slice

### Changed

- `src/getter.zig`: extraction now carries a composable orientation model through positional, names-file, and BED request resolution
- `src/getter.zig`: complement-based transforms now reject protein sequences with a clear error instead of emitting biologically invalid output
- `src/getter.zig`: complement validation now caches the last detected record type so BED / batch complement paths do not rescan the same record sample repeatedly
- `src/stats.zig`: nucleotide/protein detection now treats the full IUPAC nucleotide alphabet, including lowercase ambiguity codes and `U/u`, as nucleotide-like
- `src/indexer.zig`: fixed wrapped short-tail sequence length accounting so whole-sequence reverse extraction does not start from a trailing newline
- `bench/get/run_benchmarks.sh` and `bench/get/generate_report.py`: the main GET benchmark/report flow now includes the RC review slice by default, including RC timing and RSS comparisons against `samtools faidx -i` and `bedtools + seqtk` where applicable
- `README.md` and `bench/get/README.md`: release-surface docs updated for `--rc`, `--reverse-only`, `--complement-only`, `--annotate-rc`, and the checked-in reverse-strategy note

### Validation

- `./zig build test --summary all` (102/102 tests passed)
- `./zig build -Doptimize=ReleaseFast`
- `bash bench/get/verify_get.sh` (90/90 passed)
- `bash bench/get/verify_multi_get.sh` (22/22 passed)
- `bash bench/get/verify_bed.sh` (16/16 passed)
- `bash bench/get/verify_rc.sh` (19/19 passed)

## [0.2.7] - 2026-05-13

### Added

- `src/complement.zig`: shared IUPAC complement and reverse-complement helpers for strand-aware extraction work
- `src/bed_parser.zig`: BED line parser with 0-based half-open to 1-based inclusive conversion, comment skipping, and optional strand capture
- `z-fasta get --bed <path|->`: BED-driven extraction from files or stdin
- `z-fasta get --names <path>`: batch whole-sequence extraction from a plain text names file
- `z-fasta get --strand-aware`: strand-aware BED extraction with reverse-complement output for `-` strand rows
- `z-fasta get --honor-strand`: compatibility alias for `--strand-aware`
- `z-fasta get --summary`: stderr-only extraction summary with region count, total bases, elapsed time, and regions/sec

### Changed

- `src/getter.zig`: unified positional-region, BED, names-file, stdin, and strand-aware extraction through the same resolved-region path
- `src/main.zig`: `get` help text and CLI parsing updated for BED, names, strand, and summary flags
- `src/main.zig` and `src/getter.zig`: default BED `--chunk-size` lowered to `4096` after benchmark validation showed the 4K batch size matched or beat larger batches while using materially less memory
- `src/main.zig`: VERSION bumped to 0.2.7
- `build.zig.zon`: package version bumped to 0.2.7
- `build.zig`: standalone module tests for `src/complement.zig` and `src/bed_parser.zig` now run as part of `zig build test`
- `zig`: repo-local wrapper now points at the vendored Zig 0.16.0 toolchain instead of stale 0.15/0.14 paths

### Validation

- `zig build test --summary all` (99/99 tests passed)
- `zig build -Doptimize=ReleaseFast`
- `bench/get/verify_bed.sh` (16/16 passed)
- direct release-binary checks for `--bed`, `--bed -`, `--names`, `--strand-aware`, `--honor-strand`, and `--summary`

## [0.2.6] - 2026-04-27

### Added

- `bench/perf-recovery/run_startup.sh`: startup-floor benchmark harness for Zig 0.16 CLI entry-point overhead
- `bench/perf-recovery/README.md`: notes for the performance recovery harness

### Changed

- `src/main.zig`: VERSION bumped to 0.2.6 and CLI startup moved to `std.process.Init.Minimal`, with direct Linux stdout writes for top-level `--help` and `--version`
- `src/index_format.zig`: split index loading by use case so stats and small lookups can avoid rebuilding the full name map
- `src/getter.zig`: small multi-region batches resolve in one record scan; sequence emission now uses buffered chunk writes
- `src/stats.zig`: stats uses records-only loading and a fixed-width composition fast path for regular FASTA layouts
- `src/indexer.zig`: fixed-width records avoid recounting sequence bytes during indexing when line metrics are regular
- `README.md` and benchmark reports refreshed with v0.2.6 performance data
- `build.zig.zon`: package version bumped to 0.2.6

### Validation

- `zig build -Doptimize=ReleaseFast`
- `zig build test --summary all` (86/86 tests passed)
- `bench/index/run_tests.sh` (20/20 edge cases match samtools)
- `bench/get/verify_get.sh` (90/90 passed)
- `bench/get/verify_multi_get.sh` (22/22 passed)
- `bench/stats/verify_stats.py` (107/107 passed)

## [0.2.5] - 2026-04-27

### Changed

- Migrated the toolchain from Zig 0.14.0 to Zig 0.16.0
- `build.zig`: updated to the Zig 0.16 module-based build API for executable and test targets
- `src/main.zig`: migrated CLI entry to `std.process.Init` and `std.process.Args.Iterator`
- `src/index_format.zig`, `src/indexer.zig`, `src/getter.zig`, `src/stats.zig`: migrated file and stream handling from `std.io`/`std.fs` call sites to `std.Io`
- Collections updated for Zig 0.16 unmanaged APIs (`ArrayList`, related allocator-explicit append/deinit patterns)
- Tests updated for Zig 0.16 I/O, allocator, process spawning, and timestamp APIs
- `.github/workflows/ci.yml`: Zig CI pin updated to 0.16.0
- `README.md`: version and Zig installation/docs references updated to 0.16.0 / v0.2.5
- `build.zig.zon`: package version bumped to 0.2.5 and minimum Zig version kept at 0.16.0

### Validation

- `zig build`
- `zig build test`
- `zig build -Doptimize=ReleaseFast`

## [0.2.4] - 2026-04-02

### Added

- **Multi-region `get`**: `z-fasta get <file.fasta> <region> [region ...]` now accepts multiple regions in a single invocation, loads the index once, and streams results in CLI order
- `bench/get/verify_multi_get.sh`: byte-identical multi-region verification against `samtools faidx`, covering duplicate regions, reversed CLI order, and the sort-path used for larger region lists
- `bench/get/run_benchmarks.sh`: multi-region extraction benchmark module for 1, 10, 50, and 100 regions in a single CLI call
- `bench/get/generate_report.py`: multi-region benchmark loader and report section generation

### Changed

- `src/getter.zig` and `src/main.zig`: CLI and extraction flow updated to resolve and emit multiple requested regions per invocation without re-loading the index
- `tests/test_get.zig`: expanded coverage for multi-region parsing, ordering, duplicate requests, and output behavior
- `README.md`: `get` command usage updated to document multi-region extraction, correctness counts updated, and a v0.2.4 performance table added
- `bench/get/REPORT.md`: added the Multi-Region Extraction section and updated figures
- `src/main.zig`: VERSION bumped to 0.2.4
- `build.zig.zon`: version bumped to 0.2.4

## [0.2.3] - 2026-03-30

### Added

- **pyfaidx** added to index and get benchmark suites (`faidx --no-output` for indexing, `faidx file region` for extraction)
- **seqtk** added to get and stats benchmark suites (`seqtk subseq` via BED file for extraction, `seqtk comp` for per-sequence composition)
- **seqkit stats -a** replaces `seqkit stats` in stats benchmarks (adds N50, GC% to comparison)
- `bench/shared/install_tools.sh`: verification helper for all pinned tool versions; installs pyfaidx into `.venv` automatically
- **Messy FASTA benchmark** (`bench/shared/messy_variants/`): four derived FASTA variants (mixed_widths, crlf_endings, trailing_whitespace, all_messy) benchmarked against all indexing tools; `compatibility.csv` with per-tool exit-code results
- Messy FASTA compatibility section added to `bench/index/REPORT.md` (auto-generated from `generate_report.py`)

### Changed

- `bench/index/generate_report.py`: pyfaidx added to COLORS/TOOL_ORDER; figure paths corrected to `results/figures/`; messy compatibility section added; em-dashes removed from all generated prose
- `bench/get/generate_report.py`: seqtk and pyfaidx added; figure paths corrected; em-dashes removed
- `bench/stats/generate_report.py`: seqtk-comp added; seqkit-stats renamed to seqkit-stats-a throughout; figure paths corrected; em-dashes removed
- `README.md`: version bumped to v0.2.3; performance tables updated with pyfaidx and seqtk columns; messy FASTA correctness bullet added; roadmap updated
- `src/main.zig`: VERSION bumped to 0.2.3
- `build.zig.zon`: version bumped to 0.2.3

## [0.2.2] - 2026-03-08

### Added

- **GET benchmark suite** (`bench/get/`): hyperfine timing, RSS memory, region-size scaling, real dataset benchmarks; `--skip-real` flag; `verify_get.sh` (90 diff tests); rewritten `generate_report.py` with human-readable names, per-module figures, and auto-generated `REPORT.md`
- **STATS benchmark suite** (`bench/stats/`): index-only vs. full-scan, file-size scaling (1 MB – 1 GB), throughput CSV; `--skip-real` flag; `verify_stats.py` (107 BioPython tests); rewritten `generate_report.py` with μs index-only speedup table, seqkit in all comparisons, and auto-generated `REPORT.md`
- `bench/README.md`: GET and STATS added to Quick Start; `--skip-real` documented

### Fixed

- Memory CSV column order mismatch in `bench/stats/run_benchmarks.sh`
- Seqkit throughput not written to CSV due to incorrect `2>&1` redirect in time group
- figure folder exceptions in .gitignore to render report images in markdowns in GitHub UI

### Changed

- `README.md`: test count corrected (63/63); `--index-only` timing corrected to `< 1 ms`; performance section split into Index / Get / Stats with per-module tables and report links; tagline and description clarified (SIMD scoped to indexer, arena-allocated); roadmap reorganized

## [0.2.1] - 2026-03-07

### Changed

- Removed samtools dependency from unit tests (`test_get.zig`)
    - 17 samtools comparison integration tests removed (coverage moved to `bench/verify_get.sh`)
    - `test_get.zig` now contains only self-contained `parseRegion` unit tests and `loadIndex` tests (12 tests)
    - Total test count: 63 (19 index + 12 get + 32 stats)
- CI no longer installs samtools; `Generate test indexes` step uses only `z-fasta index`
    - Eliminated the `samtools faidx` call from the CI index generation loop

### Fixed

- CI fresh-checkout failures caused by missing `.zfi`/`.fai` index files (gitignored)
    - Added `Generate test indexes` step to CI that runs `z-fasta index` on all test FASTA files before `zig build test`

## [0.2.0] - 2026-03-06

### Added

- **`z-fasta get <file.fasta> <region>`**: O(1) byte-offset sequence extraction
    - Output byte-identical to `samtools faidx` (verified via `diff` on 90 test cases)
    - Region formats: `NAME`, `NAME:START-END`, `NAME:START-`
    - Handles Ensembl-style names with embedded colons (right-to-left parsing)
    - mmap + MADV_RANDOM for point-access extraction
    - Coordinate clamping matches samtools behavior (END > seq_len silently clamped)
- **`z-fasta stats <file.fasta>`**: assembly/proteome statistics
    - Tier 1 (index-only, `--index-only`): sequence count, total bases, min/max/mean/median, N50, L50, N90, L90, AU, duplicates -- completes in <500 μs
    - Tier 2 (full scan): branchless composition counting, nucleotide vs protein auto-detection
    - Nucleotide: A/C/G/T/N frequencies, GC content (N excluded), GC skew, soft-masked fraction
    - Protein: top 3 amino acids with full names, lowercase fraction
- **Shared index loading** (`index_format.zig`): .zfi preferred, .fai fallback with mtime/size staleness checks
- **Test suites:** 80 unit tests (19 index + 29 get + 32 stats), 90 samtools-diff tests, 107 BioPython verification tests
- `bench/verify_get.sh` -- automated byte-identical comparison against samtools
- `bench/verify_stats.py` -- automated stats verification against BioPython

### Changed

- Source split into 5 modules: `main.zig`, `indexer.zig`, `getter.zig`, `stats.zig`, `index_format.zig`
- `build.zig` updated with 3 test targets sharing a single main module
- `main.zig` is now a slim CLI dispatcher (~280 lines)

### Technical

- Zig 0.14.0
- .fai parsing: `ArenaAllocator`-owned name strings fed into `StringHashMap`
- Extraction: O(1) byte-offset formula + newline-skipping loop, wraps at 60 chars/line
- Composition: single-pass 256-element `u64` array indexed by byte value, branchless
- N50/L50/N90/L90/AU computed in one sort + one pass

## [0.1.0] - 2026-02-28

### Added

- Initial release of z-fasta FASTA indexer
- SIMD-accelerated header scanning (32-byte vector chunks)
- Memory-mapped I/O with sequential read optimization
- Output byte-identical to `samtools faidx`
- Buffered stdout for efficient multi-sequence output
- Low-memory chunked mode (`--low-mem`) for constant 4 MB usage
- Duplicate name filtering on by default (`--no-dedup` to disable)
- `--help` and `--version` flags
- ZFI binary index format alongside FAI text output
- Comprehensive benchmark suite (`bench/`)
- Edge case test generator and correctness tests

### Performance

- 6-16x faster than samtools on real datasets:
    - Human Genome (3.0 GB): 0.57s vs 9.15s (16.1x)
    - Transcriptome (972 MB): 0.20s vs 1.80s (9.0x)
    - Proteome (66 MB): 0.008s vs 0.054s (6.8x)
- Single-threaded, allocation-free scanning

### Technical

- Zig 0.14.0
- mmap + madvise(SEQUENTIAL)
- ArenaAllocator for leak-safe memory
- ~670 lines of code
