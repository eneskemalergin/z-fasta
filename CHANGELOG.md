<!-- markdownlint-disable MD024 -->
# Changelog

All notable changes to z-fasta will be documented in this file.

## [0.3.0] - Unreleased

Feature-parity release: `validate` subcommand, messy-FASTA side tables, streaming `index --low-mem`, RSS caps on `get` and `stats` composition scans, embedded `.zfi` name table, and a complete benchmark suite rebuild (hyperfine to zebrac).

### Added

- **`.zfi` source identity (`ZFID`)**: production indexes store the FASTA mtime immediately before the `ZFNM` footer. Load rejects size or embedded-mtime mismatch; legacy files without `ZFID` keep the weaker index-mtime check. `.fai` remains mtime-only.
- **Portable file mapping (`platform.zig`)**: `std.Io.File.MemoryMap` wrapper for FASTA/index views; Windows-compatible `Args.Iterator.initAllocator`; memory advice is a no-op on Windows. Windows mappings use `populate = true` to avoid a Zig 0.16 `NtCreateSection` bug with empty allocation attributes.
- **`LoadedIndex` ownership**: maps destroyed via `deinit(io)`; arena owns heap (records for `.fai`, name copies, sidecar path); name-map keys are borrows or arena copies; no separate `page_allocator` name-map path.
- **`.zfi` embedded name table**: sequence names are stored in the index at build time (`ZFNM` footer + name blob). `lookup_full_map` builds a pointer hash over the blob without touching the FASTA header scan path. **Re-index required** for existing `.zfi` files to pick up the new layout.
- **`.zfi` name footer**: tight 12-byte on-disk footer (4-byte magic + u64 length); legacy 16-byte footers remain readable.
- **`get` side-table lookup**: binary search on side tables for non-uniform records (O(log L) per base; L is line count).
- **`z-fasta validate`**: structure, alphabet, and header checks (duplicate names, invalid characters, null bytes, UTF-8 BOM, inconsistent line widths, trailing whitespace, empty sequences, missing terminal newlines, mixed line endings, long headers, schema violations). Output modes: human text; `--json` JSON Lines (`schema_version: "v1"` per event); `--json --summary` aggregate object. `--strict` promotes warnings to errors. `--fix -o <file>` streams rewritten FASTA (does not materialize a full second copy solely to write). Event list is capped at 10000 (`Summary.truncated`); CLI prints retained events then exits with a deterministic error. JSON string fields accept arbitrary FASTA header bytes (invalid UTF-8 becomes `\u00XX`) with no 256-byte message cap. Also `--fix-format-only`, `--schema uniprot`, `--schema refseq`, `--custom-alphabet`, and `--max-header-len N` (default 1024).
- **`index --low-mem` streaming `.zfi`**: bounded-RAM build path. Output bytes match mmap `index` on simple and messy fixtures.
- **Index dedup identity**: mmap and streaming share collision-safe name comparison (`NameDedup`); hash may accelerate lookup but never decides alone.
- **`tests/test_index.zig`**: v0.2 to v0.3 side-table and index compatibility coverage; `.fai` `records_only` vs `lookup_full_map` duplicate-name parity test; malformed-index and stale-identity coverage.
- **`build.zig`**: `test_validator` target for validator unit tests.
- **`complement.complementInto`**: chunked IUPAC complement into a caller buffer (shared by forward and reverse GET emit).
- **`bench/stats/verify.sh`**: L2 stats gate (BioPython oracle, `.zfi`/`.fai`/cross, layout twins messy==uniform, messy fixtures, low-mem, dedup + Duplicates policy, per-tool parity, CLI errors; 92 checks).
- **`bench/stats/run.sh` + `generate_report.py`**: L3 stats perf (full peers, four-way mode, size + seq-count scaling) and `REPORT.md` / figures. Indexes preloaded; timed work is `stats` / peers only.
- **Layout twins** (generated under `bench/stats/data/verify/layout_twins/` by `verify.sh`): same-base FASTA layouts (uniform vs mixed widths / trailing ws / blank lines / mixed CRLF) for layout-invariant stats checks.
- **`tools/noodles_wrapper` / `tools/rustbio_wrapper` `stats`**: clean-FASTA TSV comparison peers with assembly + composition fields (N50/N90/AU, GC skew, protein top-AA) so benches/verify can compare more than sequences/bases. Not messy/side-table capable.

### Changed

- **CLI unknown options**: `index`, `get`, `stats`, and `validate` reject unrecognized `-`/`--` tokens with a nonzero exit instead of treating them as a FASTA path or region.
- **CLI failure subprocess tests**: failure paths assert exact exit status, exact stderr (usage dumps must match `--help` stdout), and stdout contracts (empty on hard errors; exact `WARNING` text on validate exit 2) instead of substring-only checks.
- **Required vs optional fixtures**: missing in-repo fixtures fail with `required test fixture missing: <path>`. Required fixture directory walks must index at least one FASTA and must not swallow index errors. Generated `bench/index/edge_cases/` and downloaded `REAL_*` tests skip via `SkipZigTest` with obtain instructions (no silent pass).
- **Stable gate fixtures** (`tests/data/gates/`): duplicate names, long headers (indexable and over-`u16` reject), invalid UTF-8 headers, large-offset and zero-geometry `.fai`, plus injectable `NameDedupWith` forced-collision check proving identity is string equality.
- **GET fixture seeds**: `bench/get/run.sh` multi/BED region generation uses SHA-256-derived seeds instead of process-salted Python `hash()`, so regenerated fixtures are stable across processes and platforms.
- **Formatter and file headers**: `zig fmt --check src tests` is clean (`tests/test_index.zig` / `tests/test_get.zig` drift fixed). Missing `//!` module headers added across `src/` and `tests/` per `plan/STYLE.md` (comments only).
- **Bench harness prose**: dropped deleted-plan citations (`plan/stats-bench.md`), fixed runner names (`run.sh` / `verify.sh` instead of `run_benchmarks.sh` / `run_tests.sh` / `run_messy_benchmarks.sh`), and corrected `.zfi` preference vs `.fai` compatibility wording in suite scripts and REPORT generators.
- **Dataset download manifest**: `bench/shared/datasets.manifest` pins URL, uncompressed size, and sha256 for Genome / Transcriptome / Proteome. `download_data.sh` skips matching sizes, verifies after download, and supports `--verify` for full hash checks. Genome and Transcriptome URLs are release-pinned; Proteome uses UniProt `current_release` (refresh the pin when that file moves).
- **Bench tool verification**: `install_tools.sh` sources `tools.sh`, checks zebrac and the Tier-1/2 peers suites use (no hyperfine requirement), and uses ASCII status lines. rust-bio pin corrected to **2.2** (matches `tools/rustbio_wrapper/Cargo.toml`); unused `bench_tool_tier` / `label` / `supports_suite` helpers removed.
- **CI-required fixtures**: track `tests/data/gates/*.fai` (were blocked by `*.fai` gitignore). Mark intentional CRLF / corrupt-index fixtures as `-text` in `.gitattributes` so Git does not strip `\r` on checkout.
- **Bench shared dedup**: `get`/`stats` `verify.sh` and `download_data.sh` source `tools.sh` for tool paths; `file_size_bytes` lives in `tools.sh` (zebrac runner no longer redefines it); orphaned `messy_variants/compatibility.csv` removed.
- **Messy fixture dirs**: tiny index correctness fixtures live in `bench/index/messy_fixtures/`; proteome-derived GET/index perf fixtures live in `bench/shared/messy_perf/` (was both named `messy_variants`). Index messy zebrac cleans leftover `.fai`/`.zfi` beside those FASTA after the suite.
- **Stats source duplicates**: `Duplicates` reports source-level extras (`sum(k-1)` over repeated names) on full scans. `--index-only` prints `n/a` on a deduplicated index (never a fabricated `0`); with `--no-dedup` it reports repeats retained in the index.
- **Validator/indexer agreement fixture**: `tests/data/validator_indexer_agreement.fasta` covers variable widths, trailing whitespace, blank lines, mixed CRLF/LF, missing terminal newline, empty records, and a 1025-byte header (near validate `--max-header-len`). Unit test checks validate issue kinds, indexed non-empty names, side tables, and mmap/stream index parity.
- **`index --low-mem` side-table prefix**: when a record stays formula-uniform for two or more lines and later becomes non-uniform, streaming now reconstructs the uniform prefix rows so `.zfi` side tables match mmap (previously dropped lines after the first).
- **Default mmap `index`**: writes `.zfi` via `scanZfiIndex()` with in-memory record and side-table arrays (replaces dummy-header streaming write). `scanFastaRecords` passes sequence data and uniform-width flag to the emit callback.
- **`index --low-mem`**: default output is `{file}.zfi` (was FAI-only in v0.2.x). Shares line-metrics semantics with mmap via `ChunkParseState` and `LineMetricsBuilder`. Removed duplicate `StreamingParseState` parser.
- **`.fai` fallback loading** (`index_format.zig`): respects `LoadMode` like `.zfi` (name hash map only for `lookup_full_map`). Single-pass mmap parse; record names live in the mmap'd `.fai` via `name_offset` / `name_len`. `LoadedIndex.getRecordName` resolves names for both index sources.
- **`.fai` `lookup_full_map` loads**: build a pointer hash over mmap'd `.fai` name fields (same fast path as embedded `.zfi` names; no arena copy).
- **`get` multi-region name lookup**: load `.lookup_full_map` for any N > 1 (and BED/names); N=1 stays `.records_only`. Removes the old 2..15 record-scan cliff on large catalogs.
- **`get --names` size bound**: `--names` always loads the whole file (up to 512 MiB). Oversize inputs fail with a names-specific error; `--chunk-size` does not stream `--names`. BED keeps chunked streaming by default and the same 512 MiB cap only for `--chunk-size -1`.
- **`get` dense BED RSS**: sorted/sequential mmap batches drop FASTA pages behind the scan cursor (`MADV_DONTNEED`, 8 MiB batching) so multi-thousand-row dense runs stay near index size.
- **`get` sparse BED paths**: skip file-order sort on small catalogs, wide median byte gaps, or FASTA under 64 MiB; emit in request order. Large sparse FASTAs (Genome-scale) open a GET-only `std.Io.File` and read via `readPositionalAll` so scattered rows do not retain mmap pages. `MADV_SEQUENTIAL` only when a true sequential scan is active; end-of-batch cache drop only on the mmap path for FASTA over 256 MiB.
- **`get` uniform emit**: `emitRegionForwardUniform` / `emitRegionBackwardUniform` for fixed-width records; reverse and RC share the same chunked path.
- **`get` protein-guard sample**: `detectRecordType` samples at most 256 bases per record (was up to 100k).
- **Sequence-type classification**: shared `stats.detectType` (overflow-safe 90% IUPAC threshold). Stats uses full-file composition; GET `--rc` / `--complement-only` samples 256 bases per record; validate samples up to 100000 bases and exposes `sequence_type` / `type_bases_sampled` / `type_sample_cap` in `--json --summary`.
- **`stats` composition scan**: walks records in file-offset order and releases FASTA pages behind the cursor with batched `MADV_DONTNEED` (8 MiB stride). Peak RSS drops from ~1:1 file-size residency to near index size (Genome 3006 MB to 242 MB; Transcriptome 474 MB to 24 MB). `.fai` records are sorted by `seq_offset` before release so both index formats share the same scan path. Fixed-width records use SIMD `countCompositionSlice` on base-only slices; side-table records use separate per-line counting. Wall time stays within +10% of the pre-release baseline.
- **`stats` nucleotide/protein detection**: `detectType` now checks the full composition counts instead of a 100K-byte sample.
- **`stats` `.fai` loading**: uses `.stats_scan` mode which streams record data without retaining the full `.fai` mapping or all name strings. Shortest and longest names are fetched on demand via recorded line offsets (two seeks per stats invocation), not by rescanning the sidecar from line zero.
- **`stats` assembly metrics output**: `formatSize` renders on-disk file size in human-readable units (bytes/KB/MB/GB).
- **Benchmark suite overhaul** (from v0.2.9 hyperfine): zebrac runners `bench/index/run.sh`, `bench/get/run.sh`, and `bench/stats/run.sh` (verify to perf to report) via shared `bench/shared/zebrac_runner.sh`; rewritten `generate_report.py` / `REPORT.md` for index, GET, and stats; GET verify scripts merged into `bench/get/verify.sh` (409 checks). Stats L2 verify (`bench/stats/verify.sh`, 92 checks) plus L3 full/mode/scale report.
- **`bench/get` RC figure**: plain z-fasta hatched like seqtk (ref); `--rc` / `--complement-only` / `--reverse-only` use distinct shades of brand gold.
- **`README.md`**: benchmark commands point at `bench/index/run.sh` and the updated verification workflow.

### Fixed

- **`.fai` emit contract**: `--emit-fai` refuses non-uniform layouts instead of writing a misleading text index.
- **`--low-mem` blank-line parity**: trailing blanks (before EOF or the next header) stay uniform like mmap body trim; interior blanks still force side tables. Fixes `MissingSideTable` scan failures on end-of-record blanks.
- **`get` on messy FASTAs**: correct extraction via side tables. Previously returned garbage or failed silently on mixed-width files.
- **`.fai` on single-region GET**: no longer always builds a full name hash map (fixed large slowdown vs `.zfi` on transcriptome-scale indexes).
- **`index --low-mem`**: skip throwaway side-table work on uniform records (fixes RSS and instruction blow-up on large genomes).
- **`get` RC / complement wall**: `--rc` and `--complement-only` no longer pay a large protein-guard sample or per-byte complement on the hot path (Genome `--rc`/plain ~1.0x after fix).
- **`bench/get/run.sh`**: positional perf uses only positional fixture labels (not multi-region lists); region span metadata parsing matches `parseRegion`.

### Removed

- **Legacy benchmark tooling** (replaced by zebrac overhaul): `bench/wrappers/`, hyperfine `run_benchmarks.sh` runners, `bench/index/run_tests.sh`, split GET verify scripts, baseline snapshot scripts; stats runner/report (`bench/stats/verify_stats.py`). `bench/baselines/BASELINES.md` keeps the v2 schema for when snapshots return.
- `bench/index/results/figures/scaling_seqs.png` and `speedup.png` (legacy figures).

## [0.2.9] - 2026-06-24

Various memory safety, optimizations and re-building benchmarking framework to work with less dependencies. Also some code cleanup, standardization.

### Changed

- `src/main.zig`: write `.zfi` indexes via a temp file and rename on success
- `src/getter.zig`: cap multi-region sort buffers at 64 MiB per region and 256 MiB total; reject BED inputs over 512 MiB when `--chunk-size -1` (`--names` always uses the same size cap)
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
- **Messy FASTA benchmark** (`bench/shared/messy_perf/`): four derived FASTA variants (mixed_widths, crlf_endings, trailing_whitespace, all_messy) benchmarked against all indexing tools; `compatibility.csv` with per-tool exit-code results
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
