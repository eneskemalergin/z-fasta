<!-- markdownlint-disable MD024 -->
# Changelog

All notable changes to z-fasta will be documented in this file.

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
