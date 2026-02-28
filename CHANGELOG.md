# Changelog

All notable changes to z-fasta will be documented in this file.

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
