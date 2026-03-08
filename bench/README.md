# Benchmarking z-fasta

This directory contains all benchmarking and verification infrastructure, organized per subcommand.

## Directory Structure

```plaintext
bench/
  shared/                  # Shared infrastructure
    download_data.sh       # Downloads REAL_* files → bench/shared/data/
    data/                  # REAL_Genome.fa, REAL_Transcriptome.fa, REAL_Proteome.fasta
  index/                   # Indexer benchmarks (v0.1)
    run_benchmarks.sh      # Hyperfine benchmarks across all tools and file sizes
    run_tests.sh           # 20 edge case correctness tests vs samtools
    generate_report.py     # Produces REPORT.md + PNG figures
    REPORT.md              # Full benchmark report
    data/                  # Generated scaling files (size_*.fasta, seqs_*.fasta)
    edge_cases/            # Generated edge-case FASTAs
    results/               # Timestamped perf JSON, memory CSV, figures/
  get/                     # Getter verification (v0.2)
    verify_get.sh          # Byte-identical diff against samtools faidx
    results/               # Future hyperfine results
  stats/                   # Stats verification (v0.2)
    verify_stats.py        # BioPython verification of all stats output
    results/               # Future hyperfine results
```

## Prerequisites

The benchmark suite has its own dependencies, separate from z-fasta itself. **z-fasta has zero runtime dependencies** — but reproducing the benchmarks requires several tools.

### Required

| Tool | Purpose | Install |
| --- | --- | --- |
| [hyperfine](https://github.com/sharkdp/hyperfine) | Precise timing with warmup/cache control | `apt install hyperfine` or `cargo install hyperfine` |
| [samtools](http://www.htslib.org/) | Reference indexer (`samtools faidx`) | `apt install samtools` or `conda install samtools` |
| z-fasta | The tool being benchmarked | `zig build -Doptimize=ReleaseFast` (from project root) |

### Optional (for comparison)

| Tool | Purpose | Install |
| --- | --- | --- |
| [seqkit](https://bioinf.shenwei.me/seqkit/) | Go-based FASTA toolkit (`seqkit faidx`) | Download binary to `tools/seqkit` |
| [fastahack](https://github.com/ekg/fastahack) | C++ FASTA indexer | Build and place in `tools/fastahack-1.0.0/` |

If seqkit or fastahack are not found, benchmarks will skip them automatically.

### Report generation (Python)

Generating the Markdown report and figures requires Python 3.10+ with:

```bash
pip install matplotlib pandas tabulate
```

### Stats verification (Python)

```bash
pip install biopython
```

### Optional: Cold-cache benchmarkss

For cold-cache measurements, the benchmark script attempts `sudo sh -c 'echo 3 > /proc/sys/vm/drop_caches'`. This requires passwordless sudo for that specific command. Without it, benchmarks run in warm-cache mode (still valid, just different).

## Quick Start

```bash
# 1. Build z-fasta (from project root)
zig build -Doptimize=ReleaseFast

# 2. Download test data (~4 GB, one-time)
bash bench/shared/download_data.sh

# 3. Run indexer benchmarks
bash bench/index/run_benchmarks.sh

# 4. Run indexer correctness tests
bash bench/index/run_tests.sh

# 5. Verify get command against samtools
bash bench/get/verify_get.sh

# 6. Verify stats against BioPython
.venv/bin/python bench/stats/verify_stats.py

# 7. Generate indexer benchmark report
python3 bench/index/generate_report.py
```

## Output

- `index/REPORT.md` — Full indexer benchmark report with tables and figures (committed to repo)
- `index/results/figures/` — PNG charts (committed to repo)
- `index/results/*.csv`, `index/results/perf_*/` — Raw data (gitignored, regenerated)
- `shared/data/` — Real datasets (gitignored, downloaded by script)
