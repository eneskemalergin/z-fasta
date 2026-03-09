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
  get/                     # Getter benchmarks + verification (v0.2)
    run_benchmarks.sh      # Hyperfine benchmarks: latency, full-seq, region scaling
    verify_get.sh          # Byte-identical diff against samtools faidx
    generate_report.py     # Produces REPORT.md + PNG figures
    REPORT.md              # Full GET benchmark report
    results/               # Timestamped JSON, memory CSV, figures/
  stats/                   # Stats benchmarks + verification (v0.2)
    run_benchmarks.sh      # Hyperfine benchmarks: full/index-only, scaling, throughput
    verify_stats.py        # BioPython verification of all stats output
    generate_report.py     # Produces REPORT.md + PNG figures
    REPORT.md              # Full STATS benchmark report
    results/               # Timestamped JSON, memory/throughput CSV, figures/
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

# 5. Run GET benchmarks (includes real datasets by default)
bash bench/get/run_benchmarks.sh
# Add --skip-real to skip real dataset benchmarks (faster, ~5 min instead of ~15 min):
# bash bench/get/run_benchmarks.sh --skip-real

# 6. Run STATS benchmarks (includes real datasets by default)
bash bench/stats/run_benchmarks.sh
# Add --skip-real to skip real dataset benchmarks:
# bash bench/stats/run_benchmarks.sh --skip-real

# 7. Verify get command against samtools
bash bench/get/verify_get.sh

# 8. Verify stats against BioPython
.venv/bin/python bench/stats/verify_stats.py

# 9. Generate all reports
python3 bench/index/generate_report.py
python3 bench/get/generate_report.py
python3 bench/stats/generate_report.py
```

## Output

- `index/REPORT.md` — Full indexer benchmark report with tables and figures
- `get/REPORT.md` — Full GET benchmark report (latency, scaling, memory)
- `stats/REPORT.md` — Full STATS benchmark report (throughput, index-only speedup)
- `*/results/figures/` — PNG charts
- `*/results/*.csv`, `*/results/*/` — Raw data (gitignored, regenerated)
- `shared/data/` — Real datasets (gitignored, downloaded by script)
