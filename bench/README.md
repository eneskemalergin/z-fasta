# Benchmarking z-fasta

This directory contains all benchmarking and verification infrastructure, organized per subcommand.

## Directory Structure

```plaintext
bench/
  shared/                  # Shared infrastructure
    download_data.sh       # Downloads REAL_* files -> bench/shared/data/
    data/                  # REAL_Genome.fa, REAL_Transcriptome.fa, REAL_Proteome.fasta
  index/                   # Indexer benchmarks
    README.md              # Suite-specific workflow and outputs
    run_benchmarks.sh      # Hyperfine benchmarks across all tools and file sizes
    run_tests.sh           # 20 edge case correctness tests vs samtools
    generate_report.py     # Produces REPORT.md + PNG figures
    REPORT.md              # Full benchmark report
    data/                  # Generated scaling files (size_*.fasta, seqs_*.fasta)
    edge_cases/            # Generated edge-case FASTAs
    results/               # Timestamped perf JSON, memory CSV, figures/
  get/                     # Getter benchmarks + verification
    README.md              # Suite-specific workflow and outputs
    run_benchmarks.sh      # Hyperfine benchmarks: latency, full-seq, region scaling
    verify_get.sh          # Byte-identical diff against samtools faidx
    generate_report.py     # Produces REPORT.md + PNG figures
    REPORT.md              # Full GET benchmark report
    results/               # Timestamped JSON, memory CSV, figures/
  stats/                   # Stats benchmarks + verification
    README.md              # Suite-specific workflow and outputs
    run_benchmarks.sh      # Hyperfine benchmarks: full/index-only, scaling, throughput
    verify_stats.py        # BioPython verification of all stats output
    generate_report.py     # Produces REPORT.md + PNG figures
    REPORT.md              # Full STATS benchmark report
    results/               # Timestamped JSON, memory/throughput CSV, figures/
  wrappers/                # Tier 2 wrapper benchmarks
    bench_all.sh           # Hyperfine benchmarks for noodles-fasta, rust-bio
    sanity_check.sh        # Correctness verification against samtools
    results/               # Timestamped JSON output
```

## Tier 2 Rust Wrappers

The `bench/wrappers/` directory contains benchmarks for Rust FASTA libraries that don't have standalone CLI tools. These wrappers are built in `tools/`:

| Wrapper | Library | Indexer | Notes |
| --- | --- | --- | --- |
| `noodles_wrapper` | noodles-fasta 0.61 | Built-in (`noodles_fasta::fs::index()`) | Direct comparison to samtools |
| `rustbio_wrapper` | rust-bio 2.3 | Custom (no built-in indexer) | Must be labeled as "custom indexer" in papers |

Run Tier 2 benchmarks:

```bash
# Build wrappers first
cd tools/noodles_wrapper && cargo build --release && cd ../..
cd tools/rustbio_wrapper && cargo build --release && cd ../..

# Run cross-tool comparison
bash bench/wrappers/bench_all.sh --runs 5
```

## Prerequisites

The benchmark suite has its own dependencies, separate from z-fasta itself. **z-fasta has zero runtime dependencies**, but reproducing the benchmarks requires several tools.

Run the verification / install helper first:

```bash
bash bench/shared/install_tools.sh
```

This prints a version-pins table and installs `pyfaidx` into the project `.venv` if it is missing. The other tools (seqtk, fastahack, seqkit) ship as pre-built binaries in `tools/` and are checked but not rebuilt.

### Required

| Tool | Pinned version | Purpose | Install |
| --- | --- | --- | --- |
| [hyperfine](https://github.com/sharkdp/hyperfine) | 1.12.0 | Precise timing with warmup/cache control | `apt install hyperfine` or `cargo install hyperfine` |
| [samtools](http://www.htslib.org/) | 1.13 | Reference indexer (`samtools faidx`) | `apt install samtools` or `conda install samtools` |
| z-fasta | current | The tool being benchmarked | `./zig build -Doptimize=ReleaseFast` (from project root) |

### Comparison tools (pre-built in `tools/`)

| Tool | Pinned version | Purpose | Binary |
| --- | --- | --- | --- |
| [seqkit](https://bioinf.shenwei.me/seqkit/) | v2.13.0 | Go-based FASTA toolkit | `tools/seqkit` |
| [fastahack](https://github.com/ekg/fastahack) | 1.0.0 | C++ FASTA indexer | `tools/fastahack-1.0.0/fastahack` |
| [seqtk](https://github.com/lh3/seqtk) | 1.5-r133 | C FASTA/FASTQ toolkit | `tools/seqtk/seqtk` |
| [noodles-fasta](https://github.com/zaeleus/noodles) | 0.61 | Rust FASTA library (Tier 2) | `tools/noodles_wrapper/target/release/noodles_wrapper` |
| [rust-bio](https://github.com/rust-bio/rust-bio) | 2.3 | Rust bioinformatics library (Tier 2) | `tools/rustbio_wrapper/target/release/rustbio_wrapper` |

If a comparison tool binary is not found, the benchmark script skips it automatically.

> **Note on fastahack version:** fastahack does not expose a `--version` flag; the pinned version is tracked via the directory name `tools/fastahack-1.0.0/`.

### Python tools

| Package | Pinned version | Purpose | Install |
| --- | --- | --- | --- |
| [pyfaidx](https://github.com/mdshw5/pyfaidx) | 0.9.0.3 | Python FASTA random-access (`faidx` CLI) | `pip install pyfaidx==0.9.0.3` |
| matplotlib | latest | Figures in benchmark reports | `pip install matplotlib` |
| pandas | latest | Data processing | `pip install pandas` |
| tabulate | latest | Markdown tables | `pip install tabulate` |
| biopython | latest | Stats verification | `pip install biopython` |

`pyfaidx` is handled automatically by `install_tools.sh`. All other packages can be installed together:

```bash
pip install matplotlib pandas tabulate biopython
```

### Optional: Cold-cache benchmarks

For cold-cache measurements, the benchmark script attempts `sudo sh -c 'echo 3 > /proc/sys/vm/drop_caches'`. This requires passwordless sudo for that specific command. Without it, benchmarks run in warm-cache mode (still valid, just different).

## Quick Start

```bash
# 1. Build z-fasta from the current checkout with the repo-local Zig wrapper (from project root)
./zig build -Doptimize=ReleaseFast

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
.venv/bin/python bench/index/generate_report.py
.venv/bin/python bench/get/generate_report.py
.venv/bin/python bench/stats/generate_report.py
```

For suite-level details, see [index/README.md](index/README.md), [get/README.md](get/README.md), and [stats/README.md](stats/README.md).

## Output

- `index/REPORT.md`: Full indexer benchmark report with tables and figures
- `get/REPORT.md`: Full GET benchmark report (latency, scaling, memory)
- `stats/REPORT.md`: Full STATS benchmark report (throughput, index-only speedup)
- `*/results/figures/`: PNG charts
- `*/results/*.csv`, `*/results/*/`: Raw data (gitignored, regenerated)
- `shared/data/`: Real datasets (gitignored, downloaded by script)

## Performance baselines (local, gitignored)

The baseline scripts are an **additive layer** on top of the existing bench flow. They do not replace `generate_report.py`, `REPORT.md`, or the raw `results/` trees.

| Layer | Role | Output |
| --- | --- | --- |
| `run_benchmarks.sh` (per suite) | Produce hyperfine JSON + memory CSV | `bench/*/results/` |
| `generate_report.py` | Human-readable tables and figures | `REPORT.md`, `results/figures/` |
| `save_baseline.py` | Condensed snapshot for regression diffing | `bench/baselines/<stamp>_<rev>/` |
| `compare_baseline.py` | Diff two snapshots | stdout; exit 1 on regressions |

`compare.json` holds **timing only** (flat hyperfine means). Memory RSS and page faults stay in the per-suite `memory_*.csv` files and `REPORT.md`; they are copied into `index.json` / `get.json` / `stats.json` but are not part of the regression gate yet.

**Use one coherent bench run per snapshot.** `save_baseline.py` groups artifacts by bundle timestamp (same logic as `generate_report.py` for GET). If you run index today and get yesterday, the snapshot mixes runs and prints `bundle_warnings`. Prefer:

```bash
bash bench/run_all_and_baseline.sh
```

Faster iteration while developing:

```bash
bash bench/run_all_and_baseline.sh --skip-scaling --skip-rc --label my-run
```

Or step by step after a full suite pass:

```bash
.venv/bin/python bench/save_baseline.py --label optional-note
```

Snapshots land in `bench/baselines/<YYYYMMDD_HHMMSS>_<git-rev>/` (gitignored):

| File | Contents |
| --- | --- |
| `manifest.json` | Git rev, host/CPU, binary version, bundle timestamps, headline timings |
| `index.json` / `get.json` / `stats.json` | Condensed per-suite metrics |
| `compare.json` | Flat `suite.section.param.tool -> seconds` map for diffing |

Compare two snapshots:

```bash
.venv/bin/python bench/compare_baseline.py bench/baselines/<old> bench/baselines/<new> --z-fasta-only
```

`bench/baselines/LATEST` points at the most recent snapshot on this machine.

**Release gate:** run the full suite before tagging; compare against `LATEST` or the previous release baseline. Investigate any z-fasta regression above ~5–10%.
