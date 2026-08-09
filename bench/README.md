# Benchmarking Framework

This is the benchmarking framework for z-fasta. It is used to benchmark the performance of z-fasta and to compare it with other FASTA tools, across its functionality. It measures correctness, timings, peak RSS usage, page faults, and more via the zebrac benchmarking framework.

Note: [Zebrac](https://github.com/eneskemalergin/zebrac) is a linux-only tool that is similar to `hyperfine`, but also measures RSS, page faults, and more. For performance benches, it is the only tool used here. I added the binary under `tools/`; it is not a dependency of z-fasta, and it does not work on platforms other than Linux for now.

## Published reports

- [Index benchmark report](index/REPORT.md)
- [GET benchmark report](get/REPORT.md)
- [Stats benchmark report](stats/REPORT.md)

Each report owns its module's methods, field coverage, correctness checks, measurements, and figures.

## How to run

Each suite has one entrypoint: `bash bench/<index|get|stats>/run.sh`. That script runs correctness first, then optional zebrac perf, then writes `REPORT.md`. For suite-specific flags, pass `--help` to that script.

```bash
bash bench/shared/download_data.sh   # ~4 GB REAL_* datasets, once
bash bench/shared/install_tools.sh   # peer tools + zebrac
zig build -Doptimize=ReleaseFast

# Full suite (correctness + perf + report)
bash bench/index/run.sh
bash bench/get/run.sh
bash bench/stats/run.sh

# Correctness only
bash bench/index/run.sh --skip-benchmarks --skip-messy --skip-report
bash bench/get/run.sh --skip-benchmarks --skip-report
bash bench/stats/run.sh --skip-benchmarks --skip-report
```

Shared skips across suites (old names still work as aliases for one release cycle):

- `--skip-tests` (alias `--skip-verify`): skip correctness
- `--skip-benchmarks` (alias `--skip-perf`): skip zebrac / perf
- `--skip-report`: skip report generation
- `--skip-messy`: skip messy *perf* only. It never skips messy cases inside correctness.

Helpers live in `bench/shared/` (`tools.sh`, `zebrac_runner.sh`, `runner_common.sh`, `download_data.sh`, `install_tools.sh`, `generate_messy.py`, `generate_scaling.py`). Generated fixtures under `bench/*/data/` and `bench/shared/cache/` are gitignored (messy fixtures/perf, scaling FASTAs). Materialize with `python3 bench/shared/generate_messy.py` and `python3 bench/shared/generate_scaling.py` (`--force` after param edits). Old suite-local `bench/{index,stats}/data/size_*.fasta` / `seqs_*.fasta` trees are obsolete; safe to delete if present. What we keep in git for GitHub is each suite's `REPORT.md` plus `results/figures/*.png`. After a full regen, commit those together so the report images still render.

Get messy perf alone is heavy (lots of samples, long walls). Use `--skip-*` while you are iterating. Run the full suites before a release tag.
