# STATS Benchmarks

This suite measures `z-fasta stats` full-scan and `--index-only` modes against tools with overlapping FASTA statistics behavior: `seqkit stats -a` and `seqtk comp`.

## What It Covers

- Small checked-in test FASTAs from [tests/data](../../tests/data).
- Index-only versus full-scan timing on 10 MB, 50 MB, and 100 MB synthetic files.
- Synthetic file-size scaling from 1 MB through 1000 MB.
- Real dataset stats for Genome, Transcriptome, and Proteome.
- Peak RSS/page-fault measurements for full-scan and index-only modes.
- Throughput calculations for full-scan tools.
- Output verification against BioPython via [verify_stats.py](verify_stats.py).

## Run

From the repository root:

```bash
./zig build -Doptimize=ReleaseFast
bash bench/stats/run_benchmarks.sh
.venv/bin/python bench/stats/verify_stats.py
.venv/bin/python bench/stats/generate_report.py
```

Use `--skip-real` to avoid downloaded real datasets, or `--skip-scaling` to omit synthetic file-size scaling.

## Outputs

- [REPORT.md](REPORT.md): generated summary tables and chart links.
- `results/stats_<timestamp>/`: test-file hyperfine JSON.
- `results/indexonly_<timestamp>/`: index-only versus full-scan JSON.
- `results/scale_size_<timestamp>/`: synthetic file-size scaling JSON.
- `results/real_<timestamp>/`: real-dataset stats JSON.
- `results/memory_<timestamp>.csv`: `/usr/bin/time` RSS and page-fault measurements.
- `results/throughput_<timestamp>.csv`: full-scan throughput calculations.
- `results/figures/`: regenerated PNG charts used by the report.

## Notes

The current v0.2.6 run was collected in warm-cache mode because passwordless sudo was not available for page-cache drops.

`z-fasta stats --index-only` reads `.zfi` index data and computes length-derived metrics without scanning FASTA sequence bytes. It is intended for quick assembly QC. It has no direct equivalent in `seqkit` or `seqtk`, so full-scan comparisons should be read separately from index-only comparisons.
