# Index Benchmarks

This suite measures `z-fasta index` against FASTA indexers that produce or use FAI-style indexes: `samtools faidx`, `seqkit faidx`, `fastahack -i`, and `pyfaidx`.

## What It Covers

- Real datasets from [bench/shared/data](../shared/data): Genome, Transcriptome, and Proteome.
- Synthetic file-size scaling from 1 MB through 1000 MB.
- Synthetic sequence-count scaling from 10 through 100,000 sequences.
- Peak RSS and page-fault measurements on the largest real dataset.
- Edge-case correctness against `samtools` via [run_tests.sh](run_tests.sh).
- Messy FASTA compatibility through generated variants in [bench/shared/messy_variants](../shared/messy_variants).

## Run

From the repository root:

```bash
./zig-0.16.0/zig build -Doptimize=ReleaseFast
bash bench/index/run_benchmarks.sh
bash bench/index/run_tests.sh
.venv/bin/python bench/index/generate_report.py
```

Use `--skip-real` to skip downloaded real datasets, or `--skip-scaling` to run only the real-dataset section.

## Outputs

- [REPORT.md](REPORT.md): generated summary tables and chart links.
- `results/perf_<timestamp>/`: real-dataset hyperfine JSON.
- `results/scale_size_<timestamp>/`: synthetic file-size scaling JSON.
- `results/scale_seqs_<timestamp>/`: synthetic sequence-count scaling JSON.
- `results/memory_<timestamp>.csv`: `/usr/bin/time` RSS and page-fault measurements.
- `results/figures/`: regenerated PNG charts used by the report.

## Notes

The runner uses warm-cache mode unless passwordless sudo is available for dropping Linux page cache. The current v0.2.6 run was collected in warm-cache mode.

`z-fasta-default` uses mmap and duplicate-name detection. `z-fasta-nodedup` disables duplicate checking. `z-fasta-lowmem` uses streaming IO to trade speed for low peak RSS.
