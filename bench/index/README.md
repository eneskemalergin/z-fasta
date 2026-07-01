# Index Benchmarks

This suite measures `z-fasta index` against FASTA indexers that produce or use FAI-style indexes: `samtools faidx`, `seqkit faidx`, `fastahack -i`, `pyfaidx`, and the Tier 2 Rust wrappers for `noodles-fasta` and `rust-bio`.

## What It Covers

- Real datasets from [bench/shared/data](../shared/data): Genome, Transcriptome, and Proteome.
- Synthetic file-size scaling from 1 MB through 1000 MB.
- Two sequence-count sweeps: bounded ~50 MiB total (1k–250k records) and fixed 1024 bp per record (100k–1M).
- `zebrac` wall time, peak RSS, page faults, and hardware counters (single measurement source).
- Edge-case correctness against `samtools` (fixtures in `edge_cases/`).
- Messy FASTA: quick correctness checks (`messy_variants/`) and zebrac indexing on proteome-derived fixtures in [bench/shared/messy_variants](../shared/messy_variants/).

## Run

From the repository root:

```bash
./zig build -Doptimize=ReleaseFast
bash bench/index/run.sh
```

`run.sh` is the single entry point: correctness tests, zebrac benchmarks, messy zebrac, and `generate_report.py`.

### Options

```bash
bash bench/index/run.sh --runs 10 --warmup 2
bash bench/index/run.sh --skip-real          # scaling only (keeps manifest sections you skip)
bash bench/index/run.sh --skip-scaling       # real datasets only
bash bench/index/run.sh --skip-tests --skip-messy --skip-report   # benchmarks only
bash bench/index/run.sh --scaling-only --merge-base 20260630_231053 --runs 5 --warmup 2
bash bench/index/run.sh --allow-incomplete # local report drafts from partial data
```

`generate_report.py` refuses to overwrite the tracked report from smoke or incomplete runs unless `--allow-incomplete` is passed (via `run.sh` or directly).

## Outputs

- [REPORT.md](REPORT.md): generated summary tables and chart links.
- `results/perf_<timestamp>/`: real-dataset zebrac JSON.
- `results/scale_size_<timestamp>/`: synthetic file-size scaling zebrac JSON.
- `results/scale_seqs_budget_<timestamp>/`: bounded-bytes sequence scaling zebrac JSON.
- `results/scale_seqs_fixed_<timestamp>/`: fixed-length sequence scaling zebrac JSON.
- `results/run_<timestamp>.json`: run manifest with runner settings, tool versions, and section paths.
- `results/LATEST`: timestamp pointer for the most recent run bundle.
- `results/metadata_<timestamp>.jsonl`: suite/workload/tool metadata for the raw zebrac JSON.
- `results/tests_<timestamp>.csv`: edge-case correctness results.
- `results/messy_<timestamp>/`: zebrac messy-variant indexing.
- `results/figures/`: regenerated PNG charts used by the report.

## Notes

The zebrac runner uses warm-cache mode. It does not currently have a per-sample prepare hook, so index-file cleanup is included inside each measured command. That overhead is small and applies to every tool in the suite.

`z-fasta-default` uses mmap and duplicate-name detection. `z-fasta-nodedup` disables duplicate checking. `z-fasta-lowmem` uses streaming IO to trade speed for low peak RSS.

Peak RSS and page faults in the report come from zebrac on the same samples as wall time. There is no separate `/usr/bin/time` pass.

`rustbio-custom-index` is listed as **rust-bio** in reports. We use `tools/rustbio_wrapper` because rust-bio has no standalone index CLI.
