# GET Benchmarks

This suite measures `z-fasta get` region extraction latency and throughput against `samtools faidx`, `seqkit faidx`, `fastahack`, `seqtk subseq`, and `pyfaidx`.

## What It Covers

- Single-region extraction from synthetic files at 100 bp and 10 kbp.
- Full-sequence extraction from 1 MB, 10 MB, 50 MB, and 100 MB files.
- Region-size scaling from 100 bp through 1 Mbp.
- Real dataset extraction for Genome, Transcriptome, and Proteome.
- BED batch extraction for 100, 1K, 10K, and 100K regions in both default and strand-aware modes.
- Peak RSS/page-fault measurements for small, large, and full-sequence regions.
- Multi-region extraction in one CLI call for 1, 10, 50, and 100 regions.
- Byte-for-byte verification against `samtools faidx` via [verify_get.sh](verify_get.sh) and [verify_multi_get.sh](verify_multi_get.sh).
- BED extraction verification against `bedtools getfasta` plus default-mode `samtools faidx -r` comparisons via [verify_bed.sh](verify_bed.sh), including chunked `--bed`, stdin BED, and `--names` runs.
- Reverse / complement verification via [verify_rc.sh](verify_rc.sh), including `--rc` against `samtools faidx -i --mark-strand no`, transform annotation checks, BED `--strand-aware --rc` composition, protein rejection, and a synthetic chromosome-like full-sequence fixture.
- Focused RC lock-in benchmarking is executed inline inside [run_benchmarks.sh](run_benchmarks.sh) by default so the main GET report includes RC timing and RSS coverage without a second benchmark entrypoint.
- A checked-in reverse-path summary is captured in [RC_STRATEGY.md](RC_STRATEGY.md).

## Run

From the repository root:

```bash
./zig build -Doptimize=ReleaseFast
bash bench/get/run_benchmarks.sh
bash bench/get/verify_get.sh
bash bench/get/verify_multi_get.sh
bash bench/get/verify_bed.sh
bash bench/get/verify_rc.sh
.venv/bin/python bench/get/generate_report.py
```

Use `--skip-real` to avoid downloaded real datasets, or `--skip-scaling` to omit region-size scaling.
`run_benchmarks.sh` runs the RC orientation review inline by default and writes the report plus figures to the standard tracked GET benchmark surface. Use `--skip-rc` if you want to exclude it, or `--rc-full` if you want the heavier RC lock-in profile instead of the quick default.
`verify_bed.sh` uses `tests/data/simple.fasta` by default and generates synthetic small / medium / large / x-large BED files internally. It compares default BED output against both `bedtools getfasta` and `samtools faidx -r`, checks stranded BED output against `bedtools -s`, and covers `--names` with comments and duplicate names. It skips cleanly when `bedtools` or `samtools` is not installed.
`verify_rc.sh` stays correctness-only. It uses `samtools faidx -i --mark-strand no` as the gold reference for `--rc` when available, and adds exact-output checks for `--annotate-rc`, `--reverse-only`, `--complement-only`, BED composition, the protein error path, and a wrapped full-sequence chromosome-like fixture.

## Outputs

- [REPORT.md](REPORT.md): generated summary tables and chart links, including RC timing and RSS comparisons against `samtools` and `bedtools + seqtk` where applicable.
- [RC_STRATEGY.md](RC_STRATEGY.md): shipped reverse-path note and current measurement slice.
- `results/rc_review_<timestamp>/`: inline RC review tables, JSON, and RSS snapshots produced by `run_benchmarks.sh`.
- `results/single_<timestamp>/`: single-region hyperfine JSON.
- `results/fullseq_<timestamp>/`: full-sequence hyperfine JSON.
- `results/scale_region_<timestamp>/`: region-size scaling JSON.
- `results/real_<timestamp>/`: real-dataset extraction JSON.
- `results/bed_<timestamp>/`: BED batch extraction JSON for default and stranded modes.
- `results/multi_<timestamp>/`: multi-region extraction JSON.
- `results/memory_<timestamp>.csv`: `/usr/bin/time` RSS and page-fault measurements.
- `results/figures/`: regenerated PNG charts used by the report.

## Notes

The checked-in run was collected in warm-cache mode because passwordless sudo was not available for page-cache drops.

For indexed tools, small-region latency is dominated by process startup plus index lookup. `seqtk` has no index and scans the FASTA, so it is included as a scan-based reference rather than a direct indexed-access equivalent.
