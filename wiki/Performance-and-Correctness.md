# Performance and Correctness

I care about correct output before timing. The benchmark harness keeps correctness checks separate from measured runs, and the Wiki shows only selected results from the owning reports.

## Correctness coverage

### Index

- Samtools-compatible `.fai` comparison on representable edge cases.
- `.zfi` format golden bytes and corruption mutations.
- Read-boundary, identifier-boundary, geometry, staleness, duplicate, empty-name, and cleanup tests.
- Messy FASTA fixtures for side-table output.

### GET

- 418 correctness checks in the current development gate.
- Positional, multi-region, names, BED, stdin, strand, transform, summary, and messy-layout cases.
- Peer comparison with samtools, bedtools, and seqtk where semantics overlap.
- Exact CLI status, stdout, stderr, and valid-prefix behavior.

### Stats

- Independent exact composition oracle.
- `.zfi` and `.fai` parity.
- Layout twins and messy records.
- Exact nucleotide and protein report tests.
- Peer comparison where field sets overlap.

### Validate

- Every currently defined issue family.
- JSON escaping, including invalid UTF-8.
- Fix idempotence and output-path safety.
- Event-cap behavior and late blockers.
- Biological retrieval preservation after repair.
- Validator and indexer agreement fixtures.

## Why the hot paths are bounded

- Index reads sequence payload through a fixed 1 MiB buffer. Catalog memory grows with records, identifiers, deduplication, and side tables, not FASTA payload size.
- GET uses positional reads and fixed input/output buffers. Names and BED requests stream through bounded reusable workspaces.
- Stats uses one 256 KiB descriptor-backed read window over indexed sequence spans.
- Validate currently maps the complete FASTA, so its peak RSS can approach input size.

## Representative published results

The repository reports use an AMD Ryzen 9 3950X with warm cache. These tables are selected snapshots from those reports, not a promise for every machine or dataset.

Index `.fai` examples:

| Dataset | Size | z-fasta | samtools |
| --- | ---: | ---: | ---: |
| Human genome | ~2.9 GiB | 0.3792 s | 8.0870 s |
| Transcriptome | ~459 MiB | 0.2171 s | 1.6615 s |
| Proteome | ~13 MiB | 0.0106 s | 0.0577 s |

GET 1 kbp examples through `.zfi`:

| Dataset | z-fasta | samtools |
| --- | ---: | ---: |
| Human genome | 2.1 ms | 7.1 ms |
| Transcriptome | 5.0 ms | 315.8 ms |
| Proteome | 2.4 ms | 16.6 ms |

<details>
<summary>View the full positional GET chart</summary>

<p align="center"><img src="https://github.com/eneskemalergin/z-fasta/blob/main/bench/get/results/figures/perf_pos_wall.png?raw=1" alt="Positional GET wall time across genome, transcriptome, and proteome datasets at four region sizes" width="100%"></p>

The chart includes additional indexed peers, region sizes from 100 bp through a complete record, and repeated-run error bars. Read the labels as peer time divided by z-fasta `.zfi` time.

</details>

Complete stats examples through `.zfi`:

| Dataset | z-fasta | SeqKit `stats -a` |
| --- | ---: | ---: |
| Human genome | 2.708 s | 17.626 s |
| Transcriptome | 0.382 s | 2.409 s |
| Proteome | 12.2 ms | 58.7 ms |

<details>
<summary>View the full stats chart</summary>

<p align="center"><img src="https://github.com/eneskemalergin/z-fasta/blob/main/bench/stats/results/figures/perf_full.png?raw=1" alt="Full stats wall time, peak memory, and minor page faults across genome, transcriptome, and proteome datasets" width="100%"></p>

The stats figure shows time, peak RSS, and minor page faults together. Hatched peer lanes cover smaller field sets, so they remain context rather than equal-work comparisons.

</details>

> [!IMPORTANT]
> SeqKit reports a smaller assembly and composition field set in this comparison. Its timing is an ecosystem reference, not equal-work equivalence.

## Interpret small-region latency

Below roughly 10 kbp on the benchmark host, end-to-end GET is dominated by process startup, sidecar resolution, and output setup. The byte-offset extraction itself is not the only measured cost.

For many requests, one invocation normally beats one process per region because it loads the index once and can share reads across neighboring spans.

## Reports and methodology

The repository reports own the complete tables, peer definitions, dataset details, verification coverage, and measurement notes:

- [Index report](https://github.com/eneskemalergin/z-fasta/blob/main/bench/index/REPORT.md)
- [GET report](https://github.com/eneskemalergin/z-fasta/blob/main/bench/get/REPORT.md)
- [Stats report](https://github.com/eneskemalergin/z-fasta/blob/main/bench/stats/REPORT.md)

The Wiki intentionally keeps only representative tables and figures. The benchmark harness and its external tool setup will continue to evolve, so these pages do not promise a turnkey reproduction environment.

Validate does not yet have a published performance report.

## Related pages

- [Choosing a FASTA tool](Choosing-a-FASTA-Tool)
- [Index formats](Index-Formats)
- [For contributors](For-Contributors)
