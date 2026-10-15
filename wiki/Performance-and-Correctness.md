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

- Every finding code, with hand-counted occurrences, records, and first locations.
- Identical reports across read-buffer sizes, and vector checks that agree with the per-byte path.
- Reversible identifier bytes in JSON.
- Repair idempotence, refusals, and output-path safety.
- Biological retrieval preservation after repair.
- Validator and indexer agreement fixtures.
- Index comparison: every `pair.*` code from same-size rewrites, through `.zfi` and `.fai`.

## Why the hot paths are bounded

- Index reads sequence payload through a fixed 256 KiB buffer. Inputs above 32 MiB use compact name references and stream `.zfi` records to temporary output. Name, name-length, and side-table staging spill to disk after 1 MiB each. The `.zfi` name table is built in memory at the end, about 5 bytes per record. A non-uniform record's line table remains in memory while that record is parsed.
- GET uses positional reads and fixed input/output buffers. Names and BED requests stream through bounded reusable workspaces.
- GET resolves `.zfi` identifiers through the index's name table and checks only the records its requests use. With a `.fai`, positional requests and seekable names and BED files read the sidecar only until every requested name is found and check only the lines they keep; stdin and named FIFOs map every name and parse a line only when a request uses it.
- Stats reads the index in file order through fixed buffers and keeps 4 bytes per record for its length (8 once a record is longer than 4,294,967,295 symbols). It reads sequence spans through one 256 KiB window; an index that lists records out of file order is loaded whole for that scan.
- Validate reads through one 256 KiB buffer. Its duplicate check keeps a fingerprint and a file offset per distinct identifier and re-reads an earlier identifier only on a fingerprint match. With `--index` it reads the index in file order through fixed buffers.

## Representative published results

The selected v0.4.0 reports use an AMD Ryzen 9 3950X with warm cache, measured on 2026-10-15. These tables are snapshots from those reports; times and memory use depend on the machine and dataset.

Index `.fai` examples:

| Dataset       |     Size |  z-fasta | samtools |
| ------------- | -------: | -------: | -------: |
| Human genome  | ~2.9 GiB | 0.3279 s | 8.1025 s |
| Transcriptome | ~459 MiB | 0.1252 s | 1.7222 s |
| Proteome      |  ~13 MiB | 0.0091 s | 0.0655 s |

GET 1 kbp examples through `.zfi`:

| Dataset       | z-fasta | samtools |
| ------------- | ------: | -------: |
| Human genome  |  2.1 ms |   8.7 ms |
| Transcriptome |  2.2 ms | 292.0 ms |
| Proteome      |  2.3 ms |  18.2 ms |

<details>
<summary>View the full positional GET chart</summary>

<p align="center"><img src="https://github.com/eneskemalergin/z-fasta/blob/main/bench/get/results/figures/perf_pos_wall.png?raw=1" alt="Positional GET wall time across genome, transcriptome, and proteome datasets at four region sizes" width="100%"></p>

The chart includes additional indexed peers, region sizes from 100 bp through a complete record, and repeated-run error bars. Read the labels as peer time divided by z-fasta `.zfi` time.

</details>

Complete stats examples through `.zfi`:

| Dataset       | z-fasta | SeqKit `stats -a` |
| ------------- | ------: | ----------------: |
| Human genome  | 0.413 s |          17.619 s |
| Transcriptome | 88.9 ms |           2.405 s |
| Proteome      | 10.4 ms |           64.1 ms |

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

Against v0.3.4 and v0.3.5, timed together in one batch, v0.4.0 indexes every real input faster: the transcriptome `.fai` in 123 ms against 143 ms and 157 ms. Transcriptome `.zfi` positional GET takes 2.2 ms against 8.2 ms in the v0.3.5 report, because GET checks only the index records and side tables its requests reach, and reads a `.fai` only until every requested name is found. The [benchmark framework page](https://github.com/eneskemalergin/z-fasta/blob/main/bench/README.md#v040-measurements) lists the release comparisons.

## Reports and methodology

The repository reports own the complete tables, peer definitions, dataset details, verification coverage, and measurement notes:

- [Index report](https://github.com/eneskemalergin/z-fasta/blob/main/bench/index/REPORT.md)
- [GET report](https://github.com/eneskemalergin/z-fasta/blob/main/bench/get/REPORT.md)
- [Stats report](https://github.com/eneskemalergin/z-fasta/blob/main/bench/stats/REPORT.md)
- [Validate report](https://github.com/eneskemalergin/z-fasta/blob/main/bench/validate/REPORT.md)

The Wiki intentionally keeps only representative tables and figures. The benchmark harness and its external tool setup will continue to evolve, so these pages do not promise a turnkey reproduction environment.

## Related pages

- [Choosing a FASTA tool](Choosing-a-FASTA-Tool)
- [Index formats](Index-Formats)
- [For contributors](For-Contributors)
