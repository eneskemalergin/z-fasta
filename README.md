<!-- markdownlint-disable MD033 MD041 -->

<div align="center">
  <img src="assets/logo-readme.svg" alt="z-fasta" width="220">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="assets/tagline-dark.svg">
    <source media="(prefers-color-scheme: light)" srcset="assets/tagline-light.svg">
    <img src="assets/tagline-light.svg" alt="Lightning FASTA work without the toolchain sprawl" width="900">
  </picture>
  <p><strong>Index, extract, validate, and inspect uncompressed FASTA with one static executable built in Zig.</strong></p>
  <p>
    <a href="https://github.com/eneskemalergin/z-fasta/actions/workflows/ci.yml"><img src="https://github.com/eneskemalergin/z-fasta/actions/workflows/ci.yml/badge.svg?branch=main" alt="CI"></a>
    <a href="https://github.com/eneskemalergin/z-fasta/releases/latest"><img src="https://img.shields.io/github/v/release/eneskemalergin/z-fasta?style=flat-square" alt="Latest release"></a>
    <a href="https://github.com/eneskemalergin/z-fasta/wiki"><img src="https://img.shields.io/badge/wiki-documentation-2563eb?style=flat-square" alt="Wiki"></a>
    <a href="CHANGELOG.md"><img src="https://img.shields.io/badge/changelog-history-7c3aed?style=flat-square" alt="Changelog"></a>
    <a href="bench/"><img src="https://img.shields.io/badge/benchmarks-reports-f59e0b?style=flat-square" alt="Benchmarks"></a>
    <a href="https://ziglang.org/download/"><img src="https://img.shields.io/badge/Zig-0.16.0-F7A41D?style=flat-square&logo=zig&logoColor=white" alt="Zig 0.16.0"></a>
    <a href="LICENSE"><img src="https://img.shields.io/github/license/eneskemalergin/z-fasta?style=flat-square" alt="MIT License"></a>
  </p>
</div>

---

I built z-fasta around the four FASTA jobs I use most. One CLI, no runtime dependencies. I kept the scope small on purpose:

- `index`: Build the default `.zfi` index or emit a samtools-compatible `.fai` for representable FASTA.
- `get`: Extract records, regions, names files, or BED intervals with optional strand and orientation transforms.
- `validate`: Report structural, alphabet, header, and layout problems, with supported repair to a separate file.
- `stats`: Report indexed length, Nx, DNA, RNA, or protein statistics through a bounded sequence scan.

z-fasta accepts uncompressed FASTA and writes FASTA or reports. It does not currently handle FASTQ, gzip, or BGZF.

## A few honest tradeoffs

I would like to support Windows properly. Right now, keeping the native builds reliable takes more time than I can justify, and I would rather be honest about that than publish an executable I cannot support well. If you use Windows, you can still run z-fasta through WSL with the Linux release. I am sorry for the extra step.

FASTQ is coming as a separate project, z-fastq, built around the same priorities: speed, low memory use, and portability (Please look forward to it).

Compressed FASTA is a different tradeoff. Supporting gzip or BGZF well would mean investing real time in Zig's deflate path or writing an efficient implementation myself. I have not been able to justify that work yet, but I would be happy to revisit it if enough people need it. (Might actuall work on it over z-fastq since there working on compressed files is a lot common, when I crack a optimized implementation will move it here)

## Start

Grab a binary from [Releases](https://github.com/eneskemalergin/z-fasta/releases), or [build it from source with Zig 0.16.0](https://github.com/eneskemalergin/z-fasta/wiki/Installation).

```bash
z-fasta validate reference.fa
z-fasta index reference.fa
z-fasta get reference.fa chr1:1000-2000
z-fasta stats reference.fa
```

Work through the [five-minute example](https://github.com/eneskemalergin/z-fasta/wiki/Getting-Started), then keep the [command cheat sheet](https://github.com/eneskemalergin/z-fasta/wiki/Command-Cheat-Sheet) nearby.

## Documentation

The Wiki holds the details that would otherwise bury this page:

[Indexing](https://github.com/eneskemalergin/z-fasta/wiki/Indexing) | [Extraction](https://github.com/eneskemalergin/z-fasta/wiki/Extracting-Sequences) | [BED and strand](https://github.com/eneskemalergin/z-fasta/wiki/BED-and-Strand-Workflows) | [Validation](https://github.com/eneskemalergin/z-fasta/wiki/Validation-and-Repair) | [Statistics](https://github.com/eneskemalergin/z-fasta/wiki/Statistics)

[Index formats](https://github.com/eneskemalergin/z-fasta/wiki/Index-Formats) | [Coordinates and names](https://github.com/eneskemalergin/z-fasta/wiki/Coordinates-and-Sequence-Names) | [Limits](https://github.com/eneskemalergin/z-fasta/wiki/Limits-and-Supported-Formats) | [Recipes](https://github.com/eneskemalergin/z-fasta/wiki/Recipes) | [Troubleshooting](https://github.com/eneskemalergin/z-fasta/wiki/Troubleshooting)

## Performance

Selected warm-cache results from an AMD Ryzen 9 3950X are shown below. Lower is better. The linked reports own the complete methods, field coverage, memory results, scaling, and peer definitions.

### Index

`.fai` indexing on real datasets:

| Dataset       |     Size |      z-fasta |  noodles | rust-bio | samtools |
| ------------- | -------: | -----------: | -------: | -------: | -------: |
| Human genome  | ~2.9 GiB | **0.3714 s** | 1.3337 s | 3.6716 s | 9.0976 s |
| Transcriptome | ~459 MiB | **0.2124 s** | 0.3658 s | 0.6995 s | 1.8315 s |
| Proteome      |  ~13 MiB | **0.0105 s** | 0.0174 s | 0.0255 s | 0.0598 s |

[Index benchmark report](bench/index/REPORT.md)

### GET

One 1 kbp positional region through each indexed implementation:

| Dataset       | z-fasta `.zfi` | z-fasta `.fai` | noodles | rust-bio | samtools |
| ------------- | -------------: | -------------: | ------: | -------: | -------: |
| Human genome  |     **2.0 ms** |         2.2 ms |  2.6 ms |   2.7 ms |   3.3 ms |
| Transcriptome |     **4.8 ms** |        27.1 ms | 88.1 ms | 530.5 ms | 292.0 ms |
| Proteome      |     **2.4 ms** |         3.6 ms |  6.8 ms |  20.0 ms |  12.7 ms |

<p align="center">
  <a href="bench/get/REPORT.md"><img src="bench/get/results/figures/perf_pos_wall.png" alt="Positional GET wall time across genome, transcriptome, and proteome datasets at four region sizes" width="100%"></a>
</p>

[GET benchmark report](bench/get/REPORT.md)

### Stats

Complete z-fasta and noodles reports, with SeqKit `stats -a` as a partial ecosystem reference:

| Dataset       | z-fasta `.zfi` | z-fasta `.fai` | noodles | SeqKit `stats -a` |
| ------------- | -------------: | -------------: | ------: | ----------------: |
| Human genome  |    **2.752 s** |        2.757 s | 6.119 s |          17.777 s |
| Transcriptome |    **0.382 s** |        0.416 s | 1.107 s |           2.403 s |
| Proteome      |    **12.3 ms** |        15.5 ms | 37.0 ms |           59.6 ms |

[Stats benchmark report](bench/stats/REPORT.md) | [Benchmark framework](bench/README.md)

## License

MIT. See [LICENSE](LICENSE).

---

<p align="center"><em>Aligned life in bytes,<br>
FASTA sings through mirrored streams.<br>
Humans bloom as code.</em></p>
