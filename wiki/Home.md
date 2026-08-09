<p align="center">
  <img src="https://github.com/eneskemalergin/z-fasta/blob/main/assets/logo-readme.svg?raw=1" alt="z-fasta" width="220">
</p>

<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="https://github.com/eneskemalergin/z-fasta/blob/main/assets/tagline-dark.svg?raw=1">
    <source media="(prefers-color-scheme: light)" srcset="https://github.com/eneskemalergin/z-fasta/blob/main/assets/tagline-light.svg?raw=1">
    <img src="https://github.com/eneskemalergin/z-fasta/blob/main/assets/tagline-light.svg?raw=1" alt="Lightning FASTA work without the toolchain sprawl" width="900">
  </picture>
</p>

<p align="center"><strong>Index, extract, validate, and inspect uncompressed FASTA with one static executable.</strong></p>

<p align="center"><a href="Installation">Install</a> | <a href="Getting-Started">Five-minute start</a> | <a href="Command-Cheat-Sheet">Cheat sheet</a> | <a href="Troubleshooting">Troubleshooting</a></p>

I built z-fasta to keep four common FASTA jobs in one place: check a file, index it, pull out the records you need, and understand what is inside. When exact behavior matters, this Wiki gets specific about coordinates, sidecars, output, and failure modes.

> [!NOTE]
> This Wiki documents z-fasta v0.3.2. Run `z-fasta --version` before following version-specific behavior.

## Start with this workflow

```bash
z-fasta validate reference.fa
z-fasta index reference.fa
z-fasta get reference.fa chr1:1000-2000
z-fasta stats reference.fa
```

Four commands take you from an unfamiliar reference to a checked, indexed, queryable file. The extraction uses 1-based inclusive coordinates.

```text
reference.fa
    |
    +-> validate
    |
    +-> index -> reference.fa.zfi
                       |
                       +-> get   -> selected FASTA
                       |
                       +-> stats -> sequence report
```

## What are you here to do?

- **Try z-fasta for the first time:** [Install it](Installation), then work through the [five-minute example](Getting-Started).
- **Extract a sequence or interval:** Start with [Extracting sequences](Extracting-Sequences) for records, regions, names files, and transforms.
- **Use BED or preserve feature strand:** Go to [BED and strand workflows](BED-and-Strand-Workflows). It keeps the two coordinate systems explicit.
- **Understand or clean up an odd FASTA:** Read [Validation and repair](Validation-and-Repair) before changing the source.
- **Share an index with FAI-based tools:** Read [Index formats](Index-Formats) before choosing `--emit-fai`.
- **Recover from an error:** Match the stderr message in [Troubleshooting](Troubleshooting).

## Why use z-fasta?

- One dependency-free executable for Linux and macOS on x86_64 and arm64.
- Default `.zfi` indexes embed names, source identity, and side tables for messy records.
- Positional GET behavior is checked for byte equality with `samtools faidx` on the verified path.
- BED, names, reverse-complement, and messy-layout extraction have dedicated correctness coverage.
- Stats scans indexed sequence spans and reports nucleotide or protein composition, Nx, auN, quartiles, and extrema.
- Validate reports the documented structure, alphabet, and header problems, then rewrites supported format issues to a separate file.

## Know the boundary

I intentionally keep z-fasta focused on uncompressed FASTA. FASTQ is coming as a separate project, z-fastq, with the same focus on speed, low memory use, and portability. Supporting gzip or BGZF well would require real work on deflate performance, and I have not been able to justify that work yet. If your workflow starts outside that boundary, [Choosing a FASTA tool](Choosing-a-FASTA-Tool) will save you time.

## The reading path I recommend

1. [Install the binary](Installation).
2. [Run the first workflow](Getting-Started).
3. Keep the [command cheat sheet](Command-Cheat-Sheet) nearby.
4. Open a command page when the basic workflow needs more control.
5. Use the concept pages when index choice, coordinates, identifiers, or classification affect correctness.
