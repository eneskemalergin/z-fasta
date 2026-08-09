# Command Cheat Sheet

Keep this page nearby once you know the basic workflow. It favors commands you can paste over explanations; each command page owns the details and edge cases.

## General

```bash
z-fasta --help
z-fasta --version
```

## Index

```bash
# Default .zfi index
z-fasta index genome.fa

# Samtools-compatible .fai for uniform FASTA
z-fasta index --emit-fai genome.fa > genome.fa.fai

# Retain duplicate index records
z-fasta index --no-dedup genome.fa
```

## Get

```bash
# Complete record
z-fasta get genome.fa chr1

# 1-based inclusive region
z-fasta get genome.fa chr1:1000-2000

# Open-ended region
z-fasta get genome.fa chr1:1000-

# Several regions in request order
z-fasta get genome.fa chr1:1-100 chr2:1-100

# Literal full-record names from a file or stdin
z-fasta get genome.fa --names ids.txt
printf 'chr1\nchr2\n' | z-fasta get genome.fa --names -

# BED from a file or stdin
z-fasta get genome.fa --bed regions.bed
awk '$5 >= 100' regions.bed | z-fasta get genome.fa --bed -

# BED strand
z-fasta get genome.fa --bed regions.bed --strand-aware

# Orientation
z-fasta get genome.fa chr1:1-100 --rc
z-fasta get genome.fa chr1:1-100 --reverse-only
z-fasta get genome.fa chr1:1-100 --complement-only
z-fasta get genome.fa chr1:1-100 --rc --annotate-rc

# Successful-run summary to stderr
z-fasta get genome.fa --bed regions.bed --summary > regions.fa
```

## Validate

```bash
# Human report
z-fasta validate genome.fa

# Fail a warnings-only quality gate
z-fasta validate --strict genome.fa

# JSON Lines or one JSON summary
z-fasta validate --json genome.fa
z-fasta validate --json --summary genome.fa

# Header policies
z-fasta validate --schema refseq genome.fa
z-fasta validate --max-header-len 200 genome.fa
z-fasta validate --custom-alphabet ACGTN- aligned.fa

# Write a normalized copy
z-fasta validate --fix -o genome.clean.fa genome.fa

# Keep invalid alphabet bytes while fixing format
z-fasta validate --fix --fix-format-only -o genome.clean.fa genome.fa
```

## Stats

```bash
z-fasta stats genome.fa
```

## Coordinate reminder

- Positional GET: 1-based inclusive.
- BED input: 0-based half-open.
- GET FASTA headers: 1-based inclusive.

> [!WARNING]
> A present `.zfi` is authoritative. If it is stale or corrupt, z-fasta fails instead of falling back to a valid `.fai`.
