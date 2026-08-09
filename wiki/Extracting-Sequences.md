# Extracting Sequences

Use GET when you know the record names or intervals you want. It reads an existing index and writes complete records or 1-based inclusive regions to stdout, wrapped at 60 symbols per line.

## Choose one request source

```text
z-fasta get <file.fasta> [options] <region> [region ...]
z-fasta get <file.fasta> [options] --names <file.txt|->
z-fasta get <file.fasta> [options] --bed <file.bed|->
```

Positional regions, `--names`, and `--bed` are mutually exclusive request sources.

## Complete records

```bash
z-fasta get genome.fa chr1
```

```fasta
>chr1
ACGT...
```

The output header uses the indexed identifier, not the original description.

## Regions

```bash
z-fasta get genome.fa chr1:1000-2000
```

Regions are 1-based inclusive. `1000-2000` contains 1001 symbols.

Open-ended ranges continue through record end:

```bash
z-fasta get genome.fa chr1:1000-
```

Start must be at least 1 and no greater than sequence length. End must be at least start.

> [!NOTE]
> An end beyond sequence length is clamped for extraction, but the FASTA header preserves the requested end to match samtools-faidx behavior.

## Several regions

```bash
z-fasta get genome.fa \
  chr1:1-100 \
  chr2:500-900 \
  chrX \
  > selected.fa
```

Requests are emitted in command-line order. One invocation loads the index once, and neighboring requests may share a bounded disk read. Positional input accepts at most 1024 regions.

For a larger interval set, I recommend BED.

## Names files

```bash
z-fasta get proteins.fa --names accessions.txt > selected.fa
```

A names file is the simplest choice when you want complete records. Put one literal identifier on each line:

```text
# mitochondrial proteins
YP_003024026.1
YP_003024027.1
```

Blank lines and lines beginning with `#` are skipped. CRLF is accepted. Other whitespace is significant and is not trimmed.

Names files do not parse coordinate suffixes. A line `chr1:1-100` looks for a record whose literal identifier contains those bytes.

Read names from stdin:

```bash
printf 'chr1\nchr2\n' | z-fasta get genome.fa --names -
```

Names input streams in reusable batches and has no whole-input row limit.

## Identifiers containing colons

GET checks the rightmost colon for a valid decimal coordinate suffix. This keeps Ensembl-style identifiers usable:

```bash
z-fasta get genome.fa 'chromosome:GRCh38:1:1:248956422:1'
z-fasta get genome.fa 'chromosome:GRCh38:1:1:248956422:1:100-200'
```

An invalid or overflowing suffix makes the whole token a literal identifier.

## Reverse, complement, and reverse-complement

```bash
z-fasta get genome.fa chr1:100-200 --rc
z-fasta get genome.fa chr1:100-200 --reverse-only
z-fasta get genome.fa chr1:100-200 --complement-only
```

These flags are mutually exclusive.

- `--rc`: reverse order and complement.
- `--reverse-only`: reverse order without complementing.
- `--complement-only`: complement without changing order.

The complement map covers uppercase and lowercase IUPAC nucleotide symbols and maps `U/u` to `A/a`. Unknown bytes pass through the table, but complement-based transforms first classify up to 256 requested bases. Protein classification is rejected. Reverse-only remains available for protein.

## Annotate transformed headers

```bash
z-fasta get genome.fa chr1:100-200 --rc --annotate-rc
```

The final header receives one of these suffixes:

```text
 (reverse complement)
 (reverse)
 (complement)
```

Identity output has no annotation. `--annotate-rc` requires a global transform or strand-aware BED input.

## Summary output

```bash
z-fasta get genome.fa --names ids.txt --summary \
  > selected.fa \
  2> selected.summary.txt
```

The summary includes region count, total bases, elapsed seconds, and regions per second. Timing starts before index loading and includes request acquisition.

If any request fails, no summary is printed. Names and BED streams may already have emitted complete earlier records before a later failure.

## Output compatibility

The verified positional path matches `samtools faidx` output bytes, including headers, end clamping, request order, and 60-symbol wrapping. Names and BED behavior has separate peer and project correctness coverage.

## Related pages

- [Coordinates and sequence names](Coordinates-and-Sequence-Names)
- [BED and strand workflows](BED-and-Strand-Workflows)
- [Troubleshooting](Troubleshooting)
- [Recipes](Recipes)
