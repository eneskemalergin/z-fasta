# BED and Strand Workflows

Use BED when your intervals already live in a table, when 1024 positional arguments are not enough, or when feature strand matters. z-fasta streams the rows instead of holding the whole file in memory.

## Basic BED extraction

```bash
z-fasta get genome.fa --bed regions.bed > regions.fa
```

For this BED file:

```bed
chr1	2	8
chr2	0	4
```

GET emits requests equivalent to `chr1:3-8` and `chr2:1-4`.

> [!IMPORTANT]
> BED uses 0-based half-open intervals. Positional GET and emitted FASTA headers use 1-based inclusive intervals.

## Required fields

Each data row must be tab-separated and contain:

1. Non-empty chromosome or indexed name.
2. Unsigned decimal start.
3. Unsigned decimal end greater than start.

Columns 4 and 5 are ignored. Column 6 is used only with strand handling.

## Skipped lines

GET skips:

- empty lines;
- lines beginning with `#`;
- a line equal to `track` or beginning with `track `;
- a line equal to `browser` or beginning with `browser `.

Names such as `track1` and `browser1` remain data. The literal chromosome `track` followed by a tab also remains data.

## Strand-aware output

```bash
z-fasta get genome.fa --bed regions.bed --strand-aware > oriented.fa
```

`--honor-strand` is an alias for `--strand-aware`.

Column 6 behavior:

- `+`: forward output.
- `-`: reverse-complement output.
- `.`, an empty field, or a missing field: forward output.
- any other value: error when strand handling is enabled.

Without `--strand-aware`, column 6 is ignored, including invalid values.

## Example with strand

```bed
chr1	0	5	forward	0	+
chr1	0	5	reverse	0	-
```

```bash
z-fasta get genome.fa --bed regions.bed --strand-aware --annotate-rc
```

The plus row has no suffix. The minus row ends with ` (reverse complement)`.

## Compose strand and global orientation

BED strand is composed with the global transform. Reverse and complement each cancel when applied twice.

For a minus-strand row:

- no global flag produces reverse-complement;
- `--rc` produces identity;
- `--reverse-only` produces complement-only;
- `--complement-only` produces reverse-only.

Think of this as two steps: first orient the sequence to the BED feature, then apply the global transform you requested.

## Read BED from stdin

```bash
awk -F '\t' '$5 >= 100 && $1 !~ /^#/' raw.bed | \
  z-fasta get genome.fa --bed - --strand-aware \
  > selected.fa
```

BED streams in batches of at most 4096 requests. There is no whole-input row cap.

> [!WARNING]
> A later malformed or missing BED request can fail after earlier records reached stdout. Write pipeline output to a temporary file when partial output must never become a final artifact.

## Protein input

Minus-strand BED requires reverse-complement and is rejected when the requested record is classified as protein. If you need byte reversal for protein data, use a non-strand request with `--reverse-only` and confirm that reversal matches the biological workflow.

## BED limits

- Chromosome or record name: at most 65535 bytes.
- Interval values: unsigned `u64`, subject to coordinate validation against record length.
- Row buffering: a newline must appear within 69631 bytes.
- Request order: preserved across batches.

## Related pages

- [Extracting sequences](Extracting-Sequences)
- [Coordinates and sequence names](Coordinates-and-Sequence-Names)
- [Recipes](Recipes)
- [Troubleshooting](Troubleshooting)
