# Getting Started

You can learn the whole workflow with one tiny file. In about five minutes, you will validate it, build an index, extract a region, transform it, and inspect its statistics.

## 1. Create a FASTA

Save this as `example.fa`:

```fasta
>chr1 chromosome one
ACGTACGTACGT
>rna
AUGCRY
>protein
MPEPTIDE
```

z-fasta uses the identifier before the first space or tab. The indexed names are `chr1`, `rna`, and `protein`.

## 2. Validate it

```bash
z-fasta validate example.fa
```

Expected output:

```text
OK: no issues found
```

Validation does not require an index.

## 3. Build the default index

```bash
z-fasta index example.fa
```

Expected stderr:

```text
wrote example.fa.zfi (3 sequences)
```

You now have `example.fa.zfi`, the sidecar GET and stats will use.

## 4. Extract a complete record

```bash
z-fasta get example.fa chr1
```

```fasta
>chr1
ACGTACGTACGT
```

## 5. Extract a region

```bash
z-fasta get example.fa chr1:3-8
```

```fasta
>chr1:3-8
GTACGT
```

Positional coordinates are 1-based inclusive. The third through eighth symbols produce six output symbols.

## 6. Apply a transform

```bash
z-fasta get example.fa chr1:3-8 --rc --annotate-rc
```

```fasta
>chr1:3-8 (reverse complement)
ACGTAC
```

Complement-based transforms classify up to 256 requested bases and reject protein records.

## 7. Print statistics

```bash
z-fasta stats example.fa
```

The report brings file details, length summaries, Nx statistics, auN, and complete nucleotide or protein composition together in one place.

<details>
<summary>Why this mixed file is classified as protein</summary>

Stats classifies the complete indexed symbol population. A sample is nucleotide only when case-insensitive IUPAC nucleotide bytes are strictly more than 90 percent of all symbols. The protein record pushes this small combined file below that threshold.

</details>

## What to read next

- [Indexing](Indexing) for `.zfi`, `.fai`, duplicate identifiers, and publication safety.
- [Extracting sequences](Extracting-Sequences) for multiple regions and names files.
- [BED and strand workflows](BED-and-Strand-Workflows) for interval files.
- [Validation and repair](Validation-and-Repair) for messy FASTA.
- [Statistics](Statistics) for report definitions.
