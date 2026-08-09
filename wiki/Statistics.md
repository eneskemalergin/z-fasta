# Statistics

Run stats when you need to understand the indexed collection, not just count records. It combines index metadata with a complete bounded scan of the indexed sequence spans.

```bash
z-fasta stats genome.fa
```

An existing `.zfi` or `.fai` is required.

## Report sections

### File

- FASTA path.
- Selected sidecar path.
- FASTA size in bytes.

### Lengths

- Indexed record count.
- Total symbols.
- Shortest and longest length and identifier.
- Integer mean.
- Q1, median, and Q3.
- Range.

### Nx

- N50 and L50.
- N90 and L90.
- auN with two decimal places.

### Composition

- Nucleotide or protein type.
- Complete category counts and percentages.
- Lowercase count.
- Nucleotide GC and GC skew where defined.

## Indexed population

Stats describes records retained by the selected index. Default indexing removes later duplicate identifiers and skips empty records. `index --no-dedup` retains duplicate records, so reports from different duplicate policies are not directly comparable.

> [!IMPORTANT]
> Run `validate` when you need source-header duplicate or empty-record reporting. Stats intentionally reports the indexed population.

## Length definitions

- Ties for shortest and longest keep the first indexed record.
- Mean uses integer division.
- Even-count median uses the midpoint rounded down.
- Q1 and Q3 are medians of the lower and upper halves.
- Nx accumulation sorts lengths from longest to shortest and stops at the requested percentage of total symbols.
- auN is `sum(length * length) / total_symbols`.

The implementation uses checked `u64` totals and `u128` intermediate arithmetic for large assemblies.

## Sequence classification

Stats counts the complete indexed symbol population and classifies it as nucleotide when case-insensitive IUPAC nucleotide symbols are strictly more than 90 percent of all symbols.

The nucleotide set includes:

```text
A C G T U R Y S W K M B D H V N
```

Empty samples classify as nucleotide, though a production index contains no empty records.

## Nucleotide composition

Type labels distinguish:

- `nucleotide_t`: T present, U absent.
- `nucleotide_u`: U present, T absent.
- `nucleotide_mixed_tu`: both T and U present.
- `nucleotide`: neither T nor U present.

Fields include A, C, G, T, U, N, combined IUPAC ambiguity, invalid symbols, and lowercase.

Category percentages use `total_symbols` as denominator. GC percent uses canonical `A + C + G + T + U` as denominator. GC is `n/a` when that denominator is zero. GC skew is `(G - C) / (G + C)` and is `n/a` when `G + C` is zero.

## Protein composition

Protein output includes:

- 20 standard residues;
- B, Z, J, and X ambiguity or unknown codes;
- U selenocysteine;
- O pyrrolysine;
- stop `*`;
- invalid symbols;
- lowercase symbols.

Percentages use `total_symbols`.

## Cross-format behavior

`.zfi` and `.fai` report equal biological values when both represent the same retained uniform records. `.zfi` also supports non-uniform records through side tables.

Stats validates that counted symbols equal the index-derived total. A mismatch fails instead of printing a plausible but inconsistent report.

## Save and compare reports

```bash
z-fasta stats assembly-a.fa > assembly-a.stats.txt
z-fasta stats assembly-b.fa > assembly-b.stats.txt
diff -u assembly-a.stats.txt assembly-b.stats.txt
```

Paths and file sizes differ by design. The text report is human-readable and does not yet define a stable machine schema. If you parse it in automation, pin the z-fasta version and treat label changes as compatibility work.

## Related pages

- [Indexing](Indexing)
- [Index formats](Index-Formats)
- [Performance and correctness](Performance-and-Correctness)
- [Troubleshooting](Troubleshooting)
