# Recipes

These are copy-ready patterns I use for real jobs, not isolated flag examples. Read the linked command pages when you need to change a recipe or understand its failure boundary.

## Prepare a new reference

```bash
z-fasta validate --strict reference.fa
```

Continue only when validation exits 0:

```bash
z-fasta index reference.fa
z-fasta stats reference.fa > reference.stats.txt
```

Do not let a clean-looking report hide a nonzero status. Treat strict validation as the approval gate.

## Normalize for strict FAI tools

```bash
z-fasta validate source.fa
z-fasta validate --fix -o source.clean.fa source.fa
z-fasta validate --strict source.clean.fa
```

Emit the compatibility index only after strict validation exits 0:

```bash
z-fasta index --emit-fai source.clean.fa > source.clean.fa.fai
```

The strict post-fix validation is the quality gate.

## Extract more than 1024 intervals

Convert or store the requests as BED:

```bash
z-fasta get reference.fa --bed regions.bed > regions.fa
```

BED streams without a whole-input row cap.

## Extract complete records from a table

If the first tab-separated column contains identifiers:

```bash
cut -f1 selected.tsv | z-fasta get proteins.fa --names - > selected.fa
```

Make sure the stream contains no header row unless it begins with `#`.

## Preserve BED feature orientation

```bash
z-fasta get reference.fa \
  --bed genes.bed \
  --strand-aware \
  --annotate-rc \
  > genes.fa
```

Minus rows become reverse-complement. Plus, dot, empty, and missing strand remain forward.

## Request the opposite of feature orientation

```bash
z-fasta get reference.fa \
  --bed genes.bed \
  --strand-aware \
  --rc \
  --annotate-rc \
  > genes.opposite.fa
```

On minus rows, BED reverse-complement and global reverse-complement cancel. On plus rows, the global transform remains.

## Keep streamed output atomic

GET can emit earlier records before a later names or BED request fails. Write to a temporary path and publish only after success:

```bash
output=regions.fa
output_tmp=$(mktemp "${output}.tmp.XXXXXX")
if z-fasta get reference.fa --bed regions.bed > "$output_tmp"; then
  mv "$output_tmp" "$output"
else
  rm -f "$output_tmp"
  exit 1
fi
```

The temporary file sits beside the final path, so the successful `mv` remains on one filesystem.

## Capture data and summary separately

```bash
z-fasta get reference.fa --names ids.txt --summary \
  > selected.fa \
  2> selected.summary.txt
```

The summary appears only after successful output.

## Validate in CI

```bash
z-fasta validate \
  --strict \
  --json --summary \
  incoming.fa \
  > validation.json
```

Warnings become exit 1 while the JSON event levels remain unchanged.

## Validate RefSeq headers

```bash
z-fasta validate \
  --strict \
  --schema refseq \
  --max-header-len 200 \
  --json --summary \
  reference.fa \
  > validation.json
```

## Normalize an alignment alphabet

```bash
z-fasta validate --custom-alphabet 'ACGTN-.' aligned.fa
z-fasta validate \
  --custom-alphabet 'ACGTN-.' \
  --fix \
  -o aligned.clean.fa \
  aligned.fa
```

Use `--fix-format-only` only when invalid alphabet bytes are intentional and must remain unchanged.

## Compare two assemblies

```bash
z-fasta index assembly-a.fa
z-fasta index assembly-b.fa
z-fasta stats assembly-a.fa > assembly-a.stats.txt
z-fasta stats assembly-b.fa > assembly-b.stats.txt
diff -u assembly-a.stats.txt assembly-b.stats.txt
```

Align duplicate-index policy before comparing the biological metrics.

The raw diff also includes the FASTA path and file size. Use it for human review, not as a stable machine comparison contract.

## Rebuild after replacing a FASTA

```bash
z-fasta validate --strict reference.next.fa
```

Replace the source only after validation exits 0:

```bash
mv reference.next.fa reference.fa
z-fasta index reference.fa
z-fasta stats reference.fa > reference.stats.txt
```

Do not touch an old sidecar forward to hide staleness. Rebuild it against the new source.

## Related pages

- [Command cheat sheet](Command-Cheat-Sheet)
- [Troubleshooting](Troubleshooting)
- [Validation and repair](Validation-and-Repair)
