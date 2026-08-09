# Indexing

Index once, then GET and stats can jump to the sequence bytes they need. I recommend the default `.zfi` when z-fasta is the only consumer. Validate does not need an index.

## Build the default `.zfi`

```bash
z-fasta index genome.fa
```

The command writes `genome.fa.zfi` and reports the retained sequence count to stderr.

```text
wrote genome.fa.zfi (N sequences)
```

Use `.zfi` unless another program explicitly asks for `.fai`. It embeds identifiers, stores source size and mtime, and adds per-line side tables when fixed byte geometry is unsafe.

## Emit a compatible `.fai`

```bash
z-fasta index --emit-fai genome.fa > genome.fa.fai
```

`.fai` is written to stdout. Name the redirected file `<fasta>.fai` so z-fasta can find it when `.zfi` is absent.

> [!IMPORTANT]
> `.fai` succeeds only when every retained record has fixed line geometry. Use default `.zfi` indexing for mixed widths, interior blank lines, trailing sequence whitespace, or mixed line endings.

Before writing stdout, z-fasta scans into an exclusive temporary spool. A scan failure does not replay a partial index. Shell redirection may still create an empty destination file.

## What counts as a sequence name

Given this header:

```fasta
>chr1 chromosome one
```

The indexed identifier is `chr1`. Parsing stops at the first space, tab, CR, LF, or end of file. Header descriptions are not lookup keys.

Identifiers are byte strings and do not need to be valid UTF-8. The maximum identifier length is 65535 bytes.

## Empty records

Records with no sequence symbols are skipped:

```fasta
>empty
>kept
ACGT
```

The index contains only `kept`. Validate still reports `empty` as an empty-sequence warning.

## Duplicate identifiers

Default indexing keeps the first non-empty record for each exact identifier:

```bash
z-fasta index duplicates.fa
```

Keep every index record when stats must describe the full source population:

```bash
z-fasta index --no-dedup duplicates.fa
```

GET still resolves a duplicate identifier to the first exact match. If every record needs independent retrieval, rename duplicates in the FASTA instead of relying on index policy.

## Uniform and non-uniform records

A record is uniform when byte offsets can use one line formula. These ordinary cases stay uniform:

- equal-width sequence lines;
- a shorter final line;
- a final line without a newline;
- a final CRLF line after an LF body when base width remains valid.

These cases require a `.zfi` side table:

- an interior line with a different base width;
- a final line wider than the established width;
- trailing spaces or tabs on sequence lines;
- an interior blank line;
- an interior LF and CRLF width change.

Read [Index formats](Index-Formats) for the operational difference.

## Publication safety

The default path scans the source, checks that the source path still has the original size and mtime, writes `<fasta>.zfi.tmp`, then renames the completed file to `<fasta>.zfi`.

The `.fai` path also checks source size and mtime before replaying its spool to stdout. If the source changes during either workflow, indexing fails.

## Duplicate sidecars

You may keep both `genome.fa.zfi` and `genome.fa.fai`. GET and stats always choose `.zfi` when it exists. `.fai` is considered only when `.zfi` is absent.

## Recommended workflows

### z-fasta only

```bash
z-fasta validate genome.fa
z-fasta index genome.fa
```

### Shared with FAI-based tools

```bash
z-fasta validate --strict genome.fa
z-fasta index --emit-fai genome.fa > genome.fa.fai
```

If `.fai` emission rejects the layout, normalize a separate copy:

```bash
z-fasta validate --fix -o genome.clean.fa genome.fa
z-fasta validate --strict genome.clean.fa
z-fasta index --emit-fai genome.clean.fa > genome.clean.fa.fai
```

## Related pages

- [Index formats](Index-Formats)
- [Validation and repair](Validation-and-Repair)
- [Troubleshooting](Troubleshooting)
- [Limits and supported formats](Limits-and-Supported-Formats)
