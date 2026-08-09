# Index Formats

My default recommendation is simple: use `.zfi` for z-fasta. Reach for `.fai` when another tool needs the established five-column format.

## Default `.zfi`

Create it with:

```bash
z-fasta index genome.fa
```

`.zfi` stores:

- source FASTA size;
- source FASTA mtime in current production indexes;
- indexed record geometry;
- embedded identifiers;
- optional side tables for non-uniform records.

Uniform records use direct byte-offset arithmetic. Non-uniform records use a binary-searchable per-line table.

## Compatibility `.fai`

Create it with:

```bash
z-fasta index --emit-fai genome.fa > genome.fa.fai
```

Each FASTA record is one tab-separated line containing name, sequence length, first-base byte offset, line bases, and line bytes.

`.fai` can represent ordinary fixed-width records whose final sequence line may be shorter. It cannot represent arbitrary mixed widths, interior blank lines, or sequence-line whitespace.

## Which format is selected

GET and stats follow this order:

1. Look for `<fasta>.zfi`.
2. If `.zfi` exists, load it or fail.
3. Look for `<fasta>.fai` only when `.zfi` is absent.
4. Fail when neither usable sidecar exists.

> [!WARNING]
> A stale, corrupt, inaccessible, or unsupported `.zfi` blocks `.fai` fallback. Remove or regenerate it deliberately; z-fasta never hides an authoritative-index failure.

## Staleness and source identity

Current `.zfi` loading checks:

- sidecar mtime is not older than FASTA mtime;
- embedded source size equals FASTA size;
- embedded source mtime equals FASTA mtime.

Legacy embedded-name `.zfi` files without the source-identity trailer use the weaker sidecar-age and source-size checks.

`.fai` has no embedded source identity. z-fasta can only check that the sidecar is not older than the FASTA.

Timestamp-preserving same-size source replacement remains a residual risk. Sidecars are not full-content fingerprints.

## Messy FASTA

`.zfi` side tables support:

- mixed sequence-line widths;
- trailing spaces or tabs;
- interior blank lines;
- mixed LF and CRLF;
- a missing terminal newline.

If another program needs `.fai`, normalize a separate source copy:

```bash
z-fasta validate --fix -o genome.clean.fa genome.fa
z-fasta validate --strict genome.clean.fa
z-fasta index --emit-fai genome.clean.fa > genome.clean.fa.fai
```

## Keep both formats

You can keep both when you understand precedence:

```text
genome.fa
genome.fa.zfi
genome.fa.fai
```

z-fasta uses `.zfi`. Other FAI consumers can use `.fai`.

## Compatibility boundary

The `.zfi` layout is a versioned little-endian format with explicit magic and version bytes. Unknown versions are rejected rather than reinterpreted. `.fai` follows the five-column FASTA layout used by samtools-compatible tools.

## Related pages

- [Indexing](Indexing)
- [Limits and supported formats](Limits-and-Supported-Formats)
- [Troubleshooting](Troubleshooting)
- [Choosing a FASTA tool](Choosing-a-FASTA-Tool)
