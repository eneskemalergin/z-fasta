# Coordinates and Sequence Names

Bookmark this page if a workflow mixes command-line regions and BED. Coordinate and identifier mistakes often produce believable output, which makes the exact conventions especially important.

## Positional GET coordinates

Positional regions are 1-based inclusive:

```text
NAME:START-END
```

For this sequence:

```text
Position  1 2 3 4 5 6 7 8
Sequence  A C G T A C G T
```

`seq:3-6` returns `GTAC`.

The interval length is:

```text
END - START + 1
```

## BED coordinates

BED is 0-based half-open:

```bed
seq	2	6
```

The same bytes become positional `seq:3-6` and return `GTAC`.

Conversion:

```text
positional_start = bed_start + 1
positional_end   = bed_end
```

> [!IMPORTANT]
> Never add one to the BED end. The exclusive BED end becomes the inclusive 1-based end unchanged.

## Full and open-ended records

```bash
# Complete record
z-fasta get genome.fa chr1

# From base 1000 through record end
z-fasta get genome.fa chr1:1000-
```

A bare name and an explicit range can return the same sequence bytes but use different output headers.

## End clamping

If end exceeds sequence length, z-fasta extracts through the real end but preserves the requested end in the header.

For a 12-base record:

```bash
z-fasta get genome.fa chr1:10-20
```

```fasta
>chr1:10-20
CGT
```

Start beyond sequence length is an error.

## Indexed identifier

The indexed name begins immediately after `>` and ends at the first space, tab, CR, LF, or file end.

```fasta
>chr1 chromosome one
```

Identifier: `chr1`

Description: `chromosome one`

GET output uses the identifier and does not restore the description.

## Exact byte matching

Identifiers are compared as exact byte strings:

- case-sensitive;
- no Unicode normalization;
- no leading or trailing whitespace trimming;
- invalid UTF-8 allowed in indexes;
- maximum 65535 bytes.

Names input skips blank and `#` lines but otherwise preserves bytes.

## Names containing colons

GET parses only a valid coordinate suffix after the rightmost colon:

```text
chromosome:GRCh38:1:1:248956422:1
chromosome:GRCh38:1:1:248956422:1:100-200
```

The first token is a complete identifier. The second selects bases 100 through 200 from that identifier.

If the rightmost suffix is malformed or overflows `u64`, GET treats the whole token as a literal name.

## Names ending in coordinate syntax

A literal identifier that ends in a valid `:START-END` or `:START-` suffix is interpreted as a positional region. Use names input when you need the complete record with those literal bytes:

```bash
printf '%s\n' 'sample:1-100' | z-fasta get genome.fa --names -
```

Names input does not parse coordinate suffixes.

## Duplicate names

Default indexing keeps the first non-empty record for each name. `--no-dedup` retains all index records, but GET still resolves the first exact match.

Validate reports duplicate identifiers as errors. Rename the source records when each one must remain independently addressable.

## Empty names

Empty identifiers are representable but difficult to address. An explicit empty positional shell argument can retrieve the first one:

```bash
z-fasta get empty-name.fa ''
```

Names files cannot request an empty identifier because blank lines are skipped. BED requires a non-empty chromosome field.

## Related pages

- [Extracting sequences](Extracting-Sequences)
- [BED and strand workflows](BED-and-Strand-Workflows)
- [Indexing](Indexing)
- [FAQ](FAQ)
