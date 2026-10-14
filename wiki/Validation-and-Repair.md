# Validation and Repair

I recommend running validate before using an unfamiliar FASTA. It reads the file once, front to back, through a fixed 256 KiB buffer, and counts every problem it finds. It never changes the file and needs no index; with `--index` it also checks the index `get` would use. It is a useful input gate, but it does not prove that every possible FASTA defect is absent. When the problems are layout only, validate can write a repaired copy in the same pass.

## Validate a file

```bash
z-fasta validate genome.fa
```

Clean input prints:

```text
OK genome.fa: no findings
complete: examined 3151425851 bytes, 51662809 lines, 194 records, 3099750718 symbols; nucleotide alphabet chosen from 100000 symbols
```

A file with problems lists each finding code once, with its count, the records it affects, and its first occurrence:

```text
INVALID sample.fa: 3 errors in 2 codes, 0 warnings in 0 codes; 2 of 3 records with errors, 0 with warnings
error [fasta.duplicate_identifier] 1 record; first at line 3 column 2, record 2 "chr1": Identifier duplicates an earlier record. First defined at line 1 column 2.
error [sequence.invalid_symbol] 2 bytes in 2 records; first at line 2 column 5, record 1 "chr1": Sequence byte is outside the selected alphabet. byte=0x3f
complete: examined 35 bytes, 6 lines, 3 records, 14 symbols; protein alphabet chosen from 14 symbols
```

The report goes to stdout. Identifiers in the text report show printable ASCII as is and escape quotes, backslashes, and every other byte as `\xNN`, so a control byte in a header cannot drive the terminal. Operational failures, such as a missing file or a failed read, go to stderr in the format `--error-format` selects; see [Error diagnostics](Error-Diagnostics).

Use `-` to read standard input:

```bash
zcat genome.fa.gz | z-fasta validate -
```

From standard input, the duplicate check keeps copies of the identifiers instead of re-reading them from the file, and a duplicate's example does not name its first definition.

## What validate reads as FASTA

Validate reads FASTA the way `index` does, so the two agree on what a symbol and an identifier are.

- A line ends with LF or CRLF. Any other CR is a control byte.
- A header is a line that starts with `>`. Its identifier runs from `>` to the first space, tab, CR, or LF.
- A symbol is any sequence byte above space. Spaces and tabs in a sequence line are skipped, as `index` and `get` skip them, and reported as whitespace.
- Other bytes below 0x20 are control bytes. Tab is allowed in a header.
- A record's lines are wrapped consistently when every line with symbols but the last has as many symbols as the first; the last may be shorter.

## What validate reports

- Errors:
  - `fasta.no_records`, once per file;
  - `fasta.content_before_header`, per line;
  - `fasta.header_control_byte` and `fasta.sequence_control_byte`, per byte;
  - `fasta.empty_identifier`, per record;
  - `fasta.identifier_too_long`, per record: an identifier over 65,535 bytes, which no `.zfi` or `.fai` can hold; it takes no part in the duplicate check;
  - `fasta.duplicate_identifier`, per record after the first;
  - `sequence.invalid_symbol`, per byte.
- Warnings:
  - `fasta.empty_record` and `fasta.long_header`, per record;
  - `header.convention_mismatch`, per record, with `--schema`;
  - `format.sequence_whitespace` and `format.blank_sequence_line`, per line;
  - `format.inconsistent_wrapping`, per record;
  - `format.leading_bom`, `format.missing_terminal_newline`, and `format.mixed_line_endings`, once per file.
- With `--index`, the six `pair.*` errors described under [Compare with the index](#compare-with-the-index).

A line before the first header counts as content when it holds anything but spaces, tabs, or CR. An empty identifier (`>` or `> desc`) is an error and is not checked for duplicates. Duplicates are compared by exact bytes. Empty input reports `fasta.no_records` only.

Counts are exact and nothing is capped: a file with a million findings reports a million, with the same exit status it would have with one. Each code keeps only its first occurrence as an example; the duplicate's example names where the identifier was first defined.

## Exit codes

- 0: no findings.
- 2: warnings only.
- 1: errors, warnings with `--strict`, or a report that could not finish.

`--strict` changes the exit status, not the severities in the report.

## Sequence type and alphabet

By default (`--type auto`), validate classifies the first 100,000 symbols with the shared strict 90 percent rule, then checks every symbol against the chosen alphabet:

- nucleotide: IUPAC codes with both T and U, in either case;
- protein: the 26 letters, `*` for a stop, and `-` for a gap, in either case.

A prefix sample cannot prove the whole file is one type, so name the type when you know it. The result then does not depend on which records come first:

```bash
z-fasta validate --type dna genome.fa
z-fasta validate --type rna transcripts.fa
z-fasta validate --type protein proteins.fa
z-fasta validate --type any reads.fa
```

- `dna`: IUPAC DNA codes without U, in either case.
- `rna`: IUPAC RNA codes, U instead of T, in either case.
- `protein`: as above.
- `any`: every symbol is accepted; structure, control-byte, identifier, and layout checks still run.

Use an exact custom alphabet when other bytes are valid:

```bash
z-fasta validate --custom-alphabet 'ACGTN-.' aligned.fa
```

Matching is case-sensitive. Include lowercase letters explicitly when needed. `--custom-alphabet` cannot be combined with a `--type` other than `auto`.

## Header policies

```bash
z-fasta validate --max-header-len 200 genome.fa
z-fasta validate --schema refseq reference.fa
z-fasta validate --schema uniprot proteins.fa
```

The default header warning threshold is 1024 bytes after `>`, including the description. It is separate from the index identifier limit. `--schema` checks the identifier's prefix and pipe-separated shape; it does not certify the full UniProt or RefSeq formats.

## JSON report

```bash
z-fasta validate --json genome.fa
```

`--json` prints one JSON object on stdout: the envelope described in [Error diagnostics](Error-Diagnostics), with `validation` filled. It holds:

- `outcome`: `clean`, `warnings`, `invalid`, or `incomplete`;
- `completion` and the selected `policy`;
- `checks`: `passed`, `warnings`, `failed`, or `incomplete` for structure, identifiers, alphabet, and format;
- what was `examined`, finding totals, totals by category, and the records with errors and with warnings;
- with `--index`, an `index_match` check and category, and `index`: the `format`, the `entries` read (all of them unless the scan stopped early), how many `matched`, and the `omitted_duplicates`;
- one entry per code that occurred, and the occurrences of each byte outside the alphabet.

<details>
<summary>One code entry</summary>

```json
{
  "code": "fasta.duplicate_identifier",
  "category": "identifiers",
  "severity": "error",
  "unit": "record",
  "occurrences": "1",
  "records": "1",
  "first": {
    "offset": "13",
    "line": "3",
    "column": "2",
    "record": "2",
    "identifier": {
      "encoding": "utf8",
      "value": "chr1",
      "truncated": false,
      "original_bytes": "4"
    },
    "related": { "offset": "1", "line": "1", "column": "2", "record": "1" },
    "data": { "identifier": null }
  }
}
```

</details>

Identifiers are `encoding` (`utf8` or `base64`), `value`, `truncated`, and `original_bytes`; decoding a value gives the original bytes, so two identifiers that differ by one byte stay different. A truncated value is a prefix of at most 256 bytes, never an identifier to compare.

`completion` is `complete` when every check ran over the whole input. It is `incomplete` when the duplicate check reached its `--max-name-memory` budget: every other check still ran to the end, the identifiers check reads `incomplete`, `name_budget` names the record where it stopped, and the exit status is 1. It is `interrupted` when a read or memory failure stopped the scan: the report covers what was examined, whole-file findings are not decided, the failure goes to stderr, and the exit status is 1.

## Compare with the index

```bash
z-fasta index genome.fa
z-fasta validate --index genome.fa
z-fasta validate --index=genome.fa.fai genome.fa
z-fasta validate --index=backup/genome.fa.zfi genome.fa
```

`--index` loads the index `get` would use, `genome.fa.zfi` or else `genome.fa.fai`; `--index=PATH` loads the named file, whose `.zfi` or `.fai` suffix gives its format (any other suffix is refused). Validate then checks, in the same pass, that the index describes this FASTA. This catches an index whose recorded size and modification time still match but whose contents do not, such as one written by an older version with different rules. `index`'s own scanner runs over the same buffer, and each record it would write is matched with an index entry by exact identifier:

- `pair.record_missing_from_index`: a record with symbols has no entry;
- `pair.record_missing_from_fasta`: an entry matches no record;
- `pair.sequence_length_mismatch`: the entry's length differs;
- `pair.sequence_offset_mismatch`: the entry's offset of the first symbol differs;
- `pair.line_geometry_mismatch`: the symbols or bytes per line differ, or one side uses the fixed-width formula and the other a side table;
- `pair.side_table_mismatch`: both use a side table and the tables differ; `line` names the first line where they do.

Each is an error, counted per record. A record without symbols is never indexed, and an index that leaves out later duplicates, as a samtools `.fai` does, is not missing them. The closing line counts the matches and the duplicates left out:

```text
OK genome.fa: no findings
complete: examined 3151425851 bytes, 51662809 lines, 194 records, 3099750718 symbols; nucleotide alphabet chosen from 100000 symbols; zfi index: 194 of 194 entries matched
```

Validate checks the index with `get`'s rules. A missing or stale index, or a damaged header or name table, fails before the scan, on stderr and without a report; its code says why, for example `index.not_found` or `index.source_size_mismatch`. Each entry is checked as the scan reaches it, so a damaged entry stops validation there: the report is `interrupted` and covers what was examined, and the entry's code goes to stderr. `--index` needs a file and cannot read standard input.

## Write a repaired copy

```bash
z-fasta validate --fix -o genome.clean.fa genome.fa
z-fasta validate --fix --width 80 -o genome.clean.fa genome.fa
```

Repair writes in the same pass as validation:

- headers exactly as they are, with LF;
- every symbol in order, wrapped at `--width` symbols (default 60);
- no BOM, no blank lines, no spaces or tabs in sequence lines, and a final newline.

Duplicate identifiers, symbols outside the alphabet, and empty records are copied as they are and still reported. The report describes the input; after a published repair, the exit status leaves out the format warnings the copy no longer has.

> [!WARNING]
> Repair never overwrites the input, including an alternate path that resolves to the same existing file. It writes a separate destination through atomic replacement. Publication is refused if required output identity is unknown or an observed source size or mtime changes. Those observations do not prove byte-for-byte source identity. A refusal preserves the previous destination; replacing a separate hard-link path leaves the input inode intact.

## What repair refuses

Repair refuses, and leaves any existing destination unchanged, when the file has:

- no records;
- content before the first header, which the copy could not keep;
- a control byte in a header or a sequence line.

The report is still printed, and the refusal goes to stderr as `repair.blocked` with the first blocking finding.

## Safe normalization workflow

```bash
z-fasta validate genome.fa
z-fasta validate --fix -o genome.clean.fa genome.fa
z-fasta validate --strict genome.clean.fa
z-fasta index genome.clean.fa
z-fasta stats genome.clean.fa
```

Treat the strict second validation as the gate. A successful write only means the rewrite completed; the copy keeps any duplicates, empty records, and invalid symbols.

## Memory and speed

Validate keeps the 256 KiB buffer, the current identifier, and one entry per distinct identifier for the duplicate check: a 32-bit fingerprint and the identifier's file offset. It reads an earlier identifier back from the file only when fingerprints match.

Bound the duplicate check when a file may hold more identifiers than memory allows:

```bash
z-fasta validate --max-name-memory 268435456 proteins.fa
```

Past the budget the duplicate check stops and the report is `incomplete`; it never evicts names or guesses. `--index` reads the index in file order through three 64 KiB buffers and rebuilds the `.zfi` name table from the entries to compare with the stored one, so it adds about 5 bytes per record. Only when the index lists records in a different order from the FASTA does validate load the whole index to match the remaining records by name. On a file the check reserves about 1 MiB before its first name. Equally wrapped lines are checked in vector blocks; the first irregular byte drops to a per-byte check at that byte. The duplicate check is the only part that grows with the input: with the number of distinct identifiers, plus a copy of each identifier that repeats, never with sequence bytes.

## Related pages

- [Indexing](Indexing)
- [Limits and supported formats](Limits-and-Supported-Formats)
- [Recipes](Recipes)
- [Troubleshooting](Troubleshooting)
