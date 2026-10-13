# Validation and Repair

I recommend running validate before using an unfamiliar FASTA. It checks the issue families documented below without changing the file or requiring an index. It is a useful input gate, but it does not prove that every possible FASTA defect is absent. When the reported problems fall within the supported fix set, validate can write a separate clean copy.

## Validate a file

```bash
z-fasta validate genome.fa
```

Clean input prints:

```text
OK: no issues found
```

Findings go to stdout. Operational command, mapping and write failures go to stderr. `--error-format=json` selects coded operational reports; it leaves v1 findings and warning statuses unchanged. See [Error diagnostics](Error-Diagnostics). `--report v2` selects the [v2 report](#report-v2), with exact counts and reversible identifiers; v1 stays the default.

## What validate checks

Errors:

- no FASTA headers;
- duplicate identifiers;
- sequence bytes outside the selected alphabet;
- null bytes.

Warnings:

- UTF-8 BOM;
- inconsistent sequence-line widths;
- trailing spaces or tabs on sequence lines;
- empty sequences;
- missing terminal newline;
- mixed LF and CRLF endings;
- long headers;
- UniProt or RefSeq schema mismatch.

## Exit codes

- 0: no remaining issues.
- 1: one or more errors.
- 2: warnings only.
- 1 with `--strict`: warnings-only is promoted for automation.

Strict mode changes exit status, not event labels.

## Sequence type and alphabet

Validate samples the first 100000 sequence bytes across the file. It calls the shared strict 90 percent classifier, then applies the chosen nucleotide or protein alphabet to the complete input.

Use an exact custom alphabet when domain-specific bytes are valid:

```bash
z-fasta validate --custom-alphabet 'ACGTN-.' aligned.fa
```

Matching is case-sensitive. Include lowercase letters explicitly when needed.

## Header policies

```bash
z-fasta validate --max-header-len 200 genome.fa
z-fasta validate --schema refseq reference.fa
z-fasta validate --schema uniprot proteins.fa
```

The default header warning threshold is 1024 bytes after `>`, including description. This warning is separate from the index identifier hard limit.

## JSON Lines

```bash
z-fasta validate --json genome.fa
```

Each retained event becomes one `schema_version: v1` JSON object. Invalid UTF-8 header bytes are escaped so the output remains valid UTF-8 JSON. This legacy escaping is not a reversible encoding of arbitrary original bytes; the [v2 report](#report-v2) keeps identifiers reversible.

## One JSON summary

```bash
z-fasta validate --json --summary genome.fa
```

The object contains event-kind counts, first examples, truncation state, sequence type, sampled base count, and sample cap. `--summary` requires `--json`.

<details>
<summary>Summary shape</summary>

```json
{
  "schema_version": "v1",
  "truncated": false,
  "sequence_type": "nucleotide",
  "type_bases_sampled": 1200,
  "type_sample_cap": 100000,
  "counts": {},
  "first_examples": {}
}
```

The real `counts` object contains every validation kind, including zero counts.

</details>

## Report v2

```bash
z-fasta validate --report v2 genome.fa
z-fasta validate --report v2 --json genome.fa
```

Report v2 counts every occurrence of each finding code in the code's unit (byte, line, record, or file) and the distinct records it affects. Each code keeps its first occurrence: location, record, identifier, and the code's facts. There is no event cap, so the counts are exact and never change the exit status. v1 stays the default, with its own output, counts, and exit statuses.

```text
INVALID sample.fa: 2 error codes, 0 warning codes; 2 of 3 records with errors, 0 with warnings
error [fasta.duplicate_identifier] 1 record; first at line 3 column 2, record 2 "chr1": Identifier duplicates an earlier record. First defined at line 1 column 2.
error [sequence.invalid_symbol] 2 bytes in 2 records; first at line 2 column 5, record 1 "chr1": Sequence byte is outside the selected alphabet. byte=0x3f
complete: examined 35 bytes, 6 lines, 3 records, 14 symbols; protein alphabet chosen from 14 symbols
```

v2 reads FASTA as `index` does, so the two agree on what a symbol and an identifier are:

- A symbol is any byte above space. Spaces and tabs in a sequence line are skipped; each such line is one `format.sequence_whitespace` warning that records its whitespace byte count.
- Other bytes below 0x20 are control bytes, errors counted one per byte: `fasta.header_control_byte` and `fasta.sequence_control_byte`. Tab is allowed in headers. A line ends with LF or CRLF; any other CR is a control byte.
- An identifier runs from `>` to the first space, tab, or CR. An empty identifier (`>` or `> desc`) is an error and is not checked for duplicates.
- An empty line inside a record is a `format.blank_sequence_line` warning. A line before the first header that holds anything but spaces, tabs, or CR is `fasta.content_before_header`.
- Line widths count symbols. Every line with symbols but the last must match the first such line; the last may be shorter.
- The alphabet is chosen from the first 100,000 symbols, or set with `--custom-alphabet`; spaces, tabs, and control bytes are never alphabet symbols.
- Empty input reports `fasta.no_records` only.

| Code | Severity | Counted per |
| --- | --- | --- |
| `fasta.no_records` | error | file |
| `fasta.content_before_header` | error | line |
| `fasta.header_control_byte` | error | byte |
| `fasta.sequence_control_byte` | error | byte |
| `fasta.empty_record` | warning | record |
| `fasta.empty_identifier` | error | record |
| `fasta.duplicate_identifier` | error | record after the first |
| `fasta.long_header` | warning | record |
| `header.convention_mismatch` | warning | record (`--schema`) |
| `sequence.invalid_symbol` | error | byte |
| `format.sequence_whitespace` | warning | line |
| `format.blank_sequence_line` | warning | line |
| `format.inconsistent_wrapping` | warning | record |
| `format.leading_bom` | warning | file |
| `format.missing_terminal_newline` | warning | file |
| `format.mixed_line_endings` | warning | file |

`--json` prints one JSON object: the envelope described in [Error diagnostics](Error-Diagnostics), with `validation` filled. It holds `completion`, the selected `policy`, what was `examined`, the records with errors and with warnings, one entry per code that occurred, and the occurrences of each byte outside the alphabet. Identifiers are `encoding` (`utf8` or `base64`), `value`, `truncated`, and `original_bytes`; decoding a value gives the original bytes, and a truncated value is only a prefix, never an identifier to compare. `--summary` changes nothing under v2.

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
    "identifier": {"encoding": "utf8", "value": "chr1", "truncated": false, "original_bytes": "4"},
    "related": {"offset": "1", "line": "1", "column": "2", "record": null},
    "data": {"identifier": null}
  }
}
```

</details>

`completion` is `complete` when every check ran over the whole input. It is `interrupted` when an operational failure, such as running out of memory for identifiers, stopped the scan: the report covers the lines before the failing one, whole-file findings are not decided, the failure goes to stderr, and the exit status is 1. `incomplete` is reserved for a check that cannot finish within a budget; no current check produces it.

Exit statuses follow v1: 0 clean, 2 warnings only, 1 for errors, for warnings with `--strict`, or for a report that is not complete. `--fix` cannot be combined with `--report v2`; the rewrite keeps v1's rules and blockers.

## Write a normalized copy

```bash
z-fasta validate --fix -o genome.clean.fa genome.fa
```

Fix normalizes:

- BOM removal;
- LF line endings;
- trailing spaces and tabs on sequence lines;
- sequence wrapping;
- terminal newline.

Each record uses its modal observed non-empty line width. A frequency tie chooses the width seen first. Empty records default to width 60 but remain empty.

> [!WARNING]
> Fix never overwrites the input, including an alternate path that resolves to the same existing file. It writes a separate destination through atomic replacement. Publication is refused if required output identity is unknown or an observed source size or mtime changes. Those observations do not prove byte-for-byte source identity. A refusal preserves the previous destination; replacing a separate hard-link path leaves the input inode intact.

## What fix refuses

Fix always refuses:

- no sequences;
- duplicate identifiers;
- null bytes;
- non-null content before the first header, because the rewrite cannot preserve it.

Invalid alphabet bytes also block fix unless you choose format-only repair:

```bash
z-fasta validate \
  --fix \
  --fix-format-only \
  -o genome.clean.fa \
  genome.fa
```

Invalid bytes remain unchanged and still affect the report and exit code.

After a successful fix, the exit code ignores only the five warnings the rewrite clears: BOM, inconsistent widths, trailing sequence whitespace, missing terminal newline, and mixed line endings.

## Event cap

At most 10000 v1 event objects are retained. Scanning, kind counts, late fix blockers, and rewrite-width collection continue after that cap. Report v2 has no event cap.

Without `--fix`, truncation returns exit 1 under the existing cap policy. With `--fix`, a permitted complete rewrite can finish, and remaining findings select its status. Text mode reports capped examples without claiming scanning stopped. Machine error mode omits that prose notice; v1 JSON retains `truncated: true`.

v1 counts keep their legacy units: invalid-symbol events count distinct invalid byte values across the input, not every occurrence; null-byte events count at most once per affected line. Empty input reports the legacy no-sequences error and missing-terminal-newline warning. Report v2 counts every occurrence. These reports do not claim full header-schema validation or a fixed total memory bound.

## Safe normalization workflow

```bash
z-fasta validate genome.fa
z-fasta validate --fix -o genome.clean.fa genome.fa
z-fasta validate --strict genome.clean.fa
z-fasta index genome.clean.fa
z-fasta stats genome.clean.fa
```

Treat the strict second validation as the gate. A successful write only means the allowed rewrite completed; the new file may still contain non-format warnings or preserved invalid bytes.

## Related pages

- [Indexing](Indexing)
- [Limits and supported formats](Limits-and-Supported-Formats)
- [Recipes](Recipes)
- [Troubleshooting](Troubleshooting)
