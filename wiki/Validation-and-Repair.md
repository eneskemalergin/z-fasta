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

Findings go to stdout. Operational command, mapping and write failures go to stderr. `--error-format=json` selects coded operational reports; it leaves v1 findings and warning statuses unchanged. See [Error diagnostics](Error-Diagnostics).

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

Each retained event becomes one `schema_version: v1` JSON object. Invalid UTF-8 header bytes are escaped so the output remains valid UTF-8 JSON. This legacy escaping is not a reversible encoding of arbitrary original bytes; use v2 opaque-byte objects for operational diagnostic previews.

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

At most 10000 event objects are retained. Scanning, kind counts, late fix blockers, and rewrite-width collection continue after that cap.

Without `--fix`, truncation returns exit 1 under the existing cap policy. With `--fix`, a permitted complete rewrite can finish, and remaining findings select its status. Text mode reports capped examples without claiming scanning stopped. Machine error mode omits that prose notice; v1 JSON retains `truncated: true`.

Counts keep their legacy units: invalid-symbol events count distinct invalid byte values across the input, not every occurrence; null-byte events count at most once per affected line. Empty input reports the legacy no-sequences error and missing-terminal-newline warning. These reports do not claim full header-schema validation or a fixed total memory bound.

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
