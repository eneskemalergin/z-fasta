# Troubleshooting

My first rule is to start with the exact stderr message. Do not delete sidecars or rewrite the FASTA yet: most failures identify one specific input or index contract, and the original files are useful evidence.

## No index found

```text
error: no index found for genome.fa. Run 'z-fasta index genome.fa' first.
```

Create the default index:

```bash
z-fasta index genome.fa
```

GET and stats do not create indexes automatically.

## Index is stale

```text
error: index is stale (FASTA size or mtime does not match the index). Re-run 'z-fasta index <path>'.
```

The FASTA changed after indexing or source identity no longer matches. Rebuild:

```bash
z-fasta index genome.fa
```

If you intentionally use `.fai`, remove or relocate the stale `.zfi` before retrying because `.zfi` is authoritative.

## Corrupt index despite a valid `.fai`

```text
error: corrupt index file for: genome.fa
```

Check whether `genome.fa.zfi` exists. A present invalid `.zfi` blocks `.fai` fallback. Preserve it if you need forensic evidence, then regenerate or remove it deliberately.

> [!WARNING]
> Do not silently delete a sidecar in an automated pipeline without recording why. A corrupt authoritative index may indicate interrupted publication, storage damage, or a source/index mismatch.

## `.fai` emission rejects the FASTA

```text
error: cannot emit .fai for non-uniform sequence layout; run 'z-fasta index' (default) to write .zfi
```

Choose one:

```bash
# Keep the original layout for z-fasta
z-fasta index genome.fa

# Or normalize a separate copy for FAI tools
z-fasta validate --fix -o genome.clean.fa genome.fa
z-fasta validate --strict genome.clean.fa
z-fasta index --emit-fai genome.clean.fa > genome.clean.fa.fai
```

## Sequence not found

```text
error: sequence not found: chr1
```

Work through these checks in order:

- Identifier spelling and case.
- Header token before the first space or tab.
- Duplicate and empty-record policy.
- Names-file whitespace.
- Whether a colon suffix was interpreted as coordinates.
- Whether the chosen sidecar belongs to this FASTA.

Validate duplicate or empty records with:

```bash
z-fasta validate genome.fa
```

## Coordinate error

Common causes:

- Using BED coordinates as positional coordinates.
- Start is zero.
- Start exceeds sequence length.
- End is less than start.
- BED end equals start.

Review [Coordinates and sequence names](Coordinates-and-Sequence-Names).

## Protein reverse-complement error

```text
error: reverse complement is not defined for protein sequences: NAME (classified from up to 256 bases)
```

Use `--reverse-only` only when simple byte reversal is the intended operation. Otherwise remove the complement-based flag or correct the source classification problem.

For unusual nucleotide alphabets, validate with an explicit alphabet and inspect the requested leading sample.

## No regions provided

The names or BED source contained only blank, comment, or directive lines. Confirm the input stream:

```bash
sed -n '1,20p' ids.txt
sed -n '1,20p' regions.bed
```

Names lines beginning with `#` are skipped. BED also skips `track` and `browser` directives.

## GET wrote only a prefix

Names and BED stream in batches. Earlier valid records may reach stdout before a later request fails. The summary is suppressed on failure.

Use the temporary-output recipe in [Recipes](Recipes) when a partial file must never be published.

## Validate returned exit 2

Exit 2 means warnings only. Read stdout for the exact warnings.

Use `--strict` when warnings must fail automation:

```bash
z-fasta validate --strict genome.fa
```

## Validate stopped after 10000 events

The report is incomplete, though scanning and counts continued. Fix the visible problem class and rerun. JSON summary includes `truncated: true`.

With `--fix`, a complete rewrite can still finish, but stderr warns that retained examples were capped.

## Fix refuses the file

Fix cannot repair no-sequence, duplicate-name, or null-byte errors. Resolve those in the source.

Invalid alphabet bytes require either source correction, a correct `--custom-alphabet`, or an explicit format-only decision:

```bash
z-fasta validate --fix --fix-format-only -o output.fa input.fa
```

The invalid bytes remain and still affect exit status.

## Fix refuses the output path

z-fasta does not overwrite the input, including an alternate path spelling or symlink that resolves to the same existing file. Choose a distinct destination.

## File is empty or not FASTA

Index requires a non-empty input whose first byte is `>`. A UTF-8 BOM before the first header makes the first byte something else. Validate can report the BOM and write a clean copy:

```bash
z-fasta validate --fix -o clean.fa input.fa
z-fasta index clean.fa
```

## Unknown option

Flags are command-specific and unknown options are rejected anywhere. Confirm command placement with:

```bash
z-fasta --help
z-fasta get --help
```

Subcommand help currently prints the combined command reference.

## Still stuck?

Collect the smallest useful report:

- `z-fasta --version`;
- operating system and architecture;
- exact command;
- exit status;
- stdout and stderr separately;
- FASTA size and line-ending style;
- names of present `.zfi` and `.fai` sidecars;
- a minimized reproducible FASTA when data policy permits.

If the problem remains, [open a GitHub issue](https://github.com/eneskemalergin/z-fasta/issues) with the smallest safe reproducer. Do not publish private biological data.
