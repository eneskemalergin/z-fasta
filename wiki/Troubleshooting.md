# Troubleshooting

Start with the exact stderr code and observed facts. See [Error diagnostics](Error-Diagnostics) for text and JSON reports. Do not delete sidecars or rewrite the FASTA yet: most failures identify one specific input or index contract, and the original files are useful evidence.

## No index found

The code is `index.not_found`. GET and stats tried the named sidecars without finding an index.

Create the default index:

```bash
z-fasta index genome.fa
```

GET and stats do not create indexes automatically.

## Index is stale

The code is `index.older_than_source`, `index.source_size_mismatch` or `index.source_mtime_mismatch`. The facts name the observed timestamps or sizes.

The FASTA changed after indexing or source identity no longer matches. Rebuild:

```bash
z-fasta index genome.fa
```

If you intentionally use `.fai`, remove or relocate the stale `.zfi` before retrying because `.zfi` is authoritative.

## Corrupt index despite a valid `.fai`

Read the specific `index.*` code, such as `index.truncated`, `index.invalid_magic` or `index.name_out_of_bounds`, and its original source location.

Check whether `genome.fa.zfi` exists. A present invalid `.zfi` blocks `.fai` fallback. Preserve it if you need forensic evidence, then regenerate or remove it deliberately.

> [!WARNING]
> Do not silently delete a sidecar in an automated pipeline without recording why. A corrupt authoritative index may indicate interrupted publication, storage damage, or a source/index mismatch.

## `.fai` emission rejects the FASTA

The code is `index.nonuniform_fai`. The observed layout cannot be represented by a uniform FAI row.

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

The code is `request.sequence_not_found`. Its identifier preview and original request position identify the rejected request.

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

The code names the rejected rule:

- `request.start_before_one`: positional start is zero. Positional regions are one-based; BED coordinates used as positional coordinates often cause this.
- `request.start_after_sequence`: start exceeds the indexed sequence length.
- `request.end_before_start`: end is less than start.
- `request.bed.invalid_interval`: BED end is not greater than start.
- `request.bed.field_missing`, `request.bed.invalid_coordinate`, `request.bed.coordinate_out_of_range`, `request.bed.invalid_strand`: a malformed BED row; the report gives its line and column.

Review [Coordinates and sequence names](Coordinates-and-Sequence-Names).

## Protein reverse-complement error

The code is `transform.complement_on_protein`. The report retains the selected transform and the shared classifier's observed sample facts.

Use `--reverse-only` only when simple byte reversal is the intended operation. Otherwise remove the complement-based flag or correct the source classification problem.

For unusual nucleotide alphabets, validate with an explicit alphabet and inspect the requested leading sample.

## No regions provided

The code is `request.no_regions`. The names or BED source contained only blank, comment, or directive lines. Confirm the input stream:

```bash
sed -n '1,20p' ids.txt
sed -n '1,20p' regions.bed
```

Names lines beginning with `#` are skipped. BED also skips `track` and `browser` directives.

## GET wrote only a prefix

Names and BED stream in batches. Earlier valid records may reach stdout before a later request fails. The summary is suppressed on failure.

Use the temporary-output recipe in [Recipes](Recipes) when a partial file must never be published.

## Output closed early

The code is `io.write_failed` with `cause="BrokenPipe"`, and the exit status is 1. The reader closed stdout before z-fasta finished, for example `z-fasta get genome.fa chr1 | head`.

The data that reached the reader is correct; z-fasta reports that it could not deliver the rest. For a quick preview, request a smaller region instead of truncating the stream:

```bash
z-fasta get genome.fa chr1:1-1000
```

In a pipeline with `set -o pipefail`, an early-closing reader makes the pipeline fail.

## Missing or unreadable file

- `io.file_not_found`: the named FASTA, sidecar, names, or BED file does not exist.
- `io.access_denied`: the file exists but cannot be opened with the current permissions.
- `io.invalid_path`: the path is too long or otherwise rejected before opening; `reason` and `limit_bytes` say which.

The report's source label says which input failed, and `cause` keeps the operating-system reason.

## Input changed during work

The code is `input.changed`. Index, validate repair, and GET names or BED files compare the file's size and modification time during the work. A difference, or a file that disappeared, stops the command; index and repair do not publish their output.

Wait until the writer finishes, then rerun. The check compares size and time only; it does not prove identical bytes.

## Validate returned exit 2

Exit 2 means warnings only. Read stdout for the exact warnings.

Use `--strict` when warnings must fail automation:

```bash
z-fasta validate --strict genome.fa
```

## Validate reports an interrupted check

`INTERRUPTED` (JSON `completion: interrupted`) means a read or memory failure stopped the scan. The report covers what was examined, whole-file findings such as a missing final newline are not decided, the failure is on stderr, and the exit status is 1. Fix the cause and rerun; a partial report is never a pass.

## Fix refuses the file

`repair.blocked` names the first finding the copy could not keep: no records, content before the first header, or a control byte in a header or sequence line. The report on stdout lists every finding. Resolve those in the source; an existing destination is left unchanged.

Duplicate identifiers and symbols outside the alphabet do not block repair. The copy keeps them, and they still affect the exit status.

## Fix refuses the output path

`repair.output_is_input` rejects the same existing file, including alternate path spellings and symlinks. `repair.output_identity_unknown` refuses publication when the required comparison cannot be completed. Choose a destination whose identity can be checked.

A separate hard-link destination can be atomically replaced without modifying the input inode. Refused repair preserves an existing destination. Observed source size or mtime changes before publication produce `input.changed`; this check does not prove byte-for-byte identity.

## File is empty or not FASTA

- `fasta.no_records`: the input is empty.
- `index.initial_header_required`: the first byte is not `>`. The report's `byte` shows it; `0xef` is the start of a UTF-8 BOM.
- `index.no_indexable_records`: every record is empty. The report counts seen, empty, and duplicate records.

Index requires a non-empty input whose first byte is `>`. Validate can report a BOM and write a clean copy:

```bash
z-fasta validate --fix -o clean.fa input.fa
z-fasta index clean.fa
```

## Unknown option

The code is `usage.unknown_option`; other `usage.*` codes report missing, repeated, or conflicting arguments. Flags are command-specific and unknown options are rejected anywhere. Confirm command placement with:

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
