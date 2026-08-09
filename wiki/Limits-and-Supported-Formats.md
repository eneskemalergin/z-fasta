# Limits and Supported Formats

z-fasta has a deliberately narrow input contract. Check this page before building a pipeline around compressed data, another sequence format, or unusually large request fields.

## Supported today

- Uncompressed FASTA input.
- `.zfi` default indexes.
- Five-column FASTA `.fai` compatibility indexes.
- Linux and macOS.
- x86_64 and arm64 release archives.
- DNA, RNA, protein, lowercase, and IUPAC symbols.
- LF, CRLF, missing final newline, and non-uniform layout through `.zfi` side tables.

## Not supported

- FASTQ.
- gzip-compressed FASTA.
- BGZF FASTA.
- Native Windows releases. Use the Linux release through WSL.
- Remote URLs as FASTA paths.
- Automatic index creation during GET or stats.
- BED12 block concatenation.
- GFF or VCF interval input.
- In-place validation repair.

> [!NOTE]
> Use samtools, bedtools, SeqKit, or another suitable tool when a workflow needs compressed FASTA, FASTQ, BED12 blocks, GFF/VCF intervals, or broader format conversion. See [Choosing a FASTA tool](Choosing-a-FASTA-Tool).

FASTQ is coming as a separate project, z-fastq, built around the same priorities: speed, low memory use, and portability. Supporting gzip or BGZF well would mean investing real time in Zig's deflate path or writing an efficient implementation myself. I have not been able to justify that work yet, but I would be happy to revisit it if enough people need it.

## Input limits

- Indexed identifier: 65535 bytes.
- Positional GET regions: 1024 per invocation.
- Names file identifier: 65535 bytes.
- BED chromosome or identifier: 65535 bytes.
- BED or names line reader: 69631 bytes before newline is required.
- Retained validate events: 10000.
- Default validate header warning threshold: 1024 bytes.
- Validate type sample: first 100000 sequence bytes across the file.
- GET complement type sample: up to 256 requested bases per record.

## Streaming batches

- Names: at most 65536 active requests per batch.
- BED: at most 4096 active requests per batch.
- Unique active request-name storage: 4 MiB plus the final name that crosses that target.

These are memory bounds, not whole-input row caps. Names and BED continue with a fresh reusable batch.

## Symbol treatment

Index and GET treat bytes greater than ASCII space as sequence symbols. GET side-table extraction ignores ASCII whitespace between symbol runs.

Validate trims only spaces and tabs from the right edge of sequence lines before width and alphabet checks. Interior whitespace remains content for validation and is normally an invalid character.

## Sequence transforms

- Reverse-only: nucleotide and protein.
- Complement-only: nucleotide-classified requests.
- Reverse-complement: nucleotide-classified requests.
- Strand-aware minus BED: nucleotide-classified requests.

## Output

- GET writes FASTA only.
- GET wraps at 60 symbols per line.
- Index `.fai` writes stdout only.
- Stats writes one human-readable text report.
- Validate writes human text, JSON Lines, or one JSON summary.

## Index limits

`.zfi` record counts are stored as `u32` and cannot be zero. Names use `u16` lengths. Side-table absolute offsets use 40 bits. Other persisted offsets and lengths use `u64`, with range checks against actual files.

These format widths are implementation limits, not recommended operating targets.

## Related pages

- [Index formats](Index-Formats)
- [Validation and repair](Validation-and-Repair)
- [Troubleshooting](Troubleshooting)
