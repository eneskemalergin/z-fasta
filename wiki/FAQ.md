# FAQ

Short answers first. Follow the linked concept pages when the reason affects a pipeline decision.

## Does GET create an index automatically?

No. Run `z-fasta index <file>` first. Automatic index creation is outside the current command contract.

## Should I use `.zfi` or `.fai`?

I recommend `.zfi` for z-fasta. Use `.fai` when another tool requires the standard text format and the FASTA has fixed line geometry. See [Index formats](Index-Formats).

## Why does `.fai` work in another tool but z-fasta still fails?

A present `.zfi` is authoritative. If it is stale or corrupt, z-fasta does not fall back to `.fai`.

## Are regions 0-based or 1-based?

Positional GET is 1-based inclusive. BED is 0-based half-open. Output FASTA headers are 1-based inclusive.

## Why does the output header keep an end beyond sequence length?

GET clamps the extracted bytes but preserves the requested end in the header for samtools-faidx compatibility.

## Can sequence names contain colons?

Yes. Only a valid decimal range after the rightmost colon is treated as coordinates. Use `--names` when a literal identifier itself ends in range syntax.

## Are header descriptions preserved?

No. Index lookup and GET output use the identifier before the first space or tab.

## What happens to duplicate names?

Default indexing keeps the first non-empty record. `--no-dedup` retains all records for index population and stats, but GET still returns the first exact match. Validate reports duplicates as errors.

## What happens to empty sequences?

Validate warns. Index skips them. Stats and GET therefore operate on a population that excludes empty records.

## Can z-fasta index a file with mixed line widths?

Yes through default `.zfi` side tables. `.fai` emission rejects it.

## Can validate fix mixed widths?

Yes when no blocking error prevents repair. Each record is rewrapped to its modal observed non-empty width.

## Does fix change biological sequence bytes?

It removes spaces and tabs only from sequence-line ends, joins sequence content, and rewraps it. It does not replace invalid characters or invent bases. Use a strict post-fix validation as the gate.

## Can fix overwrite the input?

No. It requires `-o` and refuses a destination that resolves to the input.

## Does z-fasta support FASTQ?

No. FASTQ is coming as a separate project, z-fastq, built around the same priorities: speed, low memory use, and portability.

## Does z-fasta support compressed FASTA?

No. Current input is uncompressed FASTA only. I have not invested the time required to optimize deflate or write a separate decompressor for this project.

## Does z-fasta publish Windows binaries?

No. I would like to support Windows properly, but keeping native builds reliable currently takes more time than I can justify. For now, use the Linux release through Windows Subsystem for Linux. I am sorry for the extra step.

## Can I extract BED12 blocks?

No. GET uses the BED chromosome, start, end, and optional strand. It does not concatenate BED12 blocks.

## Why is reverse-complement rejected for my file?

The requested sample classified as protein. Complement-based transforms are nucleotide-only. Reverse-only is still available.

## How is sequence type detected?

Case-insensitive IUPAC nucleotide symbols must be strictly more than 90 percent of counted symbols. Stats uses the full indexed population, validate samples the first 100000 sequence bytes, and GET samples up to 256 requested bases per record for complement eligibility.

## Can names or BED come from stdin?

Yes, use `--names -` or `--bed -`.

## Is streamed output all-or-nothing?

No. GET can write earlier valid records before a later names or BED request fails. Index `.fai` emission uses a spool and does not replay a partial scan.

## Where does the GET summary go?

Stderr after successful stdout output. It is suppressed on failure.

## How do I report a problem safely?

Include version, platform, command, status, stdout, stderr, sidecar names, and a minimized non-sensitive fixture. Never attach private biological data to a public report.
