# Choosing a FASTA Tool

No FASTA tool is the right choice for every job. I start with the input format and output contract, then choose the smallest tool that fits them. z-fasta focuses on uncompressed FASTA indexing, exact random access, validation or normalization, and detailed nucleotide or protein statistics in one static binary. Other tools make better starting points for compressed data, broader interval formats, FASTQ, or language-specific APIs.

This comparison reflects official documentation reviewed on 2026-08-08. Check the linked tool documentation for current releases.

| Need | My starting point | Why |
| --- | --- | --- |
| Fast z-fasta-only random access, including messy FASTA | z-fasta `.zfi` | Embedded names, source identity, and non-uniform side tables |
| Standard `.fai` interoperability | samtools faidx or `z-fasta index --emit-fai` | Established five-column FAI ecosystem |
| BED/GFF/VCF extraction or BED12 block concatenation | bedtools getfasta | Interval-format breadth and block-aware extraction |
| Broad FASTA/Q manipulation and compressed inputs | SeqKit | Large command set, compression, FASTQ, and pipelines |
| Python random-access API and sequence objects | pyfaidx | Dictionary and slicing interface in Python |
| Minimal FASTA/Q stream transformations | seqtk | Small command-line toolkit and stream-oriented workflows |

## z-fasta and samtools faidx

Shared ground:

- first duplicate name resolves to the first sequence;
- `.fai` random access;
- positional sequence or region extraction;
- reverse-complement output;
- FASTA output to stdout.

For me, samtools is the better choice when I need BGZF, FASTQ indexing, compressed output, configurable output wrapping, or its wider HTSlib ecosystem. I reach for z-fasta when I need `.zfi` support for non-uniform FASTA, built-in validation and repair, or its indexed stats report.

Official reference: [samtools faidx manual](https://www.htslib.org/doc/samtools-faidx.html) and [faidx format manual](https://www.htslib.org/doc/faidx.html).

## z-fasta and bedtools getfasta

Shared ground:

- BED interval extraction;
- optional strand-aware reverse-complement;
- FASTA output.

bedtools wins when the input is GFF, VCF, or BED12 blocks, or when I need name-only, tabular, BED-output, full-header, or RNA-specific output controls. I keep z-fasta for one indexed workflow spanning positional, names, BED, transforms, summary timing, and messy `.zfi` records.

Official reference: [bedtools getfasta documentation](https://bedtools.readthedocs.io/en/latest/content/tools/getfasta.html).

## z-fasta and SeqKit

SeqKit covers FASTA and FASTQ, compressed inputs, many transformations, configurable identifier parsing, threads, and a large command catalog. Its documentation organizes commands by task and shows complete pipelines.

For a general sequence toolbox, I start with SeqKit. I use z-fasta when the work is specifically uncompressed FASTA indexing, indexed extraction, strict format diagnosis, normalization, or the z-fasta statistics contract.

SeqKit's `.seqkit.fai` may use full headers as keys and is not identical to ordinary samtools `.fai` behavior. Do not interchange sidecars without checking the consumer.

Official reference: [SeqKit usage guide](https://bioinf.shenwei.me/seqkit/usage/).

## z-fasta and pyfaidx

pyfaidx exposes FASTA records as Python objects with slicing, complements, reverse operations, configurable name handling, read-ahead, and a command-line interface.

If the workflow lives inside Python or needs its object model, pyfaidx is the natural fit. I use z-fasta for a standalone static command, `.zfi` side tables, validation or repair, and bounded CLI request streaming.

Coordinate syntax differs by interface. pyfaidx Python slicing is 0-based, while z-fasta positional CLI coordinates are 1-based inclusive.

Official reference: [pyfaidx documentation](https://github.com/mdshw5/pyfaidx).

## z-fasta and seqtk

seqtk is a compact FASTA/Q stream toolbox. It covers common transformations, subsequence selection, sampling, and FASTQ conversion without requiring an index for every job.

I choose seqtk for small stream-oriented transformations and z-fasta when I need reusable indexed access, explicit coordinate contracts, validation, repair, or indexed statistics.

Official reference: [seqtk documentation](https://github.com/lh3/seqtk).

## A practical rule

My practical rule is simple: do not choose by benchmark rank alone. Start with the tool whose input formats, coordinate model, output contract, and failure behavior fit the pipeline. Then measure the real workload.

> [!NOTE]
> z-fasta intentionally does not replace every FASTA/Q utility. Its narrower boundary keeps indexing, extraction, validation, and stats behavior easier to verify together.

## Related pages

- [Index formats](Index-Formats)
- [Limits and supported formats](Limits-and-Supported-Formats)
- [Performance and correctness](Performance-and-Correctness)
