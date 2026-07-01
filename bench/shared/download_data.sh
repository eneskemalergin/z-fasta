#!/usr/bin/env bash
#
# Download real datasets for benchmarking z-fasta
# Total size: ~4.1 GB uncompressed
#
# TODO(v0.3-bench): Emit a dataset manifest with source URLs, expected
# filenames, byte sizes, and checksums. The zebrac baseline layer should be
# able to record exactly which dataset revision produced a run.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="$SCRIPT_DIR/data"

mkdir -p "$DATA_DIR"
cd "$DATA_DIR"

echo "=== Downloading Real Datasets ==="
echo "Target: $DATA_DIR"
echo

# Human genome GRCh38 primary assembly (~3.0 GB)
echo "[1/3] Human Genome (GRCh38) - ~900 MB compressed, ~3.0 GB uncompressed"
if [[ -f "REAL_Genome.fa" ]]; then
    echo "  Skipping: REAL_Genome.fa already exists"
else
    curl -L -o REAL_Genome.fa.gz \
        "https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    gunzip REAL_Genome.fa.gz
    echo "  Downloaded: REAL_Genome.fa ($(du -h REAL_Genome.fa | cut -f1))"
fi

# Human transcriptome (GENCODE v46) (~260 MB compressed, ~972 MB uncompressed)
echo "[2/3] Human Transcriptome (GENCODE v46) - ~260 MB compressed, ~972 MB uncompressed"
if [[ -f "REAL_Transcriptome.fa" ]]; then
    echo "  Skipping: REAL_Transcriptome.fa already exists"
else
    curl -L -o REAL_Transcriptome.fa.gz \
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.transcripts.fa.gz"
    gunzip REAL_Transcriptome.fa.gz
    echo "  Downloaded: REAL_Transcriptome.fa ($(du -h REAL_Transcriptome.fa | cut -f1))"
fi

# Human proteome (UniProt) - canonical + isoforms (~22 MB compressed, ~66 MB uncompressed)
echo "[3/3] Human Proteome (UniProt) - ~22 MB compressed, ~66 MB uncompressed"
if [[ -f "REAL_Proteome.fasta" ]]; then
    echo "  Skipping: REAL_Proteome.fasta already exists"
else
    curl -L -o REAL_Proteome.fasta.gz \
        "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz"
    gunzip REAL_Proteome.fasta.gz
    echo "  Downloaded: REAL_Proteome.fasta ($(du -h REAL_Proteome.fasta | cut -f1))"
fi

echo
echo "=== Download Complete ==="
echo
du -h REAL_* 2>/dev/null || echo "No files downloaded"
