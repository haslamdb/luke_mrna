#!/bin/bash
# featureCounts quantification for Luke's CD4 T cell RNA-seq
# Counts reads per gene from STAR-aligned BAMs
#
# Usage: bash scripts/02_featurecounts.sh [threads]

set -euo pipefail

THREADS=${1:-24}
GTF="/bulkpool/reference_data/Mus_musculus/genes.gtf"
ALIGNED_DIR="/home/david/projects/luke_mrna/aligned"
OUT_DIR="/home/david/projects/luke_mrna/counts"

mkdir -p "$OUT_DIR"

# Collect all BAM files in metadata order
BAMS=()
for dir in "$ALIGNED_DIR"/*/; do
    bam="${dir}Aligned.sortedByCoord.out.bam"
    if [[ -f "$bam" ]]; then
        BAMS+=("$bam")
    fi
done

if [[ ${#BAMS[@]} -eq 0 ]]; then
    echo "ERROR: No BAM files found in $ALIGNED_DIR"
    exit 1
fi

echo "Found ${#BAMS[@]} BAM files"
echo "Running featureCounts..."

featureCounts \
    -T "$THREADS" \
    -p --countReadPairs \
    -s 0 \
    --extraAttributes gene_name \
    -a "$GTF" \
    -o "$OUT_DIR/gene_counts.txt" \
    "${BAMS[@]}"

# Note on -s (strandedness):
# -s 0 = unstranded
# -s 1 = stranded (sense)
# -s 2 = reversely stranded (most Illumina dUTP/TruSeq stranded)
# If results look off, re-run with -s 0 or -s 1
# STAR's ReadsPerGene.out.tab can help determine strandedness (see below)

echo ""
echo "=== featureCounts complete ==="
echo "Output: $OUT_DIR/gene_counts.txt"
echo "Summary: $OUT_DIR/gene_counts.txt.summary"
echo ""
echo "To check strandedness, compare columns 2-4 in any STAR ReadsPerGene.out.tab:"
echo "  col2 = unstranded, col3 = sense, col4 = antisense"
echo "  The column with the most counts indicates the correct -s setting"
