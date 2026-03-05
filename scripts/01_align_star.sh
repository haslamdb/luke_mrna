#!/bin/bash
# STAR alignment for Luke's CD4 T cell RNA-seq
# Paired-end, 101bp reads -> mouse genome (Ensembl)
#
# Usage: bash scripts/01_align_star.sh [threads]

set -euo pipefail

THREADS=${1:-24}
GENOME_DIR="/bulkpool/reference_data/Mus_musculus_STAR_101bp"
DATA_DIR="/bulkpool/sequence_data/rnaseq_data/luke_mRNA"
OUT_DIR="/home/david/projects/luke_mrna/aligned"
METADATA="/home/david/projects/luke_mrna/sample_metadata.csv"

# Verify index exists
if [[ ! -f "$GENOME_DIR/SA" ]]; then
    echo "ERROR: STAR index not found at $GENOME_DIR"
    exit 1
fi

# Read sample metadata (skip header)
tail -n +2 "$METADATA" | while IFS=',' read -r sample_id fastq_prefix cell_type mouse sort_date; do
    R1="${DATA_DIR}/${fastq_prefix}_R1_001.fastq.gz"
    R2="${DATA_DIR}/${fastq_prefix}_R2_001.fastq.gz"

    if [[ ! -f "$R1" ]]; then
        echo "ERROR: Missing R1 file: $R1"
        exit 1
    fi

    SAMPLE_OUT="${OUT_DIR}/${sample_id}/"
    mkdir -p "$SAMPLE_OUT"

    # Skip if already aligned
    if [[ -f "${SAMPLE_OUT}Aligned.sortedByCoord.out.bam" ]]; then
        echo "SKIP: ${sample_id} already aligned"
        continue
    fi

    echo "$(date '+%Y-%m-%d %H:%M:%S') ALIGNING: ${sample_id}"

    STAR \
        --runThreadN "$THREADS" \
        --genomeDir "$GENOME_DIR" \
        --readFilesIn "$R1" "$R2" \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix "$SAMPLE_OUT" \
        --quantMode GeneCounts \
        --outSAMattributes NH HI AS NM MD \
        --outFilterType BySJout \
        --outFilterMultimapNmax 20 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverReadLmax 0.04 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000

    # Index the BAM
    samtools index "${SAMPLE_OUT}Aligned.sortedByCoord.out.bam"

    echo "$(date '+%Y-%m-%d %H:%M:%S') DONE: ${sample_id}"
done

echo ""
echo "=== Alignment complete ==="
echo "Generating alignment summary..."

# Print summary stats from STAR logs
echo ""
printf "%-25s %12s %12s %12s\n" "Sample" "Total Reads" "Unique Map%" "Multi Map%"
printf "%-25s %12s %12s %12s\n" "-------" "-----------" "-----------" "----------"
for dir in "$OUT_DIR"/*/; do
    sample=$(basename "$dir")
    log="${dir}Log.final.out"
    if [[ -f "$log" ]]; then
        total=$(grep "Number of input reads" "$log" | awk '{print $NF}')
        uniq_pct=$(grep "Uniquely mapped reads %" "$log" | awk '{print $NF}')
        multi_pct=$(grep "% of reads mapped to multiple loci" "$log" | awk '{print $NF}')
        printf "%-25s %12s %12s %12s\n" "$sample" "$total" "$uniq_pct" "$multi_pct"
    fi
done
