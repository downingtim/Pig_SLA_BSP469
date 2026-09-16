#!/bin/bash
set -uo pipefail

mkdir -p vsearch_outputs
shopt -s nullglob

FASTQ_FILES=(ILLUMINA_FASTQS/*_1_fastp_trim.fastq.gz)
if [ ${#FASTQ_FILES[@]} -eq 0 ]; then
    FASTQ_FILES=(ILLUMINA_FASTQS/*_1_fastp_trim.fastq)
fi

if [ ${#FASTQ_FILES[@]} -eq 0 ]; then
    echo "ERROR: No trimmed FASTQ files found in ILLUMINA_FASTQS/"
    exit 1
fi

echo "Found ${#FASTQ_FILES[@]} samples to process."

for R1 in "${FASTQ_FILES[@]}"; do
    filename=$(basename "$R1")
    SAMPLE=$(echo "$filename" | sed -E 's/_1_fastp_trim\.fastq(\.gz)?$//')
    
    R2_GZ="ILLUMINA_FASTQS/${SAMPLE}_2_fastp_trim.fastq.gz"
    R2_PLAIN="ILLUMINA_FASTQS/${SAMPLE}_2_fastp_trim.fastq"
    
    if [ -f "$R2_GZ" ]; then
        R2="$R2_GZ"
    elif [ -f "$R2_PLAIN" ]; then
        R2="$R2_PLAIN"
    else
        echo "WARNING: Reverse read for $SAMPLE not found. Skipping..."
        continue
    fi

    MERGED_FASTQ="vsearch_outputs/${SAMPLE}_merged.fastq.gz"
    UNMERGED1="vsearch_outputs/${SAMPLE}_unmerged_R1.fastq.gz"
    UNMERGED2="vsearch_outputs/${SAMPLE}_unmerged_R2.fastq.gz"
    DEREP="vsearch_outputs/${SAMPLE}_derep.fasta"

    echo "========================================"
    echo "Processing sample: $SAMPLE"

    # 1. Merge paired reads with fastp
    if ! fastp -i "$R1" -I "$R2" \
         --merge \
         --merged_out "$MERGED_FASTQ" \
         --out1 "$UNMERGED1" \
         --out2 "$UNMERGED2" \
         --overlap_len_require 15 \
         --overlap_diff_percent_limit 10 \
         --html "vsearch_outputs/${SAMPLE}_fastp.html" \
         --json "vsearch_outputs/${SAMPLE}_fastp.json" \
         2> "vsearch_outputs/${SAMPLE}_fastp_merge.log"; then
        echo "ERROR: fastp failed for $SAMPLE. Check vsearch_outputs/${SAMPLE}_fastp_merge.log"
        continue
    fi

    # Clean up unmerged temporary files
    rm -f "$UNMERGED1" "$UNMERGED2"

    if [ ! -s "$MERGED_FASTQ" ]; then
        echo "WARNING: No reads merged for $SAMPLE. Skipping dereplication..."
        continue
    fi

    # 2. Dereplicate FASTQ input directly to FASTA output using fastx_uniques
    if ! vsearch --fastx_uniques "$MERGED_FASTQ" --fastaout "$DEREP" --sizeout; then
        echo "ERROR: vsearch failed for $SAMPLE."
        continue
    fi
done

echo "========================================"
echo "Merge and dereplication complete for all samples."
