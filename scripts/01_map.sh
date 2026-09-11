#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# Complete SLA Amplicon Pipeline: fastp -> BWA -> samtools -> bcftools -> consensus
# ==============================================================================

RAW_DIR="ILLUMINA_FASTQS"
REF_FASTA="sla_reference/ipd_mhc_sla_cds.fasta"
OUT_DIR="sla_analysis_output"
THREADS=8

# Create all necessary subdirectories explicitly
mkdir -p "${OUT_DIR}/merged" "${OUT_DIR}/bam" "${OUT_DIR}/vcf" "${OUT_DIR}/fasta" "${OUT_DIR}/reports"

echo "Checking SLA Reference index..."
if [ ! -f "${REF_FASTA}.bwt" ]; then
    bwa index "${REF_FASTA}"
fi
if [ ! -f "${REF_FASTA}.fai" ]; then
    samtools faidx "${REF_FASTA}"
fi

echo "Starting sample processing loop..."

# Enable nullglob to prevent empty wildcard errors
shopt -s nullglob
FASTQ_FILES=("${RAW_DIR}"/*_1_fastp_trim.fastq.gz)

if [ ${#FASTQ_FILES[@]} -eq 0 ]; then
    echo "Error: No matching FASTQ files (*_1_fastp_trim.fastq.gz) found in ${RAW_DIR}"
    exit 1
fi

for r1 in "${FASTQ_FILES[@]}"; do
    sample=$(basename "$r1" | sed 's/_1_fastp_trim\.fastq\.gz//')
    r2="${RAW_DIR}/${sample}_2_fastp_trim.fastq.gz"

    if [ ! -f "$r2" ]; then
        echo "WARNING: R2 missing for ${sample}. Skipping..."
        continue
    fi

    echo ">>> Processing Sample: ${sample}"

    # 1. Merge paired-end reads with fastp
    fastp \
        -i "$r1" \
        -I "$r2" \
        --merge \
        --merged_out "${OUT_DIR}/merged/${sample}_merged.fastq.gz" \
        --out1 "${OUT_DIR}/merged/${sample}_unmerged_R1.fastq.gz" \
        --out2 "${OUT_DIR}/merged/${sample}_unmerged_R2.fastq.gz" \
        --overlap_len_require 15 \
        --overlap_diff_percent_limit 10 \
        --thread "${THREADS}" \
        --html "${OUT_DIR}/reports/${sample}_fastp.html" \
        --json "${OUT_DIR}/reports/${sample}_fastp.json" \
        2> "${OUT_DIR}/reports/${sample}_fastp.log"

    # 2. Align merged reads to reference and sort
    bwa mem -t "${THREADS}" \
        -R "@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA" \
        "${REF_FASTA}" \
        "${OUT_DIR}/merged/${sample}_merged.fastq.gz" | \
        samtools sort -@ "${THREADS}" -O BAM -o "${OUT_DIR}/bam/${sample}_sorted.bam" -

    samtools index "${OUT_DIR}/bam/${sample}_sorted.bam"

    # 3 & 4. Call variants (diploid by default) and filter into a single compressed VCF stream
    bcftools mpileup -Ou -f "${REF_FASTA}" -d 50000 -a FORMAT/AD,FORMAT/DP --threads "${THREADS}" "${OUT_DIR}/bam/${sample}_sorted.bam" | \
    bcftools call -m -v -Ou | \
    bcftools filter -e 'QUAL < 30 || FORMAT/DP < 20' -Oz -o "${OUT_DIR}/vcf/${sample}_filtered.vcf.gz"

    bcftools index "${OUT_DIR}/vcf/${sample}_filtered.vcf.gz"

    # 5. Generate Consensus FASTA from filtered VCF
    bcftools consensus \
        -f "${REF_FASTA}" \
        "${OUT_DIR}/vcf/${sample}_filtered.vcf.gz" \
        > "${OUT_DIR}/fasta/${sample}_full_consensus.fasta"

    # 6. Identify the specific target gene and extract it
    # Use samtools coverage to find which reference sequence had the most reads mapped
    TARGET_ID=$(samtools coverage "${OUT_DIR}/bam/${sample}_sorted.bam" | \
                awk '$4 > 0 {print $0}' | \
                sort -k4,4nr | \
                head -n 1 | \
                cut -f1)

    if [ -n "$TARGET_ID" ]; then
        # Index the full consensus FASTA so we can extract specific sequences
        samtools faidx "${OUT_DIR}/fasta/${sample}_full_consensus.fasta"
        
        # Extract ONLY the target sequence that was actually amplified
        samtools faidx "${OUT_DIR}/fasta/${sample}_full_consensus.fasta" "$TARGET_ID" > "${OUT_DIR}/fasta/${sample}_SLA_target.fasta"
        echo "    Extracted target sequence: ${TARGET_ID}"
    else
        echo "    WARNING: No mapped reads found for ${sample}. Skipping extraction."
    fi

    echo "    Successfully completed: ${sample}"
done

echo "=================================================="
echo "Pipeline execution finished successfully for all samples!"
echo "Check output files in: ${OUT_DIR}/fasta/"
