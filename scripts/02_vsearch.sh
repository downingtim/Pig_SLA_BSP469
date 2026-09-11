#!/bin/bash

# Create an output directory for the merged and dereplicated files
mkdir -p vsearch_outputs

# Loop through all Forward (R1) fastp-trimmed files
for R1 in ILLUMINA_FASTQS/*_1_fastp_trim.fastq; do
    
    # Extract just the filename without the path
    filename=$(basename "$R1")
    
    # Strip the suffix to get the core sample ID (e.g., Pig-29_S56_L001)
    SAMPLE="${filename%_1_fastp_trim.fastq}"
    
    # Define the corresponding Reverse (R2) file
    R2="ILLUMINA_FASTQS/${SAMPLE}_2_fastp_trim.fastq"
    
    # Define the outputs
    MERGED="vsearch_outputs/${SAMPLE}_merged.fasta"
    DEREP="vsearch_outputs/${SAMPLE}_derep.fasta"
    
    echo "========================================"
    echo "Processing $SAMPLE"
    
    # 1. Merge Forward and Reverse reads (Note the --reverse flag)
    vsearch --fastq_mergepairs "$R1" \
            --reverse "$R2" \
            --fastaout "$MERGED"
            
    # 2. Dereplicate: Collapse identical sequences and append ';size=X' to the header
    vsearch --derep_fulllength "$MERGED" \
            --output "$DEREP" \
            --sizeout

done
echo "Merging and dereplication complete. Files saved in vsearch_outputs/"
