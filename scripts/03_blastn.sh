#!/bin/bash

# Define paths
REF_FASTA="sla_reference/ipd_mhc_sla_cds.fasta"
INPUT_DIR="vsearch_outputs"
OUTPUT_DIR="blast_outputs"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# 1. Generate the BLAST database (only if it hasn't been built yet)
# We check for the '.nhr' file, which is one of the files makeblastdb creates
if [ ! -f "${REF_FASTA}.nhr" ]; then
    echo "BLAST database not found. Building it now from ${REF_FASTA}..."
    makeblastdb -in "$REF_FASTA" -dbtype nucl -parse_seqids
    echo "BLAST database built."
else
    echo "BLAST database already exists. Proceeding to alignment."
fi

echo "========================================"

# 2. Loop through all dereplicated FASTA files
for derep_file in "${INPUT_DIR}"/*_derep.fasta; do
    
    # Extract the base filename and sample name
    filename=$(basename "$derep_file")
    SAMPLE="${filename%_derep.fasta}"
    
    BLAST_OUT="${OUTPUT_DIR}/${SAMPLE}_blast.txt"
    
    echo "Running BLASTn for $SAMPLE..."
    
    # 3. Execute BLASTn
    blastn -query "$derep_file" \
           -db "$REF_FASTA" \
           -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
           -max_target_seqs 5 \
           -num_threads 8 \
           -out "$BLAST_OUT"

done

echo "========================================"
echo "All BLAST runs complete. Results are in the ${OUTPUT_DIR}/ directory."
