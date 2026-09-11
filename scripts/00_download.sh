#!/usr/bin/env bash
set -euo pipefail

REF_DIR="sla_reference"
mkdir -p "${REF_DIR}"

TEMP_REPO="${REF_DIR}/IPDMHC_temp"
SLA_OUT="${REF_DIR}/ipd_mhc_sla_cds.fasta"

# Clean up previous failed attempts
rm -rf "${TEMP_REPO}" "${SLA_OUT}"

echo "=== Step 1: Cloning IPD-MHC Repository ==="
if ! command -v git &> /dev/null; then
    echo "Error: 'git' command not found. Please install git (e.g., sudo apt install git)."
    exit 1
fi

git clone --depth 1 https://github.com/ANHIG/IPDMHC.git "${TEMP_REPO}"

echo ""
echo "=== Step 2: Scanning All Repository Files for SLA Sequences ==="

# Recursively find all fasta/text files and extract entries with SLA in header
find "${TEMP_REPO}" -type f \( -name "*.fasta" -o -name "*.fa" -o -name "*.txt" \) | while read -r f; do
    if grep -q -i '^>.*SLA' "$f" 2>/dev/null; then
        echo "Extracting SLA entries from: $(basename "$f")"
        awk '/^>/ {keep = (tolower($0) ~ /sla/)} keep' "$f" >> "${SLA_OUT}"
    fi
done

# Remove temporary clone folder
rm -rf "${TEMP_REPO}"

echo ""
echo "=== Step 3: Verifying File Contents ==="

if [ ! -s "${SLA_OUT}" ]; then
    echo "ERROR: Extraction failed — output file is still empty."
    echo "Checking fallback source..."
    
    # Direct fallback download from EBI raw mirror
    wget -q -O "${SLA_OUT}.raw" "https://raw.githubusercontent.com/ANHIG/IPDMHC/master/MHC_nuc.fasta" || true
    if [ -s "${SLA_OUT}.raw" ]; then
        awk '/^>/ {keep = (tolower($0) ~ /sla/)} keep' "${SLA_OUT}.raw" > "${SLA_OUT}"
        rm -f "${SLA_OUT}.raw"
    fi
fi

SEQ_COUNT=$(grep -c '^>' "${SLA_OUT}" || echo 0)

if [ "${SEQ_COUNT}" -gt 0 ]; then
    echo "SUCCESS: Extracted ${SEQ_COUNT} SLA sequences into '${SLA_OUT}'"
    echo ""
    echo "Header Preview:"
    grep '^>' "${SLA_OUT}" | head -n 5
else
    echo "ERROR: Could not retrieve SLA sequences. Please check your network connection."
    exit 1
fi
