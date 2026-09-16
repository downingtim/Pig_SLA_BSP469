#!/usr/bin/env bash
set -euo pipefail

REF_DIR="sla_reference"
mkdir -p "${REF_DIR}"
SLA_OUT="${REF_DIR}/ipd_mhc_sla_cds.fasta"
RAW_OUT="${REF_DIR}/ipd_mhc_raw.fasta"

rm -rf "${SLA_OUT}" "${RAW_OUT}"

echo "=== Downloading IPD-MHC Database ==="
wget -q -O "${RAW_OUT}" "https://raw.githubusercontent.com/ANHIG/IPDMHC/master/MHC_nuc.fasta"

echo "=== Filtering Classical Loci and Formatting Headers ==="
awk '/^>/ {
    is_sla = (tolower($0) ~ /sla/)
    is_nonclass = ($0 ~ /SLA-6/ || $0 ~ /SLA-7/ || $0 ~ /SLA-8/)
    if (is_sla && !is_nonclass) {
        keep = 1
        header = $0
        gsub(/ /, "_", header)
        print header
    } else {
        keep = 0
    }
    next
}
keep { print }' "${RAW_OUT}" > "${SLA_OUT}"

rm -f "${RAW_OUT}"

SEQ_COUNT=$(grep -c '^>' "${SLA_OUT}" || echo 0)
echo "SUCCESS: Extracted ${SEQ_COUNT} classical SLA sequences into '${SLA_OUT}'"
