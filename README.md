### Pig SLA BSP469

# QC

The typical amplicon size was about 480 bp. Fastp, etc

# Read library processing and SLA genotyping

Porcine SLA reference sequence data was retrieved from the IPD-MHC database (Maccari et al 2024). This wa deduplicated and indexed with BWA v0.7.17 (REF) and SAMtools v1.18 (REF). This resulted in 1,022 sequences spanning 1,666,057 bp.

Raw reads were quality-filtered using fastp v0.20 (REF). PE reads with an overlap of at least 15 bp with a identity of 90% or higher were merged using fastp. Trimmed reads were assembled into full amplicon-length sequences using Vsearch vX (REF). Merged read pairs were dereplicated using Vsearch and were collapsed into identical consensus sequences but retaining read depth abundance data. The merged reads were aligned against the IPD-MHC SLA reference CDS using BLASTn vX (REF) with a maximum of five hits per query. To eliminate spurious alignment artifacts, PCR chimeras, and low-quality misalignments, alignment outputs were processed as follows. Firstly, hits were retained if they had a length > 400 bp and with a sequence identity of at least 98%. Secondly, to eliminate multi-mapping and ambiguous hits, the top two hits were compared based on their bitscores such that hits were assigned to the top allele only if the difference in bitscores was at least 2. Thirdly, read abundance was considered for each SLA allele to alow abundance estimation.


# References:

Maccari G, Robinson J, Barker DJ, Yates AD, Hammond JA, Marsh SGE. “The 2024 IPD-MHC database update: a compre
hensive resource for major histocompatibility complex studies.” Nucleic Acids Res. (2025), 53: D457–D461.


