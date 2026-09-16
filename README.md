### Pig SLA BSP469

# Swine Leukocyte Antigen (SLA) Amplicon Sequencing Pipeline

## Overview
This repository contains a complete bioinformatics pipeline for processing, aligning, filtering, and visualizing Swine Leukocyte Antigen (SLA) targeted amplicon sequencing data from Illumina paired-end platforms. 

Because SLA genes are highly homologous (very similar to one another), standard mapping tools often struggle to assign reads accurately. This pipeline addresses this by merging paired reads, removing identical duplicates to speed up computation, utilizing exact BLASTn local alignments, and enforcing strict thresholding to resolve multi-mapping ties and filter out sequencing noise (such as index hopping).

## Directory Structure
The workspace is organized into the following directories:

*   `DATA/`: Contains the reference databases, including the filtered classical SLA coding sequences (`ipd_mhc_sla_cds.fasta`).
*   `ILLUMINA_FASTQS/`: (Required) The expected location for your input raw `.fastq.gz` or `.fastq` sequencing files.
*   `TABLES/`: Stores the final output count matrices (e.g., `Filtered_SLA_Allele_Matrix.185.csv`).
*   `FIGURES/`: Stores all generated data visualizations (stacked bar plots and heatmaps).
*   `scripts/`: Contains the bash, python, and R scripts that make up the pipeline.

---

## Dependencies
To execute this pipeline, the following software must be installed and available in your environment (`$PATH` or module load system):
*   **QC & Trimming:** `fastqc`, `fastp`, `cutadapt`, `fastq_quality_trimmer` (FASTX-Toolkit)
*   **Merging & Dereplication:** `vsearch`
*   **Alignment:** NCBI `blast+` (`makeblastdb`, `blastn`)
*   **Python (Script 04):** `python3` with `pandas`
*   **R (Script 05):** `R` with `tidyverse` and `viridis` packages

---

## Pipeline Step-by-Step

### `00_runqc.sh` : Quality Control and Read Trimming
**What it does:** This script takes raw paired-end sequence files and runs them through a rigorous cleaning process. 
**How it works:** 
1. Generates initial QC reports via `FastQC`.
2. Uses `fastp` to trim adapter sequences, poly-G, and poly-X tails, while enforcing a minimum mean quality score of 20 and dropping reads shorter than 30bp.
3. Uses `cutadapt` to strip artificial poly-A and poly-T tails.
4. Performs a final pass with `fastq_quality_trimmer` to ensure the trailing ends of the reads are high quality.
*Note: This script takes the file prefix as an argument (e.g., `./00_runqc.sh ST-8B_S2_L001`).*

### `01_download.sh` : Reference Database Construction
**What it does:** Fetches the most up-to-date Swine MHC database and filters it for relevant targets.
**How it works:** Downloads the raw `MHC_nuc.fasta` from the IPD-MHC GitHub repository. It then uses `awk` to extract only the sequences belonging to the SLA loci, deliberately excluding non-classical loci (SLA-6, SLA-7, SLA-8) to prevent reads from being misassigned to pseudogenes or non-classical variants.

### `02_vsearch.sh` : Read Merging and Dereplication
**What it does:** Prepares the reads for alignment by combining pairs and collapsing duplicates.
**How it works:** 
1. **Merging:** It loops through `ILLUMINA_FASTQS/` and uses `fastp` to stitch the Forward (R1) and Reverse (R2) reads into a single, longer contiguous sequence. This requires at least a 15bp overlap. Longer reads are crucial for distinguishing between highly similar SLA alleles.
2. **Dereplication:** It uses `vsearch` to find all identical reads, collapses them into a single sequence, and adds a `size=X` tag to the FASTA header (where X is the number of times that read appeared). This vastly reduces the number of sequences BLAST has to process, saving massive amounts of compute time.

### `03_blastn.sh` : Sequence Alignment
**What it does:** Finds the best matching SLA allele for every unique read.
**How it works:** Builds a local BLAST database from the reference FASTA generated in step 01. It then aligns every dereplicated sequence against this database. It outputs a tabular text file (`outfmt 6`) containing the top 5 alignments for each read, recording critical metrics like percent identity, alignment length, and bitscore.

### `04_processblast.py` : Stringent Filtering and Matrix Generation
**What it does:** Translates the raw BLAST alignments into a reliable count matrix by applying biological and technical thresholds.
**How it works:**
1. **Length & Identity Filter:** Drops any alignments shorter than 185 base pairs or with less than 99.0% sequence identity.
2. **Multi-Mapper Resolution (Tie-Breaking):** If a read aligns equally well to multiple SLA alleles, it is ambiguous. The script checks the bitscore of the best hit versus the second-best hit. The top hit must win by a margin of `>= 2.0` bitscore points. If it doesn't, the read is discarded as unresolvable.
3. **Counting:** Re-expands the data using the `size=X` tags from step 02 and sums the successful reads per allele, per sample, outputting `Filtered_SLA_Allele_Matrix.185.csv`.

### `05_blast_processing.R` : Normalization and Visualization
**What it does:** Cleans up the final matrix and generates publication-ready plots.
**How it works:**
1. **Depth Filter:** Drops entirely failed samples (any sample with fewer than 200 total successfully mapped SLA reads).
2. **Noise Masking:** Calculates the relative abundance (percentage) of each allele within a sample. To filter out technical artifacts like index hopping or minor PCR cross-contamination, any allele making up less than `1.0%` of a sample's total reads is forced to `0`.
3. **Plotting:** Generates normalized stacked bar plots (Composition) and Heatmaps for the entire dataset. It also separates the data into `Pigs` (in vivo subjects) and `ST` (cohort/in vitro subjects) for focused visualizations, sorting the `ST` samples chronologically/by dominant haplotype (`14:02`).

---

## How to Execute
1. Place raw fastq files in `ILLUMINA_FASTQS/`.
2. Ensure you have executed `chmod +x scripts/*.sh`.
3. Run the scripts sequentially from the project root directory:
   ```bash
   ./scripts/00_runqc.sh <sample_prefix> # Run in a loop for all samples
   ./scripts/01_download.sh
   ./scripts/02_vsearch.sh
   ./scripts/03_blastn.sh
   python3 scripts/04_processblast.py
   Rscript scripts/05_blast_processing.R


Methods

Quality control of amplicon read libraries

FastQC v0.12.1 (www.bioinformatics.babraham.ac.uk/projects/fastqc/) was used to assess sequence quality. Adapter trimming, removal of reads with low base quality (BQ) scores (phred score <30), exclusion of ambiguous (N) bases, mismatched base pair correction in overlapping regions, a sliding-window quality filter, and cutting of poly-G and poly-X tracts at 3′ ends was completed using Fastp v0.23.4 (Chen et al 2018), including removing the first 11 bp that had lower BQ scores. Homopolymer tracts of ten bases or more were removed with Cutadapt v2.8 (Martin 2011). Low-quality bases with BQ <30 and reads with lengths <120 bp were removed using the FASTX-Toolkit v0.0.13 (http://hannonlab.cshl.edu/fastx_toolkit/). MultiQC v1.14 (Ewels et al 2016) was used to assess the effectiveness of Fastp and the FASTX-Toolkit by collating FastQC reports (Table S1).

Read alignment and SLA genotyping

Porcine SLA reference sequence data was retrieved from the IPD-MHC database (Maccari et al 2024). This was deduplicated and indexed with BWA v0.7.17 (REF) and SAMtools v1.18 (Li et al 2009). 12 non-classical MHC class 1 alleles such as SLA08517 and SLA09770 were excluded. This resulted in 1,009 sequences spanning 421,162 bp.

The median amplicon size was 480 bp. PE reads with an overlap of at least 15 bp and an identity >90% were merged using fastp. Trimmed reads were assembled into full amplicon-length sequences using Vsearch (REF). Samples with fewer than 200 reads were excluded. Merged read pairs were dereplicated using Vsearch and were collapsed into identical consensus sequences but retaining read depth abundance data. The merged reads were aligned against the IPD-MHC SLA reference CDS using BLASTn v2.9.0 (REF) with a maximum of five hits per query. To eliminate spurious alignment artifacts, PCR chimeras, and low-quality misalignments, alignment outputs were processed as follows. Firstly, hits were retained if they had a length >184 bp and with a sequence identity >99%. Secondly, to eliminate multi-mapping and ambiguous hits, the top two hits were compared based on their bitscores such that hits were assigned to the top allele only if the difference in bitscores was >2. Thirdly, read abundance was considered for each SLA allele to allow its estimation for each allele. Fourthly, matching SLA hits required a read abundance of >1% and presence in more than one sample. This retained five valid pig samples and 33 valid fly samples.




 
References

Maccari G, Robinson J, Barker DJ, Yates AD, Hammond JA, Marsh SGE. “The 2024 IPD-MHC database update: a comprehensive resource for major histocompatibility complex studies.” Nucleic Acids Res. (2025), 53: D457-D461.

Chen S, Zhou Y, Chen Y, Gu J. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics. 2018 34(17):i884–i890. doi: 10.1093/bioinformatics/bty560

Martin M. Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal. 2011;17(1):10–12.

Ewels P, et al. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics. 2016;32(19):3047–3048.

	
