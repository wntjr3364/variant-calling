# Variant Calling Pipeline

## Overview
This repository contains a Snakemake-based workflow for genomic variant calling, including quality control, alignment, and variant detection. The pipeline is designed to automatically detect paired-end read files in the specified folder and process them sequentially.

## Features
- **Automatic Input Detection**: Detects all paired-end read files in the `trim_data/` directory.
- **End-to-End Automation**: Includes all steps from quality control to variant calling.
- **Reproducibility**: Ensures consistent results with a fully automated workflow.
- **Parallel Processing**: Leverages multiple CPU cores for faster processing.

## Workflow
The pipeline consists of the following steps:
1. **Quality Control**: Run FastQC to generate quality reports for raw reads.
2. **Reference Genome Preparation**: Index the reference genome for alignment and variant calling.
3. **Read Alignment**: Align reads to the reference genome using BWA-MEM.
4. **Post-alignment Processing**: Mark duplicates and sort BAM files using Picard.
5. **Base Quality Recalibration**: Recalibrate base quality scores with GATK.
6. **Variant Calling**: Detect SNPs and InDels using GATK HaplotypeCaller.

## Dependencies
* FastQC: Quality control of raw reads.
* Trimmomatic: Quality trimming of raw reads.
* SAMtools: BAM file handling and depth calculation.
* BWA: Alignment of reads to the reference genome.
* Picard: Post-alignment processing.
* GATK: Variant calling and base recalibration.

## Installation
1. Install Snakemake and other dependencies:
   ```bash
   conda install -c bioconda snakemake bwa trimmomatic picard gatk samtools fastqc
2. Clone the repository:
    ```bash
    git clone https://github.com/username/variant-calling-pipeline.git
    cd variant-calling-pipeline
3. Ensure that your reference genome and read files are in the appropriate directories:
    * Reference genome: Place in the reference/ directory.
    * Read files: Place in the trim_data/ directory.


## Usage
1. Place your FASTQ files in the trim_data/ directory:
    * Single-end format: sample.fq.gz
    * Paired-end format: sample_1P.fq.gz and sample_2P.fq.gz
2. Update config.yaml with the reference genome and known variant sites.
    ```bash
    reference: "reference/Gmax_275_v2.0.Mt.Pltd.fa"
    known_sites: "reference/known_sites.vcf"
3. Run the workflow
    ```bash
    snakemake --cores 8
    To generate specific ouputs, specify the target:
    ```bash
    snakemake variants/SRR12345.vcf --cores 8


## Output
* Aligned BAM files: Stored in aligned/.
* Recalibrated BAM files: Stored in bqsr/.
* VCF files: Final variant calls saved in variants/.

## Directory Structure
The expected directory structure after running the pipeline is as follows:
```plaintext
.
├── Snakefile           # Workflow definition
├── config.yaml         # Configuration file
├── reference/          # Reference genome and index files
│   ├── Gmax_275_v2.0.Mt.Pltd.fa
│   ├── known_sites.vcf
├── trim_data/          # Input read files
│   ├── SRR12345_1P.fq.gz
│   ├── SRR12345_2P.fq.gz
├── qc_reports/         # FastQC reports
├── aligned/            # Aligned BAM files
├── bqsr/               # Recalibrated BAM files
├── variants/           # VCF files
├── log/                # Logs for each step
