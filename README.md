# Genomic Variant Calling Pipeline

## Overview

This Snakemake-based workflow automates the process of genomic variant calling, providing a comprehensive pipeline for processing paired-end sequencing data. The workflow includes quality control, read alignment, post-alignment processing, and variant detection.

## Prerequisites

### Software Dependencies

- [Snakemake](https://snakemake.readthedocs.io/)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
- [SAMtools](http://www.htslib.org/)
- [BWA](http://bio-bwa.sourceforge.net/)
- [Picard Tools](https://broadinstitute.github.io/picard/)
- [GATK](https://gatk.broadinstitute.org/)

### System Requirements

- Linux/Unix environment
- Conda (recommended for dependency management)
- Minimum 16 GB RAM
- Multi-core processor (recommended)

## Installation

1. Install Conda (if not already installed)
   ```bash
   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
   bash Miniconda3-latest-Linux-x86_64.sh
   ```

2. Create and activate a Conda environment
   ```bash
   conda create -n variant-calling
   conda activate variant-calling
   ```

3. Install dependencies
   ```bash
   conda install -c bioconda snakemake bwa trimmomatic picard gatk4 samtools fastqc
   ```

4. Clone the repository
   ```bash
   git clone https://github.com/username/variant-calling-pipeline.git
   cd variant-calling-pipeline
   ```

## Preparation

### Reference Genome and Known Sites

1. Place your reference genome in the `reference/` directory
   - Example: `reference/Gmax_275_v2.0.Mt.Pltd.fa`

2. Prepare known variant sites VCF file
   - Example: `reference/known_sites.vcf`

### Input Data

Place your sequencing data in the `read/` directory:
- Paired-end files: `sample_1P.fq.gz` and `sample_2P.fq.gz`
- Single-end files: `sample.fq.gz`

## Configuration

Update `config.yaml` with your specific settings:

```yaml
reference: "reference/Gmax_275_v2.0.Mt.Pltd.fa"
known_sites: "reference/known_sites.vcf"
```


## Workflow Steps

1. Quality Control (FastQC)
2. Read Trimming (Trimmomatic)
3. Reference Genome Indexing
4. Read Alignment (BWA-MEM)
5. Duplicate Marking (Picard)
6. Base Quality Score Recalibration (GATK)
7. Variant Calling (GATK)
8. Joint Variant Calling (GATK)

## Usage

### Run Full Pipeline

```bash
snakemake --cores 8
```

### Run Specific Targets

```bash
# Generate alignment for a specific sample
snakemake result/alignment/sample.bam --cores 8

# Generate variants for a specific sample
snakemake result/variants/sample.vcf --cores 8
```

## Output Directory Structure

```
variant-calling-pipeline/
├── reference/          # Reference genome files
├── read/               # Input sequencing data
├── log/
│   ├── trim/           # Trimming logs
│   ├── alignment/      # Alignment logs
│   └── variants/       # Variant calling logs
└── result/
    ├── trim/           # Trimmed reads
    ├── alignment/      # Aligned BAM files
    └── variants/       # Variant call VCF files
