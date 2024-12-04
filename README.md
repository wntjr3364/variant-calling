# Genomic Variant Calling Pipeline

## Overview
This Snakemake-based workflow provides an automated and comprehensive solution for genomic variant discovery, designed to process both single and paired-end sequencing data.


## Prerequisites

### Software Dependencies

- [Snakemake](https://snakemake.readthedocs.io/)
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
- Java 8 (JDK 1.8)

## Installation

1. Install Java 8
   ```bash
   # Ubuntu/Debian
   sudo apt-get update
   sudo apt-get install openjdk-8-jdk

   # CentOS/RHEL
   sudo yum install java-1.8.0-openjdk-devel

   # Verify Java installation
   java -version
   ```

2. Install Conda (if not already installed)
   ```bash
   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
   bash Miniconda3-latest-Linux-x86_64.sh
   ```

3. Create and activate a Conda environment
   ```bash
   conda create -n variant-calling
   conda activate variant-calling
   ```

4. Install dependencies via Conda
   ```bash
   conda install -c bioconda snakemake bwa trimmomatic picard samtools
   ```

5. Install GATK 3.7
   ```bash
   # Create program directory
   mkdir -p program
   cd program
   
   # Download GATK 3.7
   wget --no-check-certificate https://tgil.donga.ac.kr/CAPSMaker/GenomeAnalysisTK-nightly.tar.bz2
   
   # Extract the downloaded file
   tar -xjf GenomeAnalysisTK-nightly.tar.bz2
   
   cd ..
   ```

6. Clone the repository
   ```bash
   git clone https://github.com/wntjr3364/variant-calling.git
   cd variant-calling
   ```

## Preparation

### Reference Genome and Known Sites

1. Place your reference genome in the `reference/` directory
   - Example: `reference/Gmax_275_v2.0.Mt.Pltd.fa`

2. Prepare known variant sites VCF file
   - Example: `reference/known_sites.vcf`

### Input Data

Place your sequencing data in the `read/` directory:
- Paired-end files: `read/sample_1.fastq.gz` and `read/sample_2.fastq.gz`
- Single-end files: `read/sample.fastq.gz`

## Configuration

Create `config.yaml`:

```yaml
reference: "reference/Gmax_275_v2.0.Mt.Pltd.fa"
known_sites: "reference/dbsnp.vcf"
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
# Trim reads
snakemake result/trim/{sample}_trimmed.fastq.gz --cores 8

# Generate alignment
snakemake result/alignment/sample.bam --cores 8

# Generate variants
snakemake result/genotype_gvcf/combined_genotyped.vcf.gz --cores 8
```

## Output Directory Structure

```
variant-calling-pipeline/
├── snakefile           # Snakemake workflow definition
├── config.yaml         # Pipeline configuration file
├── reference/          # Reference genome files
├── read/               # Input sequencing data
├── log/                # Workflow log files
│   ├── trim/           # Trimming process logs
│   ├── alignment/      # BWA alignment logs
│   ├── read_group/     # Read group assignment logs
│   ├── mark_duplication/  # Duplicate marking logs
│   ├── indel_realignment/ # Indel realignment logs
│   ├── base_recalibration/ # Base quality recalibration logs
│   ├── haplotype_caller/  # Haplotype caller logs
│   ├── combine_gvcf/     # Combined gVCF logs
│   └── genotype_gvcf/    # Final genotyped VCF logs
└── result/             # Output files
    ├── trim/           # Trimmed read files
    ├── alignment/      # Alignment BAM files
    ├── read_group/     # Read group-annotated BAMs
    ├── mark_duplication/  # Deduplicated BAM files
    ├── indel_realignment/ # Realigned BAM files
    ├── base_recalibration/ # Recalibrated BAM files
    ├── haplotype_caller/  # Individual gVCF files
    ├── combine_gvcf/     # Combined gVCF file
    └── genotype_gvcf/    # Final genotyped VCF file
```
