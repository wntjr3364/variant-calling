# Variant Calling Snakemake Workflow

# Configurations
configfile: "config.yaml"

# Get all fastq files
ALL_FASTQ, = glob_wildcards("read/{fname}.fastq.gz")
SAMPLES = []
THREADS = 8

for fname in ALL_FASTQ:
    samplename = fname.split("_")[0]
    SAMPLES.append(samplename)

SAMPLES = list(set(SAMPLES))

# Separate rules for single-end and paired-end processing
rule all:
    input:
        # Reference indexing must complete first
        multiext("reference/" + config["reference"], ".amb", ".ann", ".bwt", ".pac", ".sa", ".dict", ".fai"),

        # Then process all samples
        #expand("result/variants/{sample}.g.vcf.gz", sample=SAMPLES),

        # Final joint calling
        "result/variants/combined_genotyped.vcf.gz"


rule reference_indexing:
    input:
        ref="reference/" + config["reference"]
    output:
        bwa_amb="reference/" + config["reference"] + ".amb",
        bwa_ann="reference/" + config["reference"] + ".ann",
        bwa_bwt="reference/" + config["reference"] + ".bwt",
        bwa_pac="reference/" + config["reference"] + ".pac",
        bwa_sa="reference/" + config["reference"] + ".sa",
        dict="reference/" + config["reference"] + ".dict",
        fai="reference/" + config["reference"] + ".fai"
    
    log:
        bwa="log/indexing_bwa.log",
        picard="log/indexing_picard.log",
        samtools="log/indexing_samtools.log"
    shell:
        """
        mkdir -p log
        
        # BWA indexing
        bwa index -a bwtsw {input.ref} 2> {log.bwa}
        
        # Create sequence dictionary
        gatk CreateSequenceDictionary -R {input.ref} -O {output.dict} 2> {log.picard}
        
        # Create FASTA index
        samtools faidx {input.ref} 2> {log.samtools}
        """

ruleorder: trimmomatic_paired_end > trimmomatic_single_end
# Trimmomatic for Single-End Reads
rule trimmomatic_single_end:
    input:
        read="read/{sample}.fastq.gz",
    output:
        trimmed="result/trim/{sample}_trimmed.fastq.gz"
    log:
        "log/trim/{sample}.log"
    params:
        threads=THREADS,
        leading=3,
        trailing=3,
        slidingwindow="4:5",
        minlen=90
    shell:
        """
        mkdir -p result/trim

        trimmomatic SE -threads {params.threads} \
        {input.read} \
        {output.trimmed} \
        LEADING:{params.leading} \
        TRAILING:{params.trailing} \
        SLIDINGWINDOW:{params.slidingwindow} \
        MINLEN:{params.minlen} \
        &> {log}
        """

# Trimmomatic for Paired-End Reads
rule trimmomatic_paired_end:
    input:
        read1="read/{sample}_1.fastq.gz",
        read2="read/{sample}_2.fastq.gz",
    output:
        paired1="result/trim/{sample}_1_trimmed.fastq.gz",
        unpaired1="result/trim/{sample}_1_untrimmed.fastq.gz",
        paired2="result/trim/{sample}_2_trimmed.fastq.gz",
        unpaired2="result/trim/{sample}_2_untrimmed.fastq.gz"
    log:
        "log/trim/{sample}.log"
    params:
        threads=THREADS,
        leading=3,
        trailing=3,
        slidingwindow="4:5",
        minlen=90
    shell:
        """
        mkdir -p result/trim

        trimmomatic PE -threads {params.threads} \
        {input.read1} {input.read2} \
        {output.paired1} {output.unpaired1} \
        {output.paired2} {output.unpaired2} \
        LEADING:{params.leading} \
        TRAILING:{params.trailing} \
        SLIDINGWINDOW:{params.slidingwindow} \
        MINLEN:{params.minlen} \
        &> {log}
        """

# BWA Alignment for Single-End Reads
rule bwa_alignment_single_end:
    input:
        ref="reference/" + config["reference"],
        reads="result/trim/{sample}_trimmed.fastq.gz",
        # Add dependencies on reference indices
        ref_indices=multiext("reference/" + config["reference"], 
            ".amb", ".ann", ".bwt", ".pac", ".sa", ".dict", ".fai")
    output:
        bam="result/alignment/{sample}.bam"
    log:
        "log/bwa/{sample}.log"
    params:
        threads=THREADS,
    shell:
        """
        mkdir -p result/alignment

        bwa mem -t {params.threads} {input.ref} {input.reads} 2> {log} | \
        samtools sort -@ 8 -o {output.bam}
        """

ruleorder: bwa_alignment_paired_end > bwa_alignment_single_end
# BWA Alignment for Paired-End Reads
rule bwa_alignment_paired_end:
    input:
        ref="reference/" + config["reference"],
        reads=["result/trim/{sample}_1_trimmed.fastq.gz", "result/trim/{sample}_2_trimmed.fastq.gz"],
        # Add dependencies on reference indices
        ref_indices=multiext("reference/" + config["reference"], 
            ".amb", ".ann", ".bwt", ".pac", ".sa", ".dict", ".fai")
    output:
        bam="result/alignment/{sample}.bam"
    log:
        "log/bwa/{sample}.log"
    params:
        threads=THREADS,
    shell:
        """
        mkdir -p result/alignment

        bwa mem -t {params.threads} {input.ref} {input.reads} 2> {log} | \
        samtools sort -@ 8 -o {output.bam}
        """

# Picard Mark Duplicates (works for both SE and PE)
rule mark_duplicates:
    input:
        bam="result/alignment/{sample}.bam",
    output:
        marked_bam="result/alignment/{sample}_marked.bam",
        metrics="result/alignment/{sample}_marked_metrics.txt"
    shell:
        """
        mkdir -p result/alignment

        gatk MarkDuplicates \
        -I {input.bam} \
        -O {output.marked_bam} \
        -M {output.metrics}
        """

# Base Quality Score Recalibration (works for both SE and PE)
rule base_recalibration:
    input:
        bam="result/alignment/{sample}_marked.bam",
        ref="reference/" + config["reference"],
        known_sites="reference/" + config["known_sites"],
    output:
        recal_table="result/alignment/{sample}_recal_data.table",
        recal_bam="result/alignment/{sample}_recal.bam"
    shell:
        """
        mkdir -p result/alignment

        gatk BaseRecalibrator \
        -I {input.bam} \
        -R {input.ref} \
        --known-sites {input.known_sites} \
        -O {output.recal_table}
        
        gatk ApplyBQSR \
        -R {input.ref} \
        -I {input.bam} \
        --bqsr-recal-file {output.recal_table} \
        -O {output.recal_bam}
        """

# HaplotypeCaller in GVCF mode (works for both SE and PE)
rule haplotype_caller_gvcf:
    input:
        bam="result/alignment/{sample}_recal.bam",
        ref="reference/" + config["reference"],
    output:
        gvcf="result/variants/{sample}.g.vcf.gz",
        gvcf_index="result/variants/{sample}.g.vcf.gz.tbi"
    shell:
        """
        mkdir -p result/variants

        gatk HaplotypeCaller \
        -R {input.ref} \
        -I {input.bam} \
        -O {output.gvcf} \
        -ERC GVCF
        """

# CombineGVCFs (combines both SE and PE samples)
rule combine_gvcfs:
    input:
        ref="reference/" + config["reference"],
        gvcfs=expand("result/variants/{sample}.g.vcf.gz", sample=SAMPLES),
    output:
        combined_gvcf="result/variants/combined.g.vcf.gz"
    params:
        variant_inputs=lambda wildcards, input: " ".join(f"-V {v}" for v in input.gvcfs)
    shell:
        """
        mkdir -p result/variants

        gatk CombineGVCFs \
        -R {input.ref} \
        {params.variant_inputs} \
        -O {output.combined_gvcf}
        """

# GenotypeGVCFs
rule genotype_gvcfs:
    input:
        combined_gvcf="result/variants/combined.g.vcf.gz",
        ref="reference/" + config["reference"],
    output:
        genotyped_vcf="result/variants/combined_genotyped.vcf.gz"
    params:
        stand_call_conf=config["variant_calling"]["stand_call_conf"],
        max_alternate_alleles=config["variant_calling"]["max_alternate_alleles"]
    shell:
        """
        mkdir -p result/variants

        gatk GenotypeGVCFs \
        -R {input.ref} \
        -V {input.combined_gvcf} \
        -O {output.genotyped_vcf} \
        --standard-min-confidence-threshold-for-calling {params.stand_call_conf} \
        --max-alternate-alleles {params.max_alternate_alleles}
        """
