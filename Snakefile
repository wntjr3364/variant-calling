# Variant Calling Snakemake Workflow

# Configurations
configfile: "config.yaml"

REFDICT = config["reference_dir"]
OUTDIR = config["output_dir"]
LOGDIR = config["log_dir"]
TEMPDIR = config["temp_dir"]

THREADS = config["threads"]
MEMORY = config["memory"]

# Sample setup
ALL_FASTQ, = glob_wildcards("read/{fname}.fastq.gz")
SAMPLES = []

for fname in ALL_FASTQ:
    samplename = fname.split("_")[0]
    SAMPLES.append(samplename)

SAMPLES = list(set(SAMPLES))

# Reference setup
REFBASE = os.path.basename(config["reference"])
REFNAME = ".".join(REFBASE.split(".")[:-1])
REFINDEXES = multiext(REFDICT + os.sep + REFNAME, ".fa.amb", ".fa.ann", ".fa.bwt", ".fa.pac", ".fa.sa", ".fa.fai", ".dict")

rule all:
    input:
        # Reference indexing must complete first
        REFINDEXES,

        # Final joint calling
        OUTDIR + os.sep + "genotype_gvcf/combined_genotyped.vcf.gz"

rule bwa_indexing:
    input:
        ref=REFDICT + os.sep + REFBASE
    output:
        bwa_amb=REFDICT + os.sep + REFBASE + ".amb",
        bwa_ann=REFDICT + os.sep + REFBASE + ".ann",
        bwa_bwt=REFDICT + os.sep + REFBASE + ".bwt",
        bwa_pac=REFDICT + os.sep + REFBASE + ".pac",
        bwa_sa=REFDICT + os.sep + REFBASE + ".sa",
    log:
        bwa=LOGDIR + os.sep + "indexing_bwa.log"
    shell:
        """
        mkdir -p log

        bwa index -a bwtsw {input.ref} 2> {log.bwa}
        """
        
rule picard_indexing:
    input:
        ref=REFDICT + os.sep + REFBASE
    output:
        dict=REFDICT + os.sep + REFNAME + ".dict"
    log:
        picard=LOGDIR + os.sep + "indexing_picard.log"
    shell:
        """
        mkdir -p log

        java -jar program/picard-2.8.3/picard.jar CreateSequenceDictionary \
        R={input.ref} \
        O={output.dict} \
        2> {log.picard}
        """
rule samtools_indexing:
    input:
        ref=REFDICT + os.sep + REFBASE
    output:
        fai=REFDICT + os.sep + REFBASE + ".fai"
    log:
        samtools=LOGDIR + os.sep + "indexing_samtools.log"
    shell:
        """
        mkdir -p log

        samtools faidx {input.ref} 2> {log.samtools}
        """

ruleorder: trimmomatic_paired_end > trimmomatic_single_end
# Trimmomatic for Single-End Reads
rule trimmomatic_single_end:
    input:
        read="read/{sample}.fastq.gz"
    output:
        trimmed=OUTDIR + os.sep + "trim/{sample}_trimmed.fastq.gz"
    log:
        LOGDIR + os.sep + "trim/{sample}.log"
    params:
        threads=THREADS,
        leading=config["trimmomatic"]["leading"],
        trailing=config["trimmomatic"]["trailing"],
        slidingwindow=config["trimmomatic"]["slidingwindow"],
        minlen=config["trimmomatic"]["minlen"]
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
        paired1=OUTDIR + os.sep + "trim/{sample}_1_trimmed.fastq.gz",
        unpaired1=OUTDIR + os.sep + "trim/{sample}_1_untrimmed.fastq.gz",
        paired2=OUTDIR + os.sep + "trim/{sample}_2_trimmed.fastq.gz",
        unpaired2=OUTDIR + os.sep + "trim/{sample}_2_untrimmed.fastq.gz"
    log:
        LOGDIR + os.sep + "trim/{sample}.log"
    params:
        threads=THREADS,
        leading=config["trimmomatic"]["leading"],
        trailing=config["trimmomatic"]["trailing"],
        slidingwindow=config["trimmomatic"]["slidingwindow"],
        minlen=config["trimmomatic"]["minlen"]
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

ruleorder: bwa_alignment_paired_end > bwa_alignment_single_end
# BWA Alignment for Single-End Reads
rule bwa_alignment_single_end:
    input:
        ref=REFDICT + os.sep + REFBASE,
        reads=OUTDIR + os.sep + "trim/{sample}_trimmed.fastq.gz",
        # Add dependencies on reference indices
        ref_indices=REFINDEXES
    output:
        bam=OUTDIR + os.sep + "alignment/{sample}.bam"
    log:
        bwa=LOGDIR + os.sep + "bwa/{sample}_bwa.log",
        samtools=LOGDIR + os.sep + "bwa/{sample}_samtools.log"
    params:
        threads=THREADS,
    shell:
        """
        mkdir -p result/alignment

        bwa mem -t {params.threads} {input.ref} {input.reads} 2> {log.bwa} | \
        samtools sort -@ {params.threads} -o {output.bam} &> {log.samtools}
        """

# BWA Alignment for Paired-End Reads
rule bwa_alignment_paired_end:
    input:
        ref=REFDICT + os.sep + config["reference"],
        reads=[OUTDIR + os.sep + "trim/{sample}_1_trimmed.fastq.gz", OUTDIR + os.sep + "trim/{sample}_2_trimmed.fastq.gz"],
        # Add dependencies on reference indices
        ref_indices=REFINDEXES
    output:
        bam=OUTDIR + os.sep + "alignment/{sample}.bam"
    log:
        bwa=LOGDIR + os.sep + "bwa/{sample}_bwa.log",
        samtools=LOGDIR + os.sep + "bwa/{sample}_samtools.log"
    params:
        threads=THREADS,
    shell:
        """
        mkdir -p result/alignment

        bwa mem -t {params.threads} {input.ref} {input.reads} 2> {log.bwa} | \
        samtools sort -@ {params.threads} -o {output.bam} &> {log.samtools}
        """

rule read_group:
    input:
        bam=OUTDIR + os.sep + "alignment/{sample}.bam"
    output:
        rg_bam=OUTDIR + os.sep + "read_group/{sample}.RGsorted.bam"
    params:
        java_mem="10g",
        tmp_dir=TEMPDIR,
        max_records="1280000"
    log:
        LOGDIR + os.sep + "read_group/{sample}.RGsorted.log"
    shell:
        """
        mkdir -p result/read_group log/read_group {params.tmp_dir}

        picard AddOrReplaceReadGroups \
        -Xmx{params.java_mem} \
        -Djava.io.tmpdir={params.tmp_dir} \
        INPUT={input.bam} \
        OUTPUT={output.rg_bam} \
        SORT_ORDER=coordinate \
        MAX_RECORDS_IN_RAM={params.max_records} \
        VALIDATION_STRINGENCY=LENIENT \
        RGID={wildcards.sample} \
        RGLB={wildcards.sample}_LIB \
        RGPL=ILLUMINA \
        RGPU=NONE \
        RGSM={wildcards.sample} \
        &> {log}
        """

rule mark_duplicates:
    input:
        bam=OUTDIR + os.sep + "read_group/{sample}.RGsorted.bam"
    output:
        marked_bam=OUTDIR + os.sep + "mark_duplication/{sample}_marked.bam",
        metrics=OUTDIR + os.sep + "mark_duplication/{sample}_marked_metrics.txt",
        bam_index=OUTDIR + os.sep + "mark_duplication/{sample}_marked.bam.bai"
    params:
        java_mem="10g",
        tmp_dir=TEMPDIR,
        max_records="1280000",
        max_handles="1024"
    log:
        dedup=LOGDIR + os.sep + "mark_duplication/{sample}.dedup.log",
        index=LOGDIR + os.sep + "mark_duplication/{sample}.index.log"
    shell:
        """
        mkdir -p result/mark_duplication log/mark_duplication {params.tmp_dir}

        picard -Xmx{params.java_mem} -Djava.io.tmpdir={params.tmp_dir} \
        MarkDuplicates \
        INPUT={input.bam} \
        OUTPUT={output.marked_bam} \
        METRICS_FILE={output.metrics} \
        REMOVE_DUPLICATES=true \
        ASSUME_SORTED=true \
        MAX_RECORDS_IN_RAM={params.max_records} \
        VALIDATION_STRINGENCY=LENIENT \
        MAX_FILE_HANDLES={params.max_handles} \
        &> {log.dedup} 

        samtools index {output.marked_bam} &> {log.index}
        """

# Indel Realignment (works for both SE and PE)
rule indel_realignment:
    input:
        bam=OUTDIR + os.sep + "mark_duplication/{sample}_marked.bam",
        ref=REFDICT + os.sep + REFBASE,
        known_sites=REFDICT + os.sep + config["known_sites"]
    output:
        realigned_bam=OUTDIR + os.sep + "indel_realignment/{sample}_realigned.bam",
        intervals=OUTDIR + os.sep + "indel_realignment/{sample}.intervals"
    params:
        java_mem=MEMORY,
        tmp_dir=TEMPDIR,
        threads=THREADS
    log:
        target=LOGDIR + os.sep + "indel_realignment/{sample}.target.log",
        realign=LOGDIR + os.sep + "indel_realignment/{sample}.realigned.log"
    shell:
        """
        mkdir -p result/indel_realignment log/indel_realignment {params.tmp_dir}

        java -Xmx{params.java_mem} -Djava.io.tmpdir={params.tmp_dir} \
        -jar program/GenomeAnalysisTK.jar \
        -T RealignerTargetCreator \
        -R {input.ref} \
        -I {input.bam} \
        -o {output.intervals} \
        -nt {params.threads} \
        &> {log.target}

        java -Xmx{params.java_mem} -Djava.io.tmpdir={params.tmp_dir} \
        -jar program/GenomeAnalysisTK.jar \
        -T IndelRealigner \
        -R {input.ref} \
        -I {input.bam} \
        -targetIntervals {output.intervals} \
        -o {output.realigned_bam} \
        &> {log.realign}
        """

# Base Quality Score Recalibration (works for both SE and PE)
rule base_recalibration:
    input:
        bam=OUTDIR + os.sep + "indel_realignment/{sample}_realigned.bam",
        ref=REFDICT + os.sep + REFBASE,
        known_sites=REFDICT + os.sep + config["known_sites"]
    output:
        recal_table=OUTDIR + os.sep + "base_recalibration/{sample}_recal_data.table",
        recal_bam=OUTDIR + os.sep + "base_recalibration/{sample}_recal.bam"
    params:
        java_mem=MEMORY,
        tmp_dir=TEMPDIR,
        threads=THREADS
    log:
        LOGDIR + os.sep + "base_recalibration/{sample}.recal.log"
    shell:
        """
        mkdir -p result/base_recalibration log/base_recalibration {params.tmp_dir}

        java -Xmx{params.java_mem} -Djava.io.tmpdir={params.tmp_dir} \
        -jar program/GenomeAnalysisTK.jar \
        -T BaseRecalibrator \
        -R {input.ref} \
        -I {input.bam} \
        -knownSites {input.known_sites} \
        -o {output.recal_table} \
        -nct {params.threads} \
        &> {log}
        
        java -Xmx{params.java_mem} -Djava.io.tmpdir={params.tmp_dir} \
        -jar program/GenomeAnalysisTK.jar \
        -T PrintReads \
        -R {input.ref} \
        -I {input.bam} \
        -BQSR {output.recal_table} \
        -o {output.recal_bam} \
        -nct {params.threads} \
        &>> {log}
        """

rule haplotype_caller_gvcf:
    input:
        bam=OUTDIR + os.sep + "base_recalibration/{sample}_recal.bam",
        ref=REFDICT + os.sep + REFBASE
    output:
        gvcf=OUTDIR + os.sep + "haplotype_caller/{sample}.g.vcf.gz",
        gvcf_index=OUTDIR + os.sep + "haplotype_caller/{sample}.g.vcf.gz.tbi"
    params:
        java_mem=MEMORY,
        tmp_dir=TEMPDIR,
        threads=THREADS
    log:
        LOGDIR + os.sep + "haplotype_caller/{sample}.log"
    shell:
        """
        mkdir -p result/haplotype_caller log/haplotype_caller {params.tmp_dir}

        java -Xmx{params.java_mem} -Djava.io.tmpdir={params.tmp_dir} \
        -jar program/GenomeAnalysisTK.jar \
        -T HaplotypeCaller \
        -R {input.ref} \
        -I {input.bam} \
        -o {output.gvcf} \
        -ERC GVCF \
        -nct {params.threads} \
        &> {log}
        """

# CombineGVCFs
rule combine_gvcf:
    input:
        ref=REFDICT + os.sep + REFBASE,
        gvcfs=expand(OUTDIR + os.sep + "haplotype_caller/{sample}.g.vcf.gz", sample=SAMPLES),
    output:
        combined_gvcf=OUTDIR + os.sep + "combine_gvcf/combined.g.vcf.gz"
    params:
        variant_inputs=lambda wildcards, input: " ".join(f"-V {v}" for v in input.gvcfs),
        java_mem=MEMORY,
        tmp_dir=TEMPDIR
    log:
        LOGDIR + os.sep + "combine_gvcf.log"
    shell:
        """
        mkdir -p result/combine_gvcf log/combine_gvcf

        java -Xmx{params.java_mem} -Djava.io.tmpdir={params.tmp_dir} \
        -jar program/GenomeAnalysisTK.jar \
        -T CombineGVCFs \
        -R {input.ref} \
        {params.variant_inputs} \
        -o {output.combined_gvcf} \
        &> {log}
        """

# GenotypeGVCFs
rule genotype_gvcf:
    input:
        combined_gvcf=OUTDIR + os.sep + "combine_gvcf/combined.g.vcf.gz",
        ref=REFDICT + os.sep + REFBASE,
        dbsnp=REFDICT + os.sep + config["known_sites"]
    output:
        genotyped_vcf=OUTDIR + os.sep + "genotype_gvcf/combined_genotyped.vcf.gz"
    params:
        java_mem="30g",
        tmp_dir=TEMPDIR,
        threads=THREADS,
        stand_call_conf=config["variant_calling"]["stand_call_conf"],
        max_alternate_alleles=config["variant_calling"]["max_alternate_alleles"]
    log:
        LOGDIR + os.sep + "genotype_gvcf/genotype.log"
    shell:
        """
        mkdir -p result/genotype_gvcf log/genotype_gvcf {params.tmp_dir}

        java -Xmx{params.java_mem} -Djava.io.tmpdir={params.tmp_dir} \
        -jar program/GenomeAnalysisTK.jar \
        -T GenotypeGVCFs \
        -R {input.ref} \
        -V {input.combined_gvcf} \
        -o {output.genotyped_vcf} \
        --dbsnp {input.dbsnp} \
        --stand_call_conf {params.stand_call_conf} \
        --max_alternate_alleles {params.max_alternate_alleles} \
        -nt {params.threads} \
        &> {log}
        """
        
