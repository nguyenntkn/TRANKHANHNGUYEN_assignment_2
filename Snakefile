# Variables
SRA="SRR1972739"
REF_ID="AF086833.2"
RESULTS_FOLDER="results"
RAW_DIR=f"{RESULTS_FOLDER}/raw"
ALIGNED_DIR=f"{RESULTS_FOLDER}/aligned"
VARIANT_DIR=f"{RESULTS_FOLDER}/variants"
ANNOTATED_DIR=f"{RESULTS_FOLDER}/annotated"
QC_DIR=f"{RESULTS_FOLDER}/qc"
SNPEFF_DIR=f"{RESULTS_FOLDER}/snpEff"
SNPEFF_DATA_DIR=f"{SNPEFF_DIR}/data/reference_db"
SNAKEMAKE_DIR=f"{RESULTS_FOLDER}/snakemake"
BUCKET="trankhanhnguyenassignment2"
S3_PREFIX="ebola"

rule all:
    input:
        f"{SNAKEMAKE_DIR}/.dirs_created",
        f"{RAW_DIR}/reference.fasta",
        f"{SNPEFF_DATA_DIR}/genes.gbk",
        f"{RAW_DIR}/{SRA}/{SRA}.sra",
        f"{RAW_DIR}/{SRA}.fastq",
        f"{QC_DIR}/{SRA}_fastqc.html",
        f"{RAW_DIR}/reference.fasta.fai",
        f"{RAW_DIR}/reference.fasta.bwt",
        f"{RAW_DIR}/reference.dict",
        f"{ALIGNED_DIR}/aligned.sam",
        f"{ALIGNED_DIR}/aligned.sorted.bam",
        f"{ALIGNED_DIR}/validation_report.txt",
        f"{ALIGNED_DIR}/dedup.bam",



rule create_dirs:
    output:
        marker = f"{SNAKEMAKE_DIR}/.dirs_created"
    shell:
        """
        mkdir -p {RAW_DIR} {ALIGNED_DIR} {VARIANT_DIR} {ANNOTATED_DIR} {QC_DIR} {SNPEFF_DATA_DIR}
        touch {output.marker}
        """

rule download_reference:
    input:
        marker = rules.create_dirs.output.marker
    output:
        reference_fasta = f"{RAW_DIR}/reference.fasta",
        reference_gbk = f"{SNPEFF_DATA_DIR}/genes.gbk",
    shell:
        """
        echo Downloading reference genome...
        efetch -db nucleotide -id {REF_ID} -format fasta > {output.reference_fasta}
        efetch -db nucleotide -id {REF_ID} -format genbank > {output.reference_gbk}
        echo Downloaded reference genome!
        """

rule download_sra:
    input:
        marker = rules.create_dirs.output.marker
    output:
        sequence_sra = f"{RAW_DIR}/{SRA}/{SRA}.sra"
    shell:
        """
        echo Downloading sequencing data...
        prefetch {SRA} -O {RAW_DIR}
        echo Downloaded sequencing data!
        """

rule extract_sequence:
    input:
        marker = rules.create_dirs.output.marker,
        sequence_sra = rules.download_sra.output.sequence_sra,
    output:
        sequence_fastq = f"{RAW_DIR}/{SRA}.fastq"
    shell:
        """
        echo Extracting sequencing data...
        fastq-dump -X 10000 {RAW_DIR}/{SRA}/{SRA}.sra -O {RAW_DIR}
        echo Extracted sequencing data!
        """

rule fastqc_raw_reads:
    input:
        marker = rules.create_dirs.output.marker,
        sequence_fastq = rules.extract_sequence.output.sequence_fastq,
    output:
        fastqc_report = f"{QC_DIR}/{SRA}_fastqc.html"
    shell:
        """
        echo Running FastQC on raw reads...
        fastqc -o {QC_DIR} {input.sequence_fastq}
        echo FastQC completed!
        """

rule index_reference:
    input:
        marker = rules.create_dirs.output.marker,
        reference_fasta = rules.download_reference.output.reference_fasta,
    output:
        reference_fasta_index = f"{RAW_DIR}/reference.fasta.fai",
    shell:
        """
        echo Indexing reference genome...
        samtools faidx {input.reference_fasta}
        echo Indexed reference genome!
        """

rule bwa_index:
    input:
        marker = rules.create_dirs.output.marker,
        reference_fasta = rules.download_reference.output.reference_fasta,
    output:
        index_bwa = f"{RAW_DIR}/reference.fasta.bwt",
    shell:
        """
        echo Building BWA index...
        bwa index {input.reference_fasta}   
        """

rule create_fasta_dict_gatk: 
    input:
        marker = rules.create_dirs.output.marker,
        reference_fasta = rules.download_reference.output.reference_fasta,
    output:
        fasta_dict = f"{RAW_DIR}/reference.dict"
    shell:
        """
        echo Creating FASTA dictionary using GATK...
        gatk CreateSequenceDictionary -R {input.reference_fasta} -O {output.fasta_dict}
        echo Created FASTA dictionary!
        """

rule read_alignment:
    input:
        marker = rules.create_dirs.output.marker,
        reference_fasta = rules.download_reference.output.reference_fasta,
        sequence_fastq = rules.extract_sequence.output.sequence_fastq,
        bwa_index = rules.bwa_index.output.index_bwa,
    output:
        aligned_sam = f"{ALIGNED_DIR}/aligned.sam"
    shell:
        """
        echo Aligning reads with read groups...
        bwa mem -R '@RG\\tID:1\\tLB:lib1\\tPL:illumina\\tPU:unit1\\tSM:sample1' {input.reference_fasta} {input.sequence_fastq} > {output.aligned_sam}
        echo Aligned reads!
        """

rule sam_to_sorted_bam:
    input:
        marker = rules.create_dirs.output.marker,
        aligned_sam = rules.read_alignment.output.aligned_sam,
    output: 
        sorted_bam = f"{ALIGNED_DIR}/aligned.sorted.bam"
    shell:
        """
        echo Converting SAM to sorted BAM...
        samtools view -b {input.aligned_sam} | samtools sort -o {output.sorted_bam}
        echo Converted SAM to sorted BAM!
        """ 

rule validate_bam:
    input:
        marker = rules.create_dirs.output.marker,
        sorted_bam = rules.sam_to_sorted_bam.output.sorted_bam,
    output:
        validation_report = f"{ALIGNED_DIR}/validation_report.txt"
    shell:
        """
        echo Validating BAM file...
        gatk ValidateSamFile -I {input.sorted_bam} -MODE SUMMARY > {output.validation_report}
        echo BAM file validation completed!
        """

rule mark_duplicates:   
    input:
        marker = rules.create_dirs.output.marker,
        sorted_bam = rules.sam_to_sorted_bam.output.sorted_bam,
    output:
        dedup_bam = f"{ALIGNED_DIR}/dedup.bam",
        metrics_txt = f"{ALIGNED_DIR}/dup_metrics.txt"
    shell:
        """
        echo Marking duplicates..
        gatk MarkDuplicates -I {input.sorted_bam} -O {output.dedup_bam} -M {output.metrics_txt}
        echo Duplicates marked!
        """