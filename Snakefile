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



