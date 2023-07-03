#!/bin/bash
# Usage:
# bash 12_run_velocyto.sh <genome_path> <skip-existing>
# genome_path: Path to the genome folder
# The appropriate genome file can be generated following the instructions on the 10x website
# https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build#mm10_2020A
# skip-existing: 1 to skip existing files, 0 to overwrite existing files

# Parse arguments. Ensure that 2 arguments are passed
if [ $# -ne 2 ]; then
    echo "Usage: bash $0 <genome_path> <skip-existing>"
    exit 1
fi

GENOME_PATH=$1
SKIP_EXISTING=$2

DATASETS_FILE=datasets.csv

# Read DATASETS_FILE and extract the following columns
# 1 - Dataset name (e.g. Cheung2018M)
# 4 - Path to sample (e.g. Cheung2018/SRR7898909/outs/filtered_feature_bc_matrix)
# 6 - Sample name (e.g. SRR7898909)

BLUEBOLD='\033[1;34m'
REDBOLD='\033[1;31m'
NC='\033[0m'

# Go through each line in the DATASETS_FILE (skip header)
while read -r line; do
    DATASET=$(echo $line | cut -d ',' -f 1)
    DATA_PATH=$(echo $line | cut -d ',' -f 4)
    SAMPLE=$(echo $line | cut -d ',' -f 6)

    if [[ $DATASET == "study_id" ]]; then
        continue
    fi

    echo -e "${BLUEBOLD}Running velocyto for ${DATASET} ${SAMPLE}${NC}"

    # If the output file exists and SKIP_EXISTING is set to 1, skip this sample otherwise delete the file
    if [ -f RNA_velocity/${DATASET}_${SAMPLE}.loom ]; then
        if [ $SKIP_EXISTING -eq 1 ]; then
            echo -e "${REDBOLD}Output file RNA_velocity/${DATASET}_${SAMPLE}.loom already exists. Skipping...${NC}"
            continue
        fi
        echo -e "${REDBOLD}Output file RNA_velocity/${DATASET}_${SAMPLE}.loom already exists. Deleting...${NC}"
        rm RNA_velocity/${DATASET}_${SAMPLE}.loom
    fi

    echo -e "Running command:\nvelocyto run --bcfile ${DATA_PATH}/cort_barcodes.tsv --sampleid ${DATASET}_${SAMPLE} --outputfolder RNA_velocity/${DATASET}_${SAMPLE} --samtools-threads 16 --samtools-memory 4096 ${DATA_PATH}/../possorted_genome_bam.bam ${GENOME_PATH}/gencode.vM23.primary_assembly.annotation.gtf.filtered"
    # Run velocyto
    velocyto run \
        --bcfile ${DATA_PATH}/cort_barcodes.tsv \
        --sampleid ${DATASET}_${SAMPLE} \
        --outputfolder RNA_velocity \
        --samtools-threads 16 \
        --samtools-memory 4096 ${DATA_PATH}/../possorted_genome_bam.bam \
        ${GENOME_PATH}/gencode.vM23.primary_assembly.annotation.gtf.filtered
done <${DATASETS_FILE}
