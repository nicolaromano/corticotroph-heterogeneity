#!/bin/bash

# Check we have exactly three parameters
if [ $# -ne 3 ]; then
    echo "Usage: $0 <SRR_id> <output_dir> <SRA_toolkit_dir>"
    exit 1
fi

SRR_ID = $1
OUTPUT_DIR = $2
SRA_TOOLKIT_DIR = $3

$SRA_TOOLKIT_DIR/bin/fastq-dump -v --split-files --gzip --outdir $OUTPUT_DIR $SRR_ID

# Now rename the files to match the expected format for CellRanger 
# (see https://kb.10xgenomics.com/hc/en-us/articles/115003802691-How-do-I-prepare-Sequence-Read-Archive-SRA-data-from-NCBI-for-Cell-Ranger-)
mv $OUTPUT_DIR/$SRR_ID_1.fastq.gz $OUTPUT_DIR/$SRRID_S1_L001_R1_001.fastq.gz
mv $OUTPUT_DIR/$SRR_ID_2.fastq.gz $OUTPUT_DIR/$SRRID_S1_L001_R2_001.fastq.gz
# If we have a _3 or a _4 file, we need to rename it as well, this would be the index
if [ -f $OUTPUT_DIR/$SRR_ID_3.fastq.gz ]; then
    mv $OUTPUT_DIR/$SRR_ID_3.fastq.gz $OUTPUT_DIR/$SRRID_S1_L001_I1_001.fastq.gz
fi

if [ -f $OUTPUT_DIR/$SRR_ID_4.fastq.gz ]; then
    mv $OUTPUT_DIR/$SRR_ID_4.fastq.gz $OUTPUT_DIR/$SRRID_S1_L001_I2_001.fastq.gz
fi
