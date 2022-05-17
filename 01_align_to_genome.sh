#!/bin/bash

# Time the script
# Note: SECONDS is a bash builtin variable that gives the number of seconds since the script was started
SECONDS=0

# Check we have four or five parameters
if [ $# -ne 4 ] && [ $# -ne 5 ]; then
    echo "Usage: $0 <fastq_dir> <id sample/output> <cellranger_dir> <genome> [<expected cells>]"
    exit 1
fi

FASTQ_DIR=$1
ID=$2
CELLRANGER_DIR=$3
GENOME=$4

if [ $# -eq 5 ]; then
    EXPECTED_CELLS=$5
else
    EXPECTED_CELLS=10000
fi

$CELLRANGER_DIR/cellranger count --id=$ID \
    --sample=$ID \
    --fastqs=$FASTQ_DIR \
    --transcriptome=$GENOME \
    --expect-cells=$EXPECTED_CELLS

# Move the output to the fastq directory
mv $ID $FASTQ_DIR

echo "Alignment completed in $(( $SECONDS / 60 )) minutes."