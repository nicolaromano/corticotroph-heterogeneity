#!/bin/bash

# Time the script
# Note: SECONDS is a bash builtin variable that gives the number of seconds since the script was started
SECONDS=0

# Check we have exactly three parameters
if [ $# -ne 3 ]; then
    echo "Usage: $0 <id> <output_dir> <cellranger_dir>"
    exit 1
fi

ID=$1
OUTPUT_DIR=$2
CELLRANGER_DIR=$3

# We use normalize=none as suggested here https://github.com/satijalab/seurat/issues/672
$CELLRANGER_DIR/cellranger aggr --id=$ID \
    --csv=$OUTPUT_DIR/aggr.csv \
    --normalize=none

# Move the output files to the output directory
mv $ID $OUTPUT_DIR

echo "Aggregation completed in $(( $SECONDS / 60 )) minutes."