#!/bin/bash


# Check we have exactly two parameters
if [ $# -ne 2 ]; then
    echo "Usage: $0 <filelist> <output_dir>"
    exit 1
fi

# Time the script
# Note: SECONDS is a bash builtin variable that gives the number of seconds since the script was started
SECONDS=0

FILELIST=$1
OUTPUT_DIR=$2

# Write the current time
echo "$(date) - Downloading data into $OUTPUT_DIR"

wget --directory-prefix $OUTPUT_DIR --input-file $FILELIST

echo "All files downloaded in $(( $SECONDS / 60 )) minutes."