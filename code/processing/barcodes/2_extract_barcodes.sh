#!/bin/bash

# index for files
GROUP=${1:-1}
PARENT_PATH=$(dirname $(greadlink -f $0))
result=${PARENT_PATH##*/}

# Go back path
PARENT_PATH=${PARENT_PATH%/*}
PARENT_PATH=${PARENT_PATH%/*}
PARENT_PATH=${PARENT_PATH%/*}

# Find data directory
FOLDER=$PARENT_PATH'/data/filtered_sequencing/barcodes'
OUT_DIR=$PARENT_PATH'/data/extracted_barcodes'

# Make directories for stored data
mkdir $OUT_DIR


# iterate through all sequencing files
for file in "$FOLDER"/*.fastq.gz; do
  if [ -f "$file" ]; then

    # extract file name
    fname=$(basename "$file")
    # extract sample name
    id=${fname%.*}
    id=${id%.*} 
    # update terminal
    echo "Extracting from "$id
    # set path to output
    OUT=$OUT_DIR"/"$id"_collapsed.txt"
    # run filter
    # Find barcodes
    gunzip -c $file | awk ' NR%4==2 {
            print $0;
        }'| sort | uniq -c | sort -bgr > $OUT
  fi
done