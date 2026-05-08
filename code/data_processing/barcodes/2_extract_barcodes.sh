#!/bin/bash

# index for files
GROUP=${1:-1}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


FOLDER=$SCRIPT_DIR'/../../../data/filtered_sequencing/barcodes'
OUT_DIR=$SCRIPT_DIR'/../../../data/extracted_barcodes'

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