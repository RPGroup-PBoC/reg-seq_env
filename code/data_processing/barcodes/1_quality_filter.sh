#!/bin/bash
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


FOLDER=$SCRIPT_DIR'/../../../data/sequencing_reads/barcodes/'
OUT_FOLDER=$SCRIPT_DIR'/../../../data/filtered_sequencing/barcodes/'

if [ -d $OUT_FOLDER ] 
then
    echo $OUT_FOLDER" exists"
else 
    mkdir $OUT_FOLDER
fi


# iterate through all sequencing files
for file in "$FOLDER"/*.fastq.gz; do
  if [ -f "$file" ]; then
    # get length to trim
    TRIM_LENGTH=$(gzip -cd "$file" | awk 'NR==2 {print length-20; exit}')
    # extract file name
    fname=$(basename "$file")   
    # extract sample name
    id=${fname%%_S*}   
    # set path to output
    OUT=$OUT_FOLDER"/"$id".fastq.gz"
    # run filter
    mamba run -n regseq2 fastp --in1 $file --out1 $OUT --trim_tail1 $TRIM_LENGTH  --verbose --disable_length_filtering --thread '6' --n_base_limit '0'
  fi
done
