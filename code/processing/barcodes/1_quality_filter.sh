#!/bin/bash
# Find working directiory
PARENT_PATH=$(dirname $(greadlink -f $0))

# Go back path
PARENT_PATH=${PARENT_PATH%/*}
PARENT_PATH=${PARENT_PATH%/*}
PARENT_PATH=${PARENT_PATH%/*}

# Data FOLDER
FOLDER=$PARENT_PATH'/data/sequencing_reads/barcodes/'
OUT_FOLDER=$PARENT_PATH'/data/filtered_sequencing/barcodes/'

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
    mamba run -n fastp fastp --in1 $file --out1 $OUT --trim_tail1 $TRIM_LENGTH  --verbose --disable_length_filtering --thread '6' --n_base_limit '0'
  fi
done
