#!/bin/bash
GROUP=${1:-1}
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


# Find read
READ_DNA=$(find $FOLDER -name $GROUP"*DNA.fastq.gz")
READ_RNA=$(find $FOLDER -name $GROUP"*RNA.fastq.gz")

# Define output file paths
OUT_DNA=$OUT_FOLDER"/"$GROUP"_DNA.fastq.gz"
OUT_RNA=$OUT_FOLDER"/"$GROUP"_RNA.fastq.gz"

# Define string to be ran on the termina
mamba run -n fastp fastp --in1 $READ_DNA --out1 $OUT_DNA --trim_tail1 '6'  --verbose --disable_length_filtering --thread '6' --n_base_limit '0'
mamba run -n fastp fastp --in1 $READ_RNA --out1 $OUT_RNA --trim_tail1 '6'  --verbose --disable_length_filtering --thread '6' --n_base_limit '0'