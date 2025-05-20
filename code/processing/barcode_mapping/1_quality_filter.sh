#!/bin/bash
GROUP=${1:-100}
# Find working directiory
PARENT_PATH=$(dirname $(greadlink -f $0))

# Go back path
PARENT_PATH=${PARENT_PATH%/*}
PARENT_PATH=${PARENT_PATH%/*}
PARENT_PATH=${PARENT_PATH%/*}

# If data is saved to other place, change the PARENT_PATH variable
#PARENT_PATH=<path_to_data>

# Data FOLDER
FOLDER=$PARENT_PATH'/data/sequencing_reads/mapping/'
OUT_FOLDER=$PARENT_PATH'/data/filtered_sequencing/mapping/'

if [ -d $OUT_FOLDER ] 
then
    echo $OUT_FOLDER" exists"
else 
    mkdir $OUT_FOLDER
fi


# Find read1 and read2 files
READ1=$(find $FOLDER -name $GROUP"*R1*.fastq.gz")
READ2=$(find $FOLDER -name $GROUP"*R2*.fastq.gz")

# Define output file paths
OUT1=$OUT_FOLDER"/"$GROUP"_R1_filtered.fastq.gz"
OUT2=$OUT_FOLDER"/"$GROUP"_R2_filtered.fastq.gz"

HTML=$OUT_FOLDER"/"$GROUP"_fastp_report.html"
JSON=$OUT_FOLDER"/"$GROUP"_fastp_report.json"

# Define string to be ran on the terminal
mamba run -n fastp fastp --in1 $READ1 --in2 $READ2 --out1 $OUT1 --out2 $OUT2 --trim_tail1 '11' --trim_tail2 '11' --verbose --disable_length_filtering --html $HTML --json $JSON --thread '6' -q '20' --n_base_limit '0' --unqualified_percent_limit '10'