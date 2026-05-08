#!/bin/bash
GROUP=${1:-100}
# Find working directiory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


FOLDER=$SCRIPT_DIR'/../../../data/sequencing_reads/mapping/230128_M05340_0408_000000000-GFKGC/'
OUT_FOLDER=$SCRIPT_DIR'/../../../data/filtered_sequencing/mapping_pl/'

if [ -d $OUT_FOLDER ] 
then
    echo $OUT_FOLDER" exists"
else 
    mkdir $OUT_FOLDER
fi


# Find read1 and read2 files
READ1=$(find $FOLDER -name "*Undetermined*R1*.fastq.gz")
READ2=$(find $FOLDER -name "*Undetermined*R2*.fastq.gz")

# Define output file paths
OUT1=$OUT_FOLDER"/"$GROUP"_R1_filtered.fastq.gz"
OUT2=$OUT_FOLDER"/"$GROUP"_R2_filtered.fastq.gz"

HTML=$OUT_FOLDER"/"$GROUP"_fastp_report.html"
JSON=$OUT_FOLDER"/"$GROUP"_fastp_report.json"

# Define string to be ran on the terminal
mamba run -n regseq2 fastp --in1 $READ1 --in2 $READ2 --out1 $OUT1 --out2 $OUT2 --trim_tail1 '6' --trim_tail2 '6' --verbose --disable_length_filtering --html $HTML --json $JSON --thread '6' -q '20' --n_base_limit '0' --unqualified_percent_limit '10'