#!/bin/bash

# Group number
GROUP=${1:-100}

# Find working directiory

PARENT_PATH=$(dirname $(greadlink -f $0))

# Go back path
PARENT_PATH=${PARENT_PATH%/*}
PARENT_PATH=${PARENT_PATH%/*}
PARENT_PATH=${PARENT_PATH%/*}

# find bbmap path
BBMAP_PATH=$(find $PARENT_PATH -name "bbmap.sh")

# set path to data
DATA_PATH=$PARENT_PATH"/data/extracted_pairs"



# Find data directory
OUT_FOLDER=$PARENT_PATH'/data/barcode_map/'
SAM_FILE=$OUT_FOLDER$GROUP'_collapsed.sam'


# Run BB
$BBMAP_PATH ref=$PARENT_PATH'/data/metadata/wt_sequences.fasta' path=$OUT_FOLDER
gunzip -c $DATA_PATH/$GROUP'_collapsed.fasta.gz' > $DATA_PATH/$GROUP'_collapsed.fasta'
$BBMAP_PATH ambiguous='best' indelfilter='0' nfilter='0' minid='0.85' trimreaddescriptions='t' in=$DATA_PATH/$GROUP'_collapsed.fasta' out=$SAM_FILE t='8' path=$OUT_FOLDER

rm -f $DATA_PATH/$group'_collapsed.fasta'

# # Start the process
START=$SECONDS

# # Input file name
echo "Assigning gene names to promoter-barcode pairs...\n"

# # Extract promoter-bc pairs and corresponding gene names
awk -v o="$OUT_FOLDER" 'BEGIN{FS="\t";OFS = ","} !(NR%500000){print NR " Promoters Processed"}; NF>10{gsub(/_/, ",", $1); print $10,$1,$3 >> (o"/"$3"_barcodes.txt")}' "$SAM_FILE"

# # terminal output message
echo "done! Output files written to " "$OUT_FOLDER\n"
# END=$SECONDS
DURATION=$(( END - START ))
echo "time elapsed: $DURATION seconds"
