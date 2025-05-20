#!/bin/bash

# Group number
group=${1:-100}

# Find working directiory

PARENT_PATH=$(dirname $(greadlink -f $0))

# Go back path
PARENT_PATH=${PARENT_PATH%/*}
PARENT_PATH=${PARENT_PATH%/*}
PARENT_PATH=${PARENT_PATH%/*}

# find bbmap path
BBMAP_PATH=$(find $PARENT_PATH -name "bbmap.sh")

# set path to data
PARENT_PATH=$PARENT_PATH"/data/extracted_pairs/mapping"



# Find data directory
processing_folder=$PARENT_PATH'/data/barcode_map/'
sam_file=$processing_folder'/'$group'_collapsed.sam'



# Test if output directory exists already
out_dir=$PARENT_PATH'/data/barcode_map/'$group'_per_gene'
if [ -d $out_dir ] 
then
    echo $out_dir" exists"
    read -p "Output directory already exists. Do you want to remove it before proceeding?.(Yy/Nn)" yn
    case $yn in
        [Yy]* ) rm -rf $out_dir && echo "Output directory reset";;
        [Nn]* ) echo "Please move Output directory manually and run script again.";;
        * ) echo "Please answer yes[Yy] or no[Nn].";;
    esac
fi
mkdir $out_dir


# Run BB
$BBMAP_PATH ref=$PARENT_PATH'/data/wt_sequences.fasta' path=$processing_folder
gunzip -c $processing_folder/$group'_collapsed.fasta.gz' > $processing_folder/$group'_collapsed.fasta'
$BBMAP_PATH ambiguous='best' indelfilter='0' nfilter='0' minid='0.85' trimreaddescriptions='t' in=$processing_folder/$group'_collapsed.fasta' out=$sam_file t='8' path=$processing_folder

# rm -f $processing_folder/$group'_collapsed.fasta'

# # Start the process
# START=$SECONDS

# # Input file name
# echo "Assigning gene names to promoter-barcode pairs...\n"

# # Extract promoter-bc pairs and corresponding gene names
# awk -v o="$out_dir" 'BEGIN{FS="\t";OFS = ","} !(NR%500000){print NR " Promoters Processed"}; NF>10{gsub(/_/, ",", $1); print $10,$1,$3 >> (o"/"$3"_barcodes.txt")}' "$sam_file"

# # terminal output message
# echo "done! Output files written to " "$out_dir\n"
# END=$SECONDS
# DURATION=$(( END - START ))
# echo "time elapsed: $DURATION seconds"
