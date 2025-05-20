#!/bin/bash

# index for files
group=${1:-100}

# Find parent path
PARENT_PATH=$(dirname $(greadlink -f $0))
result=${PARENT_PATH##*/}

# Go back path
PARENT_PATH=${PARENT_PATH%/*}
PARENT_PATH=${PARENT_PATH%/*}
PARENT_PATH=${PARENT_PATH%/*}

# Set custom path if data is stored outside github repo

# Find data directory
INFOLDER=$PARENT_PATH'/data/filtered_sequencing/mapping/'

# Make directories for stored data
mkdir $PARENT_PATH'/data/extracted_pairs/'
mkdir $PARENT_PATH'/data/extracted_pairs/mapping/'
OUT_FOLDER=$PARENT_PATH'/data/extracted_pairs/mapping/'

# Reset temp folder
rm -rf $OUT_FOLDER'/temp/'
mkdir $OUT_FOLDER'/temp/'

# Command to combine barcodes and promoters
command=paste
for i in $INFOLDER*$group*.gz; do
    command="$command <(gunzip -cd $i)"
done

echo "Splitting file...\n"
# Split file into smaller chunks
eval $command | awk 'NR%4==2 {print $0}' | parallel --pipe --block 1000M 'cat > '$OUT_FOLDER'/temp/small_chunk{#}'

echo "Sorting Chunks...\n"
# Sort each chunk separately
for X in $OUT_FOLDER/temp/small_chunk*; do 
  filename="${X##*/}"
  sort --parallel=8 < $X > $OUT_FOLDER'/temp/sorted-'$filename;
  rm $X
done

echo "Combining chunks and saving fasta file...\n"
# Combine sorted chunks and count unique occurences, store in fasta file
sort -m --parallel=8 $OUT_FOLDER/temp/sorted-small_chunk* | uniq -c | awk '{
    print ">"$3"_"$1"\n"$2;
    print "";
}' | gzip > $OUT_FOLDER'/'$group'_collapsed.fasta.gz'


rm -rf $OUT_FOLDER'/temp/'

#eval $command | awk 'NR%4==2 {print $0}' | sort --parallel='4' | uniq -c | sort --parallel='4' -bgr | awk '{
 #   print $3"\t"$2;
#}' > $group'_collapsed.txt'