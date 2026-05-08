# Go back path
which fasterq-dump
# Get the absolute path to the directory where this script resides
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


FOLDER="$SCRIPT_DIR/../../../data/sequencing_reads/barcodes/"

echo $FOLDER
# file with accession numbers
ACCESSION_FILE=$SCRIPT_DIR"/accession_numbers.txt"

COMPRESS="gzip -f"


# make sure quality scores get imported
vdb-config --simplified-quality-scores no

# Read only columns 1 and 4 (tab-delimited); skip header/blank/comment lines
awk -F'\t' 'NF>=4 && $1!="" && $1!~ /^#/ {print $1 "\t" $4}' "$ACCESSION_FILE" |
while IFS=$'\t' read -r ACC NAME; do
  # ensure .gz extension
  [[ "$NAME" =~ \.gz$ ]] || NAME="${NAME}.fastq.gz"

  echo "Downloading $ACC → $NAME"

  # dump to temp .fastq, then compress into final name
  TMP_FASTQ="${FOLDER}/${ACC}.fastq"
  fasterq-dump "$ACC" -O "$FOLDER" --split-3 
  $COMPRESS -c "$TMP_FASTQ" > "${FOLDER}/${NAME}"
  rm -f "$TMP_FASTQ"

  echo "✔ ${FOLDER}/${NAME}"
done