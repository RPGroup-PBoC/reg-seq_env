# Go back path
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


FOLDER=$SCRIPT_DIR'/../../../data/sequencing_reads/mapping/'

fastq-dump SRR33799406 -O $FOLDER --split-files