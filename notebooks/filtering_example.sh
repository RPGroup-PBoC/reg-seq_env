#!/usr/bin/env bash
READ1="../data/sequencing_reads/mapping/read1_mapping.fastq.gz"
READ2="../data/sequencing_reads/mapping/read2_mapping.fastq.gz"
OUT1="../data/filtered_sequencing/mapping/read1_mapping_filtered.fastq.gz"
OUT2="../data/filtered_sequencing/mapping/read2_mapping_filtered.fastq.gz"

mamba run -n fastp fastp        --in1 "$READ1" --in2 "$READ2"        --out1 "$OUT1" --out2 "$OUT2"        --trim_tail1 11 --trim_tail2 11        --verbose --disable_length_filtering        --thread 6 -q 20 --n_base_limit 0 --unqualified_percent_limit 10
