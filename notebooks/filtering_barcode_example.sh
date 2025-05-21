#!/usr/bin/env bash
# Find read
READ_DNA="../data/sequencing_reads/barcodes/gc-rep_DNA.fastq.gz"
READ_RNA="../data/sequencing_reads/barcodes/gc-rep_RNA.fastq.gz"

# Define output file paths
OUT_DNA="../data/filtered_sequencing/barcodes/gc-rep_DNA.fastq.gz"
OUT_RNA="../data/filtered_sequencing/barcodes/gc-rep_RNA.fastq.gz"

# Define string to be ran on the termina
mamba run -n fastp fastp --in1 $READ_DNA --out1 $OUT_DNA --trim_tail1 '6'  --verbose --disable_length_filtering --thread '6' --n_base_limit '0'
mamba run -n fastp fastp --in1 $READ_RNA --out1 $OUT_RNA --trim_tail1 '6'  --verbose --disable_length_filtering --thread '6' --n_base_limit '0'
