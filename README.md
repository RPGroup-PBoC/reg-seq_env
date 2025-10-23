# The Environment-Dependent Regulatory Landscape of the E. coli Genome

This github repository accompanies the paper ["The Environment-Dependent Regulatory Landscape of the E. coli Genome"](https://arxiv.org/abs/2505.08764). All code needed to process and analyze data, make figures and perform other computational tasks is stored in this repository.

## Software module
All code in this repository is either run in Bash or in Python. For processing of raw sequencing data we use the software module `fastp`. There is a pre-set conda environment in this repo, which can be installed by running
```
mamba create -f requirements.yaml
```

This environment includes all packages necessary to run code in this repo.

## Data
Sequencing data has been deposited in the SRA archive under the project ID PRJNA1263894. There are prewritten scripts that download the data into the correct folders. Check the `Code` section below for more details. The folder structure is pre-defined and is used by the processing codes.

## Code
Contains all code used for processing and analyzing data. Sequencing data can be downloaded using the scripts `code/processing/barcodes/0_import_data.sh` and `code/processing/barcode_mapping/0_import_data.sh`. These scripts require that the `sra-toolkit` is installed, which can be found [here](https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump) (we tested version 3.0.1). All other code is in the form of python scripts or jupyter notebooks.