# Processing of mapping data

Scripts in this folder are used to process sequencing data generated to map barcodes and promoter variants. The scripts look for paired-end sequencing data in the `data/sequencing_reads/mapping` folder. If data has not been added yet, refer to the `README` folder in said folder for instructions. 


## `1_quality_filter.sh`
Quality filtering of sequencing reads using `fastp`. Expects a Mamba environment to be present which is called `fastp`. A file to create such an environment can be found on the frontpage of this repository.
If everything is installed simply run from this folder.

```
chmod +x ./1_quality_filter.sh
./1_quality_filter.sh

```