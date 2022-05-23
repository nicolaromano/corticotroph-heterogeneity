# Corticotrophs heterogeneity

This repository contains scripts for the analysis of corticotrophs in scRNAseq pituitary datasets.

The datasets used for this study are listed in the `datasets.csv` file
The scripts have been numbered sequentially to help navigate the repository.

---

[00__download_data.sh](#00downloaddatash)

[01__align_to_genome.sh](#01aligntogenomesh)

[01b__make_cellranger_rat_reference.sh](#01bmakecellrangerratreference)

[02__aggregate_counts.sh](#02aggregatecountssh)


<a name="00downloaddatash"></a>
### `00_download_data.sh`

Downloads raw data for each dataset - Note that this takes **a lot** of space.
Usage

    00_download_data.sh <SRA_id> <output_dir> <SRA_toolkit_dir>

For example

    00_download_data.sh SRR123456 /home/pitdatasets/Example2022 /home/sratoolkit/

Note some data is deposited on ArrayExpress, which is not accessible from the NCBI SRA toolkit. A list of files to download is available in the directory and can be simply downloaded using `wget`.

<a name="01aligntogenomesh"></a>
### `01_align_to_genome.sh`

Aligns the raw data to the genome using Cell Ranger.
Usage

    01_align_to_genome.sh <fastq_dir> <id (sample/output)> <cellranger_dir> <genome> [<expected cells>]

For example

    01_align_to_genome.sh /home/pitdatasets/Example2022/ SRR12345 /home/cellranger/ /home/cellranger/genomes/refdata-gex-GRCh38-2020-A 6000

Data in this study was aligned to the mouse reference dataset 2020-A (July 7, 2020) available on [the 10x website](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest), based on mouse reference mm10 (GENCODE vM23/Ensembl 98).

<a name="01bmakecellrangerratreference"></a>
### `01b_make_cellranger_rat_reference.sh`

Creates a rat reference dataset for Cell Ranger.

<a name="02aggregatecountssh"></a>
### `02_aggregate_counts.sh`

Aggregates the counts from the aligned data using Cell Ranger.

    02_aggregate_counts.sh <id> <output_dir> <cellranger_dir>

For example

    02_aggregate_counts.sh SRR12345 /home/pitdatasets/Example2022 /home/cellranger/