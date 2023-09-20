# Corticotrophs heterogeneity

This repository contains scripts for the analysis of corticotrophs in scRNAseq pituitary datasets.

The datasets used for this study are listed in the `datasets.csv` file
The scripts have been numbered sequentially to help navigate the repository.

---

[00_download_data.sh](#00downloaddata)

[01_align_to_genome.sh](#01aligntogenome)

[01b_make_cellranger_rat_reference.sh](#01bmakecellrangerratreference)

[02_aggregate_counts.sh](#02aggregatecounts)

[03_cellranger_to_seurat_and_QC.R](#03cellrangertoseuratandqc)

[04_cell_typing.R](#04celltyping)

[05_subclustering.R](#05subclustering)

[06_compare_datasets.R](#06comparedatasets)

[07_correlation_maps.R](#07correlationmaps)

[08_export_matrices.R](#08exportmatrices)

[09_train_models.ipynb](#09trainmodels)

[10_label_transfer_prediction.py](#10labeltransfer)

[11_label_transfer_graph.R](#11labeltransfergraph)

[12_prepare_velocity_analysis.R](#11preparevelo)

[13_get_velocyto_loom.sh](#13getloom)

[14_scvelo.ipynb](#14scvelo)

---------------

<a name="00downloaddata"></a>

### `00_download_data.sh`

Downloads raw data for each dataset - Note that this takes **a lot** of space.
Usage

    00_download_data.sh <SRA_id> <output_dir> <SRA_toolkit_dir>

For example

    00_download_data.sh SRR123456 /home/pitdatasets/Example2022 /home/sratoolkit/

Note some data is deposited on ArrayExpress, which is not accessible from the NCBI SRA toolkit. A list of files to download is available in the directory and can be simply downloaded using `wget`.

<a name="01aligntogenome"></a>

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

<a name="02aggregatecounts"></a>

### `02_aggregate_counts.sh`

Aggregates the counts from the aligned data using Cell Ranger.

    02_aggregate_counts.sh <id> <output_dir> <cellranger_dir>

For example

    02_aggregate_counts.sh SRR12345 /home/pitdatasets/Example2022 /home/cellranger/

<a name="03cellrangertoseurat"></a>

### `03_cellranger_to_seurat_and_QC.R`

Imports data into R, and plots QC metrics. We use data as filtered by Cellranger to remove empty droplets then do some further filtering. This file will output RDS files containing raw counts and SCT-transformed data

<a name="04celltyping"></a>

### `04_cell_typing.R`

Clusters data, then filters them to only get POMC-expressing clusters (saved to RDS); these are further divided into melanotrophs and corticotrophs (saved to separate RDS files) by looking at the expression of Pcsk2 and Pax7.

<a name="05subclustering"></a>

### `05_subclustering.R`

Finds subclusters in corticotrophs, saves the output to RDS

<a name="06compare"></a>

### `06_compare_datasets.R`
Finds markers for each subcluster in each dataset, and compares the % of common markers between datasets. From these, it creates a graph and finds the communities of most similar clusters between datasets.

<a name="07correlationmaps"></a>

### `07_correlation_maps.R`

Plots correlation maps between the markers from each datasets on to the others. Plots intra-cluster and between-clusters metrics.

<a name="08exportmatrices"></a>
### `08_export_matrices.R`

Exports the expression data, dimension reductions, metadata and the assigned clusters to CSV files, to be used in Python for transfer learning and later for velocity analysis.

<a name="09trainmodels"></a>

### `09_train_models.ipynb`

Trains the models for transfer learning, and saves them to disk.

<a name="10labeltransfer"></a>

### `10_predict_labels.py`

Performs prediction of labels on the new dataset, and calculates saliency maps to identify the most important genes for the prediction.

<a name="11labeltransfergraph"></a
### `11_label_transfer_graph.R`

Creates a graph from the prediction results.

<a name="12preparevelo"></a>

### `12_prepare_velocity_analysis.R`

Prepares the data for velocity analysis, by creating corticotrops-specific barcode files to feed to velocyto.

<a name="13getloom"></a>

### `13_get_velocyto_loom.sh`

Runs velocyto on the data, and creates loom files containing the spliced and unspliced matrices for all datasets.

<a name="14scvelo"></a>

### `14_scvelo.ipynb`

Runs scvelo on the loom files, and plots the results.
