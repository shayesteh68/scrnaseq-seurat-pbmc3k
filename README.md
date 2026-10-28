# PBMC 3k single-cell workflow

This repository contains a Seurat pipeline for the 10x Genomics PBMC 3k example dataset. It performs QC, normalization, PCA, clustering, marker discovery, and manual cell-type labeling.

## Input and environment

Run Rscript scripts/00_download_data.R to download and unpack the filtered hg19 matrix into raw_data/filtered_gene_bc_matrices/hg19/.

Create the Conda environment with conda env create -f environment.yml and activate it with conda activate pbmc3k-seurat.

## Pipeline

From the repository root, run:

1. Rscript scripts/00_download_data.R
2. Rscript scripts/01_qc_filter.R
3. Rscript scripts/02_normalize_pca.R
4. Rscript scripts/03_clustering_umap.R
5. Rscript scripts/04_find_markers.R
6. Rscript scripts/05_cell_annotation.R

The cluster labels in step 05 are manual mappings for cluster IDs 0 through 8. The committed marker table supports the broad labels (including CCR7/TCF7/IL7R for cluster 0 and the lineage markers for monocytes, B cells, cytotoxic T/NK cells, dendritic cells, and platelets). The script requires an exact cluster-ID match and marks labels provisional if clustering parameters or package versions change.

PBMC 3k represents one donor. Cluster-marker tables describe this dataset and do not provide population-level inference. Previously published cell counts and figures are snapshots; regenerate them before citing.

The MIT license is included in LICENSE.
