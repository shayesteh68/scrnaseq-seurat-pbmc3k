# PBMC 3k single-cell workflow

This repository provides a reproducible, exploratory Seurat workflow for the public 10x Genomics PBMC 3k example. It covers input validation, cell-level QC, normalization, PCA, graph clustering, marker discovery, and provisional cell-type labels.

## Dataset and provenance

The workflow uses the official 10x Genomics filtered PBMC 3k hg19 matrix:

<https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz>

`scripts/00_download_data.R` validates the compressed archive, checks its expected MD5 (`F76D73BF7A582AAF32B798BCAE9574D3`), rejects unsafe archive paths, and verifies the required matrix, barcode, and feature files before extraction. The raw archive and extracted matrix are ignored by Git and are downloaded locally when needed.

This is an independent workflow case study using 10x Genomics’ standard PBMC3K example, which contains 2,700 cells from one healthy donor. It is not a reanalysis commissioned by, or conducted in collaboration with, a named laboratory or the data provider. See the [10x Cell Ranger PBMC3K summary](https://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_web_summary.html).

## Environment

The workflow was run on Linux x86_64 under WSL with R 4.4.3, Seurat 5.5.0, and ggplot2 4.0.3. `environment.yml` describes the main dependencies. `environment-linux-64.lock.txt` records the explicit package builds used for the verified Linux run; create that environment with:

```sh
conda create -n pbmc3k-seurat --file environment-linux-64.lock.txt
conda activate pbmc3k-seurat
```

The YAML file can be used to resolve an environment on other platforms, but package builds may differ from the Linux lock.

## Run the workflow

From the repository root, run:

```sh
Rscript scripts/00_download_data.R
Rscript scripts/01_qc_filter.R
Rscript scripts/02_normalize_pca.R
Rscript scripts/03_clustering_umap.R
Rscript scripts/04_find_markers.R
Rscript scripts/05_cell_annotation.R
```

The checked run retained 2,638 of 2,700 input cells after requiring more than 200 and fewer than 2,500 detected features and less than 5% mitochondrial reads. It selected 2,000 variable features, used the first 10 PCs for graph clustering and UMAP, and produced nine clusters at resolution 0.5. Seeds and the UMAP implementation and metric are specified in the scripts.

`FindAllMarkers` reports positive cluster markers using the Wilcoxon test, `min.pct = 0.25`, and `logfc.threshold = 0.25`. The verified run produced 3,455 markers. The annotation step checks that every expected cluster exists and that each label has at least two significant, positive markers in the marker table. Clusters with fewer than 50 cells are marked `provisional_small_cluster`; all labels remain provisional because they are manually assigned from marker signatures.

## Verified output

The checked run produced these cluster sizes and broad labels:

| Cluster | Cells | Provisional label |
| --- | ---: | --- |
| 0 | 705 | Naive CD4 T |
| 1 | 481 | CD14+ monocyte |
| 2 | 472 | Memory-like CD4 T |
| 3 | 344 | B cell |
| 4 | 274 | CD8 T |
| 5 | 162 | FCGR3A+ monocyte |
| 6 | 155 | NK cell |
| 7 | 32 | Dendritic cell, small cluster |
| 8 | 13 | Platelet, small cluster |

Figures are stored in `figures/`. The marker table, marker support evidence, and cell-type proportions are stored in `results/`; intermediate and final Seurat RDS objects are ignored because they are large generated files.

## Interpretation limits

PBMC 3k is a single-sample example, so these results do not support population-level or donor-level inference. Marker tests describe differences between the clusters in this sample; they are not differential-expression tests across biological replicates. QC cutoffs are explicit tutorial choices, and this workflow does not perform doublet detection, ambient-RNA correction, reference-based annotation, or parameter-sensitivity analysis. The 32-cell dendritic and 13-cell platelet clusters are retained but flagged as small. Re-run the workflow before citing counts or figures after changing the input, package versions, or parameters.

The repository code is released under the MIT License. The dataset remains subject to the source provider's terms.
