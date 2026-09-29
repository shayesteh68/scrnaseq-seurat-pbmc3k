suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
})

cat("==> [Step 1] Loading QC-filtered Seurat object...\n")
pbmc <- readRDS("results/qc/pbmc_qc_filtered.rds")

cat("==> [Step 2] Normalizing data (LogNormalize, scale.factor = 10000)...\n")
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

cat("==> [Step 3] Identifying top 2000 Highly Variable Features (HVGs)...\n")
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Plot top variable genes
top10 <- head(VariableFeatures(pbmc), 10)
cat("==> Top 10 Variable Genes:\n", paste(top10, collapse = ", "), "\n")

plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
ggsave("figures/02_variable_features.png", plot = plot2, width = 9, height = 5, dpi = 300)

cat("==> [Step 4] Scaling all genes across cells...\n")
all_genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all_genes)

cat("==> [Step 5] Running Linear Dimensionality Reduction (PCA)...\n")
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# PCA Plots
cat("==> [Step 6] Saving PCA DimPlot & ElbowPlot...\n")
pca_dim <- DimPlot(pbmc, reduction = "pca") + ggtitle("PBMC 3k - PCA Projection")
ggsave("figures/03_pca_dimplot.png", plot = pca_dim, width = 7, height = 5, dpi = 300)

elbow <- ElbowPlot(pbmc, ndims = 30) + ggtitle("Elbow Plot (Variance explained per PC)")
ggsave("figures/04_elbow_plot.png", plot = elbow, width = 7, height = 5, dpi = 300)

# Save intermediate object
cat("==> [Step 7] Saving PCA-completed object...\n")
saveRDS(pbmc, file = "results/qc/pbmc_pca.rds")
cat("==> PCA workflow completed successfully!\n")
