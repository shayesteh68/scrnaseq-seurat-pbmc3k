suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
})

cat("==> [Step 1] Loading QC-filtered Seurat object...\n")
pbmc <- readRDS("results/qc/pbmc_qc_filtered.rds")
set.seed(42)

cat("==> [Step 2] Normalizing data (LogNormalize, scale.factor = 10000)...\n")
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

cat("==> [Step 3] Identifying top 2000 highly variable features...\n")
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(pbmc), 10)
cat("==> Top 10 variable genes:\n", paste(top10, collapse = ", "), "\n")

plot1 <- VariableFeaturePlot(pbmc)
zero_mean_features <- sum(!is.finite(plot1$data$mean) | plot1$data$mean <= 0, na.rm = TRUE)
plot1$data <- plot1$data[is.finite(plot1$data$mean) & plot1$data$mean > 0, , drop = FALSE]
cat(sprintf("Omitting %d zero-mean features from the log-scale variability plot.\n", zero_mean_features))
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
ggsave("figures/02_variable_features.png", plot = plot2, width = 9, height = 5, dpi = 300, bg = "white")

cat("==> [Step 4] Scaling variable features across cells...\n")
pbmc <- ScaleData(pbmc, features = VariableFeatures(pbmc))

cat("==> [Step 5] Running linear dimensionality reduction (PCA)...\n")
pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc), seed.use = 42)

cat("==> [Step 6] Saving PCA DimPlot and ElbowPlot...\n")
pca_dim <- DimPlot(pbmc, reduction = "pca") + ggtitle("PBMC 3k - PCA projection")
ggsave("figures/03_pca_dimplot.png", plot = pca_dim, width = 7, height = 5, dpi = 300, bg = "white")
elbow <- ElbowPlot(pbmc, ndims = 30) + ggtitle("Elbow plot (variance explained per PC)")
ggsave("figures/04_elbow_plot.png", plot = elbow, width = 7, height = 5, dpi = 300, bg = "white")

cat("==> [Step 7] Saving PCA-completed object...\n")
saveRDS(pbmc, file = "results/qc/pbmc_pca.rds")
cat("==> PCA workflow completed successfully!\n")
