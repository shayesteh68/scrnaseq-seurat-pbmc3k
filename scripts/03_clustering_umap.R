suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
})

cat("==> [Step 1] Loading PCA-completed Seurat object...\n")
pbmc <- readRDS("results/qc/pbmc_pca.rds")
set.seed(42)

cat("==> [Step 2] Constructing SNN Graph (dims = 1:10)...\n")
pbmc <- FindNeighbors(pbmc, dims = 1:10)

cat("==> [Step 3] Identifying cell clusters (Louvain algorithm, resolution = 0.5)...\n")
pbmc <- FindClusters(pbmc, resolution = 0.5, algorithm = 1, random.seed = 42)

cat("==> Number of cells per cluster:\n")
print(table(Idents(pbmc)))

cat("==> [Step 4] Running UMAP (Non-linear Dimensionality Reduction)...\n")
pbmc <- RunUMAP(pbmc, dims = 1:10, seed.use = 42, umap.method = "uwot", metric = "cosine")

cat("==> [Step 5] Saving UMAP Clustering Plot...\n")
umap_plot <- DimPlot(pbmc, reduction = "umap", label = TRUE, repel = TRUE, label.size = 3.5, pt.size = 0.5) +
  ggtitle("PBMC 3k - Unsupervised Cell Clusters (UMAP)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))

ggsave("figures/05_umap_clusters.png", plot = umap_plot, width = 8, height = 6, dpi = 300, bg = "white")

cat("==> [Step 6] Saving clustered Seurat object...\n")
saveRDS(pbmc, file = "results/clustering/pbmc_clustered.rds")
cat("==> Clustering and UMAP completed successfully!\n")
