suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
})

cat("==> [Step 1] Loading clustered Seurat object...\n")
pbmc <- readRDS("results/clustering/pbmc_clustered.rds")

# In Seurat v5, join layers if needed before DE
pbmc <- JoinLayers(pbmc)

cat("==> [Step 2] Finding differentially expressed markers for all clusters...\n")
cat("    (Testing genes detected in >= 25% of cells, log2FC threshold = 0.25)...\n")

pbmc_markers <- FindAllMarkers(
  pbmc,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

cat("==> Total markers identified:", nrow(pbmc_markers), "\n")

# Top 5 markers per cluster sorted by avg_log2FC
top5 <- pbmc_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

cat("==> Top 5 markers per cluster:\n")
print(top5[, c("cluster", "gene", "avg_log2FC", "p_val_adj")], n = 50)

# Save marker table
write.csv(pbmc_markers, file = "results/markers/all_cluster_markers.csv", row.names = FALSE)
cat("==> Saved marker table to 'results/markers/all_cluster_markers.csv'\n")

# Top 2 markers per cluster for DotPlot
top2_genes <- pbmc_markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC) %>%
  pull(gene) %>%
  unique()

cat("==> [Step 3] Generating publication DotPlot and Heatmap...\n")
dot_plot <- DotPlot(pbmc, features = top2_genes) + 
  RotatedAxis() +
  ggtitle("PBMC 3k - Top Cluster Specific Markers") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("figures/06_markers_dotplot.png", plot = dot_plot, width = 10, height = 5, dpi = 300)

# Heatmap for top 5 markers
top5_genes <- top5 %>% pull(gene) %>% unique()
heatmap_plot <- DoHeatmap(pbmc, features = top5_genes) + 
  theme(text = element_text(size = 9))

ggsave("figures/07_markers_heatmap.png", plot = heatmap_plot, width = 10, height = 8, dpi = 300)

cat("==> Biomarker discovery workflow completed successfully!\n")
