suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
})

cat("==> [Step 1] Loading clustered Seurat object...\n")
pbmc <- readRDS("results/clustering/pbmc_clustered.rds")

cat("==> [Step 2] Assigning biological cell type labels...\n")
new_cluster_ids <- c(
  "0" = "Naive CD4+ T",
  "1" = "CD14+ Mono",
  "2" = "Memory CD4+ T",
  "3" = "B",
  "4" = "CD8+ T",
  "5" = "FCGR3A+ Mono",
  "6" = "NK",
  "7" = "DC",
  "8" = "Platelet"
)

observed_clusters <- unique(as.character(Idents(pbmc)))
if (!setequal(observed_clusters, names(new_cluster_ids))) {
  missing_clusters <- setdiff(observed_clusters, names(new_cluster_ids))
  unexpected_map <- setdiff(names(new_cluster_ids), observed_clusters)
  stop(
    "Cluster IDs differ from the reviewed marker map. Missing labels: ",
    paste(missing_clusters, collapse = ", "),
    "; unused map entries: ", paste(unexpected_map, collapse = ", ")
  )
}
dir.create("figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results", recursive = TRUE, showWarnings = FALSE)
pbmc <- RenameIdents(pbmc, new_cluster_ids)
pbmc$cell_type <- Idents(pbmc)
pbmc$annotation_status <- "provisional_manual_marker_reviewed"

cat("==> Cell count summary per annotated cell type:\n")
print(table(pbmc$cell_type))

cat("==> [Step 3] Generating annotated UMAP plot...\n")
umap_annotated <- DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.6, repel = TRUE) +
  ggtitle("PBMC 3k - Cell Type Annotation (Seurat v5)") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position = "right"
  )

ggsave("figures/08_umap_annotated_celltypes.png", plot = umap_annotated, width = 8, height = 6, dpi = 300)

cat("==> [Step 4] Generating Cell Type Proportion Barplot...\n")
cell_counts <- as.data.frame(table(pbmc$cell_type))
colnames(cell_counts) <- c("CellType", "Count")
cell_counts$Percentage <- (cell_counts$Count / sum(cell_counts$Count)) * 100

prop_plot <- ggplot(cell_counts, aes(x = reorder(CellType, -Percentage), y = Percentage, fill = CellType)) +
  geom_col(show.legend = FALSE, color = "black", width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", Percentage)), vjust = -0.5, size = 3.5, fontface = "bold") +
  theme_classic() +
  labs(
    title = "PBMC 3k - Cell Type Proportions",
    x = "Cell Type",
    y = "Percentage (%)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color = "black")
  ) +
  ylim(0, max(cell_counts$Percentage) + 5)

ggsave("figures/09_cell_type_proportions.png", plot = prop_plot, width = 8, height = 5, dpi = 300)

cat("==> [Step 5] Saving final annotated Seurat object...\n")
saveRDS(pbmc, file = "results/final_annotated_pbmc.rds")

cat("==> Full scRNA-seq pipeline executed successfully!\n")
