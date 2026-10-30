suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
})

cat("==> [1/5] Loading clustered object and reviewed markers...\n")
pbmc <- readRDS("results/clustering/pbmc_clustered.rds")
marker_path <- "results/markers/all_cluster_markers.csv"
if (!file.exists(marker_path)) stop("Run scripts/04_find_markers.R before annotation.")
marker_results <- read.csv(marker_path, stringsAsFactors = FALSE)

new_cluster_ids <- c(
  "0" = "Naive CD4 T",
  "1" = "CD14+ monocyte",
  "2" = "Memory-like CD4 T",
  "3" = "B cell",
  "4" = "CD8 T",
  "5" = "FCGR3A+ monocyte",
  "6" = "NK cell",
  "7" = "Dendritic cell",
  "8" = "Platelet"
)
expected_markers <- list(
  "0" = c("CCR7", "LEF1", "TCF7", "IL7R"),
  "1" = c("CD14", "LYZ", "S100A8", "S100A9", "FCN1"),
  "2" = c("CD3D", "CD3E", "IL7R", "CD40LG", "AQP3"),
  "3" = c("MS4A1", "CD79A", "CD79B", "CD74"),
  "4" = c("CD8A", "CCL5", "GZMK", "CD3D", "NKG7"),
  "5" = c("FCGR3A", "MS4A7", "LST1", "IFITM3", "LILRB1"),
  "6" = c("GNLY", "PRF1", "NKG7", "GZMB", "KLRD1"),
  "7" = c("FCER1A", "CLEC10A", "CD1C", "CD74"),
  "8" = c("ITGA2B", "GP9", "PF4", "PPBP", "GNG11")
)
observed_clusters <- sort(unique(as.character(Idents(pbmc))))
if (!setequal(observed_clusters, names(new_cluster_ids))) {
  missing_clusters <- setdiff(observed_clusters, names(new_cluster_ids))
  unused_map <- setdiff(names(new_cluster_ids), observed_clusters)
  stop(
    "Cluster IDs differ from the reviewed marker map. Missing labels: ",
    paste(missing_clusters, collapse = ", "),
    "; unused map entries: ", paste(unused_map, collapse = ", ")
  )
}

cluster_sizes <- table(factor(as.character(Idents(pbmc)), levels = names(new_cluster_ids)))
small_cluster_cutoff <- 50L
evidence_rows <- lapply(names(new_cluster_ids), function(cluster_id) {
  cluster_table <- marker_results[
    as.character(marker_results$cluster) == cluster_id &
      !is.na(marker_results$p_val_adj) & marker_results$p_val_adj < 0.05 &
      !is.na(marker_results$avg_log2FC) & marker_results$avg_log2FC > 0,
    ,
    drop = FALSE
  ]
  observed <- unique(cluster_table$gene)
  supported <- intersect(expected_markers[[cluster_id]], observed)
  data.frame(
    cluster_id = cluster_id,
    cell_type = unname(new_cluster_ids[[cluster_id]]),
    cell_count = as.integer(cluster_sizes[[cluster_id]]),
    expected_markers = paste(expected_markers[[cluster_id]], collapse = ";"),
    observed_markers = paste(supported, collapse = ";"),
    evidence_count = length(supported),
    annotation_status = if (length(supported) < 2L) {
      "insufficient_marker_evidence"
    } else if (as.integer(cluster_sizes[[cluster_id]]) < small_cluster_cutoff) {
      "provisional_small_cluster"
    } else {
      "provisional_marker_supported"
    },
    stringsAsFactors = FALSE
  )
})
annotation_evidence <- do.call(rbind, evidence_rows)
if (any(annotation_evidence$evidence_count < 2L)) {
  failed <- annotation_evidence$cluster_id[annotation_evidence$evidence_count < 2L]
  stop(
    "At least two expected significant markers are required per cluster. Review cluster(s): ",
    paste(failed, collapse = ", ")
  )
}

cat("==> Marker support per cluster:\n")
print(annotation_evidence[, c("cluster_id", "cell_type", "cell_count", "observed_markers", "annotation_status")], row.names = FALSE)

dir.create("figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/annotation", recursive = TRUE, showWarnings = FALSE)
write.csv(annotation_evidence, "results/annotation/cluster_annotation_evidence.csv", row.names = FALSE)

pbmc <- RenameIdents(pbmc, new_cluster_ids)
pbmc$cell_type <- Idents(pbmc)
status_by_type <- setNames(annotation_evidence$annotation_status, annotation_evidence$cell_type)
pbmc$annotation_status <- unname(status_by_type[as.character(pbmc$cell_type)])

cat("==> [2/5] Cell counts per annotated cell type:\n")
cell_counts <- as.data.frame(table(pbmc$cell_type), stringsAsFactors = FALSE)
colnames(cell_counts) <- c("cell_type", "cell_count")
cell_counts$percentage <- 100 * cell_counts$cell_count / sum(cell_counts$cell_count)
cell_counts$cell_type <- as.character(cell_counts$cell_type)
write.csv(cell_counts, "results/annotation/cell_type_proportions.csv", row.names = FALSE)
print(cell_counts, row.names = FALSE)

cat("==> [3/5] Generating annotated UMAP...\n")
umap_annotated <- DimPlot(pbmc, reduction = "umap", group.by = "cell_type", label = FALSE, pt.size = 0.6) +
  ggtitle("PBMC 3k: provisional marker-reviewed labels") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position = "right"
  )
ggsave("figures/08_umap_annotated_celltypes.png", plot = umap_annotated, width = 9, height = 6.5, dpi = 300, bg = "white")

cat("==> [4/5] Generating cell-type proportions...\n")
prop_plot <- ggplot(cell_counts, aes(x = reorder(cell_type, -percentage), y = percentage, fill = cell_type)) +
  geom_col(show.legend = FALSE, color = "black", width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", percentage)), vjust = -0.5, size = 3.5, fontface = "bold") +
  theme_classic() +
  labs(
    title = "PBMC 3k: cell-type proportions",
    x = "Cell type",
    y = "Cells (%)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color = "black")
  ) +
  coord_cartesian(ylim = c(0, max(cell_counts$percentage) + 5), clip = "off")
ggsave("figures/09_cell_type_proportions.png", plot = prop_plot, width = 9, height = 5.5, dpi = 300, bg = "white")

cat("==> [5/5] Saving final annotated Seurat object...\n")
saveRDS(pbmc, file = "results/final_annotated_pbmc.rds")
cat("==> Annotation completed; labels remain provisional for this single-donor dataset.\n")
