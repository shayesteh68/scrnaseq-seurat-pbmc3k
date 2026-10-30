suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(Matrix)
})

cat("==> [Step 1] Loading raw 10x Genomics matrix...\n")
data_dir <- "raw_data/filtered_gene_bc_matrices/hg19/"
if (!dir.exists(data_dir)) stop("Input folder missing: ", data_dir, ". Run scripts/00_download_data.R first.")
for (dir in c("figures", "results/qc", "results/clustering", "results/markers")) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
}
pbmc_data <- Read10X(data.dir = data_dir)
pbmc <- CreateSeuratObject(counts = pbmc_data, project = "PBMC3k", min.cells = 3, min.features = 200)

cat("==> Initial object summary:\n")
print(pbmc)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

qc_vln <- VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
ggsave("figures/01_qc_prefilter_vlnplot.png", qc_vln, width = 10, height = 5, dpi = 300, bg = "white")
qc_scatter <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ggsave("figures/01_qc_scatter.png", qc_scatter, width = 10, height = 5, dpi = 300, bg = "white")

pbmc_filtered <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
saveRDS(pbmc_filtered, "results/qc/pbmc_qc_filtered.rds")
cat("==> QC complete. Retained cells:", ncol(pbmc_filtered), "\n")
