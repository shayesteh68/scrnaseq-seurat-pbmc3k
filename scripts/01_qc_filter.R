suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(Matrix)
})

cat("==> [Step 1] Loading raw 10x Genomics matrix...\n")
data_dir <- "raw_data/filtered_gene_bc_matrices/hg19/"
pbmc_data <- Read10X(data.dir = data_dir)

# Initialize Seurat object
pbmc <- CreateSeuratObject(
  counts = pbmc_data,
  project = "PBMC3k",
  min.cells = 3,
  min.features = 200
)

cat("==> Initial Seurat Object Summary:\n")
print(pbmc)

cat("==> [Step 2] Calculating mitochondrial QC metrics...\n")
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Save pre-filter QC plots
cat("==> [Step 3] Generating pre-filter QC Violin Plots...\n")
qc_vln <- VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
ggsave("figures/01_qc_prefilter_vlnplot.png", plot = qc_vln, width = 10, height = 5, dpi = 300)

qc_scatter1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
qc_scatter2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
qc_scatter <- qc_scatter1 + qc_scatter2
ggsave("figures/01_qc_scatter.png", plot = qc_scatter, width = 10, height = 5, dpi = 300)

cat("==> [Step 4] Applying QC Filters (nFeature: 200-2500, percent.mt < 5%)...\n")
pbmc_filtered <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

cat("==> Filtered Seurat Object Summary:\n")
print(pbmc_filtered)

# Save filtered object
cat("==> [Step 5] Saving clean Seurat object...\n")
saveRDS(pbmc_filtered, file = "results/qc/pbmc_qc_filtered.rds")
cat("==> QC & Filtering completed successfully! Output saved to results/qc/pbmc_qc_filtered.rds\n")
