#!/usr/bin/env Rscript

raw_dir <- "raw_data"
archive <- file.path(raw_dir, "pbmc3k_filtered_gene_bc_matrices.tar.gz")
matrix_root <- file.path(raw_dir, "filtered_gene_bc_matrices")
dest <- file.path(matrix_root, "hg19")
url <- "https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
options(timeout = max(600, getOption("timeout")))

validate_gzip <- function(path) {
  if (!file.exists(path) || file.info(path)$size <= 0) return(FALSE)
  con <- gzfile(path, open = "rb")
  on.exit(close(con), add = TRUE)
  repeat {
    chunk <- readBin(con, what = "raw", n = 1024L * 1024L)
    if (!length(chunk)) break
  }
  TRUE
}

validate_matrix_dir <- function(path) {
  has_one <- function(stems) any(file.exists(file.path(path, stems)))
  has_one(c("matrix.mtx", "matrix.mtx.gz")) &&
    has_one(c("barcodes.tsv", "barcodes.tsv.gz")) &&
    has_one(c("genes.tsv", "genes.tsv.gz", "features.tsv", "features.tsv.gz"))
}

prepare_data <- function() {
  if (dir.exists(matrix_root)) {
    if (!validate_matrix_dir(dest)) {
      stop("An incomplete matrix directory already exists; inspect it before retrying: ", dest)
    }
    message("PBMC3k matrix is already present and its required files were found.")
    return(invisible(dest))
  }

  if (file.exists(archive)) {
    if (!validate_gzip(archive)) {
      stop("Existing archive is not a complete gzip file; preserve and inspect it: ", archive)
    }
  } else {
    tmp_archive <- paste0(archive, ".part")
    on.exit(if (file.exists(tmp_archive)) unlink(tmp_archive), add = TRUE)
    download.file(url, tmp_archive, mode = "wb", method = "libcurl", quiet = FALSE)
    if (!validate_gzip(tmp_archive)) stop("Downloaded archive failed full gzip integrity validation.")
    if (!file.rename(tmp_archive, archive)) stop("Could not publish verified archive to ", archive)
  }

  archive_files <- utils::untar(archive, list = TRUE)
  if (!length(archive_files) || any(grepl("(^|/)\\.\\.(/|$)", archive_files))) {
    stop("Archive is empty or contains unsafe relative paths.")
  }
  expected <- c(
    "(^|/)filtered_gene_bc_matrices/hg19/matrix\\.mtx(\\.gz)?$",
    "(^|/)filtered_gene_bc_matrices/hg19/barcodes\\.tsv(\\.gz)?$",
    "(^|/)filtered_gene_bc_matrices/hg19/(genes|features)\\.tsv(\\.gz)?$"
  )
  if (!all(vapply(expected, function(pattern) any(grepl(pattern, archive_files)), logical(1)))) {
    stop("Archive does not contain the expected hg19 matrix, barcode, and feature files.")
  }

  stage <- tempfile(pattern = "pbmc3k-extract-", tmpdir = raw_dir)
  if (!dir.create(stage)) stop("Could not create a temporary extraction directory.")
  on.exit(unlink(stage, recursive = TRUE, force = TRUE), add = TRUE)
  utils::untar(archive, exdir = stage)
  staged_matrix_root <- file.path(stage, "filtered_gene_bc_matrices")
  staged_dest <- file.path(staged_matrix_root, "hg19")
  if (!validate_matrix_dir(staged_dest)) {
    stop("Extracted archive is missing required hg19 matrix files.")
  }
  if (!file.rename(staged_matrix_root, matrix_root)) {
    stop("Could not publish the validated hg19 matrix directory.")
  }
  message("PBMC3k matrix downloaded, validated, and installed in ", dest)
  invisible(dest)
}

prepare_data()
