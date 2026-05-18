#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(utils)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
    if (!(flag %in% args)) {
        return(default)
    }
    index <- match(flag, args)
    if (index == length(args)) {
        stop(sprintf("Missing value for %s", flag), call. = FALSE)
    }
    args[[index + 1]]
}

output_dir <- get_arg("--output-dir", file.path("data", "raw", "dlpfc"))
sample_ids_arg <- get_arg("--sample-ids", "151507,151508,151673,151674")
sample_ids <- unlist(strsplit(sample_ids_arg, ",", fixed = TRUE))
sample_ids <- trimws(sample_ids)
sample_ids <- sample_ids[nzchar(sample_ids)]
r_lib <- get_arg("--r-lib", file.path(".r-lib"))

dir.create(r_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(normalizePath(r_lib), .libPaths()))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

required_packages <- c("spatialLIBD", "SingleCellExperiment", "SpatialExperiment", "Matrix")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
    BiocManager::install(missing_packages, ask = FALSE, update = FALSE)
}

suppressPackageStartupMessages({
    library(spatialLIBD)
    library(SingleCellExperiment)
    library(Matrix)
})

message("Downloading DLPFC SpatialExperiment from spatialLIBD::fetch_data('spe')")
spe <- spatialLIBD::fetch_data(type = "spe")

if ("sample_id" %in% names(colData(spe))) {
    keep <- colData(spe)$sample_id %in% sample_ids
    spe <- spe[, keep]
} else {
    stop("sample_id column is missing from fetched DLPFC object", call. = FALSE)
}

if ("in_tissue" %in% names(colData(spe))) {
    spe <- spe[, spe$in_tissue == 1]
}

counts_name <- if ("counts" %in% assayNames(spe)) "counts" else assayNames(spe)[[1]]
counts <- assay(spe, counts_name)
if (!inherits(counts, "sparseMatrix")) {
    counts <- Matrix(counts, sparse = TRUE)
}

output_dir <- normalizePath(output_dir, mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

Matrix::writeMM(counts, file.path(output_dir, "counts.mtx"))
write.csv(as.data.frame(colData(spe)), file.path(output_dir, "obs.csv"), quote = TRUE)
write.csv(as.data.frame(rowData(spe)), file.path(output_dir, "var.csv"), quote = TRUE)
write.csv(as.data.frame(spatialCoords(spe)), file.path(output_dir, "spatial.csv"), quote = TRUE)
writeLines(colnames(spe), file.path(output_dir, "obs_names.txt"))
writeLines(rownames(spe), file.path(output_dir, "var_names.txt"))
writeLines(sample_ids, file.path(output_dir, "selected_sample_ids.txt"))

message(sprintf("Wrote DLPFC export to %s", output_dir))
