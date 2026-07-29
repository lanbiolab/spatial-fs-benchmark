#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(Matrix)
    library(matrixStats)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7) {
    stop("Expected: input_dir output_path method max_cells_per_slice max_genes n_threads seed")
}

input_dir <- args[[1]]
output_path <- args[[2]]
method <- tolower(args[[3]])
max_cells_per_slice <- as.integer(args[[4]])
max_genes <- as.integer(args[[5]])
n_threads <- as.integer(args[[6]])
seed <- as.integer(args[[7]])

matrix <- readMM(file.path(input_dir, "matrix.mtx"))
genes <- readLines(file.path(input_dir, "genes.tsv"), warn = FALSE)
observations <- read.delim(file.path(input_dir, "observations.tsv"), check.names = FALSE)
rownames(matrix) <- genes

row_variances_sparse <- function(x) {
    means <- Matrix::rowMeans(x)
    means_sq <- Matrix::rowMeans(x ^ 2)
    pmax(means_sq - means ^ 2, 0)
}

percentile_scores <- function(statistic, decreasing = TRUE) {
    finite <- is.finite(statistic)
    scores <- rep(0, length(statistic))
    if (any(finite)) {
        ranks <- rank(statistic[finite], ties.method = "average", na.last = "keep")
        if (decreasing) {
            scores[finite] <- (sum(finite) - ranks + 1) / sum(finite)
        } else {
            scores[finite] <- ranks / sum(finite)
        }
    }
    scores
}

slice_ids <- sort(unique(as.character(observations$slice)))
slice_scores <- list()
for (slice_number in seq_along(slice_ids)) {
    slice_id <- slice_ids[[slice_number]]
    indices <- which(as.character(observations$slice) == slice_id)
    if (max_cells_per_slice > 0 && length(indices) > max_cells_per_slice) {
        set.seed(seed + slice_number * 1009L)
        indices <- sort(sample(indices, max_cells_per_slice, replace = FALSE))
    }
    minimum_cells <- if (method == "nnsvg") 65L else 10L
    if (length(indices) < minimum_cells) {
        next
    }

    slice_matrix <- matrix[, indices, drop = FALSE]
    coords <- as.matrix(observations[indices, c("x", "y"), drop = FALSE])
    scores <- rep(0, length(genes))
    if (method == "sparkx") {
        suppressPackageStartupMessages(library(SPARK))
        candidates <- which(!grepl("^(MT-|mt-)", genes))
        if (max_genes > 0 && length(candidates) > max_genes) {
            variances <- row_variances_sparse(slice_matrix[candidates, , drop = FALSE])
            candidates <- candidates[order(variances, decreasing = TRUE)[seq_len(max_genes)]]
        }
        result <- SPARK::sparkx(
            slice_matrix[candidates, , drop = FALSE],
            coords,
            numCores = n_threads,
            option = "mixture",
            verbose = FALSE
        )
        pvalues <- result$res_mtest$combinedPval
        names(pvalues) <- rownames(result$res_mtest)
        matched <- match(names(pvalues), genes)
        valid <- !is.na(matched) & is.finite(pvalues)
        scores[matched[valid]] <- percentile_scores(
            -log10(pmax(pvalues[valid], .Machine$double.xmin)),
            decreasing = FALSE
        )
    } else if (method == "nnsvg") {
        suppressPackageStartupMessages({
            library(nnSVG)
            library(scuttle)
            library(SpatialExperiment)
        })
        colnames(slice_matrix) <- paste0("spot", seq_len(ncol(slice_matrix)))
        spe <- SpatialExperiment::SpatialExperiment(
            assays = list(counts = slice_matrix),
            rowData = S4Vectors::DataFrame(gene_name = genes),
            spatialCoords = coords
        )
        spe <- nnSVG::filter_genes(
            spe,
            filter_genes_ncounts = 3,
            filter_genes_pcspots = 0.5,
            filter_mito = TRUE
        )
        zero_spots <- Matrix::colSums(SummarizedExperiment::assay(spe, "counts")) == 0
        if (any(zero_spots)) {
            spe <- spe[, !zero_spots]
        }
        if (nrow(spe) < 2 || ncol(spe) < minimum_cells) {
            next
        }
        spe <- scuttle::computeLibraryFactors(spe)
        spe <- scuttle::logNormCounts(spe)
        if (max_genes > 0 && nrow(spe) > max_genes) {
            logcounts <- SummarizedExperiment::assay(spe, "logcounts")
            variances <- row_variances_sparse(logcounts)
            keep <- order(variances, decreasing = TRUE)[seq_len(max_genes)]
            spe <- spe[keep, ]
        }
        set.seed(seed + slice_number * 1009L)
        result <- nnSVG::nnSVG(
            spe,
            assay_name = "logcounts",
            n_neighbors = 10L,
            n_threads = n_threads,
            verbose = FALSE
        )
        result_table <- as.data.frame(SummarizedExperiment::rowData(result))
        result_genes <- rownames(result_table)
        matched <- match(result_genes, genes)
        valid <- !is.na(matched) & is.finite(result_table$rank)
        scores[matched[valid]] <- percentile_scores(-result_table$rank[valid], decreasing = FALSE)
    } else {
        stop("Unsupported method: ", method)
    }
    slice_scores[[slice_id]] <- scores
}

if (length(slice_scores) == 0) {
    stop("No slice produced a valid ", method, " ranking")
}

score_matrix <- do.call(cbind, slice_scores)
aggregate_score <- rowMeans(score_matrix, na.rm = TRUE)
output <- data.frame(Feature = genes, Score = aggregate_score)
output <- output[order(output$Score, decreasing = TRUE), , drop = FALSE]
write.table(output, output_path, sep = "\t", quote = FALSE, row.names = FALSE)
