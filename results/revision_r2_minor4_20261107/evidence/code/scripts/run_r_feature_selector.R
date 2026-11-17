#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    options(stringsAsFactors = FALSE)
})

future_maxsize_env <- Sys.getenv("R_FUTURE_GLOBALS_MAXSIZE", unset = NA_character_)
if (!is.na(future_maxsize_env) && nzchar(future_maxsize_env)) {
    future_maxsize_num <- suppressWarnings(as.numeric(future_maxsize_env))
    if (is.finite(future_maxsize_num) && future_maxsize_num > 0) {
        options(future.globals.maxSize = future_maxsize_num)
    }
}

parse_args <- function() {
    raw <- commandArgs(trailingOnly = TRUE)
    out <- list()
    i <- 1
    while (i <= length(raw)) {
        key <- raw[[i]]
        if (!startsWith(key, "--")) {
            stop("Unexpected positional argument: ", key, call. = FALSE)
        }
        if (i == length(raw)) {
            stop("Missing value for ", key, call. = FALSE)
        }
        out[[substring(key, 3)]] <- raw[[i + 1]]
        i <- i + 2
    }
    required <- c("input-dir", "output", "method", "n-features")
    missing <- required[!required %in% names(out)]
    if (length(missing) > 0) {
        stop("Missing required arguments: ", paste(missing, collapse = ", "), call. = FALSE)
    }
    out[["n-features"]] <- as.numeric(out[["n-features"]])
    out
}

ensure_namespace <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        stop("Missing required R package: ", pkg, call. = FALSE)
    }
}

read_input_sce <- function(input_dir) {
    ensure_namespace("Matrix")
    ensure_namespace("SingleCellExperiment")
    ensure_namespace("S4Vectors")
    input_dir <- normalizePath(input_dir, mustWork = TRUE)
    counts <- Matrix::readMM(file.path(input_dir, "counts.mtx"))
    counts <- methods::as(counts, "dgCMatrix")
    storage.mode(counts@x) <- "double"
    genes <- read.delim(
        file.path(input_dir, "genes.tsv"),
        sep = "\t",
        header = FALSE,
        stringsAsFactors = FALSE
    )[[1]]
    cells <- read.delim(
        file.path(input_dir, "cells.tsv"),
        sep = "\t",
        header = FALSE,
        stringsAsFactors = FALSE
    )[[1]]
    obs <- read.delim(
        file.path(input_dir, "obs.tsv"),
        sep = "\t",
        header = TRUE,
        row.names = 1,
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    rownames(counts) <- as.character(genes)
    colnames(counts) <- as.character(cells)
    if (!identical(colnames(counts), rownames(obs))) {
        obs <- obs[colnames(counts), , drop = FALSE]
    }
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts),
        colData = S4Vectors::DataFrame(obs)
    )
    sce
}

sanitize_sce_counts <- function(sce, drop_zero_cells = TRUE) {
    ensure_namespace("Matrix")
    counts_mat <- SingleCellExperiment::counts(sce)
    if (inherits(counts_mat, "sparseMatrix")) {
        counts_mat <- methods::as(counts_mat, "dgCMatrix")
        storage.mode(counts_mat@x) <- "double"
    } else {
        counts_mat <- as.matrix(counts_mat)
        storage.mode(counts_mat) <- "double"
    }
    min_count <- suppressWarnings(min(counts_mat, na.rm = TRUE))
    if (is.finite(min_count) && min_count < 0) {
        counts_mat <- counts_mat - min_count
    }
    if (drop_zero_cells) {
        libsizes <- Matrix::colSums(counts_mat)
        keep <- is.finite(libsizes) & (libsizes > 0)
        if (!all(keep)) {
            sce <- sce[, keep, drop = FALSE]
            counts_mat <- counts_mat[, keep, drop = FALSE]
        }
    }
    SingleCellExperiment::counts(sce) <- counts_mat
    sce
}

to_seurat <- function(sce) {
    ensure_namespace("Seurat")
    ensure_namespace("SeuratObject")
    SeuratObject::as.Seurat(sce, counts = "counts", data = NULL)
}

get_seurat_data <- function(seurat, layer = "data") {
    if ("LayerData" %in% getNamespaceExports("SeuratObject")) {
        return(SeuratObject::LayerData(seurat, assay = Seurat::DefaultAssay(seurat), layer = layer))
    }
    SeuratObject::GetAssayData(seurat, assay = Seurat::DefaultAssay(seurat), slot = layer)
}

rank_output <- function(features, scores = NULL) {
    features <- as.character(features)
    keep <- !is.na(features) & nzchar(features)
    features <- features[keep]
    if (is.null(scores)) {
        scores <- rev(seq_along(features))
    } else {
        scores <- as.numeric(scores)[keep]
    }
    output <- data.frame(Feature = features, Score = scores)
    output <- output[!duplicated(output$Feature), , drop = FALSE]
    if (nrow(output) == 0) {
        stop("No features were selected.", call. = FALSE)
    }
    output[order(output$Score, decreasing = TRUE), , drop = FALSE]
}

first_existing <- function(x, candidates) {
    hits <- candidates[candidates %in% colnames(x)]
    if (length(hits) == 0) {
        return(NULL)
    }
    hits[[1]]
}

rank_hvf_info <- function(seurat, n_features) {
    hvf <- tryCatch(
        Seurat::HVFInfo(seurat),
        error = function(...) NULL
    )
    if (is.null(hvf) || nrow(hvf) == 0) {
        return(NULL)
    }
    hvf$Feature <- rownames(hvf)
    score_col <- first_existing(
        hvf,
        c(
            "variance.standardized",
            "dispersion.scaled",
            "dispersion",
            "residual_variance",
            "vst.variance.standardized",
            "mvp.dispersion",
            "mvp.dispersion.scaled"
        )
    )
    scores <- if (!is.null(score_col)) hvf[[score_col]] else rev(seq_len(nrow(hvf)))
    hvf <- hvf[order(scores, decreasing = TRUE), , drop = FALSE]
    hvf <- head(hvf, n_features)
    rank_output(hvf$Feature, head(scores, nrow(hvf)))
}

select_seurat <- function(sce, n_features, method) {
    if (identical(method, "seurat_sct")) {
        sce <- sanitize_sce_counts(sce)
    }
    seurat <- to_seurat(sce)
    method <- switch(
        method,
        seurat_vst = "vst",
        seurat_mvp = "mean.var.plot",
        seurat_disp = "dispersion",
        seurat_sct = "sctransform",
        stop("Unsupported Seurat method: ", method, call. = FALSE)
    )
    if (method == "sctransform") {
        ensure_namespace("future")
        Sys.setenv(R_FUTURE_GLOBALS_MAXSIZE = as.character(8 * 1024^3))
        old_plan <- future::plan()
        old_maxsize <- getOption("future.globals.maxSize")
        old_maxsize_num <- if (is.null(old_maxsize)) 0 else as.numeric(old_maxsize)
        on.exit({
            future::plan(old_plan)
            options(future.globals.maxSize = old_maxsize)
        }, add = TRUE)
        future::plan(future::sequential)
        options(future.globals.maxSize = max(8 * 1024^3, old_maxsize_num))
        counts <- get_seurat_data(seurat, layer = "counts")
        libsizes <- Matrix::colSums(counts)
        keep <- is.finite(libsizes) & (libsizes > 0)
        if (!all(keep)) {
            seurat <- subset(seurat, cells = colnames(seurat)[keep])
        }
        seurat <- Seurat::SCTransform(
            seurat,
            assay = Seurat::DefaultAssay(seurat),
            variable.features.n = n_features,
            ncells = min(5000, ncol(seurat)),
            do.correct.umi = FALSE,
            conserve.memory = TRUE,
            verbose = FALSE
        )
    } else {
        seurat <- Seurat::NormalizeData(seurat, verbose = FALSE)
        seurat <- Seurat::FindVariableFeatures(
            seurat,
            selection.method = method,
            nfeatures = n_features,
            verbose = FALSE
        )
    }
    selected <- Seurat::VariableFeatures(seurat)
    if (length(selected) > 0) {
        return(rank_output(selected))
    }
    hvf_ranked <- rank_hvf_info(seurat, n_features)
    if (!is.null(hvf_ranked)) {
        return(hvf_ranked)
    }
    stop("No features were selected.", call. = FALSE)
}

select_scsegindex <- function(sce, n_features) {
    ensure_namespace("scMerge")
    ensure_namespace("scuttle")
    ensure_namespace("scPNMF")
    suppressPackageStartupMessages(library(scPNMF))
    sce <- sanitize_sce_counts(sce)
    counts <- SingleCellExperiment::counts(sce)
    logcounts <- scuttle::normalizeCounts(counts)
    result <- scMerge::scSEGIndex(logcounts)
    result$Feature <- rownames(result)
    result <- result[order(result$segIdx, decreasing = TRUE), , drop = FALSE]
    result <- head(result, n_features)
    rank_output(result$Feature, result$segIdx)
}

select_nbumi <- function(sce, n_features) {
    ensure_namespace("M3Drop")
    sce <- sanitize_sce_counts(sce)
    counts <- SingleCellExperiment::counts(sce)
    counts_dense <- as.matrix(counts)
    min_count <- suppressWarnings(min(counts_dense, na.rm = TRUE))
    if (is.finite(min_count) && min_count < 0) {
        counts_dense <- counts_dense - min_count
    }
    counts_dense <- round(pmax(counts_dense, 0))
    count_mat <- M3Drop::NBumiConvertData(counts_dense, is.counts = TRUE)
    danb_fit <- M3Drop::NBumiFitModel(count_mat)
    result <- M3Drop::NBumiFeatureSelectionCombinedDrop(
        danb_fit,
        qval.thres = 0.01,
        suppress.plot = TRUE
    )
    if (nrow(result) < n_features) {
        result <- M3Drop::NBumiFeatureSelectionCombinedDrop(
            danb_fit,
            ntop = max(500, n_features),
            suppress.plot = TRUE
        )
    }
    result$Feature <- rownames(result)
    score_col <- first_existing(result, c("q.value", "q.value.combined", "p.value", "combined.p.value"))
    if (!is.null(score_col)) {
        scores <- -log10(pmax(result[[score_col]], 1e-300))
    } else {
        scores <- rev(seq_len(nrow(result)))
    }
    result <- result[order(scores, decreasing = TRUE), , drop = FALSE]
    result <- head(result, n_features)
    rank_output(result$Feature, head(scores, nrow(result)))
}

select_osca <- function(sce, n_features) {
    ensure_namespace("batchelor")
    ensure_namespace("scran")
    sce <- sanitize_sce_counts(sce)
    batch <- if ("Batch" %in% colnames(SummarizedExperiment::colData(sce))) sce$Batch else factor("batch")
    sce <- batchelor::multiBatchNorm(sce, batch = batch)
    feature_stats <- scran::modelGeneVar(sce, block = batch)
    top_hvgs <- scran::getTopHVGs(feature_stats, n = n_features)
    result <- feature_stats[top_hvgs, , drop = FALSE]
    score_col <- first_existing(as.data.frame(result), c("bio", "total", "ratio"))
    scores <- if (!is.null(score_col)) result[[score_col]] else rev(seq_len(length(top_hvgs)))
    rank_output(top_hvgs, scores)
}

select_dubstepr <- function(sce, n_features) {
    ensure_namespace("DUBStepR")
    seurat <- to_seurat(sce)
    seurat <- Seurat::NormalizeData(seurat, verbose = FALSE)
    result_list <- DUBStepR::DUBStepR(
        get_seurat_data(seurat, layer = "data"),
        optimise.features = FALSE
    )
    corr_info <- result_list$corr.info
    selected <- if (!is.null(result_list$optimal.feature.genes) && length(result_list$optimal.feature.genes) > 0) {
        result_list$optimal.feature.genes
    } else {
        rownames(corr_info)
    }
    if (length(selected) == 0) {
        stop("DUBStepR returned no selected features.", call. = FALSE)
    }
    corr_info <- corr_info[selected, , drop = FALSE]
    corr_info$Feature <- rownames(corr_info)
    numeric_cols <- colnames(corr_info)[vapply(corr_info, is.numeric, logical(1))]
    scores <- if (length(numeric_cols) > 0) corr_info[[numeric_cols[[1]]]] else rev(seq_len(nrow(corr_info)))
    corr_info <- corr_info[order(scores, decreasing = TRUE), , drop = FALSE]
    corr_info <- head(corr_info, n_features)
    rank_output(corr_info$Feature, head(scores, nrow(corr_info)))
}

select_scry <- function(sce, n_features) {
    ensure_namespace("scry")
    sce <- scry::devianceFeatureSelection(sce, nkeep = n_features)
    output <- as.data.frame(SummarizedExperiment::rowData(sce))
    output$Feature <- rownames(sce)
    score_col <- first_existing(output, c("binomial_deviance", "deviance", "residual_variance"))
    scores <- if (!is.null(score_col)) output[[score_col]] else rev(seq_len(nrow(output)))
    output <- output[order(scores, decreasing = TRUE), , drop = FALSE]
    output <- head(output, n_features)
    rank_output(output$Feature, head(scores, nrow(output)))
}

select_singlecellhaystack <- function(sce, n_features) {
    ensure_namespace("singleCellHaystack")
    seurat <- to_seurat(sce)
    seurat <- Seurat::NormalizeData(seurat, verbose = FALSE)
    seurat <- Seurat::FindVariableFeatures(seurat, verbose = FALSE)
    seurat <- Seurat::ScaleData(seurat, verbose = FALSE)
    seurat <- Seurat::RunPCA(seurat, features = Seurat::VariableFeatures(seurat), verbose = FALSE)
    coords <- Seurat::Embeddings(seurat, reduction = "pca")
    detection <- as.matrix(get_seurat_data(seurat, layer = "data") > 0)
    results <- singleCellHaystack::haystack(
        coords,
        detection = detection,
        method = "highD",
        grid.points = min(25, max(5, nrow(coords) - 1))
    )
    top_results <- singleCellHaystack::show_result_haystack(results, n = n_features)
    top_results$Feature <- rownames(top_results)
    score_col <- first_existing(top_results, c("D_KL", "logpval", "pval_random"))
    if (!is.null(score_col) && score_col == "pval_random") {
        scores <- -log10(pmax(top_results[[score_col]], 1e-300))
    } else if (!is.null(score_col)) {
        scores <- top_results[[score_col]]
    } else {
        scores <- rev(seq_len(nrow(top_results)))
    }
    rank_output(top_results$Feature, scores)
}

select_brennecke <- function(sce, n_features) {
    ensure_namespace("scuttle")
    ensure_namespace("scran")
    ensure_namespace("Matrix")
    sce <- sanitize_sce_counts(sce)
    counts <- SingleCellExperiment::counts(sce)
    keep_genes <- Matrix::rowSums(counts > 0) > 1 & Matrix::rowSums(counts) > 0
    if (sum(keep_genes) >= min(100, n_features)) {
        sce <- sce[keep_genes, , drop = FALSE]
    }
    sce <- scuttle::logNormCounts(sce)
    feature_stats <- scran::modelGeneCV2(sce)
    top_hvgs <- scran::getTopHVGs(feature_stats, var.field = "ratio", n = n_features)
    result <- feature_stats[top_hvgs, , drop = FALSE]
    scores <- result$ratio
    rank_output(top_hvgs, scores)
}

select_scpnmf <- function(sce, n_features) {
    ensure_namespace("scuttle")
    ensure_namespace("scPNMF")
    sce <- sanitize_sce_counts(sce)
    sce <- scuttle::logNormCounts(sce)
    logcounts <- as.matrix(SingleCellExperiment::logcounts(sce))
    pnmf <- scPNMF::PNMFfun(
        logcounts,
        K = 20,
        method = "EucDist",
        tol = 1e-4,
        maxIter = 1000,
        verboseN = FALSE
    )
    w_select <- scPNMF::basisSelect(
        W = pnmf$Weight,
        S = pnmf$Score,
        X = logcounts,
        toTest = TRUE,
        toAnnotate = FALSE,
        mc.cores = 1
    )
    genes <- scPNMF::getInfoGene(
        w_select,
        M = n_features,
        by_basis = FALSE,
        return_trunW = TRUE,
        dim_use = NULL
    )
    rank_output(genes$InfoGene)
}

run_selector <- function(method, sce, n_features) {
    switch(
        method,
        seurat_vst = select_seurat(sce, n_features, method),
        seurat_mvp = select_seurat(sce, n_features, method),
        seurat_disp = select_seurat(sce, n_features, method),
        seurat_sct = select_seurat(sce, n_features, method),
        scsegindex = select_scsegindex(sce, n_features),
        nbumi = select_nbumi(sce, n_features),
        osca = select_osca(sce, n_features),
        dubstepr = select_dubstepr(sce, n_features),
        scry = select_scry(sce, n_features),
        singleCellHaystack = select_singlecellhaystack(sce, n_features),
        Brennecke = select_brennecke(sce, n_features),
        scPNMF = select_scpnmf(sce, n_features),
        stop("Unsupported R selector method: ", method, call. = FALSE)
    )
}

main <- function() {
    cli_args <- parse_args()
    sce <- read_input_sce(cli_args[["input-dir"]])
    output <- run_selector(cli_args[["method"]], sce, cli_args[["n-features"]])
    write.table(
        output,
        file = cli_args[["output"]],
        quote = FALSE,
        sep = "\t",
        row.names = FALSE
    )
}

main()
