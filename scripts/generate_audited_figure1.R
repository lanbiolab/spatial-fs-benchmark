#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(patchwork)
    library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[[1]] else file.path("results", "atlas_style_final_all_slice", "data")
output_dir <- if (length(args) >= 2) args[[2]] else file.path("results", "atlas_style_audit", "figures")
base_dir <- normalizePath(file.path(data_dir, ".."), mustWork = FALSE)
ranges_path <- file.path(base_dir, "output", "baseline-ranges.tsv")

source(file.path("analysis", "atlas_style", "R", "spatial_summarisation.R"))
source(file.path("analysis", "atlas_style", "R", "spatial_plotting.R"))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

metrics <- readr::read_tsv(file.path(data_dir, "benchmark.tsv"), show_col_types = FALSE)
methods_meta_all <- readr::read_tsv(file.path(data_dir, "methods-metadata.tsv"), show_col_types = FALSE)
metric_ranges <- readr::read_tsv(ranges_path, show_col_types = FALSE)

metrics_summary <- summarise_spatial_metrics(
    metrics,
    metric_ranges,
    type_weights = c(
        "Integration" = 1 / 3,
        "Clustering" = 1 / 3,
        "Alignment" = 1 / 3,
        "SliceRepresentation" = 0
    )
)

figure1 <- patchwork::wrap_plots(
    plot_overview_heatmap(metrics_summary, methods_meta_all),
    plot_overview_ranking(metrics_summary, methods_meta_all),
    ncol = 2,
    widths = c(2, 1)
)

save_figure_files(figure1, file.path(output_dir, "figure1_overview_audited"), width = 8.2, height = 9.2)
