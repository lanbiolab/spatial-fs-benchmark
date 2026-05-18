#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[[1]] else file.path("results", "atlas_style", "data")
output_dir <- if (length(args) >= 2) args[[2]] else file.path("results", "atlas_style", "figures")
base_dir <- normalizePath(file.path(data_dir, ".."), mustWork = FALSE)
ranges_path <- file.path(base_dir, "output", "baseline-ranges.tsv")

source(file.path("analysis", "atlas_style", "R", "spatial_summarisation.R"))
source(file.path("analysis", "atlas_style", "R", "spatial_plotting.R"))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

metrics <- readr::read_tsv(file.path(data_dir, "benchmark.tsv"), show_col_types = FALSE)
datasets_meta <- readr::read_tsv(file.path(data_dir, "datasets-metadata.tsv"), show_col_types = FALSE)
methods_meta <- readr::read_tsv(file.path(data_dir, "methods-metadata.tsv"), show_col_types = FALSE)
methods_meta <- dplyr::filter(methods_meta, .data$Kind == "selector-family")
methods_meta_all <- readr::read_tsv(file.path(data_dir, "methods-metadata.tsv"), show_col_types = FALSE)
metric_ranges <- readr::read_tsv(ranges_path, show_col_types = FALSE)

metrics_scaled <- scale_spatial_metrics(metrics, metric_ranges)
metrics_summary <- summarise_spatial_metrics(metrics, metric_ranges)

figure1 <- patchwork::wrap_plots(
    plot_overview_heatmap(metrics_summary, methods_meta_all),
    plot_overview_ranking(metrics_summary, methods_meta_all),
    ncol = 2,
    widths = c(2, 1)
)
save_figure_files(figure1, file.path(output_dir, "figure1_overview"), width = 8.2, height = 9.2)

figure1b <- plot_setting_heatmap(metrics_summary, methods_meta_all)
save_figure_files(figure1b, file.path(output_dir, "figure1b_setting_heatmap"), width = 8.2, height = 9.2)

figure2 <- plot_task_panels(metrics_scaled, methods_meta_all)
save_figure_files(figure2, file.path(output_dir, "figure2_task_panels"), width = 8.2, height = 10.8)

figure3 <- plot_num_features_summary(metrics_summary, methods_meta, datasets_meta)
save_figure_files(figure3, file.path(output_dir, "figure3_num_features"), width = 8, height = 8)

figure4 <- plot_stability_summary(metrics_summary, datasets_meta, methods_meta_all)
save_figure_files(figure4, file.path(output_dir, "figure4_stability"), width = 8.2, height = 7.2)
