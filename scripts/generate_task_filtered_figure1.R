#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[[1]] else file.path("results", "atlas_style_final_all_slice", "data")
output_dir <- if (length(args) >= 2) args[[2]] else file.path("results", "atlas_style_task_filtered", "figures")
base_dir <- normalizePath(file.path(data_dir, ".."), mustWork = FALSE)
ranges_path <- file.path(base_dir, "output", "baseline-ranges.tsv")
audit_dir <- file.path(dirname(output_dir), "data")

source(file.path("analysis", "atlas_style", "R", "spatial_summarisation.R"))
source(file.path("analysis", "atlas_style", "R", "spatial_plotting.R"))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

metrics <- readr::read_tsv(file.path(data_dir, "benchmark.tsv"), show_col_types = FALSE)
methods_meta <- readr::read_tsv(file.path(data_dir, "methods-metadata.tsv"), show_col_types = FALSE)
metric_ranges <- readr::read_tsv(ranges_path, show_col_types = FALSE)

# Task-specific dataset filtering based on current applicability audit.
integration_datasets <- c("DLPFC", "MouseBrainSerialSections", "STOmicsVisium5Samples", "SagittalAtlas")
clustering_datasets <- sort(unique(metrics$Dataset))
alignment_datasets <- c("DLPFC", "MouseBrainSerialSections", "STOmicsVisium5Samples", "SagittalAtlas")
slice_repr_datasets <- c("DLPFC")

metrics_filtered <- metrics |>
    dplyr::filter(
        (.data$Type == "Integration" & .data$Dataset %in% integration_datasets) |
        (.data$Type == "Clustering" & .data$Dataset %in% clustering_datasets) |
        (.data$Type == "Alignment" & .data$Dataset %in% alignment_datasets) |
        (.data$Type == "SliceRepresentation" & .data$Dataset %in% slice_repr_datasets)
    )

coverage_map <- tibble::tribble(
    ~Type, ~Datasets,
    "Integration", paste(integration_datasets, collapse = ", "),
    "Clustering", paste(clustering_datasets, collapse = ", "),
    "Alignment", paste(alignment_datasets, collapse = ", "),
    "SliceRepresentation", paste(slice_repr_datasets, collapse = ", ")
)

readr::write_tsv(metrics_filtered, file.path(audit_dir, "benchmark_task_filtered.tsv"))
readr::write_tsv(coverage_map, file.path(audit_dir, "task_dataset_map.tsv"))

metrics_summary <- summarise_spatial_metrics(metrics_filtered, metric_ranges)

ranking_export <- metrics_summary |>
    dplyr::group_by(.data$Dataset, .data$IntegrationLabel) |>
    dplyr::mutate(
        RankOverall = rank(-.data$Overall),
        RankIntegration = rank(-.data$Integration),
        RankClustering = rank(-.data$Clustering),
        RankAlignment = rank(-.data$Alignment),
        RankSliceRepresentation = rank(-.data$SliceRepresentation)
    ) |>
    dplyr::ungroup() |>
    dplyr::group_by(.data$Method) |>
    dplyr::summarise(
        MeanRankOverall = mean(.data$RankOverall, na.rm = TRUE),
        MeanRankIntegration = mean(.data$RankIntegration, na.rm = TRUE),
        MeanRankClustering = mean(.data$RankClustering, na.rm = TRUE),
        MeanRankAlignment = mean(.data$RankAlignment, na.rm = TRUE),
        MeanRankSliceRepresentation = mean(.data$RankSliceRepresentation, na.rm = TRUE),
        .groups = "drop"
    ) |>
    dplyr::mutate(MethodLabel = methods_meta$Name[match(.data$Method, methods_meta$Method)]) |>
    dplyr::arrange(.data$MeanRankOverall)

readr::write_tsv(ranking_export, file.path(audit_dir, "figure1_task_filtered_ranks.tsv"))

figure1 <- patchwork::wrap_plots(
    plot_overview_heatmap(metrics_summary, methods_meta),
    plot_overview_ranking(metrics_summary, methods_meta),
    ncol = 2,
    widths = c(2, 1)
)

save_figure_files(figure1, file.path(output_dir, "figure1_overview_task_filtered"), width = 8.2, height = 9.2)
