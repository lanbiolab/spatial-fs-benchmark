#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[[1]] else file.path("results", "current_rank", "data")
output_dir <- if (length(args) >= 2) args[[2]] else file.path("results", "current_rank", "figures")
base_dir <- normalizePath(file.path(data_dir, ".."), mustWork = FALSE)
ranges_path <- file.path(base_dir, "output", "baseline-ranges.tsv")

source(file.path("external", "atlas-feature-selection-benchmark", "analysis", "R", "plotting.R"))
source(file.path("analysis", "atlas_style", "R", "spatial_summarisation.R"))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

metrics <- readr::read_tsv(file.path(data_dir, "benchmark.tsv"), show_col_types = FALSE) |>
    dplyr::filter(.data$Type %in% c("Integration", "Clustering", "Alignment"))
methods_meta <- readr::read_tsv(file.path(data_dir, "methods-metadata.tsv"), show_col_types = FALSE)
metric_ranges <- readr::read_tsv(ranges_path, show_col_types = FALSE) |>
    dplyr::filter(.data$Type %in% c("Integration", "Clustering", "Alignment"))

metrics_summary <- summarise_spatial_metrics(
    metrics,
    metric_ranges,
    type_weights = c("Integration" = 1/3, "Clustering" = 1/3, "Alignment" = 1/3),
    require_types_for_overall = c("Integration", "Clustering", "Alignment")
)

compute_method_scores_current <- function(metrics_summary, methods_meta) {
    metrics_summary |>
        dplyr::group_by(.data$Method) |>
        dplyr::summarise(
            GlobalScoreOverall = mean(.data$Overall, na.rm = TRUE),
            GlobalScoreIntegration = mean(.data$Integration, na.rm = TRUE),
            GlobalScoreClustering = mean(.data$Clustering, na.rm = TRUE),
            GlobalScoreAlignment = mean(.data$Alignment, na.rm = TRUE),
            .groups = "drop"
        ) |>
        dplyr::mutate(
            GlobalRankOverall = rank(-.data$GlobalScoreOverall),
            GlobalRankIntegration = rank(-.data$GlobalScoreIntegration),
            GlobalRankClustering = rank(-.data$GlobalScoreClustering),
            GlobalRankAlignment = rank(-.data$GlobalScoreAlignment),
            MethodLabel = methods_meta$Name[match(.data$Method, methods_meta$Method)]
        ) |>
        dplyr::arrange(.data$GlobalRankOverall, dplyr::desc(.data$GlobalScoreOverall))
}

compute_method_order_current <- function(metrics_summary, methods_meta) {
    compute_method_scores_current(metrics_summary, methods_meta)$MethodLabel
}

plot_overview_heatmap_current <- function(metrics_summary, methods_meta) {
    method_order <- compute_method_order_current(metrics_summary, methods_meta)

    overview <- metrics_summary |>
        dplyr::group_by(.data$Method, .data$IntegrationLabel) |>
        dplyr::summarise(
            Integration = mean(.data$Integration, na.rm = TRUE),
            Clustering = mean(.data$Clustering, na.rm = TRUE),
            Alignment = mean(.data$Alignment, na.rm = TRUE),
            Overall = mean(.data$Overall, na.rm = TRUE),
            .groups = "drop"
        ) |>
        tidyr::pivot_longer(
            cols = c("Integration", "Clustering", "Alignment", "Overall"),
            names_to = "Type",
            values_to = "Value"
        ) |>
        dplyr::mutate(
            Value = dplyr::if_else(is.nan(.data$Value), NA_real_, .data$Value),
            Column = ifelse(.data$Type == "Overall", "Overall", paste(.data$IntegrationLabel, .data$Type, sep = "\n")),
            MethodLabel = methods_meta$Name[match(.data$Method, methods_meta$Method)]
        )

    overview |>
        dplyr::mutate(
            MethodLabel = factor(.data$MethodLabel, levels = rev(method_order)),
            Column = factor(.data$Column, levels = unique(.data$Column))
        ) |>
        ggplot2::ggplot(ggplot2::aes(x = .data$Column, y = .data$MethodLabel, fill = .data$Value)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_gradient2(
            low = "#762a83",
            mid = "white",
            high = "#1b7837",
            midpoint = 0.5,
            na.value = "grey92"
        ) +
        ggplot2::labs(x = NULL, y = NULL, fill = "Mean scaled\nvalue") +
        theme_features_pub() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
            panel.grid = ggplot2::element_blank()
        )
}

plot_overview_ranking_current <- function(metrics_summary, methods_meta) {
    ranking <- compute_method_scores_current(metrics_summary, methods_meta) |>
        dplyr::select(
            .data$Method,
            .data$GlobalRankOverall,
            .data$GlobalRankIntegration,
            .data$GlobalRankClustering,
            .data$GlobalRankAlignment
        ) |>
        tidyr::pivot_longer(
            cols = c("GlobalRankOverall", "GlobalRankIntegration", "GlobalRankClustering", "GlobalRankAlignment"),
            names_to = "Key",
            values_to = "Rank"
        ) |>
        dplyr::mutate(
            Type = dplyr::recode(
                .data$Key,
                "GlobalRankOverall" = "Overall",
                "GlobalRankIntegration" = "Integration",
                "GlobalRankClustering" = "Clustering",
                "GlobalRankAlignment" = "Alignment"
            )
        ) |>
        dplyr::mutate(
            MethodLabel = methods_meta$Name[match(.data$Method, methods_meta$Method)],
            Type = factor(.data$Type, levels = c("Overall", "Integration", "Clustering", "Alignment"))
        )

    min_rank <- min(ranking$Rank, na.rm = TRUE)
    max_rank <- max(ranking$Rank, na.rm = TRUE)
    mid_rank <- round((min_rank + max_rank) / 2)

    ranking |>
        dplyr::mutate(MethodLabel = factor(.data$MethodLabel, levels = rev(compute_method_order_current(metrics_summary, methods_meta)))) |>
        ggplot2::ggplot(ggplot2::aes(x = .data$Type, y = .data$MethodLabel, colour = .data$Type)) +
        ggplot2::geom_point(
            ggplot2::aes(alpha = .data$Rank),
            shape = "square",
            size = 3.2,
            na.rm = TRUE
        ) +
        ggplot2::scale_colour_manual(
            values = c(
                "Overall" = "#f781bf",
                "Integration" = "#e41a1c",
                "Clustering" = "#377eb8",
                "Alignment" = "#4daf4a"
            )
        ) +
        ggplot2::scale_alpha_continuous(
            range = c(1, 0.2),
            breaks = c(min_rank, mid_rank, max_rank),
            labels = c(
                sprintf("%.0f (best)", min_rank),
                sprintf("%.0f", mid_rank),
                sprintf("%.0f (worst)", max_rank)
            ),
            name = "Global rank"
        ) +
        ggplot2::labs(x = NULL, y = NULL, colour = NULL) +
        theme_features_pub() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
            panel.grid = ggplot2::element_blank(),
            legend.position = "right"
        )
}

ranking_export <- compute_method_scores_current(metrics_summary, methods_meta) |>
    dplyr::select(
        .data$Method,
        .data$MethodLabel,
        .data$GlobalScoreOverall,
        .data$GlobalScoreIntegration,
        .data$GlobalScoreClustering,
        .data$GlobalScoreAlignment,
        .data$GlobalRankOverall,
        .data$GlobalRankIntegration,
        .data$GlobalRankClustering,
        .data$GlobalRankAlignment
    )

readr::write_tsv(ranking_export, file.path(output_dir, "current_rank_table.tsv"))

figure1 <- patchwork::wrap_plots(
    plot_overview_heatmap_current(metrics_summary, methods_meta),
    plot_overview_ranking_current(metrics_summary, methods_meta),
    ncol = 2,
    widths = c(2, 1)
)

save_figure_files(figure1, file.path(output_dir, "figure1_overview_current"), width = 8.2, height = 9.2)
