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
output_dir <- if (length(args) >= 2) args[[2]] else file.path("results", "integration_benchmark", "figures")
base_dir <- normalizePath(file.path(data_dir, ".."), mustWork = FALSE)
ranges_path <- file.path(base_dir, "output", "baseline-ranges.tsv")

source(file.path("external", "atlas-feature-selection-benchmark", "analysis", "R", "plotting.R"))
source(file.path("analysis", "atlas_style", "R", "spatial_summarisation.R"))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

types_palette <- c(
    "Overall" = "#f781bf",
    "Integration" = "#e41a1c",
    "Clustering" = "#377eb8",
    "Alignment" = "#4daf4a"
)

theme_features_integration <- theme_features_pub() +
    theme(
        legend.position = "bottom",
        legend.title.position = "top",
        panel.grid = element_blank(),
        panel.spacing.x = unit(0.02, "cm"),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.ticks.y = element_blank(),
        legend.key.width = unit(0.35, "cm"),
        legend.key.spacing.x = unit(0.03, "cm"),
        legend.box = "vertical"
    )

metrics <- readr::read_tsv(file.path(data_dir, "benchmark.tsv"), show_col_types = FALSE) |>
    dplyr::filter(.data$Type %in% c("Integration", "Clustering", "Alignment"))

metric_ranges <- readr::read_tsv(ranges_path, show_col_types = FALSE) |>
    dplyr::filter(.data$Type %in% c("Integration", "Clustering", "Alignment"))

integrations_meta <- readr::read_tsv(file.path(data_dir, "integrations-metadata.tsv"), show_col_types = FALSE)
integration_names <- integrations_meta$Name
names(integration_names) <- integrations_meta$Integration

metrics_summary <- summarise_spatial_metrics(
    metrics,
    metric_ranges,
    type_weights = c("Integration" = 1/3, "Clustering" = 1/3, "Alignment" = 1/3),
    require_types_for_overall = c("Integration", "Clustering", "Alignment")
)

task_long <- metrics_summary |>
    dplyr::select(
        .data$Dataset,
        .data$Method,
        IntegrationMethod = .data$IntegrationMethod,
        .data$IntegrationLabel,
        .data$Overall,
        .data$Integration,
        .data$Clustering,
        .data$Alignment
    ) |>
    dplyr::rename(
        ScoreIntegration = .data$Integration,
        ScoreClustering = .data$Clustering,
        ScoreAlignment = .data$Alignment
    ) |>
    tidyr::pivot_longer(
        cols = c("Overall", "ScoreIntegration", "ScoreClustering", "ScoreAlignment"),
        names_to = "Type",
        values_to = "Value"
    ) |>
    dplyr::mutate(
        Type = dplyr::recode(
            .data$Type,
            "Overall" = "Overall",
            "ScoreIntegration" = "Integration",
            "ScoreClustering" = "Clustering",
            "ScoreAlignment" = "Alignment"
        )
    )

integration_order <- task_long |>
    dplyr::group_by(.data$IntegrationLabel) |>
    dplyr::summarise(GlobalScore = mean(.data$Value, na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(.data$GlobalScore)) |>
    dplyr::pull(.data$IntegrationLabel)

integration_scores <- task_long |>
    dplyr::group_by(.data$IntegrationLabel, .data$Type) |>
    dplyr::summarise(
        Mean = mean(.data$Value, na.rm = TRUE),
        SD = sd(.data$Value, na.rm = TRUE),
        .groups = "drop"
    ) |>
    dplyr::mutate(
        IntegrationLabel = factor(.data$IntegrationLabel, levels = rev(integration_order)),
        Type = factor(.data$Type, levels = c("Overall", "Integration", "Clustering", "Alignment"))
    )

scvi_reference <- task_long |>
    dplyr::filter(.data$IntegrationLabel == "scVI") |>
    dplyr::select(.data$Dataset, .data$Method, .data$Type, scVIValue = .data$Value)

integration_score_diffs <- task_long |>
    dplyr::left_join(scvi_reference, by = c("Dataset", "Method", "Type")) |>
    dplyr::mutate(Difference = .data$Value - .data$scVIValue) |>
    dplyr::group_by(.data$IntegrationLabel, .data$Type) |>
    dplyr::summarise(
        Mean = mean(.data$Difference, na.rm = TRUE),
        SD = sd(.data$Difference, na.rm = TRUE),
        .groups = "drop"
    ) |>
    dplyr::mutate(
        IntegrationLabel = factor(.data$IntegrationLabel, levels = rev(integration_order)),
        Type = factor(.data$Type, levels = c("Overall", "Integration", "Clustering", "Alignment"))
    )

rank_long <- task_long |>
    dplyr::group_by(.data$Dataset, .data$Method, .data$Type) |>
    dplyr::mutate(Rank = rank(-.data$Value, ties.method = "average")) |>
    dplyr::ungroup()

integration_ranks <- rank_long |>
    dplyr::group_by(.data$IntegrationLabel, .data$Type) |>
    dplyr::summarise(
        MeanRank = mean(.data$Rank, na.rm = TRUE),
        SDRank = sd(.data$Rank, na.rm = TRUE),
        .groups = "drop"
    ) |>
    dplyr::mutate(
        IntegrationLabel = factor(.data$IntegrationLabel, levels = rev(integration_order)),
        Type = factor(.data$Type, levels = c("Overall", "Integration", "Clustering", "Alignment"))
    )

scvi_rank_reference <- rank_long |>
    dplyr::filter(.data$IntegrationLabel == "scVI") |>
    dplyr::select(.data$Dataset, .data$Method, .data$Type, scVIRank = .data$Rank)

integration_rank_diffs <- rank_long |>
    dplyr::left_join(scvi_rank_reference, by = c("Dataset", "Method", "Type")) |>
    dplyr::mutate(Difference = .data$Rank - .data$scVIRank) |>
    dplyr::group_by(.data$IntegrationLabel, .data$Type) |>
    dplyr::summarise(
        Mean = mean(.data$Difference, na.rm = TRUE),
        SD = sd(.data$Difference, na.rm = TRUE),
        .groups = "drop"
    ) |>
    dplyr::mutate(
        IntegrationLabel = factor(.data$IntegrationLabel, levels = rev(integration_order)),
        Type = factor(.data$Type, levels = c("Overall", "Integration", "Clustering", "Alignment"))
    )

integration_scores_plot <- ggplot(
    integration_scores,
    aes(x = .data$Type, y = .data$IntegrationLabel, colour = .data$Mean, size = .data$SD)
) +
    geom_point(shape = "square") +
    scale_colour_viridis_c(
        option = "magma",
        limits = c(0, max(integration_scores$Mean, na.rm = TRUE))
    ) +
    scale_size_continuous(
        trans = "reverse",
        range = c(0.2, 2.4),
        breaks = scales::breaks_extended(n = 5)
    ) +
    labs(
        title = "Mean scores",
        colour = "Mean score",
        size = "Standard deviation"
    ) +
    guides(
        colour = guide_colourbar(order = 1),
        size = guide_legend(theme = theme(legend.text.position = "bottom"), order = 2)
    ) +
    theme_features_integration

integration_scores_diffs_plot <- ggplot(
    integration_score_diffs,
    aes(x = .data$Type, y = .data$IntegrationLabel, colour = .data$Mean, size = .data$SD)
) +
    geom_point(shape = "square") +
    colorspace::scale_colour_continuous_diverging(
        palette = "PurpleGreen"
    ) +
    scale_size_continuous(
        trans = "reverse",
        range = c(0.2, 2.4),
        breaks = scales::breaks_extended(n = 5)
    ) +
    labs(
        title = "Difference in\nmean scores",
        colour = "Difference to scVI",
        size = "Standard deviation"
    ) +
    guides(
        colour = guide_colourbar(order = 1),
        size = guide_legend(theme = theme(legend.text.position = "bottom"), order = 2)
    ) +
    theme_features_integration +
    theme(axis.text.y = element_blank())

integration_ranks_plot <- ggplot(
    integration_ranks,
    aes(x = .data$Type, y = .data$IntegrationLabel, colour = .data$Type)
) +
    geom_point(
        aes(alpha = .data$MeanRank, size = .data$SDRank),
        shape = "square"
    ) +
    scale_colour_manual(values = types_palette, guide = "none") +
    scale_alpha_continuous(
        range = c(1, 0.2),
        breaks = c(1, 2, 3, 4),
        labels = c("1 (best)", "2", "3", "4 (worst)"),
        name = "Mean rank"
    ) +
    scale_size_continuous(
        trans = "reverse",
        limits = c(max(integration_ranks$SDRank, na.rm = TRUE), 0),
        range = c(0.2, 2.4),
        breaks = scales::breaks_extended(n = 5)
    ) +
    labs(
        title = "Mean ranks",
        size = "Standard deviation"
    ) +
    guides(
        alpha = guide_legend(theme = theme(legend.text.position = "bottom"), nrow = 1, order = 1),
        size = guide_legend(theme = theme(legend.text.position = "bottom"), order = 2)
    ) +
    theme_features_integration +
    theme(axis.text.y = element_blank())

integration_diffs_ranks_plot <- ggplot(
    integration_rank_diffs,
    aes(x = .data$Type, y = .data$IntegrationLabel, colour = .data$Mean, size = .data$SD)
) +
    geom_point(shape = "square") +
    colorspace::scale_colour_continuous_diverging(
        palette = "Tropic",
        rev = TRUE
    ) +
    scale_size_continuous(
        trans = "reverse",
        range = c(0.2, 2.4),
        breaks = scales::breaks_extended(n = 5)
    ) +
    labs(
        title = "Difference in\nmean ranks",
        colour = "Difference to scVI",
        size = "Standard deviation"
    ) +
    guides(
        colour = guide_colourbar(order = 1),
        size = guide_legend(theme = theme(legend.text.position = "bottom"), nrow = 1, order = 2)
    ) +
    theme_features_integration +
    theme(axis.text.y = element_blank())

integration_figure <- wrap_plots(
    integration_scores_plot,
    integration_scores_diffs_plot,
    integration_ranks_plot,
    integration_diffs_ranks_plot,
    nrow = 1,
    widths = c(3, 2, 3, 2)
)

readr::write_tsv(integration_scores, file.path(output_dir, "integration_scores_summary.tsv"))
readr::write_tsv(integration_score_diffs, file.path(output_dir, "integration_score_diffs_summary.tsv"))
readr::write_tsv(integration_ranks, file.path(output_dir, "integration_ranks_summary.tsv"))
readr::write_tsv(integration_rank_diffs, file.path(output_dir, "integration_rank_diffs_summary.tsv"))

save_figure_files(integration_figure, file.path(output_dir, "figure_integration_benchmark"), width = 8.2, height = 4.8)
