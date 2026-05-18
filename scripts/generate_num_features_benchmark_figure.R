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
output_dir <- if (length(args) >= 2) args[[2]] else file.path("results", "num_features_benchmark", "figures")
base_dir <- normalizePath(file.path(data_dir, ".."), mustWork = FALSE)
ranges_path <- file.path(base_dir, "output", "baseline-ranges.tsv")

source(file.path("external", "atlas-feature-selection-benchmark", "analysis", "R", "plotting.R"))
source(file.path("analysis", "atlas_style", "R", "spatial_summarisation.R"))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

types_palette <- c(
    "Integration" = "#e41a1c",
    "Clustering" = "#377eb8",
    "Alignment" = "#4daf4a",
    "Overall" = "#f781bf"
)

metrics <- readr::read_tsv(file.path(data_dir, "num-features.tsv"), show_col_types = FALSE) |>
    dplyr::filter(.data$Type %in% c("Integration", "Clustering", "Alignment")) |>
    dplyr::filter(as.character(.data$SelFeatures) != "all")

metric_ranges <- readr::read_tsv(ranges_path, show_col_types = FALSE) |>
    dplyr::filter(.data$Type %in% c("Integration", "Clustering", "Alignment"))

methods_meta <- readr::read_tsv(file.path(data_dir, "methods-metadata.tsv"), show_col_types = FALSE) |>
    dplyr::filter(.data$Kind == "selector-family")

datasets_meta <- readr::read_tsv(file.path(data_dir, "datasets-metadata.tsv"), show_col_types = FALSE)

method_names <- methods_meta$Name
names(method_names) <- methods_meta$Method

dataset_names <- datasets_meta$Name
names(dataset_names) <- datasets_meta$Dataset

metrics_summary <- summarise_spatial_metrics(
    metrics,
    metric_ranges,
    type_weights = c("Integration" = 1/3, "Clustering" = 1/3, "Alignment" = 1/3),
    require_types_for_overall = c("Integration", "Clustering", "Alignment")
)

sel_levels <- metrics_summary |>
    dplyr::distinct(.data$SelFeatures) |>
    dplyr::mutate(
        SelFeaturesText = as.character(.data$SelFeatures),
        SelFeaturesOrder = dplyr::if_else(
            .data$SelFeaturesText == "all",
            Inf,
            suppressWarnings(as.numeric(.data$SelFeaturesText))
        )
    ) |>
    dplyr::arrange(.data$SelFeaturesOrder) |>
    dplyr::pull(.data$SelFeaturesText)

metrics_summary_plotting <- metrics_summary |>
    dplyr::transmute(
        Dataset = dplyr::recode(.data$Dataset, !!!dataset_names, .default = .data$Dataset),
        Method = dplyr::recode(.data$MethodBase, !!!method_names, .default = .data$MethodBase),
        IntegrationMethod = .data$IntegrationLabel,
        SelFeatures = factor(as.character(.data$SelFeatures), levels = sel_levels),
        IntegrationScore = .data$Integration,
        ClusteringScore = .data$Clustering,
        AlignmentScore = .data$Alignment,
        Overall = .data$Overall
    ) |>
    tidyr::pivot_longer(
        cols = c("IntegrationScore", "ClusteringScore", "AlignmentScore", "Overall"),
        names_to = "Type",
        values_to = "Value"
    ) |>
    dplyr::mutate(
        Type = dplyr::recode(
            .data$Type,
            "IntegrationScore" = "Integration",
            "ClusteringScore" = "Clustering",
            "AlignmentScore" = "Alignment",
            "Overall" = "Overall"
        ),
        Type = factor(.data$Type, levels = c("Integration", "Clustering", "Alignment", "Overall"))
    ) |>
    dplyr::filter(!is.na(.data$Value))

metrics_summary_plotting <- metrics_summary_plotting |>
    dplyr::group_by(.data$Dataset, .data$Method, .data$IntegrationMethod, .data$Type) |>
    dplyr::mutate(StandardValue = as.vector(scale(.data$Value))) |>
    dplyr::ungroup()

overall_means <- metrics_summary_plotting |>
    dplyr::group_by(.data$Type, .data$SelFeatures) |>
    dplyr::summarise(
        Value = mean(.data$Value, na.rm = TRUE),
        StandardValue = mean(.data$StandardValue, na.rm = TRUE),
        .groups = "drop"
    )

dataset_means <- metrics_summary_plotting |>
    dplyr::group_by(.data$Dataset, .data$Type, .data$SelFeatures) |>
    dplyr::summarise(
        Value = mean(.data$Value, na.rm = TRUE),
        SD = sd(.data$StandardValue, na.rm = TRUE),
        StandardValue = mean(.data$StandardValue, na.rm = TRUE),
        .groups = "drop"
    ) |>
    dplyr::mutate(SD = dplyr::if_else(is.na(.data$SD), 0, .data$SD))

method_means <- metrics_summary_plotting |>
    dplyr::group_by(.data$Method, .data$Type, .data$SelFeatures) |>
    dplyr::summarise(
        Value = mean(.data$Value, na.rm = TRUE),
        SD = sd(.data$StandardValue, na.rm = TRUE),
        StandardValue = mean(.data$StandardValue, na.rm = TRUE),
        .groups = "drop"
    ) |>
    dplyr::mutate(SD = dplyr::if_else(is.na(.data$SD), 0, .data$SD))

method_order <- method_means |>
    dplyr::filter(.data$Type == "Overall") |>
    dplyr::group_by(.data$Method) |>
    dplyr::summarise(Score = mean(.data$Value, na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(.data$Score)) |>
    dplyr::pull(.data$Method)

dataset_order <- datasets_meta |>
    dplyr::mutate(Name = dplyr::recode(.data$Dataset, !!!dataset_names, .default = .data$Dataset)) |>
    dplyr::pull(.data$Name)

mean_limits <- c(
    min(c(method_means$StandardValue, dataset_means$StandardValue), na.rm = TRUE),
    max(c(method_means$StandardValue, dataset_means$StandardValue), na.rm = TRUE)
)

sd_limits <- c(
    max(c(dataset_means$SD, method_means$SD), na.rm = TRUE),
    0
)

type_levels <- c("Integration", "Clustering", "Alignment", "Overall")

dataset_means_complete <- tidyr::expand_grid(
    Dataset = dataset_order,
    Type = factor(type_levels, levels = type_levels),
    SelFeatures = factor(sel_levels, levels = sel_levels)
) |>
    dplyr::left_join(dataset_means, by = c("Dataset", "Type", "SelFeatures")) |>
    dplyr::mutate(Missing = is.na(.data$StandardValue))

method_means_complete <- tidyr::expand_grid(
    Method = method_order,
    Type = factor(type_levels, levels = type_levels),
    SelFeatures = factor(sel_levels, levels = sel_levels)
) |>
    dplyr::left_join(method_means, by = c("Method", "Type", "SelFeatures")) |>
    dplyr::mutate(Missing = is.na(.data$StandardValue))

theme_features_heatmap <- theme_features_pub() +
    theme(
        legend.position = "bottom",
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank()
    )

overall_lineplot <- ggplot(
    metrics_summary_plotting,
    aes(x = .data$SelFeatures, y = .data$StandardValue, colour = .data$Type, fill = .data$Type)
) +
    geom_hline(yintercept = 0, colour = "red", linewidth = 0.3) +
    geom_point(
        shape = 21,
        stroke = 0,
        size = 0.55,
        alpha = 0.18,
        position = position_jitter(width = 0.10, height = 0)
    ) +
    geom_line(
        data = overall_means,
        aes(group = .data$Type),
        linewidth = 0.9
    ) +
    geom_point(
        data = overall_means,
        shape = 23,
        size = 1.5,
        stroke = 0.55,
        fill = "white"
    ) +
    scale_color_manual(values = types_palette, guide = "none") +
    scale_fill_manual(values = types_palette, guide = "none") +
    facet_grid(. ~ .data$Type) +
    labs(y = "Standardised value") +
    theme_features_pub() +
    theme(
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(margin = margin(r = 0)),
        axis.title.y.left = element_text(vjust = 0.15, margin = margin(r = -58)),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 2)
    )

datasets_heatmap <- ggplot(
    dataset_means_complete |>
        dplyr::mutate(Dataset = factor(.data$Dataset, levels = dataset_order)),
    aes(x = .data$SelFeatures, y = .data$Dataset)
) +
    geom_point(
        data = \(x) dplyr::filter(x, .data$Missing),
        shape = 15,
        size = 1.45,
        colour = "#d9d9d9"
    ) +
    geom_point(
        data = \(x) dplyr::filter(x, !.data$Missing),
        aes(colour = .data$StandardValue, size = .data$SD),
        shape = 15
    ) +
    scale_colour_gradient2(
        low = "#762a83",
        mid = "white",
        high = "#1b7837",
        midpoint = 0,
        limits = mean_limits,
        na.value = "#d9d9d9"
    ) +
    scale_size_continuous(
        trans = "reverse",
        limits = sd_limits,
        range = c(0.15, 2.4)
    ) +
    facet_grid(. ~ .data$Type) +
    labs(
        title = "Datasets",
        x = "Number of selected features",
        colour = "Mean standardised value",
        size = "Standard deviation of\nstandardised values"
    ) +
    guides(
        colour = guide_colourbar(order = 1),
        size = guide_legend(theme = theme(legend.text.position = "bottom"), order = 2)
    ) +
    theme_features_heatmap +
    theme(
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 5.2, face = "bold", colour = "black")
    )

methods_heatmap <- ggplot(
    method_means_complete |>
        dplyr::mutate(Method = factor(.data$Method, levels = method_order)),
    aes(x = .data$SelFeatures, y = .data$Method)
) +
    geom_point(
        data = \(x) dplyr::filter(x, .data$Missing),
        shape = 15,
        size = 1.45,
        colour = "#d9d9d9"
    ) +
    geom_point(
        data = \(x) dplyr::filter(x, !.data$Missing),
        aes(colour = .data$StandardValue, size = .data$SD),
        shape = 15
    ) +
    geom_point(
        data = data.frame(SelFeatures = factor(sel_levels[[1]], levels = sel_levels), Method = factor(method_order[[1]], levels = method_order), MissingLegend = "NA"),
        aes(x = .data$SelFeatures, y = .data$Method, shape = .data$MissingLegend),
        inherit.aes = FALSE,
        alpha = 0,
        colour = "#d9d9d9",
        size = 3,
        show.legend = TRUE
    ) +
    scale_colour_gradient2(
        low = "#762a83",
        mid = "white",
        high = "#1b7837",
        midpoint = 0,
        limits = mean_limits,
        na.value = "#d9d9d9"
    ) +
    scale_shape_manual(
        values = c("NA" = 15),
        name = NULL,
        guide = guide_legend(
            order = 3,
            override.aes = list(alpha = 1, colour = "#d9d9d9", size = 4)
        )
    ) +
    scale_size_continuous(
        trans = "reverse",
        limits = sd_limits,
        range = c(0.15, 2.4)
    ) +
    facet_grid(. ~ .data$Type) +
    labs(
        title = "Methods",
        x = "Number of selected features",
        colour = "Mean standardised value",
        size = NULL
    ) +
    guides(
        colour = guide_colourbar(order = 1),
        size = "none"
    ) +
    theme_features_heatmap +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 4.2, face = "bold", colour = "black")
    )

summary_plot_main <- wrap_plots(
    overall_lineplot,
    datasets_heatmap,
    methods_heatmap,
    ncol = 1,
    heights = c(1, 1, 1.7),
    guides = "collect"
) &
    theme(
        legend.position = "bottom",
        legend.title.position = "top",
        legend.box = "horizontal"
    )
summary_plot <- summary_plot_main

readr::write_tsv(overall_means, file.path(output_dir, "overall_feature_number_means.tsv"))
readr::write_tsv(dataset_means, file.path(output_dir, "dataset_feature_number_summary.tsv"))
readr::write_tsv(method_means, file.path(output_dir, "method_feature_number_summary.tsv"))

save_figure_files(summary_plot, file.path(output_dir, "figure_num_features_benchmark"), width = 7.8, height = 6.8)
