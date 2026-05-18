#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(patchwork)
    library(purrr)
    library(colorspace)
})

source(file.path("external", "atlas-feature-selection-benchmark", "analysis", "R", "plotting.R"))

args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[[1]] else file.path("results", "current_rank", "data")
output_dir <- if (length(args) >= 2) args[[2]] else file.path("results", "extended_data_fig1", "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

metric_levels <- c(
    "bASW", "dASW", "iLISI", "dLISI", "ILL", "GC",
    "ARI", "NMI", "CHAOS", "PAS", "Silhouette",
    "Accuracy", "Ratio"
)

metric_labels <- c(
    "ari" = "ARI",
    "nmi" = "NMI",
    "silhouette" = "Silhouette"
)

type_levels <- c("Integration", "Clustering", "Alignment")
feature_key <- c(
    feat_genes = "Genes",
    feat_spots = "Spots",
    feat_slices = "Slices",
    feat_labels = "Labels",
    feat_spots_per_slice = "Spots/slice",
    feat_platform = "Platform"
)
feature_levels <- unname(feature_key)

safe_cor <- function(x, y) {
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) <= 2) return(NA_real_)
    suppressWarnings(stats::cor(x[ok], y[ok], use = "pairwise.complete.obs", method = "pearson"))
}

safe_sd <- function(x) {
    if (sum(is.finite(x)) <= 1) return(NA_real_)
    stats::sd(x, na.rm = TRUE)
}

metrics <- read_tsv(file.path(data_dir, "benchmark.tsv"), show_col_types = FALSE) |>
    filter(.data$Type %in% type_levels) |>
    mutate(
        MetricDisplay = recode(.data$MetricName, !!!metric_labels, .default = .data$MetricName),
        MetricDisplay = factor(.data$MetricDisplay, levels = metric_levels),
        Type = factor(.data$Type, levels = type_levels)
    )

datasets_meta <- read_tsv(file.path(data_dir, "datasets-metadata.tsv"), show_col_types = FALSE) |>
    transmute(
        Dataset,
        feat_genes = log10(.data$Features),
        feat_spots = log10(.data$Spots),
        feat_slices = log10(.data$NSlices),
        feat_labels = if_else(.data$Labels > 0, log10(.data$Labels), NA_real_),
        feat_spots_per_slice = if_else(.data$SpotsPerSlice > 0, log10(.data$SpotsPerSlice), NA_real_),
        feat_platform = if_else(.data$Platform == "Slide-seq", 1, 0)
    )

metrics_aug <- metrics |>
    left_join(datasets_meta, by = "Dataset")

feature_corr_long <- bind_rows(
    lapply(names(feature_key), function(feat) {
        metrics_aug |>
            mutate(FeatureValue = .data[[feat]]) |>
            group_by(.data$Method, .data$IntegrationLabel, .data$Type, .data$MetricDisplay) |>
            summarise(
                Corr = safe_cor(.data$Value, .data$FeatureValue),
                .groups = "drop"
            ) |>
            mutate(Feature = feature_key[[feat]])
    })
) |>
    filter(is.finite(.data$Corr))

feature_corr_summary <- feature_corr_long |>
    group_by(.data$MetricDisplay, .data$Type, .data$Feature) |>
    summarise(
        MeanCorr = mean(.data$Corr, na.rm = TRUE),
        SDCorr = safe_sd(.data$Corr),
        .groups = "drop"
    ) |>
    mutate(
        Feature = factor(.data$Feature, levels = feature_levels),
        MetricDisplay = factor(.data$MetricDisplay, levels = rev(metric_levels)),
        Type = factor(.data$Type, levels = type_levels)
    )

metric_types <- metrics_aug |>
    select(.data$MetricDisplay, .data$Type) |>
    distinct()

metric_corr_summary <- metrics_aug |>
    select(.data$Dataset, .data$IntegrationLabel, .data$Method, .data$MetricDisplay, .data$Value) |>
    pivot_wider(names_from = .data$MetricDisplay, values_from = .data$Value) |>
    group_by(.data$Dataset, .data$IntegrationLabel) |>
    group_split() |>
    map_dfr(function(df) {
        num_df <- df |>
            ungroup() |>
            select(-.data$Dataset, -.data$IntegrationLabel, -.data$Method)
        cmat <- suppressWarnings(stats::cor(num_df, use = "pairwise.complete.obs", method = "pearson"))
        out <- as.data.frame(as.table(cmat), stringsAsFactors = FALSE)
        names(out) <- c("RowMetric", "ColMetric", "Corr")
        out
    }) |>
    filter(is.finite(.data$Corr)) |>
    group_by(.data$RowMetric, .data$ColMetric) |>
    summarise(
        Mean = mean(.data$Corr, na.rm = TRUE),
        SD = safe_sd(.data$Corr),
        .groups = "drop"
    ) |>
    left_join(metric_types, by = c("RowMetric" = "MetricDisplay")) |>
    rename(RowType = .data$Type) |>
    left_join(metric_types, by = c("ColMetric" = "MetricDisplay")) |>
    rename(ColType = .data$Type) |>
    mutate(
        RowMetric = factor(.data$RowMetric, levels = rev(metric_levels)),
        ColMetric = factor(.data$ColMetric, levels = metric_levels),
        RowType = factor(.data$RowType, levels = type_levels),
        ColType = factor(.data$ColType, levels = type_levels)
    )

theme_ext <- theme_features_pub() +
    theme(
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 4.6, face = "bold", colour = "black"),
        axis.text.y = element_text(size = 4.8, face = "bold", colour = "black"),
        panel.grid = element_blank(),
        strip.text = element_text(size = 4.8, face = "bold", colour = "white"),
        strip.background = element_rect(fill = "black", colour = "black"),
        panel.spacing = unit(0.025, "cm"),
        legend.position = "bottom",
        legend.title.position = "top",
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.35, "cm"),
        legend.text = element_text(size = 5.5),
        legend.title = element_text(size = 6.5, face = "bold"),
        plot.margin = margin(0, 0, 0, 0)
    )

plot_tech_mean <- ggplot(feature_corr_summary) +
    geom_point(
        aes(x = .data$Feature, y = .data$MetricDisplay, colour = .data$MeanCorr),
        shape = "square", size = 3.2
    ) +
    colorspace::scale_colour_continuous_diverging(
        palette = "Tropic",
        rev = TRUE,
        limits = c(-1, 1)
    ) +
    facet_grid(.data$Type ~ ., scales = "free_y", space = "free_y") +
    labs(colour = "Mean") +
    theme_ext +
    theme(
        strip.text.y = element_blank()
    )

plot_tech_sd <- ggplot(feature_corr_summary) +
    geom_point(
        aes(x = .data$Feature, y = .data$MetricDisplay, colour = .data$SDCorr),
        shape = "square", size = 3.2
    ) +
    scale_colour_viridis_c(option = "cividis", limits = c(0, max(feature_corr_summary$SDCorr, na.rm = TRUE))) +
    facet_grid(.data$Type ~ ., scales = "free_y", space = "free_y") +
    labs(colour = "SD") +
    theme_ext +
    theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.y = element_blank()
    )

plot_metric_mean <- ggplot(metric_corr_summary) +
    geom_point(
        aes(x = .data$ColMetric, y = .data$RowMetric, colour = .data$Mean),
        shape = "square", size = 3.2
    ) +
    colorspace::scale_colour_continuous_diverging(
        palette = "Purple-Green",
        limits = c(-1, 1)
    ) +
    facet_grid(.data$RowType ~ .data$ColType, scales = "free", space = "free") +
    labs(colour = "Mean") +
    theme_ext +
    theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
    )

plot_metric_sd <- ggplot(metric_corr_summary) +
    geom_point(
        aes(x = .data$ColMetric, y = .data$RowMetric, colour = .data$SD),
        shape = "square", size = 3.2
    ) +
    scale_colour_viridis_c(option = "cividis", limits = c(0, max(metric_corr_summary$SD, na.rm = TRUE))) +
    facet_grid(.data$RowType ~ .data$ColType, scales = "free", space = "free") +
    labs(colour = "SD") +
    theme_ext +
    theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
    )

panel_a <- wrap_plots(
    plot_tech_mean,
    plot_tech_sd,
    nrow = 1,
    widths = c(1, 0.95),
    guides = "collect"
) &
    theme(
        legend.position = "bottom",
        plot.margin = margin(0, 0, 0, 0)
    )

panel_b <- wrap_plots(
    plot_metric_mean,
    plot_metric_sd,
    nrow = 1,
    widths = c(1, 0.95),
    guides = "collect"
) &
    theme(
        legend.position = "bottom",
        plot.margin = margin(0, 0, 0, 0)
    )

fig <- wrap_plots(
    panel_a,
    panel_b,
    nrow = 1,
    widths = c(1.0, 1.9)
) &
    theme(
        plot.margin = margin(0, 0, 0, 0)
    )

write_tsv(feature_corr_summary, file.path(output_dir, "extended_fig1_feature_corr_summary.tsv"))
write_tsv(metric_corr_summary, file.path(output_dir, "extended_fig1_metric_corr_summary.tsv"))

save_figure_files(fig, file.path(output_dir, "extended_data_fig1"), width = 8.3, height = 5.15)
