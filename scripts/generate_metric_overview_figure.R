#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(patchwork)
    library(ggridges)
})

args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[[1]] else file.path("results", "current_rank", "data")
output_dir <- if (length(args) >= 2) args[[2]] else file.path("results", "metric_overview", "figures")

source(file.path("external", "atlas-feature-selection-benchmark", "analysis", "R", "plotting.R"))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

frozen_settings_path <- file.path(data_dir, "setting_metric_seed_summary.tsv")
use_frozen_scores <- file.exists(frozen_settings_path)

if (use_frozen_scores) {
    frozen_ranges <- readr::read_tsv(
        file.path(data_dir, "frozen_metric_ranges.tsv"),
        show_col_types = FALSE
    )

    metrics <- readr::read_tsv(frozen_settings_path, show_col_types = FALSE) |>
        dplyr::left_join(
            frozen_ranges |>
                dplyr::select(
                    .data$dataset, .data$task, .data$metric_name,
                    .data$Lower, .data$Upper, .data$FrozenScaleDenominator,
                    .data$RangeStatus
                ),
            by = c("dataset", "task", "metric_name")
        ) |>
        dplyr::filter(.data$RangeStatus == "ok") |>
        dplyr::transmute(
            Dataset = .data$dataset,
            MethodBase = .data$fs_method,
            Method = paste0(
                .data$fs_method, "-N",
                dplyr::if_else(.data$n_features == "all", "all", .data$n_features)
            ),
            IntegrationLabel = dplyr::recode(
                .data$integration_method,
                "scvi" = "scVI",
                "cellcharter" = "CellCharter",
                .default = .data$integration_method
            ),
            SelFeatures = .data$n_features,
            Type = dplyr::recode(
                .data$task,
                "integration_eval" = "Integration",
                "clustering_eval" = "Clustering",
                "alignment_eval" = "Alignment"
            ),
            MetricName = .data$metric_name,
            Value = (.data$OrientedMean - .data$Lower) / .data$FrozenScaleDenominator
        )
} else {
    metrics <- readr::read_tsv(file.path(data_dir, "benchmark.tsv"), show_col_types = FALSE) |>
        dplyr::filter(.data$Type %in% c("Integration", "Clustering", "Alignment"))
}

datasets_meta <- readr::read_tsv(
    file.path("results", "current_rank", "data", "datasets-metadata.tsv"),
    show_col_types = FALSE
) |>
    dplyr::filter(.data$Dataset %in% unique(metrics$Dataset)) |>
    dplyr::mutate(
        LogSpots = log10(.data$Spots),
        LogFeatures = log10(.data$Features),
        LogLabels = dplyr::if_else(.data$Labels > 0, log10(.data$Labels), NA_real_),
        LogSpotsPerSlice = dplyr::if_else(.data$SpotsPerSlice > 0, log10(.data$SpotsPerSlice), NA_real_),
        IsSlideSeq = ifelse(.data$Platform == "Slide-seq", 1, 0)
    ) |>
    dplyr::select(
        .data$Dataset,
        .data$Platform,
        .data$NSlices,
        .data$Spots,
        .data$Features,
        .data$Labels,
        .data$SpotsPerSlice,
        .data$LogSpots,
        .data$LogFeatures,
        .data$LogLabels,
        .data$LogSpotsPerSlice,
        .data$IsSlideSeq
    )

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
type_colors <- c(
    "Integration" = "#6a3d9a",
    "Clustering" = "#e84a8a",
    "Alignment" = "#f0b429"
)

metrics_aug <- metrics |>
    dplyr::mutate(
        MetricDisplay = dplyr::recode(.data$MetricName, !!!metric_labels, .default = .data$MetricName),
        MetricDisplay = factor(.data$MetricDisplay, levels = rev(metric_levels)),
        Type = factor(.data$Type, levels = type_levels),
        SelFeaturesNum = suppressWarnings(as.numeric(.data$SelFeatures))
    ) |>
    dplyr::left_join(datasets_meta, by = "Dataset") |>
    dplyr::mutate(
        SelFeaturesNum = dplyr::if_else(.data$MethodBase == "all_features", .data$Features, .data$SelFeaturesNum),
        SelFeaturesNum = dplyr::if_else(.data$SelFeaturesNum <= .data$Features, .data$SelFeaturesNum, NA_real_)
    )

has_random_method <- any(metrics_aug$MethodBase == "random", na.rm = TRUE)

baseline_method_bases <- c(
    "all_features",
    "random",
    "scanpy_cell_ranger",
    "scanpy_cell_ranger_batch",
    "scsegindex"
)

metrics_random <- if (has_random_method) {
    metrics_aug |>
        dplyr::filter(.data$MethodBase == "all_features" | .data$MethodBase == "random")
} else {
    metrics_aug
}

metrics_nonrandom <- if (has_random_method) {
    metrics_aug |>
        dplyr::filter(.data$MethodBase == "all_features" | .data$MethodBase != "random")
} else {
    metrics_aug
}

if (use_frozen_scores) {
    metrics_aug <- metrics_aug |>
        dplyr::mutate(ScaledValue = .data$Value)
} else {
    baseline_ranges <- metrics_aug |>
        dplyr::filter(.data$MethodBase %in% baseline_method_bases) |>
        dplyr::filter(!is.na(.data$Value)) |>
        dplyr::group_by(.data$Dataset, .data$MetricDisplay, .data$Type) |>
        dplyr::summarise(
            Lower = min(.data$Value),
            Upper = max(.data$Value),
            .groups = "drop"
        )

    metrics_aug <- metrics_aug |>
        dplyr::left_join(baseline_ranges, by = c("Dataset", "MetricDisplay", "Type")) |>
        dplyr::mutate(
            ScaledValue = dplyr::if_else(
                is.finite(.data$Lower) & is.finite(.data$Upper) & (.data$Upper - .data$Lower) > 0,
                (.data$Value - .data$Lower) / (.data$Upper - .data$Lower),
                NA_real_
            )
        )
}

metrics_random <- if (has_random_method) {
    metrics_aug |>
        dplyr::filter(.data$MethodBase == "all_features" | .data$MethodBase == "random")
} else {
    metrics_aug
}

metrics_nonrandom <- if (has_random_method) {
    metrics_aug |>
        dplyr::filter(.data$MethodBase == "all_features" | .data$MethodBase != "random")
} else {
    metrics_aug
}

observed_value_data <- metrics_aug |>
    dplyr::filter(!is.na(.data$ScaledValue)) |>
    dplyr::group_by(.data$MetricDisplay, .data$Type) |>
    dplyr::mutate(MeanValue = mean(.data$ScaledValue, na.rm = TRUE)) |>
    dplyr::ungroup()

obs_range <- ggplot2::ggplot(
    metrics_aug |>
        dplyr::filter(!is.na(.data$ScaledValue)) |>
        dplyr::group_by(.data$MetricDisplay, .data$Type) |>
        dplyr::mutate(MeanValue = mean(.data$ScaledValue, na.rm = TRUE)) |>
        dplyr::ungroup(),
    ggplot2::aes(x = .data$ScaledValue, y = .data$MetricDisplay, fill = .data$MeanValue, group = .data$MetricDisplay)
) +
    ggridges::geom_density_ridges(
        ggplot2::aes(height = after_stat(ndensity)),
        scale = 0.9,
        rel_min_height = 0.01,
        colour = "white",
        linewidth = 0.3,
        alpha = 1,
        jittered_points = FALSE,
        bandwidth = 0.05,
        from = 0,
        to = 1,
        quantile_lines = TRUE,
        quantiles = 2,
        quantile_fun = function(x, probs) stats::quantile(x, probs = 0.5, na.rm = TRUE),
        vline_size = 0.25,
        vline_color = "#f7f7f7",
        vline_linetype = "solid"
    ) +
    ggplot2::scale_fill_viridis_c(
        option = "magma",
        limits = c(0, 1),
        guide = "none"
    ) +
    ggplot2::scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    ggplot2::facet_grid(
        rows = ggplot2::vars(.data$Type),
        scales = "free_y",
        space = "free_y"
    ) +
    ggplot2::labs(title = "Observed\nrange", x = NULL, y = NULL) +
    theme_features_pub() +
    ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 8.2, lineheight = 0.92, margin = ggplot2::margin(b = 4)),
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(fill = NA, colour = "grey20", linewidth = 0.35),
        strip.text.y = ggplot2::element_blank(),
        strip.background = ggplot2::element_blank(),
        panel.spacing = ggplot2::unit(0.08, "cm")
    )

feature_cor <- metrics_nonrandom |>
    dplyr::filter(is.finite(.data$SelFeaturesNum)) |>
    dplyr::group_by(.data$MetricDisplay, .data$Type, .data$Dataset, .data$IntegrationLabel) |>
    dplyr::summarise(
        Corr = suppressWarnings(stats::cor(.data$ScaledValue, .data$SelFeaturesNum, method = "spearman", use = "pairwise.complete.obs")),
        .groups = "drop"
    ) |>
    dplyr::filter(is.finite(.data$Corr))

feature_cor_data <- feature_cor |>
    dplyr::group_by(.data$MetricDisplay, .data$Type) |>
    dplyr::mutate(Mean = mean(.data$Corr, na.rm = TRUE)) |>
    dplyr::ungroup()

feature_cor_plot <- ggplot2::ggplot(
    feature_cor_data,
    ggplot2::aes(x = .data$Corr, y = .data$MetricDisplay, fill = .data$Mean, group = .data$MetricDisplay)
) +
    ggplot2::geom_vline(xintercept = c(-0.5, 0, 0.5), colour = "grey80", linewidth = 0.25) +
    ggridges::geom_density_ridges(
        ggplot2::aes(height = after_stat(ndensity)),
        scale = 0.9,
        rel_min_height = 0.01,
        colour = "white",
        linewidth = 0.3,
        alpha = 1,
        jittered_points = FALSE,
        bandwidth = 0.1,
        from = -1,
        to = 1,
        quantile_lines = TRUE,
        quantiles = 2,
        quantile_fun = function(x, probs) stats::quantile(x, probs = 0.5, na.rm = TRUE),
        vline_size = 0.28,
        vline_color = "#f7f7f7",
        vline_linetype = "solid"
    ) +
    colorspace::scale_fill_continuous_divergingx(
        palette = "Zissou 1",
        limits = c(-1, 1),
        guide = "none"
    ) +
    ggplot2::scale_x_continuous(limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1)) +
    ggplot2::facet_grid(
        rows = ggplot2::vars(.data$Type),
        scales = "free_y",
        space = "free_y"
    ) +
    ggplot2::labs(title = "Correlation\nwith number\nof features", x = NULL, y = NULL, fill = NULL) +
    theme_features_pub() +
    ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 8.2, lineheight = 0.92, margin = ggplot2::margin(b = 4)),
        legend.position = "none",
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(fill = NA, colour = "grey20", linewidth = 0.35),
        strip.text.y = ggplot2::element_blank(),
        strip.background = ggplot2::element_rect(fill = "#2d2d2d", colour = "#2d2d2d", linewidth = 0.4),
        panel.spacing = ggplot2::unit(0.08, "cm")
    )

dataset_feature_cor <- dplyr::bind_rows(
    metrics_aug |>
        dplyr::group_by(.data$MetricDisplay, .data$Type) |>
        dplyr::summarise(
            Feature = "NSlices",
            Corr = suppressWarnings(stats::cor(.data$ScaledValue, .data$NSlices, method = "spearman", use = "pairwise.complete.obs")),
            .groups = "drop"
        ),
    metrics_aug |>
        dplyr::group_by(.data$MetricDisplay, .data$Type) |>
        dplyr::summarise(
            Feature = "log10(Spots)",
            Corr = suppressWarnings(stats::cor(.data$ScaledValue, .data$LogSpots, method = "spearman", use = "pairwise.complete.obs")),
            .groups = "drop"
        ),
    metrics_aug |>
        dplyr::group_by(.data$MetricDisplay, .data$Type) |>
        dplyr::summarise(
            Feature = "log10(Genes)",
            Corr = suppressWarnings(stats::cor(.data$ScaledValue, .data$LogFeatures, method = "spearman", use = "pairwise.complete.obs")),
            .groups = "drop"
        ),
    metrics_aug |>
        dplyr::group_by(.data$MetricDisplay, .data$Type) |>
        dplyr::summarise(
            Feature = "log10(Labels)",
            Corr = suppressWarnings(stats::cor(.data$ScaledValue, .data$LogLabels, method = "spearman", use = "pairwise.complete.obs")),
            .groups = "drop"
        ),
    metrics_aug |>
        dplyr::group_by(.data$MetricDisplay, .data$Type) |>
        dplyr::summarise(
            Feature = "log10(SpotsPerSlice)",
            Corr = suppressWarnings(stats::cor(.data$ScaledValue, .data$LogSpotsPerSlice, method = "spearman", use = "pairwise.complete.obs")),
            .groups = "drop"
        ),
    metrics_aug |>
        dplyr::group_by(.data$MetricDisplay, .data$Type) |>
        dplyr::summarise(
            Feature = "Slide-seq",
            Corr = suppressWarnings(stats::cor(.data$ScaledValue, .data$IsSlideSeq, method = "pearson", use = "pairwise.complete.obs")),
            .groups = "drop"
        )
    ) |>
    dplyr::mutate(
        Feature = dplyr::recode(
            .data$Feature,
            "NSlices" = "slices",
            "log10(Spots)" = "spots",
            "log10(Genes)" = "genes",
            "log10(Labels)" = "labels",
            "log10(SpotsPerSlice)" = "spots/slice",
            "Slide-seq" = "platform"
        ),
        Feature = factor(.data$Feature, levels = c("slices", "spots", "genes", "labels", "spots/slice", "platform")),
        MetricDisplay = factor(.data$MetricDisplay, levels = rev(metric_levels)),
        Type = factor(.data$Type, levels = type_levels)
    )

dataset_feature_plot <- ggplot2::ggplot(
    dataset_feature_cor,
    ggplot2::aes(x = .data$Feature, y = .data$MetricDisplay, colour = .data$Corr)
) +
    ggplot2::geom_tile(fill = "#f4f4f4", colour = "#e9e9e9", linewidth = 0.16) +
    ggplot2::geom_point(shape = 15, size = 3.0, alpha = 0.98) +
    ggplot2::scale_colour_gradient2(
        low = "#b2182b",
        mid = "#f7f7f7",
        high = "#2166ac",
        space = "Lab",
        midpoint = 0,
        limits = c(-1, 1),
        na.value = "#cccccc",
        guide = ggplot2::guide_colourbar(
            order = 2,
            direction = "horizontal",
            title.position = "top",
            barwidth = grid::unit(2.2, "cm"),
            barheight = grid::unit(0.26, "cm")
        )
    ) +
    ggplot2::facet_grid(
        rows = ggplot2::vars(.data$Type),
        scales = "free_y",
        space = "free_y"
    ) +
    ggplot2::labs(title = "Correlation with\ndataset features", x = NULL, y = NULL, colour = "Mean dataset-feature\ncorrelation") +
    theme_features_pub() +
    ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 8.2, lineheight = 0.92, margin = ggplot2::margin(b = 4)),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(size = 7, angle = 90, hjust = 0.5, vjust = 0.5),
        panel.grid = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(fill = NA, colour = "grey20", linewidth = 0.35),
        strip.text.y = ggplot2::element_blank(),
        strip.background = ggplot2::element_rect(fill = "#2d2d2d", colour = "#2d2d2d", linewidth = 0.4),
        panel.spacing = ggplot2::unit(0.08, "cm")
    )

metric_matrix <- metrics_random |>
    dplyr::select(
        .data$Dataset,
        .data$Method,
        .data$IntegrationLabel,
        .data$SelFeatures,
        .data$MetricDisplay,
        .data$ScaledValue
    ) |>
    dplyr::group_by(.data$Dataset, .data$Method, .data$IntegrationLabel, .data$MetricDisplay) |>
    dplyr::summarise(Value = mean(.data$ScaledValue, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(
        names_from = "MetricDisplay",
        values_from = "Value"
    )

metric_info <- metrics_random |>
    dplyr::distinct(.data$MetricDisplay, .data$Type) |>
    dplyr::arrange(.data$Type, .data$MetricDisplay) |>
    dplyr::mutate(
        MetricDisplay = factor(.data$MetricDisplay, levels = metric_levels),
        Type = factor(.data$Type, levels = type_levels)
    ) |>
    dplyr::arrange(.data$Type, .data$MetricDisplay)

wide_cols <- setdiff(names(metric_matrix), c("Dataset", "Method", "IntegrationLabel"))
metric_col_map <- metric_info |>
    dplyr::mutate(ColName = as.character(.data$MetricDisplay)) |>
    dplyr::filter(.data$ColName %in% wide_cols)

metric_corr_raw <- metric_matrix |>
    dplyr::group_by(.data$Dataset, .data$IntegrationLabel) |>
    dplyr::group_split() |>
    lapply(function(.df) {
        vals <- .df[, wide_cols, drop = FALSE]
        if (nrow(vals) < 2) {
            return(NULL)
        }
        out <- list()
        for (i in seq_len(nrow(metric_col_map))) {
            for (j in seq_len(nrow(metric_col_map))) {
                xi <- vals[[metric_col_map$ColName[[i]]]]
                xj <- vals[[metric_col_map$ColName[[j]]]]
                corr <- suppressWarnings(stats::cor(xi, xj, method = "spearman", use = "pairwise.complete.obs"))
                out[[length(out) + 1]] <- data.frame(
                    RowMetric = metric_col_map$MetricDisplay[[i]],
                    RowType = metric_col_map$Type[[i]],
                    ColMetric = metric_col_map$MetricDisplay[[j]],
                    ColType = metric_col_map$Type[[j]],
                    Corr = corr
                )
            }
        }
        dplyr::bind_rows(out)
    }) |>
    dplyr::bind_rows()

metric_corr <- metric_corr_raw |>
    dplyr::group_by(.data$RowMetric, .data$RowType, .data$ColMetric, .data$ColType) |>
    dplyr::summarise(
        CorrSD = stats::sd(.data$Corr, na.rm = TRUE),
        Corr = mean(.data$Corr, na.rm = TRUE),
        .groups = "drop"
    ) |>
    dplyr::mutate(
        RowMetric = factor(.data$RowMetric, levels = rev(metric_levels)),
        ColMetric = factor(.data$ColMetric, levels = metric_levels),
        RowType = factor(.data$RowType, levels = type_levels),
        ColType = factor(.data$ColType, levels = type_levels),
        CorrSD = dplyr::if_else(is.na(.data$CorrSD), 0, .data$CorrSD)
    )

metric_corr_plot <- ggplot2::ggplot(
    metric_corr,
    ggplot2::aes(x = .data$ColMetric, y = .data$RowMetric)
) +
    ggplot2::geom_tile(fill = "#f4f4f4", colour = "#e9e9e9", linewidth = 0.16) +
    ggplot2::geom_point(
        ggplot2::aes(colour = .data$Corr, size = .data$CorrSD),
        shape = 15,
        alpha = 0.98
    ) +
    colorspace::scale_colour_continuous_diverging(
        palette = "Purple-Green",
        rev = TRUE,
        limits = c(-1, 1),
        guide = ggplot2::guide_colourbar(
            order = 4,
            direction = "horizontal",
            title.position = "top",
            barwidth = grid::unit(2.2, "cm"),
            barheight = grid::unit(0.24, "cm")
        )
    ) +
    ggplot2::scale_size_continuous(
        trans = "reverse",
        range = c(0.7, 3.2),
        guide = ggplot2::guide_legend(
            order = 3,
            nrow = 1,
            byrow = TRUE,
            title.position = "top",
            override.aes = list(colour = "#9f9f9f", alpha = 1, shape = 15)
        )
    ) +
    ggplot2::facet_grid(
        rows = ggplot2::vars(.data$RowType),
        cols = ggplot2::vars(.data$ColType),
        scales = "free",
        space = "free"
    ) +
    ggplot2::labs(
        title = "Correlations between metrics",
        x = NULL,
        y = NULL,
        colour = "Mean metric\ncorrelation",
        size = "Correlation\ns.d."
    ) +
    theme_features_pub() +
    ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 10),
        axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.spacing = grid::unit(0.12, "cm"),
        strip.text.x = ggplot2::element_text(colour = "white", face = "bold", size = 8, margin = ggplot2::margin(3, 0, 3, 0)),
        strip.text.y = ggplot2::element_text(angle = 270, colour = "white", face = "bold", size = 7.5, margin = ggplot2::margin(0, 4, 0, 4)),
        strip.background = ggplot2::element_rect(fill = "black", colour = "black", linewidth = 0.55),
        panel.border = ggplot2::element_rect(fill = NA, colour = "grey15", linewidth = 0.35)
    )

readr::write_tsv(
    feature_cor |>
        dplyr::rename(FeatureCorrelation = .data$Corr),
    file.path(output_dir, "metric_feature_number_correlations.tsv")
)
readr::write_tsv(
    dataset_feature_cor,
    file.path(output_dir, "metric_dataset_feature_correlations.tsv")
)
readr::write_tsv(
    metric_corr |>
        dplyr::mutate(RowMetric = as.character(.data$RowMetric), ColMetric = as.character(.data$ColMetric)),
    file.path(output_dir, "metric_metric_correlations.tsv")
)

legend_palette <- function(n, low, mid = "#ffffff", high) {
    grDevices::colorRampPalette(c(low, mid, high))(n)
}

legend_gradients <- dplyr::bind_rows(
    data.frame(
        section = "feature",
        x = seq(0.22, 0.78, length.out = 180),
        fill_hex = legend_palette(180, "#c47d00", "#f5f5f5", "#1a6fa8")
    ),
    data.frame(
        section = "metric",
        x = seq(0.10, 0.90, length.out = 180),
        fill_hex = legend_palette(180, "#8a5bb0", "#ffffff", "#4d9a56")
    )
)

legend_base <- ggplot2::theme_void() +
    ggplot2::theme(
        plot.margin = ggplot2::margin(0, 1.5, 0, 1.5, unit = "mm")
    )

value_gradient <- data.frame(
    x = seq(0.14, 0.86, length.out = 180),
    fill_hex = viridisLite::magma(180)
)
legend_mean_value <- ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0.00, y = 0.95, hjust = 0, label = "Mean range", size = 2.5, fontface = "bold") +
    ggplot2::geom_tile(data = value_gradient, ggplot2::aes(x = .data$x, y = 0.56, fill = .data$fill_hex), width = 0.005, height = 0.22, inherit.aes = FALSE) +
    ggplot2::annotate("segment", x = c(0.14, 0.32, 0.50, 0.68, 0.86), xend = c(0.14, 0.32, 0.50, 0.68, 0.86), y = 0.45, yend = 0.42, colour = "grey40", linewidth = 0.25) +
    ggplot2::annotate("text", x = c(0.14, 0.32, 0.50, 0.68, 0.86), y = 0.30, label = c("0.0", "0.25", "0.5", "0.75", "1.0"), size = 2.0, fontface = "bold") +
    ggplot2::scale_fill_identity() +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    legend_base

feature_gradient <- legend_gradients[legend_gradients$section == "feature", ]
feature_number_gradient <- data.frame(
    x = seq(0.14, 0.86, length.out = 180),
    fill_hex = colorspace::divergingx_hcl(180, palette = "Zissou 1")
)
dataset_feature_gradient <- data.frame(
    x = seq(0.14, 0.86, length.out = 180),
    fill_hex = scales::div_gradient_pal(low = "#b2182b", mid = "#f7f7f7", high = "#2166ac", space = "Lab")(seq(0, 1, length.out = 180))
)
legend_feature_number_corr <- ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0.00, y = 0.95, hjust = 0, label = "Mean feature-number correlation", size = 2.5, fontface = "bold") +
    ggplot2::geom_tile(data = feature_number_gradient, ggplot2::aes(x = .data$x, y = 0.56, fill = .data$fill_hex), width = 0.004, height = 0.22, inherit.aes = FALSE) +
    ggplot2::annotate("segment", x = c(0.14, 0.32, 0.50, 0.68, 0.86), xend = c(0.14, 0.32, 0.50, 0.68, 0.86), y = 0.45, yend = 0.42, colour = "grey40", linewidth = 0.25) +
    ggplot2::annotate("text", x = c(0.14, 0.32, 0.50, 0.68, 0.86), y = 0.30, label = c("-1.0", "-0.5", "0.0", "0.5", "1.0"), size = 2.0, fontface = "bold") +
    ggplot2::scale_fill_identity() +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    legend_base

legend_dataset_feature_corr <- ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0.00, y = 0.95, hjust = 0, label = "Mean dataset-feature correlation", size = 2.5, fontface = "bold") +
    ggplot2::geom_tile(data = dataset_feature_gradient, ggplot2::aes(x = .data$x, y = 0.56, fill = .data$fill_hex), width = 0.004, height = 0.22, inherit.aes = FALSE) +
    ggplot2::annotate("segment", x = c(0.14, 0.32, 0.50, 0.68, 0.86), xend = c(0.14, 0.32, 0.50, 0.68, 0.86), y = 0.45, yend = 0.42, colour = "grey40", linewidth = 0.25) +
    ggplot2::annotate("text", x = c(0.14, 0.32, 0.50, 0.68, 0.86), y = 0.30, label = c("-1.0", "-0.5", "0.0", "0.5", "1.0"), size = 2.0, fontface = "bold") +
    ggplot2::annotate("rect", xmin = 0.84, xmax = 0.91, ymin = 0.45, ymax = 0.67, fill = "#cccccc", colour = "#aaaaaa", linewidth = 0.25) +
    ggplot2::annotate("text", x = 0.94, y = 0.56, hjust = 0, label = "NA", size = 2.1, fontface = "bold") +
    ggplot2::scale_fill_identity() +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    legend_base

legend_sd <- ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0.00, y = 0.95, hjust = 0, label = "Correlation s.d.", size = 2.5, fontface = "bold") +
    ggplot2::annotate("point", x = c(0.10, 0.28, 0.50, 0.72, 0.90), y = rep(0.56, 5),
        size = c(4.5, 3.5, 2.6, 1.8, 1.1), shape = 15, colour = "#b0b0b0") +
    ggplot2::annotate("text", x = c(0.10, 0.28, 0.50, 0.72, 0.90), y = 0.30, label = c("0.0", "0.25", "0.5", "0.75", "1.0"), size = 2.0, fontface = "bold") +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    legend_base

metric_gradient <- data.frame(
    x = seq(0.10, 0.90, length.out = 180),
    fill_hex = colorspace::diverging_hcl(180, palette = "Purple-Green", rev = TRUE)
)
legend_metric_corr <- ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0.00, y = 0.95, hjust = 0, label = "Mean metric correlation", size = 2.5, fontface = "bold") +
    ggplot2::geom_tile(data = metric_gradient, ggplot2::aes(x = .data$x, y = 0.56, fill = .data$fill_hex), width = 0.0048, height = 0.22, inherit.aes = FALSE) +
    ggplot2::annotate("segment", x = c(0.10, 0.30, 0.50, 0.70, 0.90), xend = c(0.10, 0.30, 0.50, 0.70, 0.90), y = 0.45, yend = 0.42, colour = "grey40", linewidth = 0.25) +
    ggplot2::annotate("text", x = c(0.10, 0.30, 0.50, 0.70, 0.90), y = 0.30, label = c("-1.0", "-0.5", "0.0", "0.5", "1.0"), size = 2.0, fontface = "bold") +
    ggplot2::scale_fill_identity() +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    legend_base

legend_strip <- patchwork::wrap_plots(
    legend_mean_value,
    legend_feature_number_corr,
    legend_dataset_feature_corr,
    legend_sd,
    legend_metric_corr,
    nrow = 1,
    widths = c(1.95, 2.05, 2.20, 1.45, 2.05)
)

figure <- patchwork::wrap_plots(
    obs_range,
    feature_cor_plot,
    dataset_feature_plot,
    metric_corr_plot,
    ncol = 4,
    widths = c(1.35, 1.2, 1.2, 4.2)
) &
    ggplot2::theme(
        legend.position = "none",
        text = ggplot2::element_text(family = "Arial", face = "plain"),
        axis.text = ggplot2::element_text(size = 6.2),
        axis.title = ggplot2::element_text(size = 7.0),
        strip.text = ggplot2::element_text(face = "bold", size = 8.2),
        plot.title = ggplot2::element_text(face = "bold", size = 9.2)
    )

figure <- figure / legend_strip + patchwork::plot_layout(heights = c(1, 0.145))

png_path <- file.path(output_dir, "figure_metric_overview.png")
pdf_path <- file.path(output_dir, "figure_metric_overview.pdf")
saveRDS(figure, file.path(output_dir, "figure_metric_overview.rds"))
ggplot2::ggsave(png_path, figure, width = 11.2, height = 7.5, dpi = 600)
register_arial_pdf_font()
ggplot2::ggsave(
    pdf_path, figure,
    width = 11.2, height = 7.5,
    device = grDevices::pdf, family = "Arial", useDingbats = FALSE
)
