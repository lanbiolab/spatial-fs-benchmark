#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(args) >= 1) args[[1]] else file.path("results", "num_features_benchmark", "figures")
output_dir <- if (length(args) >= 2) args[[2]] else file.path("results", "extended_data_fig2", "figures")

source(file.path("external", "atlas-feature-selection-benchmark", "analysis", "R", "plotting.R"))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

dataset_summary <- readr::read_tsv(
    file.path(input_dir, "dataset_feature_number_summary.tsv"),
    show_col_types = FALSE
)

method_summary <- readr::read_tsv(
    file.path(input_dir, "method_feature_number_summary.tsv"),
    show_col_types = FALSE
)

sel_levels <- c("100", "200", "500", "1000", "2000", "5000", "10000")
type_levels <- c("Integration", "Clustering", "Alignment", "Overall")
dataset_order <- c(
    "DLPFC",
    "MouseBrainSerialSections",
    "STOmicsVisium5Samples",
    "STOmics0212",
    "STOmics0218",
    "STOmics0224",
    "E8p5Embryo",
    "E9p5Embryo"
)

# Match the original Extended Data Fig. 2 philosophy: show a compact set of
# representative methods rather than the full method pool.
method_order <- c(
    "Statistic variance",
    "Statistic mean",
    "Seurat-VST",
    "Seurat-scTransform",
    "Seurat-Dispersion",
    "scanpy-SeuratV3",
    "scanpy-Seurat",
    "scanpy-Pearson",
    "scanpy-CellRanger"
)

dataset_complete <- tidyr::expand_grid(
    Dataset = dataset_order,
    Type = factor(type_levels, levels = type_levels),
    SelFeatures = factor(sel_levels, levels = sel_levels)
) |>
    dplyr::left_join(
        dataset_summary |>
            dplyr::mutate(
                Type = factor(.data$Type, levels = type_levels),
                SelFeatures = factor(as.character(.data$SelFeatures), levels = sel_levels)
            ),
        by = c("Dataset", "Type", "SelFeatures")
    ) |>
    dplyr::mutate(Missing = is.na(.data$StandardValue))

method_complete <- tidyr::expand_grid(
    Method = method_order,
    Type = factor(type_levels, levels = type_levels),
    SelFeatures = factor(sel_levels, levels = sel_levels)
) |>
    dplyr::left_join(
        method_summary |>
            dplyr::filter(.data$Method %in% method_order) |>
            dplyr::mutate(
                Method = factor(.data$Method, levels = method_order),
                Type = factor(.data$Type, levels = type_levels),
                SelFeatures = factor(as.character(.data$SelFeatures), levels = sel_levels)
            ),
        by = c("Method", "Type", "SelFeatures")
    ) |>
    dplyr::mutate(Missing = is.na(.data$StandardValue))

mean_limits <- range(
    c(dataset_complete$StandardValue, method_complete$StandardValue),
    na.rm = TRUE
)
sd_limits <- range(
    c(dataset_complete$SD, method_complete$SD),
    na.rm = TRUE
)

theme_ext2 <- theme_features_pub() +
    theme(
        legend.position = "bottom",
        legend.title.position = "top",
        panel.grid = element_blank(),
        panel.spacing = unit(0.03, "cm"),
        plot.margin = margin(0, 0, 0, 0),
        strip.background = element_rect(fill = "black", colour = "black", linewidth = 0.3),
        strip.text = element_text(size = 5.2, colour = "white", margin = margin(2.5, 2.5, 2.5, 2.5)),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 5.0),
        legend.key.width = unit(0.72, "cm"),
        legend.key.height = unit(0.18, "cm"),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5.3)
    )

plot_heat <- function(df, y_var, y_levels, value_col, title = NULL, is_bottom = FALSE, legend_title = NULL) {
    ggplot(df, aes(x = .data$SelFeatures, y = .data[[y_var]])) +
        geom_point(
            data = \(x) dplyr::filter(x, .data$Missing),
            shape = 15,
            size = 2.85,
            colour = "#d9d9d9"
        ) +
        geom_point(
            data = \(x) dplyr::filter(x, !.data$Missing),
            aes(colour = .data[[value_col]]),
            shape = 15,
            size = 2.85
        ) +
        {
            if (value_col == "StandardValue") {
                colorspace::scale_colour_continuous_diverging(
                    palette = "Purple-Green",
                    limits = mean_limits,
                    name = legend_title %||% "Mean standardised value",
                    na.value = "#d9d9d9"
                )
            } else {
                scale_colour_viridis_c(
                    option = "cividis",
                    limits = sd_limits,
                    name = legend_title %||% "Standard deviation of\nstandardised values",
                    na.value = "#d9d9d9"
                )
            }
        } +
        facet_grid(. ~ .data$Type) +
        labs(title = title) +
        theme_ext2 +
        theme(
            plot.title = element_text(size = 6.3, face = "plain", hjust = 0, margin = margin(b = 2)),
            axis.text.x = if (is_bottom) element_text(angle = 90, hjust = 1, vjust = 0.5, size = 4.6) else element_blank(),
            axis.ticks.x = if (is_bottom) element_line(linewidth = 0.2) else element_blank(),
            axis.text.y = element_text(size = 5.0)
        )
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

datasets_means_plot <- plot_heat(
    dataset_complete |> dplyr::mutate(Dataset = factor(.data$Dataset, levels = rev(dataset_order))),
    "Dataset",
    dataset_order,
    "StandardValue",
    title = "Datasets",
    is_bottom = FALSE,
    legend_title = "Mean standardised value"
)

datasets_sds_plot <- plot_heat(
    dataset_complete |> dplyr::mutate(Dataset = factor(.data$Dataset, levels = rev(dataset_order))),
    "Dataset",
    dataset_order,
    "SD",
    title = NULL,
    is_bottom = FALSE,
    legend_title = "Standard deviation of\nstandardised values"
)

methods_means_plot <- plot_heat(
    method_complete |> dplyr::mutate(Method = factor(.data$Method, levels = method_order)),
    "Method",
    method_order,
    "StandardValue",
    title = "Methods",
    is_bottom = FALSE,
    legend_title = "Mean standardised value"
)

methods_sds_plot <- plot_heat(
    method_complete |> dplyr::mutate(Method = factor(.data$Method, levels = method_order)),
    "Method",
    method_order,
    "SD",
    title = NULL,
    is_bottom = TRUE,
    legend_title = "Standard deviation of\nstandardised values"
)

figure <- wrap_plots(
    datasets_means_plot,
    datasets_sds_plot,
    methods_means_plot,
    methods_sds_plot,
    ncol = 1,
    heights = c(1, 1, 1, 1),
    guides = "collect"
) &
    theme(
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.title.position = "top"
    )

readr::write_tsv(dataset_complete, file.path(output_dir, "extended_fig2_dataset_summary.tsv"))
readr::write_tsv(method_complete, file.path(output_dir, "extended_fig2_method_summary.tsv"))

save_figure_files(
    figure,
    file.path(output_dir, "extended_data_fig2"),
    width = 8.2,
    height = 6.2
)
