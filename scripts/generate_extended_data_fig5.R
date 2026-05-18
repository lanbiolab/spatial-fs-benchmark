#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(ggplot2)
    library(patchwork)
    library(cowplot)
    library(grid)
})

source(file.path("external", "atlas-feature-selection-benchmark", "analysis", "R", "plotting.R"))

args <- commandArgs(trailingOnly = TRUE)
summary_path <- if (length(args) >= 1) args[[1]] else file.path("results", "extended_data_fig5", "figures", "extended_fig5_markers_summary.tsv")
celltype_path <- if (length(args) >= 2) args[[2]] else file.path("results", "extended_data_fig5", "figures", "extended_fig5_celltype_markers.tsv")
methods_path <- if (length(args) >= 3) args[[3]] else file.path("results", "fig5a_spatial_lineages", "figures", "fig5a_lineage_methods.tsv")
output_dir <- if (length(args) >= 4) args[[4]] else file.path("results", "extended_data_fig5", "figures")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

methods_tbl <- read_tsv(methods_path, show_col_types = FALSE) |>
    filter(.data$MethodBase != "random") |>
    mutate(DisplayName = gsub(" \\(N=[0-9]+\\)$", "", .data$Name))
method_levels <- rev(methods_tbl$DisplayName)

lineage_levels <- c("Immune", "Epithelial")
dataset_levels <- c("Full", "Immune", "Epithelial")
label_levels <- c(
    "Langerhans Cell", "T Cell",
    "Keratinocyte", "Merkel Cell", "Stem Cell", "Progenitor Cell"
)

summary_df <- read_tsv(summary_path, show_col_types = FALSE) |>
    filter(.data$Lineage %in% lineage_levels) |>
    mutate(
        Method = gsub(" \\(N=[0-9]+\\)$", "", .data$Method),
        Method = factor(.data$Method, levels = method_levels),
        Dataset = factor(.data$Dataset, levels = dataset_levels),
        Lineage = factor(.data$Lineage, levels = lineage_levels)
    ) |>
    mutate(
        Column = factor(
            paste(.data$Lineage, .data$Dataset, sep = "::"),
            levels = c(
                "Immune::Full", "Immune::Immune", "Immune::Epithelial",
                "Epithelial::Full", "Epithelial::Immune", "Epithelial::Epithelial"
            )
        )
    )

celltype_df <- read_tsv(celltype_path, show_col_types = FALSE) |>
    filter(.data$Lineage %in% lineage_levels) |>
    mutate(
        Method = gsub(" \\(N=[0-9]+\\)$", "", .data$Method),
        Method = factor(.data$Method, levels = method_levels),
        Dataset = factor(.data$Dataset, levels = dataset_levels),
        Label = factor(.data$Label, levels = label_levels),
        Lineage = factor(.data$Lineage, levels = lineage_levels)
    )

make_group_header <- function(labels, n_per_group) {
    n_groups <- length(labels)
    xmin <- seq(0, n_groups - 1) / n_groups
    xmax <- seq(1, n_groups) / n_groups
    ggplot() +
        annotate("rect", xmin = xmin, xmax = xmax, ymin = 0, ymax = 1, fill = "black", colour = "white", linewidth = 0.22) +
        annotate("text", x = (xmin + xmax) / 2, y = 0.5, label = labels, colour = "white", size = 3.0, fontface = "bold") +
        coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE, clip = "off") +
        theme_void() +
        theme(plot.margin = margin(0, 0, 0, 0))
}

make_title_bar <- function(text) {
    ggplot() +
        annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1, fill = "black", colour = "black") +
        annotate("text", x = 0.5, y = 0.5, label = text, colour = "white", size = 3.1, fontface = "bold") +
        coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE, clip = "off") +
        theme_void() +
        theme(plot.margin = margin(0, 0, 0, 0))
}

base_theme <- theme_features_pub() +
    theme(
        panel.grid = element_blank(),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, unit = "mm"),
        axis.title = element_blank(),
        legend.position = "bottom",
        legend.title.position = "top",
        legend.text = element_text(size = 5.2, face = "bold", colour = "black"),
        legend.title = element_text(size = 5.3, face = "bold", colour = "black"),
        legend.key.width = unit(0.95, "cm"),
        legend.key.height = unit(0.13, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, 0, 0)
    )

a_mean_heat <- ggplot(summary_df, aes(x = .data$Column, y = .data$Method, fill = .data$MeanProp)) +
    geom_tile(width = 0.78, height = 0.72, colour = "white", linewidth = 0.18) +
    geom_vline(xintercept = 3.5, colour = "grey20", linewidth = 0.32) +
    scale_fill_viridis_c(option = "C", limits = c(0, 1), name = "Mean proportion of markers") +
    base_theme +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 4.5, face = "bold", colour = "black"),
        axis.text.y = element_text(size = 5.1, face = "bold", colour = "black")
    )

a_sd_heat <- ggplot(summary_df, aes(x = .data$Column, y = .data$Method, fill = .data$SDProp)) +
    geom_tile(width = 0.78, height = 0.72, colour = "white", linewidth = 0.18) +
    geom_vline(xintercept = 3.5, colour = "grey20", linewidth = 0.32) +
    scale_fill_viridis_c(option = "cividis", limits = c(0, max(summary_df$SDProp, na.rm = TRUE)), name = "SD proportion of markers") +
    base_theme +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 4.5, face = "bold", colour = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
    )

a_mean <- wrap_plots(
    make_group_header(lineage_levels, 3),
    a_mean_heat + theme(legend.position = "none"),
    ncol = 1,
    heights = c(0.08, 1)
)

a_sd <- wrap_plots(
    make_group_header(lineage_levels, 3),
    a_sd_heat + theme(legend.position = "none"),
    ncol = 1,
    heights = c(0.08, 1)
)

a_panel <- wrap_plots(a_mean, a_sd, nrow = 1, widths = c(1, 1))

make_dataset_panel <- function(ds) {
    dat <- filter(celltype_df, .data$Dataset == ds)
    heat <- ggplot(dat, aes(x = .data$Label, y = .data$Method, fill = .data$Value)) +
        geom_tile(width = 0.82, height = 0.72, colour = "white", linewidth = 0.18) +
        geom_vline(xintercept = 2.5, colour = "grey20", linewidth = 0.32) +
        scale_fill_viridis_c(option = "plasma", limits = c(0, 1), name = "Proportion of marker genes") +
        base_theme +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, lineheight = 0.8, size = 4.4, face = "bold", colour = "black"),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank()
        )
    wrap_plots(
        make_title_bar(ds),
        heat + theme(legend.position = "none"),
        ncol = 1,
        heights = c(0.08, 1)
    )
}

b_panel <- wrap_plots(
    make_dataset_panel("Full"),
    make_dataset_panel("Immune"),
    make_dataset_panel("Epithelial"),
    nrow = 1,
    widths = c(1, 1, 1)
)

legend_mean <- cowplot::get_legend(a_mean_heat)
legend_sd <- cowplot::get_legend(a_sd_heat)
legend_prop <- cowplot::get_legend(
    ggplot(celltype_df, aes(x = .data$Label, y = .data$Method, fill = .data$Value)) +
        geom_tile() +
        scale_fill_viridis_c(option = "plasma", limits = c(0, 1), name = "Proportion of marker genes") +
        base_theme +
        theme_void()
)

legend_row <- wrap_plots(
    wrap_plots(legend_mean, legend_sd, nrow = 1, widths = c(1, 1)),
    legend_prop,
    ncol = 1,
    heights = c(1, 1)
)

figure <- wrap_plots(
    a_panel,
    b_panel,
    legend_row,
    ncol = 1,
    heights = c(1.0, 1.12, 0.24)
) &
    theme(plot.margin = margin(1.5, 1.5, 1.5, 1.5, unit = "mm"))

save_figure_files(
    figure,
    file.path(output_dir, "extended_data_fig5"),
    width = 8.3,
    height = 6.4
)
