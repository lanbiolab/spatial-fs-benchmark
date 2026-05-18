#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(ggplot2)
    library(patchwork)
    library(colorspace)
})

source(file.path("external", "atlas-feature-selection-benchmark", "analysis", "R", "plotting.R"))

args <- commandArgs(trailingOnly = TRUE)
methods_path <- if (length(args) >= 1) args[[1]] else file.path("results", "fig5a_spatial_lineages", "figures", "fig5a_lineage_methods.tsv")
scores_path <- if (length(args) >= 2) args[[2]] else file.path("results", "fig5d_spatial_lineages", "figures", "fig5d_lineage_celltype_overlap.tsv")
diffs_path <- if (length(args) >= 3) args[[3]] else file.path("results", "fig5d_spatial_lineages", "figures", "fig5d_lineage_celltype_overlap_diff.tsv")
output_dir <- if (length(args) >= 4) args[[4]] else file.path("results", "fig5d_spatial_lineages", "figures")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

methods_tbl <- read_tsv(methods_path, show_col_types = FALSE) |>
    mutate(DisplayName = gsub(" \\(N=[0-9]+\\)$", "", .data$Name))
method_levels <- rev(methods_tbl$DisplayName)

label_levels <- c("Langerhans Cell", "T Cell", "Keratinocyte", "Merkel Cell", "Stem Cell", "Progenitor Cell")
diff_label_levels <- c("Langerhans Cell", "T Cell", "Keratinocyte", "Merkel Cell", "Stem Cell", "Progenitor Cell")

scores <- read_tsv(scores_path, show_col_types = FALSE) |>
    mutate(
        Method = gsub(" \\(N=[0-9]+\\)$", "", .data$Method),
        Method = factor(.data$Method, levels = method_levels),
        Dataset = factor(.data$Dataset, levels = c("Full", "Immune", "Epithelial")),
        Label = factor(.data$Label, levels = label_levels)
    )

diffs <- read_tsv(diffs_path, show_col_types = FALSE) |>
    mutate(
        Method = gsub(" \\(N=[0-9]+\\)$", "", .data$Method),
        Method = factor(.data$Method, levels = method_levels),
        Dataset = factor(.data$Dataset, levels = c("Immune", "Epithelial")),
        Label = factor(.data$Label, levels = diff_label_levels)
    )

left_plot <- ggplot(scores, aes(x = .data$Label, y = .data$Method, fill = .data$Value)) +
    geom_tile() +
    facet_grid(. ~ .data$Dataset, scales = "free_x", space = "free_x") +
    scale_fill_viridis_c(option = "plasma", limits = c(0, 1), name = "Marker\noverlap") +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0), drop = FALSE) +
    labs(x = NULL, y = NULL) +
    theme_features_pub() +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, lineheight = 0.8, face = "bold", colour = "black", size = 5),
        axis.text.y = element_text(face = "bold", colour = "black", size = 5.6),
        panel.grid = element_blank(),
        plot.margin = margin(0.05, 0.05, 0.05, 0.05, "cm")
    )

right_plot <- ggplot(diffs, aes(x = .data$Label, y = .data$Method, fill = .data$Difference)) +
    geom_tile() +
    facet_grid(. ~ .data$Dataset, scales = "free_x", space = "free_x") +
    colorspace::scale_fill_continuous_diverging(
        palette = "Purple-Green",
        limits = c(-1, 1),
        name = "Difference"
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0), drop = FALSE) +
    labs(x = NULL, y = NULL) +
    theme_features_pub() +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, lineheight = 0.8, face = "bold", colour = "black", size = 5),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        plot.margin = margin(0.05, 0.05, 0.05, 0.02, "cm")
    )

fig <- wrap_plots(
    left_plot,
    right_plot,
    nrow = 1,
    widths = c(1, 0.42),
    guides = "collect"
) &
    theme(
        legend.position = "bottom",
        legend.title.position = "top",
        legend.box = "horizontal",
        legend.box.just = "left"
    )

save_figure_files(fig, file.path(output_dir, "figure_fig5d_spatial_lineages"), width = 8.2, height = 5.7)
