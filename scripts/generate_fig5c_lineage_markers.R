#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(ggplot2)
})

source(file.path("external", "atlas-feature-selection-benchmark", "analysis", "R", "plotting.R"))

args <- commandArgs(trailingOnly = TRUE)
markers_path <- if (length(args) >= 1) args[[1]] else file.path("results", "fig5c_spatial_lineages", "figures", "fig5c_lineage_marker_overlap.tsv")
methods_path <- if (length(args) >= 2) args[[2]] else file.path("results", "fig5a_spatial_lineages", "figures", "fig5a_lineage_methods.tsv")
output_dir <- if (length(args) >= 3) args[[3]] else file.path("results", "fig5c_spatial_lineages", "figures")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

methods_tbl <- read_tsv(methods_path, show_col_types = FALSE)
method_levels <- rev(methods_tbl$Name)

markers <- read_tsv(markers_path, show_col_types = FALSE) |>
    mutate(
        Method = factor(.data$Method, levels = method_levels),
        Dataset = factor(.data$Dataset, levels = c("Full", "Immune", "Epithelial")),
        Lineage = factor(.data$Lineage, levels = c("Epithelial", "Immune"))
    )

fig <- ggplot(markers, aes(x = .data$Dataset, y = .data$Method)) +
    geom_point(
        aes(fill = .data$MeanProp, size = .data$SDProp),
        shape = 22, colour = "white", stroke = 0.2,
        na.rm = TRUE
    ) +
    facet_wrap(~ .data$Lineage, nrow = 1) +
    scale_fill_viridis_c(option = "C", limits = c(0, 1), name = "Mean proportion") +
    scale_size_continuous(
        trans = "reverse",
        range = c(0.7, 4.2),
        name = "s.d. proportion"
    ) +
    labs(x = NULL, y = NULL) +
    theme_features_pub() +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "bold", colour = "black"),
        axis.text.y = element_text(face = "bold", colour = "black", size = 6),
        strip.background = element_rect(fill = "black", colour = "black"),
        strip.text = element_text(colour = "white", face = "bold"),
        panel.grid = element_blank(),
        legend.position = "bottom",
        plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm")
    )

save_figure_files(fig, file.path(output_dir, "figure_fig5c_spatial_lineages"), width = 4.8, height = 7.8)
