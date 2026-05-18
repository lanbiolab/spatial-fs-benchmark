#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(ggplot2)
    library(patchwork)
    library(grid)
})

source(file.path("external", "atlas-feature-selection-benchmark", "analysis", "R", "plotting.R"))

args <- commandArgs(trailingOnly = TRUE)
methods_path <- if (length(args) >= 1) args[[1]] else file.path("results", "fig5a_spatial_lineages", "figures", "fig5a_lineage_methods.tsv")
a_path <- if (length(args) >= 2) args[[2]] else file.path("results", "fig5a_spatial_lineages", "figures", "fig5a_lineage_ranks.tsv")
b_path <- if (length(args) >= 3) args[[3]] else file.path("results", "fig5b_spatial_lineages", "figures", "fig5b_lineage_overlap.tsv")
c_path <- if (length(args) >= 4) args[[4]] else file.path("results", "fig5c_spatial_lineages", "figures", "fig5c_lineage_marker_overlap.tsv")
output_dir <- if (length(args) >= 5) args[[5]] else file.path("results", "fig5_spatial_lineages", "figures")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

types_palette <- c(
    "Overall" = "#f781bf",
    "Integration" = "#e41a1c",
    "Clustering" = "#377eb8",
    "Alignment" = "#4daf4a"
)

methods_tbl <- read_tsv(methods_path, show_col_types = FALSE)
methods_tbl <- methods_tbl |>
    mutate(DisplayName = gsub(" \\(N=[0-9]+\\)$", "", .data$Name))
method_levels <- rev(methods_tbl$DisplayName)

a_df <- read_tsv(a_path, show_col_types = FALSE) |>
    mutate(
        MethodLabel = gsub(" \\(N=[0-9]+\\)$", "", .data$MethodLabel),
        MethodLabel = factor(.data$MethodLabel, levels = method_levels),
        Dataset = factor(.data$Dataset, levels = c("Full", "Immune", "Epithelial")),
        Type = factor(.data$Type, levels = c("Overall", "Integration", "Clustering", "Alignment"))
    )

b_df <- read_tsv(b_path, show_col_types = FALSE) |>
    mutate(
        Method = gsub(" \\(N=[0-9]+\\)$", "", .data$Method),
        Method = factor(.data$Method, levels = method_levels),
        Combination = factor(.data$Combination, levels = c("Full - Immune", "Full - Epithelial", "Immune - Epithelial"))
    )

c_df <- read_tsv(c_path, show_col_types = FALSE) |>
    mutate(
        Method = gsub(" \\(N=[0-9]+\\)$", "", .data$Method),
        Method = factor(.data$Method, levels = method_levels),
        Dataset = factor(.data$Dataset, levels = c("Full", "Immune", "Epithelial")),
        Lineage = factor(.data$Lineage, levels = c("Epithelial", "Immune"))
    )

plot_a <- ggplot(a_df, aes(x = .data$Dataset, y = .data$MethodLabel)) +
    geom_tile(aes(fill = .data$Type, alpha = .data$Rank), colour = "white", linewidth = 0.35) +
    facet_wrap(~ .data$Type, nrow = 1) +
    scale_fill_manual(values = types_palette, guide = "none") +
    scale_alpha_continuous(
        limits = c(max(a_df$Rank, na.rm = TRUE), 1),
        range = c(0.22, 1.0),
        trans = "reverse",
        breaks = c(5, 10, 15, 20, 25),
        name = "Rank",
        guide = guide_legend(order = 1)
    ) +
    labs(x = NULL, y = NULL) +
    theme_features_pub() +
    theme(
        strip.background = element_rect(fill = "black", colour = "black"),
        strip.text = element_text(colour = "white", face = "bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "bold", colour = "black"),
        axis.text.y = element_text(face = "bold", colour = "black", size = 5.8),
        panel.grid = element_blank(),
        plot.margin = margin(0.08, 0.08, 0.06, 0.08, "cm")
    )

plot_b <- ggplot(b_df, aes(x = .data$Combination, y = .data$Method, fill = .data$Jaccard)) +
    geom_tile(colour = "white", linewidth = 0.35) +
    scale_fill_viridis_c(
        option = "C",
        limits = c(0, 1),
        name = "Jaccard\nindex",
        guide = guide_colourbar(order = 2)
    ) +
    labs(x = NULL, y = NULL) +
    theme_features_pub() +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "bold", colour = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        plot.margin = margin(0.08, 0.04, 0.06, 0.02, "cm")
    )

plot_c <- ggplot(c_df, aes(x = .data$Dataset, y = .data$Method)) +
    geom_point(
        aes(colour = .data$MeanProp, size = .data$SDProp),
        shape = 15, stroke = 0,
        na.rm = TRUE
    ) +
    facet_wrap(~ .data$Lineage, nrow = 1) +
    scale_colour_viridis_c(
        option = "C",
        limits = c(0, 1),
        name = "Mean proportion",
        guide = guide_colourbar(
            order = 4,
            theme = theme(legend.margin = margin(0, 0, 0, 8))
        )
    ) +
    scale_size_continuous(
        trans = "reverse",
        range = c(0.7, 4.0),
        name = "s.d.",
        guide = guide_legend(order = 3)
    ) +
    labs(x = NULL, y = NULL) +
    theme_features_pub() +
    theme(
        strip.background = element_rect(fill = "black", colour = "black"),
        strip.text = element_text(colour = "white", face = "bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "bold", colour = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        plot.margin = margin(0.08, 0.08, 0.06, 0.02, "cm")
    )
top_fig <- wrap_plots(
    plot_a, plot_b, plot_c,
    nrow = 1,
    widths = c(1, 0.22, 0.56),
    guides = "collect"
) &
    theme(
        legend.position = "bottom",
        legend.title.position = "top",
        legend.box = "horizontal",
        legend.box.just = "left",
        legend.margin = margin(0, 0, 0, 0, "cm"),
        legend.key.height = unit(0.22, "cm"),
        legend.key.width = unit(0.42, "cm"),
        legend.spacing.x = unit(0.10, "cm"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7, face = "bold")
    )

save_figure_files(top_fig, file.path(output_dir, "figure_fig5abc_spatial_lineages"), width = 7.0, height = 5.75)
