#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(patchwork)
    library(readr)
    library(scales)
    library(tidyr)
})

source("external/atlas-feature-selection-benchmark/analysis/R/plotting.R")

args <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(args) >= 1) args[[1]] else "results/submission_validations_v1"
output_dir <- if (length(args) >= 2) args[[2]] else file.path(input_dir, "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

method_labels <- c(
    "random" = "Random",
    "scanpy_seurat_v3_batch" = "Scanpy Seurat v3, batch",
    "scanpy_cell_ranger_batch" = "Scanpy Cell Ranger, batch",
    "seurat_vst" = "Seurat VST",
    "triku" = "triku",
    "singleCellHaystack" = "singleCellHaystack",
    "scsegindex" = "scSEGIndex",
    "morans_i" = "Moran's I",
    "sparkx" = "SPARK-X",
    "nnsvg" = "nnSVG",
    "spatialde" = "SpatialDE",
    "somde" = "SOMDE",
    "wilcoxon" = "Wilcoxon",
    "hybrid_hvg_svg_union" = "Balanced HVG-SVG union",
    "hybrid_hvg_svg_intersection" = "Exact HVG-SVG intersection"
)

base_theme <- theme_features_pub() +
    theme(
        axis.text = element_text(colour = "black"),
        axis.line = element_line(linewidth = 0.25, colour = "black"),
        axis.ticks = element_line(linewidth = 0.25, colour = "black"),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "black", colour = "black", linewidth = 0.3),
        strip.text = element_text(size = 6.0, face = "bold", colour = "white"),
        plot.tag = element_text(size = 9, face = "bold"),
        plot.tag.position = c(0, 1),
        plot.margin = margin(2, 3, 2, 2)
    )

save_figure <- function(plot, filename, height_mm) {
    width_in <- 183 / 25.4
    height_in <- height_mm / 25.4
    base <- file.path(output_dir, filename)
    svglite::svglite(paste0(base, ".svg"), width = width_in, height = height_in)
    print(plot)
    dev.off()
    register_arial_pdf_font()
    grDevices::pdf(
        paste0(base, ".pdf"), width = width_in, height = height_in,
        family = "Arial", useDingbats = FALSE, bg = "white"
    )
    print(plot)
    dev.off()
    ragg::agg_tiff(
        paste0(base, ".tiff"), width = width_in, height = height_in,
        units = "in", res = 600, compression = "lzw"
    )
    print(plot)
    dev.off()
}

# Supplemental Figure S7: cross-statistic held-out validation.
held_summary <- read_tsv(
    file.path(input_dir, "heldout_cross_reference_summary.tsv"), show_col_types = FALSE
) |>
    filter(.data$method != "all_features") |>
    mutate(Method = recode(.data$method, !!!method_labels))
held_values <- read_tsv(
    file.path(input_dir, "heldout_cross_reference_values.tsv"), show_col_types = FALSE
) |>
    filter(.data$method != "all_features") |>
    mutate(Method = recode(.data$method, !!!method_labels))
held_concordance <- read_tsv(
    file.path(input_dir, "heldout_cross_reference_concordance.tsv"), show_col_types = FALSE
)

method_order <- held_summary |>
    filter(.data$Metric == "MeanReferencePercentile") |>
    group_by(.data$Method) |>
    summarise(Mean = mean(.data$Mean), .groups = "drop") |>
    arrange(.data$Mean) |>
    pull(.data$Method)
held_summary <- held_summary |>
    mutate(Method = factor(.data$Method, levels = method_order))
held_values <- held_values |>
    mutate(Method = factor(.data$Method, levels = method_order))

heatmap_panel <- function(metric, strip_label, legend_title, panel_tag, show_y = TRUE) {
    frame <- filter(held_summary, .data$Metric == metric) |>
        mutate(
            Panel = strip_label,
            TextColour = if_else(.data$Mean < 0.45, "white", "black")
        )
    ggplot(frame, aes(x = .data$Reference, y = .data$Method)) +
        geom_tile(aes(fill = .data$Mean), width = 0.86, height = 0.86, colour = "white", linewidth = 0.25) +
        geom_text(aes(label = sprintf("%.2f", .data$Mean), colour = .data$TextColour), size = 1.65) +
        scale_colour_identity() +
        scale_fill_viridis_c(option = "magma", limits = c(0, 1), name = legend_title) +
        facet_wrap(~ .data$Panel, nrow = 1) +
        labs(tag = panel_tag, x = NULL, y = NULL) +
        base_theme +
        theme(
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.text.y = if (show_y) element_text(size = 5.1) else element_blank(),
            legend.position = "right",
            legend.key.height = unit(0.78, "cm"),
            legend.key.width = unit(0.22, "cm")
        )
}

rank_frame <- held_summary |>
    filter(.data$Metric %in% c("Jaccard", "MeanReferencePercentile", "RankSpearman")) |>
    group_by(.data$Reference, .data$Metric) |>
    mutate(Rank = rank(-.data$Mean, ties.method = "average")) |>
    ungroup() |>
    select(Method, Reference, Metric, Rank) |>
    pivot_wider(names_from = Reference, values_from = Rank)
rho_labels <- held_concordance |>
    filter(.data$Metric %in% c("Jaccard", "MeanReferencePercentile", "RankSpearman")) |>
    transmute(
        Metric = .data$Metric,
        Label = sprintf("rho = %.2f", .data$SpearmanRho),
        x = 2,
        y = 12.5
    )
panel_c <- ggplot(rank_frame, aes(x = .data[["Moran's I"]], y = .data$nnSVG)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.35, colour = "grey65") +
    geom_point(shape = 21, fill = "#f781bf", colour = "white", stroke = 0.3, size = 1.8) +
    geom_text(data = rho_labels, aes(x = .data$x, y = .data$y, label = .data$Label), inherit.aes = FALSE, hjust = 0, size = 2.1) +
    facet_wrap(~Metric, nrow = 1, labeller = as_labeller(c(
        "Jaccard" = "Top-N Jaccard",
        "MeanReferencePercentile" = "Mean percentile",
        "RankSpearman" = "Rank agreement"
    ))) +
    scale_x_continuous(limits = c(1, 13), breaks = c(1, 5, 9, 13)) +
    scale_y_continuous(limits = c(1, 13), breaks = c(1, 5, 9, 13)) +
    labs(tag = "c", x = "Rank using Moran's I reference", y = "Rank using nnSVG reference") +
    base_theme

association <- held_concordance |>
    filter(grepl("vs ARI", .data$Metric)) |>
    mutate(
        Reference = if_else(grepl("nnSVG", .data$Metric), "nnSVG", "Moran's I"),
        Measure = case_when(
            grepl("^Jaccard", .data$Metric) ~ "Top-N Jaccard",
            grepl("^MeanReference", .data$Metric) ~ "Mean reference percentile",
            TRUE ~ "Full-ranking agreement"
        )
    )
panel_d <- ggplot(association, aes(x = .data$SpearmanRho, y = .data$Measure, colour = .data$Reference, shape = .data$Reference)) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey75") +
    geom_errorbar(
        aes(xmin = .data$Lower95, xmax = .data$Upper95),
        orientation = "y", width = 0.16, linewidth = 0.5,
        position = position_dodge(width = 0.38)
    ) +
    geom_point(size = 1.9, position = position_dodge(width = 0.38)) +
    scale_colour_manual(values = c("Moran's I" = "#377eb8", "nnSVG" = "#4daf4a"), name = "Held-out reference") +
    scale_shape_manual(values = c("Moran's I" = 16, "nnSVG" = 17), name = "Held-out reference") +
    scale_x_continuous(limits = c(-0.25, 0.72), breaks = c(-0.2, 0, 0.2, 0.4, 0.6)) +
    labs(tag = "d", x = "Fold-bootstrap Spearman correlation with ARI", y = NULL) +
    base_theme +
    theme(legend.position = "bottom", legend.title.position = "top")

panel_a <- heatmap_panel("MeanReferencePercentile", "Held-out reference percentile", "Mean percentile", "a")
panel_b <- heatmap_panel("Jaccard", "Held-out top-N overlap", "Mean Jaccard", "b", show_y = FALSE)
figure_s7 <- (panel_a | panel_b) / (panel_c | panel_d) +
    plot_layout(heights = c(1.28, 0.82))
save_figure(figure_s7, "Supplemental_Figure_S7", 170)

# Supplemental Figure S8: hybrid controls.
hybrid_comparison <- read_tsv(
    file.path(input_dir, "hybrid_comparison_task_scores.tsv"), show_col_types = FALSE
) |>
    mutate(
        Method = recode(.data$fs_method, !!!method_labels),
        Integrator = recode(.data$integration_method, "scvi" = "scVI", "cellcharter" = "CellCharter")
    )
hybrid_effective <- read_tsv(
    file.path(input_dir, "hybrid_effective_features.tsv"), show_col_types = FALSE
) |>
    mutate(
        Method = recode(.data$fs_method, !!!method_labels),
        Dataset = recode(.data$dataset, "MouseBrainSerialSections" = "Mouse Brain")
    )
hybrid_diff <- read_tsv(
    file.path(input_dir, "hybrid_differences_from_best_parent.tsv"), show_col_types = FALSE
) |>
    mutate(
        Method = recode(.data$fs_method, !!!method_labels),
        Integrator = recode(.data$integration_method, "scvi" = "scVI", "cellcharter" = "CellCharter"),
        Dataset = recode(.data$dataset, "MouseBrainSerialSections" = "Mouse Brain")
    )
hybrid_global <- read_tsv(
    file.path(input_dir, "hybrid_global_rank_summary.tsv"), show_col_types = FALSE
) |>
    filter(.data$fs_method %in% c(
        "hybrid_hvg_svg_union", "hybrid_hvg_svg_intersection",
        "scanpy_seurat_v3_batch", "morans_i", "triku", "somde", "seurat_vst"
    )) |>
    mutate(
        Method = recode(.data$fs_method, !!!method_labels),
        Integrator = recode(.data$integration_method, "scvi" = "scVI", "cellcharter" = "CellCharter")
    )

hybrid_palette <- c(
    "Balanced HVG-SVG union" = "#f781bf",
    "Exact HVG-SVG intersection" = "#377eb8",
    "Scanpy Seurat v3, batch" = "#984ea3",
    "Moran's I" = "#4daf4a",
    "triku" = "#666666",
    "SOMDE" = "#a6cee3",
    "Seurat VST" = "#ff7f00"
)

panel_a8 <- ggplot(hybrid_effective, aes(x = .data$effective_n_features, y = reorder(.data$Dataset, .data$effective_n_features), colour = .data$Method)) +
    geom_vline(xintercept = 2000, linetype = "dashed", linewidth = 0.35, colour = "grey65") +
    geom_point(size = 2.2) +
    scale_colour_manual(
        values = hybrid_palette[c("Balanced HVG-SVG union", "Exact HVG-SVG intersection")],
        name = NULL
    ) +
    scale_x_continuous(limits = c(0, 2100), breaks = c(0, 500, 1000, 1500, 2000)) +
    labs(tag = "a", x = "Effective number of selected genes", y = NULL) +
    base_theme +
    theme(legend.position = "bottom", legend.title.position = "top")

method_order8 <- c(
    "Balanced HVG-SVG union", "Exact HVG-SVG intersection",
    "Moran's I", "Scanpy Seurat v3, batch", "SOMDE", "Seurat VST", "triku"
)
hybrid_comparison <- hybrid_comparison |>
    mutate(Method = factor(.data$Method, levels = rev(method_order8)))
panel_b8 <- ggplot(hybrid_comparison, aes(x = .data$CoreOverallMean, y = .data$Method, colour = .data$Method)) +
    geom_point(alpha = 0.42, size = 1.15, position = position_jitter(height = 0.10, width = 0)) +
    stat_summary(fun = mean, geom = "point", shape = 23, fill = "white", size = 2.4, stroke = 0.5) +
    facet_wrap(~Integrator, nrow = 1) +
    scale_colour_manual(values = hybrid_palette, guide = "none") +
    labs(tag = "b", x = "CoreOverall score", y = NULL) +
    base_theme

panel_c8 <- ggplot(hybrid_diff, aes(x = .data$DifferenceFromBestParent, y = reorder(.data$Dataset, .data$DifferenceFromBestParent), colour = .data$Method)) +
    geom_vline(xintercept = 0, linewidth = 0.35, colour = "grey45") +
    geom_segment(aes(x = 0, xend = .data$DifferenceFromBestParent, yend = reorder(.data$Dataset, .data$DifferenceFromBestParent)), linewidth = 0.5) +
    geom_point(size = 1.8) +
    facet_wrap(~Integrator, nrow = 1) +
    scale_colour_manual(
        values = hybrid_palette[c("Balanced HVG-SVG union", "Exact HVG-SVG intersection")],
        name = NULL,
        guide = "none"
    ) +
    labs(tag = "c", x = "CoreOverall difference from best parent", y = NULL) +
    base_theme +
    theme(legend.position = "none")

hybrid_global <- hybrid_global |>
    mutate(Method = factor(.data$Method, levels = rev(method_order8)))
panel_d8 <- ggplot(hybrid_global, aes(x = .data$GlobalRankWithHybrids, y = .data$Method, colour = .data$Method)) +
    geom_segment(aes(x = 1, xend = .data$GlobalRankWithHybrids, yend = .data$Method), colour = "grey80", linewidth = 0.45) +
    geom_point(size = 2) +
    facet_wrap(~Integrator, nrow = 1) +
    scale_colour_manual(values = hybrid_palette, guide = "none") +
    scale_x_continuous(limits = c(1, 31), breaks = c(1, 5, 10, 20, 30)) +
    labs(tag = "d", x = "Global rank (1 = best)", y = NULL) +
    base_theme

figure_s8 <- (panel_a8 | panel_b8) / (panel_c8 | panel_d8) +
    plot_layout(heights = c(0.95, 1.05), guides = "collect") &
    theme(legend.position = "bottom")
save_figure(figure_s8, "Supplemental_Figure_S8", 158)

write_tsv(held_summary, file.path(output_dir, "Supplemental_Figure_S7a_b_summary.tsv"), na = "NA")
write_tsv(rank_frame, file.path(output_dir, "Supplemental_Figure_S7c_rank_concordance.tsv"), na = "NA")
write_tsv(association, file.path(output_dir, "Supplemental_Figure_S7d_ari_association.tsv"), na = "NA")
write_tsv(hybrid_effective, file.path(output_dir, "Supplemental_Figure_S8a_effective_features.tsv"), na = "NA")
write_tsv(hybrid_comparison, file.path(output_dir, "Supplemental_Figure_S8b_scores.tsv"), na = "NA")
write_tsv(hybrid_diff, file.path(output_dir, "Supplemental_Figure_S8c_parent_differences.tsv"), na = "NA")
write_tsv(hybrid_global, file.path(output_dir, "Supplemental_Figure_S8d_global_ranks.tsv"), na = "NA")
