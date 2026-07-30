#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(patchwork)
    library(readr)
})

source(file.path("external", "atlas-feature-selection-benchmark", "analysis", "R", "plotting.R"))

args <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(args) >= 1) args[[1]] else "results/submission_validations_v1"
output_dir <- if (length(args) >= 2) args[[2]] else "results/reference_aligned_v2/fig4e_hybrid"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

method_labels <- c(
    "hybrid_hvg_svg_union" = "Balanced HVG-SVG union",
    "hybrid_hvg_svg_intersection" = "Exact HVG-SVG intersection"
)
method_order <- unname(method_labels)
hybrid_palette <- c(
    "Balanced HVG-SVG union" = "#f781bf",
    "Exact HVG-SVG intersection" = "#377eb8"
)
dataset_labels <- c(
    "MouseBrainSerialSections" = "Mouse Brain",
    "STOmics0212" = "STOmics-0212",
    "STOmics0218" = "STOmics-0218",
    "STOmics0224" = "STOmics-0224",
    "E8p5Embryo" = "E8.5 Embryo",
    "E9p5Embryo" = "E9.5 Embryo",
    "DLPFC" = "DLPFC"
)

base_theme <- theme_features_pub() +
    theme(
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.30),
        strip.background = element_rect(fill = "black", colour = "black", linewidth = 0.30),
        strip.text = element_text(colour = "white", face = "bold", size = 5.8),
        axis.text = element_text(colour = "black"),
        axis.title = element_text(size = 5.8),
        plot.margin = margin(0.5, 1.0, 0.5, 0.5, "mm")
    )

effective <- read_tsv(
    file.path(input_dir, "hybrid_effective_features.tsv"), show_col_types = FALSE
) |>
    mutate(
        Method = factor(recode(.data$fs_method, !!!method_labels), levels = method_order),
        Dataset = recode(.data$dataset, !!!dataset_labels)
    )
dataset_order <- effective |>
    filter(.data$Method == "Exact HVG-SVG intersection") |>
    arrange(.data$effective_n_features) |>
    pull(.data$Dataset)
effective <- effective |>
    mutate(Dataset = factor(.data$Dataset, levels = dataset_order))

panel_size <- ggplot(
    effective,
    aes(x = .data$effective_n_features, y = .data$Dataset, colour = .data$Method)
) +
    geom_vline(xintercept = 2000, linetype = "dashed", linewidth = 0.30, colour = "grey60") +
    geom_point(size = 1.8) +
    facet_wrap(~ factor("Effective feature count"), nrow = 1) +
    scale_colour_manual(values = hybrid_palette, name = "Hybrid construction") +
    scale_x_continuous(limits = c(0, 2100), breaks = c(0, 1000, 2000)) +
    labs(x = "Number of selected genes", y = NULL) +
    base_theme +
    theme(
        axis.text.y = element_text(size = 4.65),
        legend.position = "bottom",
        legend.title.position = "top",
        legend.title = element_text(face = "bold"),
        legend.key.size = unit(0.26, "cm")
    )

differences <- read_tsv(
    file.path(input_dir, "hybrid_differences_from_best_parent.tsv"), show_col_types = FALSE
) |>
    mutate(
        Method = factor(recode(.data$fs_method, !!!method_labels), levels = method_order),
        Integrator = factor(
            recode(.data$integration_method, "cellcharter" = "CellCharter", "scvi" = "scVI"),
            levels = c("CellCharter", "scVI")
        ),
        Dataset = factor(recode(.data$dataset, !!!dataset_labels), levels = dataset_order)
    )

panel_difference <- ggplot(
    differences,
    aes(
        x = .data$DifferenceFromBestParent, y = .data$Dataset,
        colour = .data$Method, group = .data$Method
    )
) +
    geom_vline(xintercept = 0, linewidth = 0.30, colour = "grey45") +
    geom_segment(
        aes(x = 0, xend = .data$DifferenceFromBestParent, yend = .data$Dataset),
        linewidth = 0.40, alpha = 0.72,
        position = position_dodge(width = 0.34)
    ) +
    geom_point(size = 1.55, position = position_dodge(width = 0.34)) +
    facet_wrap(~ .data$Integrator, nrow = 1) +
    scale_colour_manual(values = hybrid_palette, guide = "none") +
    scale_x_continuous(limits = c(-0.17, 0.12), breaks = c(-0.15, 0, 0.1)) +
    labs(x = "CoreOverall difference from better parent", y = NULL) +
    base_theme +
    theme(axis.text.y = element_text(size = 4.55))

global_ranks <- read_tsv(
    file.path(input_dir, "hybrid_global_rank_summary.tsv"), show_col_types = FALSE
) |>
    filter(.data$fs_method %in% names(method_labels)) |>
    mutate(
        Method = factor(recode(.data$fs_method, !!!method_labels), levels = rev(method_order)),
        Integrator = factor(
            recode(.data$integration_method, "cellcharter" = "CellCharter", "scvi" = "scVI"),
            levels = c("CellCharter", "scVI")
        )
    )

panel_rank <- ggplot(
    global_ranks,
    aes(x = .data$GlobalRankWithHybrids, y = .data$Method, colour = .data$Method)
) +
    geom_segment(
        aes(x = 1, xend = .data$GlobalRankWithHybrids, yend = .data$Method),
        colour = "grey80", linewidth = 0.42
    ) +
    geom_point(size = 1.9) +
    geom_text(
        aes(label = sprintf("%.1f", .data$GlobalRankWithHybrids)),
        nudge_x = 1.1, hjust = 0, size = 1.75, colour = "black"
    ) +
    facet_wrap(~ .data$Integrator, nrow = 1) +
    scale_colour_manual(values = hybrid_palette, guide = "none") +
    scale_x_continuous(limits = c(1, 31), breaks = c(1, 10, 20, 30)) +
    labs(x = "Global rank (1 = best)", y = NULL) +
    base_theme +
    theme(axis.text.y = element_text(size = 4.45))

figure <- panel_size + panel_difference + panel_rank +
    plot_layout(widths = c(0.28, 0.44, 0.28), guides = "collect") &
    theme(
        legend.position = "bottom",
        legend.title.position = "top",
        legend.box = "horizontal",
        legend.box.just = "left",
        legend.text = element_text(size = 5.0),
        legend.title = element_text(size = 5.6, face = "bold")
    )

saveRDS(figure, file.path(output_dir, "figure_fig4e_hybrid_panel.rds"))
ggsave(
    file.path(output_dir, "figure_fig4e_hybrid_panel.png"), figure,
    width = 8.2, height = 2.45, dpi = 300,
    device = ragg::agg_png, bg = "white"
)
register_arial_pdf_font()
ggsave(
    file.path(output_dir, "figure_fig4e_hybrid_panel.pdf"), figure,
    width = 8.2, height = 2.45,
    device = grDevices::pdf, family = "Arial", useDingbats = FALSE, bg = "white"
)

write_tsv(effective, file.path(output_dir, "fig4e_hybrid_effective_features.tsv"), na = "NA")
write_tsv(differences, file.path(output_dir, "fig4e_hybrid_parent_differences.tsv"), na = "NA")
write_tsv(global_ranks, file.path(output_dir, "fig4e_hybrid_global_ranks.tsv"), na = "NA")

message("Hybrid-control panel written to ", output_dir)
