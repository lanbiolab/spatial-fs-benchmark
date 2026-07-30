#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(patchwork)
    library(readr)
})

source("external/atlas-feature-selection-benchmark/analysis/R/plotting.R")

args <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(args) >= 1) args[[1]] else "results/submission_robustness_v1"
output_dir <- if (length(args) >= 2) args[[2]] else "results/submission_robustness_v1/figure"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

method_labels <- c(
    "scanpy_seurat_v3" = "Scanpy Seurat v3",
    "scanpy_seurat_v3_batch" = "Scanpy Seurat v3, batch",
    "scanpy_cell_ranger" = "Scanpy Cell Ranger",
    "scanpy_cell_ranger_batch" = "Scanpy Cell Ranger, batch",
    "scanpy_seurat" = "Scanpy Seurat",
    "scanpy_seurat_batch" = "Scanpy Seurat, batch",
    "scanpy_pearson" = "Scanpy Pearson",
    "scanpy_pearson_batch" = "Scanpy Pearson, batch",
    "seurat_vst" = "Seurat VST",
    "seurat_mvp" = "Seurat MVP",
    "seurat_disp" = "Seurat Dispersion",
    "seurat_sct" = "Seurat sctransform",
    "singleCellHaystack" = "singleCellHaystack",
    "statistic_mean" = "Mean expression",
    "statistic_variance" = "Expression variance",
    "morans_i" = "Moran's I",
    "spatialde" = "SpatialDE",
    "sparkx" = "SPARK-X",
    "nnsvg" = "nnSVG",
    "somde" = "SOMDE",
    "Brennecke" = "Brennecke",
    "scPNMF" = "scPNMF",
    "nbumi" = "NBumi",
    "dubstepr" = "DUBStepR",
    "osca" = "OSCA",
    "scry" = "scry",
    "anticor" = "Anticor",
    "hotspot" = "Hotspot",
    "triku" = "triku"
)
scenario_labels <- c(
    "observed_range" = "Observed\nrange",
    "rank_percentile" = "Percentile\nscaling",
    "z_score" = "z-score\nscaling",
    "integration_weight_0.4" = "Integration\nweight 0.4",
    "integration_weight_0.6" = "Integration\nweight 0.6",
    "labelled_datasets_only" = "Annotated\ndatasets"
)
scenario_order <- names(scenario_labels)

ranks <- read_tsv(file.path(input_dir, "scaling_weight_sensitivity_ranks.tsv"), show_col_types = FALSE) |>
    mutate(
        Method = recode(.data$fs_method, !!!method_labels),
        Scenario = factor(.data$Scenario, levels = scenario_order, labels = unname(scenario_labels))
    )
bootstrap <- read_tsv(file.path(input_dir, "dataset_bootstrap_rank_intervals.tsv"), show_col_types = FALSE) |>
    mutate(Method = recode(.data$fs_method, !!!method_labels))
concordance <- read_tsv(file.path(input_dir, "scaling_weight_concordance.tsv"), show_col_types = FALSE) |>
    filter(.data$NScenarios == 6)

scvi_order <- ranks |>
    filter(.data$integration_method == "scvi", .data$Scenario == scenario_labels[["observed_range"]]) |>
    arrange(desc(.data$GlobalRank)) |>
    pull(.data$Method)
ranks <- ranks |> mutate(Method = factor(.data$Method, levels = scvi_order))

base_theme <- theme_features_pub() +
    theme(
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_text(colour = "black"),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "black", colour = "black", linewidth = 0.3),
        strip.text = element_text(colour = "white", face = "bold", size = 6.0),
        plot.subtitle = element_text(size = 5.8, colour = "#4A4A4A", margin = margin(b = 1.5)),
        plot.tag = element_text(size = 9, face = "bold"),
        plot.tag.position = c(0, 1),
        plot.margin = margin(2, 3, 2, 2)
    )

heatmap_panel <- function(integrator, show_y, panel_tag) {
    w <- concordance |> filter(.data$integration_method == integrator) |> pull(.data$KendallsW)
    frame <- ranks |>
        filter(.data$integration_method == integrator) |>
        mutate(Integrator = if_else(integrator == "scvi", "scVI", "CellCharter"))
    ggplot(frame, aes(x = .data$Scenario, y = .data$Method)) +
        geom_tile(aes(fill = .data$GlobalRank), width = 0.90, height = 0.88, colour = "white", linewidth = 0.16) +
        scale_fill_viridis_c(
            option = "magma",
            direction = -1,
            limits = c(1, 29),
            breaks = c(1, 10, 20, 29),
            name = "Global rank\n(1 = best)"
        ) +
        facet_wrap(~ .data$Integrator, nrow = 1) +
        scale_x_discrete(position = "top") +
        labs(tag = panel_tag, subtitle = sprintf("Kendall's W = %.3f", w)) +
        base_theme +
        theme(
            axis.text.x = element_text(size = 5.5, lineheight = 0.9, angle = 0, hjust = 0.5),
            axis.text.y = if (show_y) element_text(size = 5.5) else element_blank(),
            axis.ticks.y = if (show_y) element_line(linewidth = 0.25) else element_blank(),
            legend.position = if (integrator == "cellcharter") "right" else "none",
            legend.key.height = unit(0.78, "cm"),
            legend.key.width = unit(0.22, "cm")
        )
}

interval_panel <- function(integrator, panel_tag) {
    top_methods <- ranks |>
        filter(.data$integration_method == integrator, .data$Scenario == scenario_labels[["observed_range"]]) |>
        arrange(.data$GlobalRank) |>
        slice_head(n = 10) |>
        pull(.data$Method)
    frame <- bootstrap |>
        filter(.data$integration_method == integrator, .data$Method %in% top_methods) |>
        mutate(
            Method = factor(.data$Method, levels = rev(top_methods)),
            Integrator = if_else(integrator == "scvi", "scVI", "CellCharter")
        )
    ggplot(frame, aes(y = .data$Method)) +
        geom_segment(
            aes(x = .data$BootstrapGlobalRankLower95, xend = .data$BootstrapGlobalRankUpper95, yend = .data$Method),
            linewidth = 0.45,
            colour = "grey70",
            lineend = "round"
        ) +
        geom_segment(
            aes(x = .data$BootstrapGlobalRankLowerQuartile, xend = .data$BootstrapGlobalRankUpperQuartile, yend = .data$Method),
            linewidth = 1.6,
            colour = "#f781bf",
            lineend = "round"
        ) +
        geom_point(aes(x = .data$ObservedGlobalRank), shape = 23, size = 1.9, stroke = 0.35, fill = "white", colour = "black") +
        facet_wrap(~ .data$Integrator, nrow = 1) +
        scale_x_continuous(breaks = c(1, 5, 10, 15, 20, 25, 29), limits = c(0.5, 29.5)) +
        labs(tag = panel_tag, x = "Global rank (1 = best)", y = NULL) +
        base_theme +
        theme(
            axis.line.x = element_line(linewidth = 0.3),
            axis.ticks.x = element_line(linewidth = 0.3),
            axis.text.y = element_text(size = 5.7),
            axis.title.x = element_text(size = 6.2, margin = margin(t = 3))
        )
}

panel_a <- heatmap_panel("scvi", TRUE, "a")
panel_b <- heatmap_panel("cellcharter", FALSE, "b")
panel_c <- interval_panel("scvi", "c")
panel_d <- interval_panel("cellcharter", "d")

figure <- (panel_a | panel_b) / (panel_c | panel_d) +
    plot_layout(heights = c(1.52, 0.82))

width_mm <- 183
height_mm <- 190
width_in <- width_mm / 25.4
height_in <- height_mm / 25.4
base <- file.path(output_dir, "Supplemental_Figure_S6")

svglite::svglite(paste0(base, ".svg"), width = width_in, height = height_in)
print(figure)
dev.off()
register_arial_pdf_font()
grDevices::pdf(
    paste0(base, ".pdf"), width = width_in, height = height_in,
    family = "Arial", useDingbats = FALSE, bg = "white"
)
print(figure)
dev.off()
ragg::agg_tiff(paste0(base, ".tiff"), width = width_in, height = height_in, units = "in", res = 600, compression = "lzw")
print(figure)
dev.off()

write_tsv(ranks, file.path(output_dir, "Supplemental_Figure_S6a_b_rank_sensitivity.tsv"), na = "NA")
write_tsv(bootstrap, file.path(output_dir, "Supplemental_Figure_S6c_d_bootstrap_intervals.tsv"), na = "NA")
