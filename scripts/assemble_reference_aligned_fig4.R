#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ggplot2)
    library(patchwork)
})

source(file.path("external", "atlas-feature-selection-benchmark", "analysis", "R", "plotting.R"))

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) args[[1]] else file.path("results", "reference_aligned_v2")
output_dir <- if (length(args) >= 2) args[[2]] else file.path(root, "figure4")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

panel_a <- readRDS(file.path(root, "fig4a", "figure_fig4a_spatial_benchmark.rds"))
panel_b <- readRDS(file.path(root, "fig4bcd", "figure_fig4b_panel_b.rds"))
panel_c <- readRDS(file.path(root, "fig4bcd", "figure_fig4c_panel_c.rds"))
panel_d <- readRDS(file.path(root, "fig4bcd", "figure_fig4d_panel_d.rds"))
panel_e <- readRDS(file.path(root, "fig4e_hybrid", "figure_fig4e_hybrid_panel.rds"))

panel_a <- wrap_elements(full = panel_a)
panel_b <- wrap_elements(full = panel_b)
panel_c <- wrap_elements(full = panel_c)
panel_d <- wrap_elements(full = panel_d)
panel_e <- wrap_elements(full = panel_e)

middle <- wrap_plots(panel_b, panel_c, panel_d, nrow = 1, widths = c(1.20, 0.78, 0.92))

figure <- wrap_plots(
    panel_a,
    middle,
    panel_e,
    ncol = 1,
    heights = c(1.65, 1.0, 0.68)
) +
    plot_annotation(tag_levels = "a") &
    theme(
        plot.tag = element_text(face = "bold", size = 10),
        plot.margin = margin(2, 2, 2, 2)
    )

saveRDS(figure, file.path(output_dir, "figure4_complete_reference_aligned.rds"))
save_figure_files(
    figure,
    file.path(output_dir, "figure4_complete_reference_aligned"),
    width = 8.2,
    height = 10.2
)
ggsave(
    file.path(output_dir, "figure4_complete_reference_aligned.svg"), figure,
    width = 8.2, height = 10.2,
    device = svglite::svglite, bg = "white"
)
ggsave(
    file.path(output_dir, "figure4_complete_reference_aligned.tiff"), figure,
    width = 8.2, height = 10.2, dpi = 600,
    device = ragg::agg_tiff, compression = "lzw", bg = "white"
)
