#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ggplot2)
    library(ggplotify)
    library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
panel_a_rds <- if (length(args) >= 1) args[[1]] else file.path(
    "results", "figure1_reference_dag_v5", "Figure1_reference_dag_v5.rds"
)
panel_b_rds <- if (length(args) >= 2) args[[2]] else file.path(
    "results", "reference_aligned_v2", "metric_overview", "figure_metric_overview.rds"
)
output_dir <- if (length(args) >= 3) args[[3]] else file.path(
    "results", "figure1_reference_dag_v5", "complete"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(panel_a_rds)) stop("Missing panel a RDS: ", panel_a_rds)
if (!file.exists(panel_b_rds)) stop("Missing panel b RDS: ", panel_b_rds)

panel_a <- readRDS(panel_a_rds)
panel_b <- ggplotify::as.ggplot(readRDS(panel_b_rds)) +
    labs(tag = "b") +
    theme(
        plot.tag = element_text(family = "Arial", face = "bold", size = 8.5, colour = "#26333F"),
        plot.tag.position = c(0.003, 0.998),
        plot.margin = margin(0, 0, 0, 0)
    )

figure <- wrap_plots(
    wrap_elements(full = panel_a),
    wrap_elements(full = panel_b),
    ncol = 1,
    heights = c(0.35, 0.65)
) & theme(plot.margin = margin(1, 1, 1, 1, "mm"))

stem <- file.path(output_dir, "Figure1")
width_mm <- 183
height_mm <- 185

ggsave(paste0(stem, ".svg"), figure, width = width_mm, height = height_mm, units = "mm",
       device = svglite::svglite, bg = "white")
svg <- readLines(paste0(stem, ".svg"), warn = FALSE)
svg <- gsub(" textLength='[^']+' lengthAdjust='spacingAndGlyphs'", "", svg)
writeLines(svg, paste0(stem, ".svg"), useBytes = TRUE)

ggsave(paste0(stem, ".pdf"), figure, width = width_mm, height = height_mm, units = "mm",
       device = grDevices::cairo_pdf, family = "Arial", bg = "white")
ggsave(paste0(stem, ".tiff"), figure, width = width_mm, height = height_mm, units = "mm",
       dpi = 600, device = ragg::agg_tiff, compression = "lzw", bg = "white")
ggsave(paste0(stem, ".png"), figure, width = width_mm, height = height_mm, units = "mm",
       dpi = 300, device = ragg::agg_png, bg = "white")

message("Complete reference-aligned Figure 1 written to ", output_dir)
