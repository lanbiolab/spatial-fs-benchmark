#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) >= 1) args[[1]] else file.path("results", "metric_overview", "figures")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cols <- c(
    overall = "#f06aa7",
    integration = "#6a3d9a",
    clustering = "#e84a8a",
    alignment = "#f0b429"
)

p <- ggplot() +
    annotate("text", x = 1.2, y = 2.55, label = "Overall", colour = cols["overall"],
             fontface = "italic", family = "sans", size = 11.5, hjust = 0) +
    annotate("text", x = 3.75, y = 2.55, label = "=", colour = "black",
             family = "sans", size = 9.8) +
    annotate("text", x = 5.20, y = 2.95, label = "1", colour = "black",
             family = "sans", size = 8.8) +
    annotate("segment", x = 4.82, xend = 5.58, y = 2.63, yend = 2.63,
             linewidth = 0.65, colour = "black") +
    annotate("text", x = 5.20, y = 2.17, label = "3", colour = "black",
             family = "sans", size = 8.8) +
    annotate("text", x = 6.08, y = 2.55, label = "×", colour = "black",
             family = "sans", size = 8.6) +
    annotate("text", x = 6.72, y = 2.55, label = "(", colour = "black",
             family = "sans", size = 10.5) +
    annotate("text", x = 8.35, y = 2.55, label = "Integration", colour = cols["integration"],
             fontface = "italic", family = "sans", size = 9.6) +
    annotate("segment", x = 7.12, xend = 9.58, y = 2.28, yend = 2.28,
             linewidth = 0.7, colour = "black") +
    annotate("text", x = 10.18, y = 2.55, label = "+", colour = "black",
             family = "sans", size = 8.2) +
    annotate("text", x = 12.10, y = 2.55, label = "Clustering", colour = cols["clustering"],
             fontface = "italic", family = "sans", size = 9.6) +
    annotate("segment", x = 10.78, xend = 13.38, y = 2.28, yend = 2.28,
             linewidth = 0.7, colour = "black") +
    annotate("text", x = 14.18, y = 2.55, label = "+", colour = "black",
             family = "sans", size = 8.2) +
    annotate("text", x = 16.20, y = 2.55, label = "Alignment", colour = cols["alignment"],
             fontface = "italic", family = "sans", size = 9.6) +
    annotate("segment", x = 14.78, xend = 17.55, y = 2.28, yend = 2.28,
             linewidth = 0.7, colour = "black") +
    annotate("text", x = 18.10, y = 2.55, label = ")", colour = "black",
             family = "sans", size = 10.5) +
    annotate("segment", x = 6.92, xend = 17.92, y = 1.72, yend = 1.72,
             linewidth = 0.5, colour = "black") +
    annotate("segment", x = 6.92, xend = 6.92, y = 1.72, yend = 2.00,
             linewidth = 0.5, colour = "black") +
    annotate("segment", x = 17.92, xend = 17.92, y = 1.72, yend = 2.00,
             linewidth = 0.5, colour = "black") +
    annotate("segment", x = 12.42, xend = 12.42, y = 1.72, yend = 1.46,
             linewidth = 0.5, colour = "black") +
    annotate("text", x = 12.42, y = 0.88, label = "Task means", colour = "black",
             family = "sans", size = 8.8) +
    coord_cartesian(xlim = c(0.4, 18.7), ylim = c(0.45, 3.2), clip = "off") +
    theme_void() +
    theme(
        plot.margin = margin(8, 10, 8, 10),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA)
    )

png_path <- file.path(output_dir, "figure_overall_formula.png")
pdf_path <- file.path(output_dir, "figure_overall_formula.pdf")

ggsave(png_path, p, width = 14.5, height = 2.4, dpi = 600, bg = "white")
ggsave(pdf_path, p, width = 14.5, height = 2.4, device = grDevices::cairo_pdf, bg = "white")

cat("Saved:", png_path, "\n")
