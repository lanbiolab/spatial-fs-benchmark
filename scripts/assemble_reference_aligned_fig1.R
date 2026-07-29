#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ggplot2)
    library(ggplotify)
    library(patchwork)
    library(png)
})

source(file.path("external", "atlas-feature-selection-benchmark", "analysis", "R", "plotting.R"))

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) args[[1]] else file.path("results", "reference_aligned_v2")
output_dir <- if (length(args) >= 2) args[[2]] else file.path(root, "figure1")
workflow_pdf <- if (length(args) >= 3) args[[3]] else file.path(
    "manuscript_genome_research_special_issue", "Benchmark_fig_paper", "fig1.pdf"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(workflow_pdf)) {
    stop("Existing workflow figure not found: ", workflow_pdf)
}
if (Sys.which("pdftocairo") == "") {
    stop("pdftocairo is required to reuse the existing vector workflow panel")
}

# The manuscript Figure 1 contains the approved workflow above the old metric
# overview. Render only that workflow region, leaving its geometry untouched.
workflow_stem <- file.path(output_dir, "workflow_panel_a_source")
status <- system2(
    "pdftocairo",
    c(
        "-png", "-singlefile", "-r", "600",
        "-x", "300", "-y", "0", "-W", "7394", "-H", "3100",
        shQuote(workflow_pdf), shQuote(workflow_stem)
    )
)
if (!identical(status, 0L)) {
    stop("Failed to render the approved workflow panel from ", workflow_pdf)
}

workflow_png <- paste0(workflow_stem, ".png")
workflow_raster <- png::readPNG(workflow_png, native = TRUE)
panel_a <- ggplotify::as.ggplot(
    grid::rasterGrob(
        workflow_raster,
        width = grid::unit(1, "npc"),
        height = grid::unit(1, "npc"),
        interpolate = TRUE
    )
) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))

panel_b <- readRDS(file.path(root, "metric_overview", "figure_metric_overview.rds"))
panel_b_tagged <- ggplotify::as.ggplot(panel_b) +
    labs(tag = "b") +
    theme(
        plot.tag = element_text(face = "bold", size = 10, colour = "black"),
        plot.tag.position = c(0.005, 0.995),
        plot.margin = margin(0, 0, 0, 0)
    )

figure <- wrap_plots(
    wrap_elements(full = panel_a),
    wrap_elements(full = panel_b_tagged),
    ncol = 1,
    heights = c(0.38, 0.62)
) &
    theme(
        plot.margin = margin(1, 1, 1, 1)
    )

saveRDS(figure, file.path(output_dir, "figure1_complete_reference_aligned.rds"))
ggsave(
    file.path(output_dir, "figure1_complete_reference_aligned.png"),
    figure,
    width = 11.2,
    height = 11.6,
    dpi = 300,
    device = ragg::agg_png,
    bg = "white"
)
register_arial_pdf_font()
ggsave(
    file.path(output_dir, "figure1_complete_reference_aligned.pdf"),
    figure,
    width = 11.2,
    height = 11.6,
    device = grDevices::pdf,
    family = "Arial",
    useDingbats = FALSE,
    bg = "white"
)

cat("Saved assembled Figure 1 to", output_dir, "\n")
