#!/usr/bin/env Rscript

# Figure 1 combines the workflow framework with the full metric-characterization
# panel. Both panels are drawn at their intended physical size on one device so
# that grid text is never distorted by a post-hoc patchwork rescale.

suppressPackageStartupMessages({
    library(ggplot2)
    library(ggplotify)
    library(grid)
    library(svglite)
    library(ragg)
})

args <- commandArgs(trailingOnly = TRUE)
dag_rds <- if (length(args) >= 1) args[[1]] else file.path(
    "results", "figure1_reference_dag_v5", "Figure1_reference_dag_v5.rds"
)
metric_rds <- if (length(args) >= 2) args[[2]] else file.path(
    "results", "reference_aligned_v2", "metric_overview", "figure_metric_overview.rds"
)
output_dir <- if (length(args) >= 3) args[[3]] else file.path(
    "results", "figure1_complete_v6"
)

if (!file.exists(dag_rds)) stop("Missing workflow panel: ", dag_rds)
if (!file.exists(metric_rds)) stop("Missing metric panel: ", metric_rds)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

replace_font <- function(x, family = "Helvetica") {
    if (inherits(x, "grob")) {
        if (!is.null(x$gp)) x$gp$fontfamily <- family
        if (inherits(x, "text") && length(x$label) == 1) {
            if (identical(x$label, "Mean feature-number correlation") ||
                identical(x$label, "Mean dataset-feature correlation")) {
                x$gp$fontsize <- 5.0
            }
            if (identical(x$label, "Alignment")) x$gp$fontsize <- 5.0
        }
        if (length(x$children)) {
            x$children <- lapply(x$children, replace_font, family = family)
        }
        if (length(x$grobs)) {
            x$grobs <- lapply(x$grobs, replace_font, family = family)
        }
    }
    x
}

dag <- replace_font(ggplotify::as.grob(readRDS(dag_rds)))

# The original metric overview was exported at two-column width.  In Figure 1
# it occupies a single full-width panel, so only its header hierarchy is made
# more compact; all tiles, distributions, legends and data encodings are kept.
metric_plot <- readRDS(metric_rds) &
    theme(
        plot.title = element_text(size = 6.4, lineheight = 0.90,
                                  margin = margin(b = 1.5, unit = "mm")),
        strip.text = element_text(size = 5.6)
    )
metrics <- replace_font(ggplotify::as.grob(metric_plot))

# The lower panel's native aspect ratio is 11.2 / 7.5.  Keeping it intact is
# essential because the small correlation tiles and legends encode data.
width_mm <- 183
height_mm <- 204
margin_mm <- 4
dag_h_mm <- 70
gap_mm <- 5
metric_h_mm <- height_mm - 2 * margin_mm - dag_h_mm - gap_mm

draw_figure <- function() {
    grid.newpage()
    pushViewport(viewport(gp = gpar(fontfamily = "Helvetica")))

    # Panel a: native 183 x 70-mm workflow graph.
    pushViewport(viewport(
        x = unit(0.5, "npc"),
        y = unit(1, "npc") - unit(margin_mm, "mm") - unit(dag_h_mm / 2, "mm"),
        width = unit(1, "npc") - unit(2 * margin_mm, "mm"),
        height = unit(dag_h_mm, "mm"),
        just = c("center", "center"), clip = "off"
    ))
    grid.draw(dag)
    upViewport()
    grid.text("a", x = unit(margin_mm, "mm"),
              y = unit(1, "npc") - unit(margin_mm - 0.6, "mm"),
              just = c("left", "top"),
              gp = gpar(fontfamily = "Helvetica", fontsize = 9, fontface = "bold", col = "#26333F"))

    # Panel b: complete metric overview, kept at its intended physical height.
    pushViewport(viewport(
        x = unit(0.5, "npc"),
        y = unit(margin_mm, "mm") + unit(metric_h_mm / 2, "mm"),
        # Use the full page width for the dense metric matrix so its rightmost
        # Alignment strip and individual legend titles are never clipped.
        width = unit(1, "npc"),
        height = unit(metric_h_mm, "mm"),
        just = c("center", "center"), clip = "off"
    ))
    grid.draw(metrics)
    upViewport()
    grid.text("b", x = unit(margin_mm, "mm"),
              y = unit(margin_mm, "mm") + unit(metric_h_mm, "mm") + unit(1.2, "mm"),
              just = c("left", "bottom"),
              gp = gpar(fontfamily = "Helvetica", fontsize = 9, fontface = "bold", col = "#26333F"))
    upViewport()
}

stem <- file.path(output_dir, "Figure1")
inches <- c(width_mm, height_mm) / 25.4

svglite::svglite(paste0(stem, ".svg"), width = inches[[1]], height = inches[[2]], bg = "white")
draw_figure()
dev.off()

# Remove SVG text-length constraints, which can cause character-spacing artefacts
# in some TeX/PDF preview paths.
svg <- readLines(paste0(stem, ".svg"), warn = FALSE)
svg <- gsub(" textLength='[^']+' lengthAdjust='spacingAndGlyphs'", "", svg)
writeLines(svg, paste0(stem, ".svg"), useBytes = TRUE)

grDevices::pdf(paste0(stem, ".pdf"), width = inches[[1]], height = inches[[2]],
               family = "Helvetica", useDingbats = FALSE, bg = "white")
draw_figure()
dev.off()

ragg::agg_png(paste0(stem, ".png"), width = width_mm, height = height_mm,
              units = "mm", res = 300, background = "white")
draw_figure()
dev.off()

ragg::agg_tiff(paste0(stem, ".tiff"), width = width_mm, height = height_mm,
               units = "mm", res = 600, compression = "lzw", background = "white")
draw_figure()
dev.off()

ragg::agg_png(file.path(output_dir, "Figure1_qa.png"), width = width_mm,
              height = height_mm, units = "mm", res = 180, background = "white")
draw_figure()
dev.off()

saveRDS(list(dag = dag, metrics = metrics, dimensions_mm = c(width_mm, height_mm)),
        file.path(output_dir, "Figure1_panels.rds"))

message("Complete Figure 1 written to ", output_dir)
