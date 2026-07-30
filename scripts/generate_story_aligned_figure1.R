#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ggplot2)
    library(ggplotify)
    library(grid)
    library(patchwork)
})

source(file.path("external", "atlas-feature-selection-benchmark", "analysis", "R", "plotting.R"))

args <- commandArgs(trailingOnly = TRUE)
metric_rds <- if (length(args) >= 1) {
    args[[1]]
} else {
    file.path("results", "metric_overview", "figures", "figure_metric_overview.rds")
}
output_dir <- if (length(args) >= 2) {
    args[[2]]
} else {
    file.path("results", "main_figures", "figure1_story_aligned")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(metric_rds)) {
    stop("Metric-overview RDS not found: ", metric_rds)
}

palette <- c(
    navy = "#2F6DB0",
    teal = "#238B83",
    coral = "#D95F4B",
    indigo = "#6857A6",
    gold = "#D9A62E",
    green = "#4E9A65",
    ink = "#27313D",
    mid = "#667085",
    line = "#CAD2DD",
    pale_navy = "#EEF4FB",
    pale_teal = "#ECF7F5",
    pale_coral = "#FCEFEB",
    pale_indigo = "#F2EFF9",
    pale_gold = "#FCF6E7",
    pale_gray = "#F5F7FA"
)

grobs <- list()
diagram_text_scale <- 1.50
add_grob <- function(g) {
    grobs[[length(grobs) + 1]] <<- g
}

rr <- function(x, y, w, h, fill = "white", col = palette[["line"]], lwd = 0.8, radius = 1.8) {
    roundrectGrob(
        x = unit(x, "npc"), y = unit(y, "npc"),
        width = unit(w, "npc"), height = unit(h, "npc"),
        r = unit(radius, "mm"),
        gp = gpar(fill = fill, col = col, lwd = lwd)
    )
}

txt <- function(label, x, y, size = 6, col = palette[["ink"]], face = "plain",
                just = "centre", rot = 0, lineheight = 0.92) {
    textGrob(
        label, x = unit(x, "npc"), y = unit(y, "npc"),
        just = just, rot = rot,
        gp = gpar(
            fontfamily = "Helvetica", fontsize = size * diagram_text_scale, col = col,
            fontface = face, lineheight = lineheight
        )
    )
}

seg <- function(x0, y0, x1, y1, col = palette[["mid"]], lwd = 0.8,
                arrow_end = FALSE, lty = 1) {
    segmentsGrob(
        x0 = unit(x0, "npc"), y0 = unit(y0, "npc"),
        x1 = unit(x1, "npc"), y1 = unit(y1, "npc"),
        arrow = if (arrow_end) arrow(type = "closed", length = unit(1.4, "mm")) else NULL,
        gp = gpar(col = col, lwd = lwd, lty = lty, lineend = "round")
    )
}

pt <- function(x, y, fill, size = 1.5, col = "white", lwd = 0.35, pch = 21) {
    pointsGrob(
        x = unit(x, "npc"), y = unit(y, "npc"), pch = pch,
        size = unit(size, "mm"),
        gp = gpar(fill = fill, col = col, lwd = lwd)
    )
}

poly <- function(x, y, fill = "white", col = palette[["line"]], lwd = 0.8) {
    polygonGrob(
        x = unit(x, "npc"), y = unit(y, "npc"),
        gp = gpar(fill = fill, col = col, lwd = lwd, linejoin = "round")
    )
}

pill <- function(label, x, y, w, fill, col, size = 5.0) {
    add_grob(rr(x, y, w, 0.057, fill = fill, col = col, lwd = 0.75, radius = 1.7))
    add_grob(txt(label, x, y, size = size, col = col, face = "bold"))
}

# Panel background and tag.
add_grob(rectGrob(gp = gpar(fill = "white", col = NA)))
add_grob(txt("a", 0.007, 0.975, size = 9, face = "bold", just = c("left", "top")))

# Four principal columns.
col1 <- list(x = 0.025, y = 0.11, w = 0.155, h = 0.80)
col2 <- list(x = 0.215, y = 0.11, w = 0.205, h = 0.80)
col3 <- list(x = 0.455, y = 0.11, w = 0.340, h = 0.80)
col4 <- list(x = 0.830, y = 0.11, w = 0.150, h = 0.80)

add_grob(rr(col1$x + col1$w / 2, col1$y + col1$h / 2, col1$w, col1$h,
            fill = "white", col = palette[["navy"]], lwd = 1.15, radius = 2.6))
add_grob(rr(col2$x + col2$w / 2, col2$y + col2$h / 2, col2$w, col2$h,
            fill = "white", col = palette[["teal"]], lwd = 1.15, radius = 2.6))
add_grob(rr(col4$x + col4$w / 2, col4$y + col4$h / 2, col4$w, col4$h,
            fill = "white", col = palette[["gold"]], lwd = 1.15, radius = 2.6))

add_grob(txt("Spatial inputs", col1$x + col1$w / 2, 0.865, size = 7.0,
             col = palette[["navy"]], face = "bold"))
add_grob(txt("Feature-selection objectives", col2$x + col2$w / 2, 0.865, size = 7.0,
             col = palette[["teal"]], face = "bold"))
add_grob(txt("Two-track evaluation", col3$x + col3$w / 2, 0.920, size = 7.0,
             col = palette[["ink"]], face = "bold"))
add_grob(txt("Objective-aware\nguidance", col4$x + col4$w / 2, 0.852, size = 6.8,
             col = palette[["gold"]], face = "bold"))

# Column 1: two spatial slices and context labels.
slice_shape <- list(
    x = c(-0.060, -0.026, 0.018, 0.057, 0.068, 0.039, -0.012, -0.057),
    y = c(0.010, 0.047, 0.054, 0.029, -0.012, -0.048, -0.055, -0.028)
)
for (i in seq_len(2)) {
    cy <- c(0.690, 0.545)[i]
    add_grob(poly(col1$x + col1$w / 2 + slice_shape$x, cy + slice_shape$y,
                  fill = palette[["pale_navy"]], col = "#AFC3DA", lwd = 0.7))
    xs <- seq(col1$x + 0.045, col1$x + 0.115, length.out = 5)
    ys <- seq(cy - 0.028, cy + 0.028, length.out = 4)
    grid_xy <- expand.grid(x = xs, y = ys)
    keep <- with(grid_xy, abs(x - (col1$x + col1$w / 2)) / 0.075 + abs(y - cy) / 0.055 < 1.25)
    grid_xy <- grid_xy[keep, , drop = FALSE]
    fills <- if (i == 1) {
        rep(c(palette[["navy"]], "#6FA8DC", palette[["teal"]]), length.out = nrow(grid_xy))
    } else {
        rep(c(palette[["green"]], palette[["navy"]], "#86B96B"), length.out = nrow(grid_xy))
    }
    for (j in seq_len(nrow(grid_xy))) {
        add_grob(pt(grid_xy$x[j], grid_xy$y[j], fills[j], size = 1.1, col = NA))
    }
}
add_grob(seg(col1$x + 0.025, 0.455, col1$x + col1$w - 0.025, 0.455,
             col = "#C8D5E5", lwd = 0.6))

input_rows <- list(
    list(y = 0.390, col = palette[["navy"]], head = "7 public datasets", sub = "multi-slice tissues"),
    list(y = 0.290, col = palette[["teal"]], head = "Biological contexts", sub = "full and lineage subsets"),
    list(y = 0.190, col = palette[["coral"]], head = "Validation designs", sub = "simulation and\nheld-out slices")
)
for (row in input_rows) {
    add_grob(rr(col1$x + 0.030, row$y, 0.020, 0.036, fill = row$col, col = row$col, radius = 1.1))
    add_grob(txt(row$head, col1$x + 0.048, row$y + 0.012, size = 5.35,
                 face = "bold", just = "left"))
    add_grob(txt(row$sub, col1$x + 0.048, row$y - 0.018, size = 4.75,
                 col = palette[["mid"]], just = "left"))
}

# Column 2: objective cards.
objective_cards <- list(
    list(y = 0.690, fill = palette[["pale_navy"]], col = palette[["navy"]],
         title = "Expression-driven / HVG", sub = "variability and local\nexpression structure", icon = "bars"),
    list(y = 0.500, fill = palette[["pale_teal"]], col = palette[["teal"]],
         title = "Spatially informed / SVG", sub = "coordinate-dependent\nexpression patterns", icon = "spatial"),
    list(y = 0.310, fill = palette[["pale_gray"]], col = palette[["mid"]],
         title = "Controls and oracle", sub = "random, all features,\nstable genes and labels", icon = "control")
)
for (card in objective_cards) {
    add_grob(rr(col2$x + col2$w / 2, card$y, col2$w - 0.030, 0.148,
                fill = card$fill, col = card$col, lwd = 0.75, radius = 2.0))
    icon_x <- col2$x + 0.038
    if (card$icon == "bars") {
        heights <- c(0.025, 0.048, 0.035, 0.060)
        for (k in seq_along(heights)) {
            add_grob(rectGrob(
                x = unit(icon_x - 0.016 + k * 0.010, "npc"),
                y = unit(card$y - 0.026 + heights[k] / 2, "npc"),
                width = unit(0.007, "npc"), height = unit(heights[k], "npc"),
                gp = gpar(fill = card$col, col = NA)
            ))
        }
    } else if (card$icon == "spatial") {
        locs <- expand.grid(x = c(icon_x - 0.016, icon_x, icon_x + 0.016),
                            y = c(card$y - 0.025, card$y, card$y + 0.025))
        fills <- c("#A8D8D2", card$col, "#A8D8D2", card$col, card$col,
                   "#A8D8D2", "#A8D8D2", card$col, "#A8D8D2")
        for (k in seq_len(nrow(locs))) {
            add_grob(pt(locs$x[k], locs$y[k], fills[k], size = 1.5, col = "white"))
        }
    } else {
        add_grob(seg(icon_x - 0.020, card$y - 0.025, icon_x + 0.020, card$y + 0.025,
                     col = card$col, lwd = 0.8, lty = 2))
        add_grob(seg(icon_x - 0.020, card$y + 0.025, icon_x + 0.020, card$y - 0.025,
                     col = card$col, lwd = 0.8, lty = 2))
        add_grob(pt(icon_x, card$y, "white", size = 2.2, col = card$col, lwd = 0.8))
    }
    add_grob(txt(card$title, col2$x + 0.074, card$y + 0.025, size = 5.65,
                 col = card$col, face = "bold", just = "left"))
    add_grob(txt(card$sub, col2$x + 0.074, card$y - 0.024, size = 4.55,
                 col = palette[["mid"]], just = "left"))
}
add_grob(txt("Feature budget: 100 to 10,000 genes", col2$x + col2$w / 2, 0.180,
             size = 5.0, col = palette[["mid"]], face = "bold"))

# Arrows into and out of feature selection.
add_grob(seg(col1$x + col1$w + 0.006, 0.510, col2$x - 0.009, 0.510,
             col = palette[["navy"]], lwd = 1.25, arrow_end = TRUE))
add_grob(seg(col2$x + col2$w + 0.006, 0.510, col3$x - 0.022, 0.510,
             col = palette[["ink"]], lwd = 1.0))
add_grob(seg(col3$x - 0.022, 0.510, col3$x - 0.022, 0.690,
             col = palette[["coral"]], lwd = 0.95))
add_grob(seg(col3$x - 0.022, 0.510, col3$x - 0.022, 0.325,
             col = palette[["indigo"]], lwd = 0.95))
add_grob(seg(col3$x - 0.022, 0.690, col3$x - 0.006, 0.690,
             col = palette[["coral"]], lwd = 0.95, arrow_end = TRUE))
add_grob(seg(col3$x - 0.022, 0.325, col3$x - 0.006, 0.325,
             col = palette[["indigo"]], lwd = 0.95, arrow_end = TRUE))

# Column 3 top branch: spatial signal evidence.
top_y <- 0.535
top_h <- 0.335
add_grob(rr(col3$x + col3$w / 2, top_y + top_h / 2, col3$w, top_h,
            fill = palette[["pale_coral"]], col = palette[["coral"]], lwd = 1.0, radius = 2.4))
add_grob(txt("Spatial signal evidence", col3$x + 0.018, 0.830, size = 6.35,
             col = palette[["coral"]], face = "bold", just = "left"))
add_grob(txt("Does the selector recover reproducible spatial structure?", col3$x + 0.018, 0.792,
             size = 4.65, col = palette[["mid"]], just = "left"))

top_cards <- list(
    list(x = col3$x + 0.060, title = "Simulation truth", sub = "known spatial genes", kind = "truth"),
    list(x = col3$x + 0.170, title = "Held-out slices", sub = "no reselection", kind = "holdout"),
    list(x = col3$x + 0.280, title = "Cross-statistic", sub = "Moran's I / nnSVG", kind = "cross")
)
for (card in top_cards) {
    add_grob(rr(card$x, 0.675, 0.095, 0.170, fill = "white", col = "#E7B9AF", lwd = 0.65, radius = 1.8))
    if (card$kind == "truth") {
        locs <- expand.grid(x = card$x + c(-0.023, 0, 0.023), y = 0.705 + c(-0.020, 0, 0.020))
        fills <- ifelse(locs$x > card$x - 0.005 & locs$y > 0.695, palette[["coral"]], "#F3C8BF")
        for (k in seq_len(nrow(locs))) add_grob(pt(locs$x[k], locs$y[k], fills[k], size = 1.2, col = "white"))
    } else if (card$kind == "holdout") {
        add_grob(rr(card$x - 0.014, 0.705, 0.038, 0.050, fill = "white", col = palette[["coral"]], lwd = 0.7, radius = 0.8))
        add_grob(rr(card$x + 0.019, 0.690, 0.038, 0.050, fill = "#F5D8D1", col = palette[["coral"]], lwd = 0.7, radius = 0.8))
        add_grob(seg(card$x - 0.002, 0.705, card$x + 0.006, 0.697, col = palette[["coral"]], lwd = 0.7, arrow_end = TRUE))
    } else {
        add_grob(seg(card$x - 0.027, 0.687, card$x - 0.008, 0.722, col = palette[["navy"]], lwd = 1.0))
        add_grob(seg(card$x - 0.008, 0.722, card$x + 0.010, 0.696, col = palette[["navy"]], lwd = 1.0))
        add_grob(seg(card$x + 0.010, 0.696, card$x + 0.029, 0.730, col = palette[["green"]], lwd = 1.0))
    }
    add_grob(txt(card$title, card$x, 0.630, size = 4.95, col = palette[["coral"]], face = "bold"))
    add_grob(txt(card$sub, card$x, 0.593, size = 4.35, col = palette[["mid"]]))
}
add_grob(txt("Spatial recovery is evaluated separately from downstream performance",
             col3$x + col3$w / 2, 0.553, size = 4.65, col = palette[["coral"]], face = "bold"))

# Column 3 bottom branch: downstream representation utility.
bot_y <- 0.130
bot_h <- 0.335
add_grob(rr(col3$x + col3$w / 2, bot_y + bot_h / 2, col3$w, bot_h,
            fill = palette[["pale_indigo"]], col = palette[["indigo"]], lwd = 1.0, radius = 2.4))
add_grob(txt("Downstream representation utility", col3$x + 0.018, 0.426, size = 6.25,
             col = palette[["indigo"]], face = "bold", just = "left"))
add_grob(txt("Does the selected representation support the intended workflow?", col3$x + 0.018, 0.389,
             size = 4.65, col = palette[["mid"]], just = "left"))

pill("scVI", col3$x + 0.060, 0.330, 0.070, "white", palette[["navy"]])
pill("CellCharter", col3$x + 0.170, 0.330, 0.105, "white", palette[["teal"]], size = 4.7)
pill("GraphST", col3$x + 0.280, 0.330, 0.085, "white", palette[["indigo"]])

task_cards <- list(
    list(x = col3$x + 0.060, title = "Integration", col = palette[["coral"]]),
    list(x = col3$x + 0.170, title = "Clustering", col = palette[["navy"]]),
    list(x = col3$x + 0.280, title = "Alignment", col = palette[["green"]])
)
for (card in task_cards) {
    add_grob(rr(card$x, 0.230, 0.095, 0.095, fill = "white", col = "#CFC8E4", lwd = 0.65, radius = 1.6))
    add_grob(pt(card$x - 0.020, 0.250, card$col, size = 1.35, col = "white"))
    add_grob(pt(card$x, 0.265, card$col, size = 1.35, col = "white"))
    add_grob(pt(card$x + 0.020, 0.242, card$col, size = 1.35, col = "white"))
    add_grob(txt(card$title, card$x, 0.194, size = 4.85, col = card$col, face = "bold"))
}
add_grob(txt("13 metrics; unavailable endpoints remain missing, not zero",
             col3$x + col3$w / 2, 0.148, size = 4.65, col = palette[["indigo"]], face = "bold"))

# Arrows from both evidence tracks into guidance.
add_grob(seg(col3$x + col3$w + 0.006, 0.690, col4$x - 0.012, 0.690,
             col = palette[["coral"]], lwd = 0.95, arrow_end = TRUE))
add_grob(seg(col3$x + col3$w + 0.006, 0.325, col4$x - 0.012, 0.325,
             col = palette[["indigo"]], lwd = 0.95, arrow_end = TRUE))

# Column 4: recommendation cards and closing statement.
guidance_cards <- list(
    list(y = 0.690, col = palette[["coral"]], fill = palette[["pale_coral"]],
         title = "Spatial discovery", sub = "SVG + truth or\nheld-out evidence"),
    list(y = 0.505, col = palette[["indigo"]], fill = palette[["pale_indigo"]],
         title = "Integration /\nclustering", sub = "HVG candidates +\nintended model"),
    list(y = 0.320, col = palette[["teal"]], fill = palette[["pale_teal"]],
         title = "Context and\nbudget", sub = "500-2,000 genes +\nsensitivity checks")
)
for (card in guidance_cards) {
    add_grob(rr(col4$x + col4$w / 2, card$y, col4$w - 0.025, 0.145,
                fill = card$fill, col = card$col, lwd = 0.75, radius = 2.0))
    add_grob(rr(col4$x + 0.025, card$y + 0.028, 0.022, 0.040,
                fill = card$col, col = card$col, lwd = 0.5, radius = 1.0))
    add_grob(txt(card$title, col4$x + 0.044, card$y + 0.028, size = 5.15,
                 col = card$col, face = "bold", just = "left"))
    add_grob(txt(card$sub, col4$x + 0.044, card$y - 0.030, size = 4.45,
                 col = palette[["mid"]], just = "left"))
}
add_grob(rr(col4$x + col4$w / 2, 0.175, col4$w - 0.025, 0.105,
            fill = palette[["pale_gold"]], col = palette[["gold"]], lwd = 0.9, radius = 2.0))
add_grob(txt("No universal winner", col4$x + col4$w / 2, 0.190, size = 5.65,
             col = palette[["gold"]], face = "bold"))
add_grob(txt("Match selector to objective\nand workflow", col4$x + col4$w / 2, 0.148, size = 4.65,
             col = palette[["ink"]], face = "bold"))

panel_a_grob <- do.call(grobTree, grobs)
panel_a <- as.ggplot(panel_a_grob) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))

metric_panel <- readRDS(metric_rds)
panel_b <- as.ggplot(metric_panel) +
    labs(tag = "b", title = "Downstream metric characterization") +
    theme(
        plot.tag = element_text(family = "Helvetica", face = "bold", size = 9, colour = "black"),
        plot.tag.position = c(0.002, 0.998),
        plot.title = element_text(
            family = "Helvetica", face = "bold", size = 8,
            hjust = 0.5, margin = margin(b = 1)
        ),
        plot.margin = margin(0, 0, 0, 0)
    )

figure <- wrap_plots(
    wrap_elements(full = panel_a),
    wrap_elements(full = panel_b),
    ncol = 1,
    heights = c(0.40, 0.60)
) &
    theme(plot.margin = margin(1, 1, 1, 1))

stem <- file.path(output_dir, "Figure1_story_aligned")
width_mm <- 284.5
height_mm <- 294.6

ggsave(
    paste0(stem, ".svg"), figure,
    width = width_mm, height = height_mm, units = "mm",
    device = svglite::svglite, bg = "white"
)
svg_text <- readLines(paste0(stem, ".svg"), warn = FALSE)
svg_text <- gsub(" textLength='[^']+' lengthAdjust='spacingAndGlyphs'", "", svg_text)
svg_text <- gsub(
    'font-family: "Liberation Sans";',
    "font-family: Arial, Helvetica, sans-serif;",
    svg_text, fixed = TRUE
)
writeLines(svg_text, paste0(stem, ".svg"), useBytes = TRUE)

register_arial_pdf_font()
ggsave(
    paste0(stem, ".pdf"), figure,
    width = width_mm, height = height_mm, units = "mm",
    device = grDevices::pdf, family = "Arial",
    useDingbats = FALSE, bg = "white"
)
ggsave(
    paste0(stem, ".tiff"), figure,
    width = width_mm, height = height_mm, units = "mm", dpi = 600,
    device = ragg::agg_tiff, compression = "lzw", bg = "white"
)
ggsave(
    paste0(stem, ".png"), figure,
    width = width_mm, height = height_mm, units = "mm", dpi = 300,
    device = ragg::agg_png, bg = "white"
)
saveRDS(figure, paste0(stem, ".rds"))

message("Story-aligned Figure 1 written to ", output_dir)
