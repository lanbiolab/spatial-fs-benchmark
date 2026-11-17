#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ggplot2)
    library(ggplotify)
    library(grid)
})

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) >= 1) args[[1]] else file.path("results", "figure1_flowchart_review")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

pal <- c(
    blue = "#356FAE",
    teal = "#1F766F",
    coral = "#B64B3C",
    indigo = "#6B5AA5",
    green = "#347A52",
    gold = "#8F6400",
    gold_accent = "#D0A020",
    ink = "#263241",
    slate = "#667487",
    line = "#C5CFDA",
    blue_bg = "#F0F5FA",
    teal_bg = "#EFF7F5",
    coral_bg = "#FDF2EF",
    indigo_bg = "#F4F2F9",
    gold_bg = "#FCF7EA",
    neutral_bg = "#F6F8FA"
)

grobs <- list()
add <- function(x) grobs[[length(grobs) + 1L]] <<- x

rr <- function(x, y, w, h, fill = "white", col = pal[["line"]], lwd = 0.70, r = 2.0) {
    roundrectGrob(
        x = unit(x, "npc"), y = unit(y, "npc"),
        width = unit(w, "npc"), height = unit(h, "npc"), r = unit(r, "mm"),
        gp = gpar(fill = fill, col = col, lwd = lwd)
    )
}

tx <- function(text, x, y, size = 6.0, col = pal[["ink"]], face = "plain",
               just = "centre", lineheight = 0.92) {
    textGrob(
        text, x = unit(x, "npc"), y = unit(y, "npc"), just = just,
        gp = gpar(
            fontfamily = "Arial", fontsize = size, col = col,
            fontface = face, lineheight = lineheight
        )
    )
}

seg <- function(x0, y0, x1, y1, col = pal[["slate"]], lwd = 0.75,
                arrow_end = FALSE, lty = 1) {
    segmentsGrob(
        x0 = unit(x0, "npc"), y0 = unit(y0, "npc"),
        x1 = unit(x1, "npc"), y1 = unit(y1, "npc"),
        arrow = if (arrow_end) arrow(type = "closed", length = unit(1.35, "mm")) else NULL,
        gp = gpar(col = col, lwd = lwd, lty = lty, lineend = "round")
    )
}

pt <- function(x, y, fill, size = 1.0, col = "white", lwd = 0.25, pch = 21) {
    pointsGrob(
        x = unit(x, "npc"), y = unit(y, "npc"), pch = pch,
        size = unit(size, "mm"), gp = gpar(fill = fill, col = col, lwd = lwd)
    )
}

poly <- function(x, y, fill = "white", col = pal[["line"]], lwd = 0.70, lty = 1) {
    polygonGrob(
        x = unit(x, "npc"), y = unit(y, "npc"),
        gp = gpar(fill = fill, col = col, lwd = lwd, lty = lty, linejoin = "round")
    )
}

tile <- function(x, y, w, h, fill, col = NA, lwd = 0.15) {
    rectGrob(
        x = unit(x, "npc"), y = unit(y, "npc"),
        width = unit(w, "npc"), height = unit(h, "npc"),
        gp = gpar(fill = fill, col = col, lwd = lwd)
    )
}

tissue <- function(cx, cy, sx = 1, sy = 1, fill = pal[["blue_bg"]], col = "#AFC3DA") {
    px <- c(-0.060, -0.034, 0.004, 0.046, 0.061, 0.036, -0.012, -0.057) * sx
    py <- c(0.005, 0.038, 0.047, 0.029, -0.008, -0.039, -0.046, -0.025) * sy
    add(poly(cx + px, cy + py, fill = fill, col = col, lwd = 0.62))
}

matrix_icon <- function(cx, cy, ncol = 7, nrow = 4, w = 0.007, h = 0.014,
                        colours = c("#EEF3F8", "#BFD2E5", "#75A4D0", "#356FAE"),
                        selected = integer()) {
    xs <- cx + (seq_len(ncol) - (ncol + 1) / 2) * w * 1.18
    ys <- cy + (seq_len(nrow) - (nrow + 1) / 2) * h * 1.15
    for (j in seq_len(ncol)) {
        for (i in seq_len(nrow)) {
            idx <- 1L + ((i * 2L + j * 3L + i * j) %% length(colours))
            add(tile(xs[j], ys[i], w, h, colours[idx], col = "white"))
        }
        if (j %in% selected) {
            add(rr(xs[j], cy, w * 1.30, h * nrow * 1.30,
                   fill = NA, col = pal[["gold_accent"]], lwd = 0.62, r = 0.35))
        }
    }
}

stage_label <- function(number, label, x, accent) {
    add(tx(number, x, 0.935, size = 4.8, col = accent, face = "bold", just = "left"))
    add(seg(x, 0.905, x + 0.065, 0.905, col = accent, lwd = 1.45))
    add(tx(label, x, 0.875, size = 5.5, col = pal[["ink"]], face = "bold", just = "left"))
}

add(rectGrob(gp = gpar(fill = "white", col = NA)))

stage_label("01", "INPUT", 0.030, pal[["blue"]])
stage_label("02", "SELECT", 0.205, pal[["teal"]])
stage_label("03", "EVALUATE", 0.445, pal[["slate"]])
stage_label("04", "GUIDANCE", 0.785, pal[["gold"]])

# Input process node ---------------------------------------------------------
add(rr(0.095, 0.505, 0.140, 0.535, fill = "white", col = "#AAB8C7", lwd = 0.78, r = 2.3))
add(tx("Spatial data", 0.095, 0.720, size = 6.6, col = pal[["blue"]], face = "bold"))

tissue(0.095, 0.615, sx = 0.75, sy = 0.72)
xy <- expand.grid(x = 0.095 + seq(-0.032, 0.032, length.out = 5),
                  y = 0.615 + seq(-0.016, 0.016, length.out = 3))
cols <- rep(c(pal[["blue"]], "#75A9D6", pal[["teal"]]), length.out = nrow(xy))
for (i in seq_len(nrow(xy))) add(pt(xy$x[i], xy$y[i], cols[i], size = 0.85, col = NA))
add(seg(0.095, 0.550, 0.095, 0.510, col = pal[["blue"]], lwd = 0.70, arrow_end = TRUE))
matrix_icon(0.095, 0.430, ncol = 8, nrow = 5, w = 0.008, h = 0.015)
add(tx("Expression × space", 0.095, 0.340, size = 4.9, col = pal[["blue"]]))

# Feature-selection process node -------------------------------------------
add(rr(0.285, 0.505, 0.180, 0.535, fill = "white", col = "#AAB8C7", lwd = 0.78, r = 2.3))
add(tx("Feature selection", 0.285, 0.720, size = 6.6, col = pal[["teal"]], face = "bold"))

routes <- list(
    list(y = 0.610, label = "Expression-driven", fill = pal[["blue_bg"]], col = pal[["blue"]], kind = "matrix"),
    list(y = 0.505, label = "Spatially informed", fill = pal[["teal_bg"]], col = pal[["teal"]], kind = "spatial"),
    list(y = 0.400, label = "Controls", fill = pal[["neutral_bg"]], col = pal[["slate"]], kind = "control")
)
for (route in routes) {
    add(rr(0.285, route$y, 0.150, 0.082, fill = route$fill, col = "#C5CFDA", lwd = 0.55, r = 1.6))
    ix <- 0.235
    if (route$kind == "matrix") {
        matrix_icon(ix, route$y, ncol = 4, nrow = 3, w = 0.005, h = 0.012, selected = c(2, 4))
    } else if (route$kind == "spatial") {
        tissue(ix, route$y, sx = 0.22, sy = 0.28, fill = "white", col = "#93C9C2")
        for (j in seq_len(5)) {
            add(pt(ix + c(-0.010, 0, 0.010, -0.006, 0.008)[j],
                   route$y + c(-0.010, -0.006, 0, 0.009, 0.011)[j],
                   if (j >= 4) pal[["teal"]] else "#AED9D4", size = 0.72, col = NA))
        }
    } else {
        for (j in seq_len(6)) {
            add(rr(ix - 0.015 + (j - 1) * 0.006, route$y + 0.008, 0.004, 0.012,
                   fill = "#D3DAE2", col = pal[["slate"]], lwd = 0.35, r = 0.25))
            add(rr(ix - 0.015 + (j - 1) * 0.006, route$y - 0.009, 0.004, 0.012,
                   fill = if (j %in% c(2, 5)) "#8C99A9" else "white",
                   col = pal[["slate"]], lwd = 0.35, r = 0.25))
        }
    }
    add(tx(route$label, 0.260, route$y, size = 5.15, col = route$col, just = "left"))
}

# Main flow into the evidence fork.
add(seg(0.167, 0.505, 0.192, 0.505, col = pal[["slate"]], lwd = 0.90, arrow_end = TRUE))
add(seg(0.377, 0.505, 0.400, 0.505, col = pal[["slate"]], lwd = 0.90, arrow_end = TRUE))
add(pt(0.410, 0.505, pal[["ink"]], size = 1.55, col = "white", lwd = 0.35))

# Two-track evaluation ------------------------------------------------------
add(rr(0.575, 0.680, 0.270, 0.225, fill = pal[["coral_bg"]], col = "#E3A99F", lwd = 0.68, r = 2.0))
add(tx("Spatial signal recovery", 0.458, 0.755, size = 6.0, col = pal[["coral"]], face = "bold", just = "left"))
add(seg(0.548, 0.620, 0.548, 0.715, col = "#E9C5BE", lwd = 0.55))

# Simulation truth icon.
tissue(0.500, 0.674, sx = 0.24, sy = 0.28, fill = "white", col = "#E8B9B0")
tissue(0.532, 0.674, sx = 0.24, sy = 0.28, fill = "white", col = "#E8B9B0")
for (side in c(0.500, 0.532)) {
    xy <- expand.grid(x = side + c(-0.010, 0, 0.010), y = 0.674 + c(-0.008, 0, 0.008))
    for (i in seq_len(nrow(xy))) {
        hot <- xy$x[i] > side - 0.002 && xy$y[i] > 0.670
        add(pt(xy$x[i], xy$y[i], if (hot) pal[["coral"]] else "#F0C0B8", size = 0.65, col = NA))
    }
}
add(seg(0.512, 0.674, 0.520, 0.674, col = pal[["coral"]], lwd = 0.60, arrow_end = TRUE))
add(tx("Simulation\ntruth", 0.516, 0.615, size = 4.30, col = pal[["coral"]]))

# Held-out icon.
for (k in 0:2) {
    xx <- 0.590 + k * 0.030
    add(rr(xx, 0.677 - k * 0.002, 0.040, 0.045,
           fill = if (k < 2) "#F7D9D3" else "white", col = pal[["coral"]], lwd = 0.55, r = 0.8))
    if (k == 2) add(seg(xx - 0.016, 0.677, xx + 0.016, 0.677, col = pal[["coral"]], lwd = 0.45, lty = 2))
}
add(seg(0.612, 0.638, 0.641, 0.638, col = pal[["coral"]], lwd = 0.60, arrow_end = TRUE))
add(tx("Held-out\nvalidation", 0.620, 0.615, size = 4.30, col = pal[["coral"]]))

add(rr(0.575, 0.330, 0.270, 0.225, fill = pal[["indigo_bg"]], col = "#BAB1D4", lwd = 0.68, r = 2.0))
add(tx("Downstream utility", 0.458, 0.405, size = 6.0, col = pal[["indigo"]], face = "bold", just = "left"))

task_x <- c(0.495, 0.575, 0.655)
for (xx in task_x) add(rr(xx, 0.318, 0.070, 0.130, fill = "white", col = "#D3CEE5", lwd = 0.52, r = 1.3))

# Integration: separation to mixing.
ix <- task_x[1]
for (i in seq_len(8)) {
    add(pt(ix - 0.022 + c(-0.008, 0, 0.008, -0.006)[1 + ((i - 1) %% 4)],
           0.340 + c(0.020, 0.008, -0.010, -0.023)[1 + ((i - 1) %% 4)],
           if (i <= 4) pal[["coral"]] else pal[["blue"]], size = 0.70, col = NA))
}
add(seg(ix - 0.004, 0.340, ix + 0.008, 0.340, col = pal[["indigo"]], lwd = 0.55, arrow_end = TRUE))
for (i in seq_len(6)) {
    add(pt(ix + 0.016 + c(-0.006, 0.003, 0.010)[1 + ((i - 1) %% 3)],
           0.340 + c(0.018, 0, -0.018)[1 + ((i - 1) %% 3)],
           if (i %% 2) pal[["coral"]] else pal[["blue"]], size = 0.70, col = NA))
}
add(tx("Integration", ix, 0.282, size = 4.45, col = pal[["blue"]]))

# Clustering: coherent domains.
cx <- task_x[2]
tissue(cx, 0.340, sx = 0.43, sy = 0.42, fill = "#FBFCFD", col = "#BFC9D6")
xy <- expand.grid(x = cx + c(-0.020, -0.007, 0.007, 0.020), y = 0.340 + c(-0.015, 0, 0.015))
for (i in seq_len(nrow(xy))) {
    cc <- if (xy$x[i] < cx - 0.003) pal[["blue"]] else if (xy$y[i] > 0.340) pal[["coral"]] else pal[["green"]]
    add(pt(xy$x[i], xy$y[i], cc, size = 0.66, col = NA))
}
add(tx("Clustering", cx, 0.282, size = 4.45, col = pal[["coral"]]))

# Alignment: slice correspondence.
ax <- task_x[3]
tissue(ax - 0.013, 0.350, sx = 0.28, sy = 0.29, fill = pal[["blue_bg"]], col = pal[["blue"]])
tissue(ax + 0.013, 0.327, sx = 0.28, sy = 0.29, fill = pal[["teal_bg"]], col = pal[["teal"]])
for (k in -1:1) add(seg(ax - 0.020 + k * 0.008, 0.350 + k * 0.003,
                        ax + 0.006 + k * 0.008, 0.327 + k * 0.003,
                        col = "#9AA7B7", lwd = 0.38, lty = 2))
add(tx("Alignment", ax, 0.282, size = 4.45, col = pal[["green"]]))

# Fork and merge connectors.
add(seg(0.410, 0.505, 0.410, 0.680, col = pal[["coral"]], lwd = 0.72))
add(seg(0.410, 0.680, 0.434, 0.680, col = pal[["coral"]], lwd = 0.72, arrow_end = TRUE))
add(seg(0.410, 0.505, 0.410, 0.330, col = pal[["indigo"]], lwd = 0.72))
add(seg(0.410, 0.330, 0.434, 0.330, col = pal[["indigo"]], lwd = 0.72, arrow_end = TRUE))

add(seg(0.712, 0.680, 0.742, 0.680, col = pal[["coral"]], lwd = 0.72))
add(seg(0.742, 0.680, 0.742, 0.505, col = pal[["coral"]], lwd = 0.72))
add(seg(0.712, 0.330, 0.742, 0.330, col = pal[["indigo"]], lwd = 0.72))
add(seg(0.742, 0.330, 0.742, 0.505, col = pal[["indigo"]], lwd = 0.72))
add(pt(0.742, 0.505, pal[["ink"]], size = 1.45, col = "white", lwd = 0.30))

# Objective decision and guidance ------------------------------------------
add(seg(0.750, 0.505, 0.765, 0.505, col = pal[["slate"]], lwd = 0.86, arrow_end = TRUE))
add(poly(c(0.800, 0.845, 0.800, 0.755), c(0.575, 0.505, 0.435, 0.505),
         fill = pal[["gold_bg"]], col = pal[["gold"]], lwd = 0.70))
add(tx("Analysis\ngoal", 0.800, 0.505, size = 5.3, col = pal[["gold"]], face = "bold"))

# Goal branches.
add(seg(0.845, 0.505, 0.862, 0.505, col = pal[["gold"]], lwd = 0.64))
add(seg(0.862, 0.505, 0.862, 0.650, col = pal[["gold"]], lwd = 0.64))
add(seg(0.862, 0.650, 0.872, 0.650, col = pal[["gold"]], lwd = 0.64, arrow_end = TRUE))
add(seg(0.862, 0.505, 0.862, 0.360, col = pal[["gold"]], lwd = 0.64))
add(seg(0.862, 0.360, 0.872, 0.360, col = pal[["gold"]], lwd = 0.64, arrow_end = TRUE))

add(rr(0.925, 0.650, 0.105, 0.145, fill = pal[["coral_bg"]], col = "#D9988D", lwd = 0.64, r = 1.7))
for (i in seq_len(5)) add(pt(0.925 + c(-0.010, 0, 0.010, -0.005, 0.007)[i],
                              0.681 + c(-0.009, -0.008, -0.003, 0.007, 0.011)[i],
                              if (i >= 4) pal[["coral"]] else "#EAB7AE", size = 0.78, col = NA))
add(tx("Spatial\ndiscovery", 0.925, 0.625, size = 5.05, col = pal[["coral"]]))

add(rr(0.925, 0.360, 0.105, 0.145, fill = pal[["indigo_bg"]], col = "#AAA0CC", lwd = 0.64, r = 1.7))
add(pt(0.915, 0.391, pal[["blue"]], size = 0.78, col = NA))
add(pt(0.925, 0.370, pal[["coral"]], size = 0.78, col = NA))
add(pt(0.935, 0.390, pal[["green"]], size = 0.78, col = NA))
add(seg(0.916, 0.389, 0.924, 0.372, col = pal[["indigo"]], lwd = 0.42))
add(seg(0.926, 0.372, 0.934, 0.388, col = pal[["indigo"]], lwd = 0.42))
add(tx("Integrated\nanalysis", 0.925, 0.335, size = 5.05, col = pal[["indigo"]]))

# Both goals converge on a task-aware method choice.
add(seg(0.978, 0.650, 0.990, 0.650, col = pal[["gold"]], lwd = 0.60))
add(seg(0.978, 0.360, 0.990, 0.360, col = pal[["gold"]], lwd = 0.60))
add(seg(0.990, 0.650, 0.990, 0.205, col = pal[["gold"]], lwd = 0.60))
add(seg(0.990, 0.205, 0.925, 0.205, col = pal[["gold"]], lwd = 0.60))
add(seg(0.925, 0.205, 0.925, 0.183, col = pal[["gold"]], lwd = 0.68, arrow_end = TRUE))
add(rr(0.925, 0.125, 0.120, 0.105, fill = pal[["gold_bg"]], col = pal[["gold"]], lwd = 0.70, r = 1.7))
add(tx("Objective-aware\nmethod choice", 0.925, 0.125, size = 4.8, col = pal[["gold"]], face = "bold"))

figure_grob <- do.call(grobTree, grobs)
figure <- as.ggplot(figure_grob) +
    theme_void() +
    theme(plot.margin = margin(1.0, 1.0, 1.0, 1.0, "mm"))

stem <- file.path(output_dir, "Figure1_flowchart_review")
width_mm <- 183
height_mm <- 92

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

message("Figure 1 flowchart review written to ", output_dir)
