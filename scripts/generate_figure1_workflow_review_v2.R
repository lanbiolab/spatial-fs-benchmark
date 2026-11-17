#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ggplot2)
    library(ggplotify)
    library(grid)
})

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) >= 1) args[[1]] else file.path("results", "figure1_workflow_review_v3")
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

g <- list()
add <- function(x) g[[length(g) + 1L]] <<- x

rr <- function(x, y, w, h, fill = "white", col = pal[["line"]], lwd = 0.75, r = 2.0) {
    roundrectGrob(
        x = unit(x, "npc"), y = unit(y, "npc"),
        width = unit(w, "npc"), height = unit(h, "npc"), r = unit(r, "mm"),
        gp = gpar(fill = fill, col = col, lwd = lwd)
    )
}

tx <- function(text, x, y, size = 6.2, col = pal[["ink"]], face = "plain",
               just = "centre", lineheight = 0.92) {
    textGrob(
        text, x = unit(x, "npc"), y = unit(y, "npc"), just = just,
        gp = gpar(
            fontfamily = "Arial", fontsize = size, fontface = face,
            col = col, lineheight = lineheight
        )
    )
}

ln <- function(x0, y0, x1, y1, col = pal[["slate"]], lwd = 0.8,
               arrow_end = FALSE, lty = 1) {
    segmentsGrob(
        x0 = unit(x0, "npc"), y0 = unit(y0, "npc"),
        x1 = unit(x1, "npc"), y1 = unit(y1, "npc"),
        arrow = if (arrow_end) arrow(type = "closed", length = unit(1.45, "mm")) else NULL,
        gp = gpar(col = col, lwd = lwd, lty = lty, lineend = "round")
    )
}

pt <- function(x, y, fill, size = 1.15, col = "white", lwd = 0.3, pch = 21) {
    pointsGrob(
        x = unit(x, "npc"), y = unit(y, "npc"), pch = pch,
        size = unit(size, "mm"), gp = gpar(fill = fill, col = col, lwd = lwd)
    )
}

poly <- function(x, y, fill = "white", col = pal[["line"]], lwd = 0.7, lty = 1) {
    polygonGrob(
        x = unit(x, "npc"), y = unit(y, "npc"),
        gp = gpar(fill = fill, col = col, lwd = lwd, lty = lty, linejoin = "round")
    )
}

tile <- function(x, y, w, h, fill, col = NA, lwd = 0.2) {
    rectGrob(
        x = unit(x, "npc"), y = unit(y, "npc"),
        width = unit(w, "npc"), height = unit(h, "npc"),
        gp = gpar(fill = fill, col = col, lwd = lwd)
    )
}

stage <- function(x, w, accent, title, number, fill = "#FDFEFE") {
    add(rr(x + w / 2, 0.50, w, 0.84, fill = fill, col = "#BFC9D5", lwd = 0.86, r = 2.5))
    add(tx(number, x + 0.018, 0.858, size = 5.3, col = accent, face = "bold", just = "left"))
    add(ln(x + 0.046, 0.838, x + 0.046, 0.878, col = accent, lwd = 1.7))
    add(tx(title, x + 0.057, 0.858, size = 7.0, col = pal[["ink"]], face = "bold", just = "left"))
    add(ln(x + 0.017, 0.807, x + w - 0.017, 0.807, col = "#E1E6EC", lwd = 0.55))
}

tissue_shape <- function(cx, cy, sx = 1, sy = 1, fill = pal[["blue_bg"]], col = "#AFC3DA") {
    px <- c(-0.060, -0.034, 0.004, 0.046, 0.061, 0.036, -0.012, -0.057) * sx
    py <- c(0.005, 0.038, 0.047, 0.029, -0.008, -0.039, -0.046, -0.025) * sy
    add(poly(cx + px, cy + py, fill = fill, col = col, lwd = 0.72))
}

spot_grid <- function(cx, cy, xs, ys, fills, size = 1.0, keep = NULL) {
    coords <- expand.grid(x = cx + xs, y = cy + ys)
    if (!is.null(keep)) coords <- coords[keep(coords), , drop = FALSE]
    cols <- rep(fills, length.out = nrow(coords))
    for (i in seq_len(nrow(coords))) add(pt(coords$x[i], coords$y[i], cols[i], size = size, col = NA))
}

mini_matrix <- function(cx, cy, ncol = 8, nrow = 5, w = 0.008, h = 0.015,
                        palette = c("#EDF2F7", "#B9D0E8", "#5E93C9", "#2465A5"),
                        selected = integer()) {
    xvals <- cx + (seq_len(ncol) - (ncol + 1) / 2) * (w * 1.18)
    yvals <- cy + (seq_len(nrow) - (nrow + 1) / 2) * (h * 1.14)
    for (j in seq_len(ncol)) {
        for (i in seq_len(nrow)) {
            idx <- 1L + ((i * 2L + j * 3L + i * j) %% length(palette))
            add(tile(xvals[j], yvals[i], w, h, palette[idx], col = "white", lwd = 0.18))
        }
        if (j %in% selected) {
            add(rr(xvals[j], cy, w * 1.30, h * nrow * 1.28,
                   fill = NA, col = pal[["gold_accent"]], lwd = 0.75, r = 0.4))
        }
    }
}

add(rectGrob(gp = gpar(fill = "white", col = NA)))

input <- c(x = 0.022, w = 0.164)
select <- c(x = 0.222, w = 0.198)
eval <- c(x = 0.456, w = 0.318)
guide <- c(x = 0.805, w = 0.173)

stage(input[["x"]], input[["w"]], pal[["blue"]], "Spatial data", "01")
stage(select[["x"]], select[["w"]], pal[["teal"]], "Feature selection", "02")
stage(eval[["x"]], eval[["w"]], pal[["slate"]], "Evaluation", "03")
stage(guide[["x"]], guide[["w"]], pal[["gold"]], "Guidance", "04")

# Stage 1: paired spatial maps and their expression matrix.
tissue_shape(input[["x"]] + input[["w"]] / 2, 0.660, sx = 0.95, sy = 0.90)
spot_grid(
    input[["x"]] + input[["w"]] / 2, 0.660,
    xs = seq(-0.040, 0.040, length.out = 5), ys = seq(-0.022, 0.022, length.out = 3),
    fills = c(pal[["blue"]], "#75A9DA", pal[["teal"]]), size = 1.05
)
tissue_shape(input[["x"]] + input[["w"]] / 2, 0.505, sx = 0.95, sy = 0.90)
spot_grid(
    input[["x"]] + input[["w"]] / 2, 0.505,
    xs = seq(-0.040, 0.040, length.out = 5), ys = seq(-0.022, 0.022, length.out = 3),
    fills = c(pal[["green"]], pal[["blue"]], "#85B96D"), size = 1.05
)
add(ln(input[["x"]] + input[["w"]] / 2, 0.450,
       input[["x"]] + input[["w"]] / 2, 0.400,
       col = pal[["blue"]], lwd = 0.85, arrow_end = TRUE))
mini_matrix(input[["x"]] + input[["w"]] / 2, 0.300, ncol = 9, nrow = 5,
            w = 0.009, h = 0.016,
            palette = c("#F0F4F8", "#C4D7EA", "#77A6D2", "#2F6DB0"))
add(tx("spot × gene matrix", input[["x"]] + input[["w"]] / 2, 0.215, size = 5.0,
       col = pal[["blue"]], face = "plain"))

# Stage 2: distinct mechanisms rather than generic category blocks.
selector_cards <- list(
    list(y = 0.670, fill = pal[["blue_bg"]], col = pal[["blue"]], border = "#A8C1DA", label = "Expression-driven"),
    list(y = 0.455, fill = pal[["teal_bg"]], col = pal[["teal"]], border = "#A8D1CB", label = "Spatially informed"),
    list(y = 0.240, fill = pal[["neutral_bg"]], col = pal[["slate"]], border = "#BBC5D0", label = "Controls")
)
for (card in selector_cards) {
    add(rr(select[["x"]] + select[["w"]] / 2, card$y, select[["w"]] - 0.030, 0.145,
           fill = card$fill, col = card$border, lwd = 0.67, r = 1.8))
    add(tx(card$label, select[["x"]] + 0.079, card$y,
           size = 5.5, col = card$col, face = "plain", just = "left"))
}

# Expression-based icon: a small matrix with high-variance columns selected.
mini_matrix(select[["x"]] + 0.048, 0.670, ncol = 5, nrow = 4, w = 0.006, h = 0.014,
            palette = c("#E8F0F8", "#AFC9E2", "#5B91C7", "#2F6DB0"), selected = c(2, 5))

# Spatial icon: a coherent hotspot rather than an undifferentiated dot grid.
tissue_shape(select[["x"]] + 0.048, 0.455, sx = 0.38, sy = 0.48,
             fill = "white", col = "#9BCBC5")
spatial_xy <- expand.grid(x = select[["x"]] + 0.048 + c(-0.020, 0, 0.020),
                          y = 0.455 + c(-0.018, 0, 0.018))
for (i in seq_len(nrow(spatial_xy))) {
    d <- sqrt((spatial_xy$x[i] - (select[["x"]] + 0.060))^2 + (spatial_xy$y[i] - 0.465)^2)
    add(pt(spatial_xy$x[i], spatial_xy$y[i], if (d < 0.024) pal[["teal"]] else "#B7DDD8",
           size = 1.15, col = NA))
}

# Controls icon: all features above, fixed reference subset below.
for (i in seq_len(6)) {
    xx <- select[["x"]] + 0.025 + (i - 1) * 0.009
    add(rr(xx, 0.253, 0.006, 0.018,
           fill = "#D3DAE2", col = pal[["slate"]], lwd = 0.42, r = 0.35))
    add(rr(xx, 0.227, 0.006, 0.018,
           fill = if (i %in% c(2, 5)) "#8C99A9" else "white",
           col = pal[["slate"]], lwd = 0.42, r = 0.35))
}

# Stage-to-stage arrows.
add(ln(input[["x"]] + input[["w"]] + 0.007, 0.500,
       select[["x"]] - 0.010, 0.500, col = pal[["slate"]], lwd = 1.05, arrow_end = TRUE))
add(ln(select[["x"]] + select[["w"]] + 0.007, 0.500,
       eval[["x"]] - 0.010, 0.500, col = pal[["slate"]], lwd = 1.05, arrow_end = TRUE))

# Stage 3: two evidence tracks.
add(rr(eval[["x"]] + eval[["w"]] / 2, 0.640, eval[["w"]] - 0.022, 0.290,
       fill = pal[["coral_bg"]], col = "#E5A99E", lwd = 0.76, r = 2.0))
add(tx("Spatial signal recovery", eval[["x"]] + 0.020, 0.742,
       size = 6.5, col = pal[["coral"]], face = "bold", just = "left"))

# Simulation truth: true and recovered spatial patterns side by side.
sim_x <- eval[["x"]] + 0.092
add(rr(sim_x, 0.625, 0.124, 0.165, fill = "white", col = "#E9B9B0", lwd = 0.70, r = 1.7))
tissue_shape(sim_x - 0.027, 0.650, sx = 0.28, sy = 0.34, fill = "#FFF9F7", col = "#E9B9B0")
tissue_shape(sim_x + 0.027, 0.650, sx = 0.28, sy = 0.34, fill = "#FFF9F7", col = "#E9B9B0")
for (side in c(-0.027, 0.027)) {
    xy <- expand.grid(x = sim_x + side + c(-0.012, 0, 0.012), y = 0.650 + c(-0.011, 0, 0.011))
    for (i in seq_len(nrow(xy))) {
        hot <- xy$x[i] > sim_x + side - 0.002 && xy$y[i] > 0.645
        add(pt(xy$x[i], xy$y[i], if (hot) pal[["coral"]] else "#F4C8C0", size = 0.78, col = NA))
    }
}
add(ln(sim_x - 0.007, 0.650, sim_x + 0.007, 0.650,
       col = pal[["coral"]], lwd = 0.70, arrow_end = TRUE))
add(tx("Simulation truth", sim_x, 0.568, size = 5.3, col = pal[["coral"]], face = "plain"))

# Held-out validation: train slices transfer to an untouched slice.
hold_x <- eval[["x"]] + 0.236
add(rr(hold_x, 0.625, 0.124, 0.165, fill = "white", col = "#E9B9B0", lwd = 0.70, r = 1.7))
for (k in 0:2) {
    xx <- hold_x - 0.033 + k * 0.031
    add(rr(xx, 0.650 - k * 0.003, 0.041, 0.050,
           fill = if (k < 2) "#F8D9D3" else "white",
           col = pal[["coral"]], lwd = 0.70, r = 1.0))
    if (k == 2) {
        add(ln(xx - 0.017, 0.672, xx + 0.017, 0.672, col = pal[["coral"]], lwd = 0.55, lty = 2))
        add(ln(xx - 0.017, 0.628, xx + 0.017, 0.628, col = pal[["coral"]], lwd = 0.55, lty = 2))
    }
}
add(ln(hold_x - 0.004, 0.615, hold_x + 0.027, 0.615,
       col = pal[["coral"]], lwd = 0.75, arrow_end = TRUE))
add(tx("Held-out validation", hold_x, 0.568, size = 5.2, col = pal[["coral"]], face = "plain"))

add(rr(eval[["x"]] + eval[["w"]] / 2, 0.315, eval[["w"]] - 0.022, 0.290,
       fill = pal[["indigo_bg"]], col = "#B9B0D5", lwd = 0.76, r = 2.0))
add(tx("Downstream utility", eval[["x"]] + 0.020, 0.422,
       size = 6.5, col = pal[["indigo"]], face = "bold", just = "left"))

# Integration: separated batches become mixed while structure remains visible.
task_x <- c(eval[["x"]] + 0.058, eval[["x"]] + 0.159, eval[["x"]] + 0.258)
for (xx in task_x) add(rr(xx, 0.300, 0.089, 0.165, fill = "white", col = "#CEC6E4", lwd = 0.70, r = 1.7))
ix <- task_x[1]
left_pts <- list(
    c(ix - 0.030, 0.330, pal[["coral"]]), c(ix - 0.024, 0.345, pal[["coral"]]),
    c(ix - 0.017, 0.326, pal[["coral"]]), c(ix - 0.030, 0.298, pal[["blue"]]),
    c(ix - 0.021, 0.286, pal[["blue"]]), c(ix - 0.013, 0.303, pal[["blue"]])
)
for (p in left_pts) add(pt(p[1], p[2], p[3], size = 0.95, col = NA))
add(ln(ix - 0.006, 0.316, ix + 0.008, 0.316, col = pal[["indigo"]], lwd = 0.70, arrow_end = TRUE))
mixed <- data.frame(
    x = ix + c(0.016, 0.024, 0.032, 0.018, 0.029, 0.036),
    y = c(0.334, 0.322, 0.341, 0.296, 0.304, 0.286),
    col = c(pal[["coral"]], pal[["blue"]], pal[["coral"]], pal[["blue"]], pal[["coral"]], pal[["blue"]])
)
for (i in seq_len(nrow(mixed))) add(pt(mixed$x[i], mixed$y[i], mixed$col[i], size = 0.95, col = NA))
add(tx("Integration", ix, 0.240, size = 5.2, col = pal[["blue"]], face = "plain"))

# Clustering: coherent spatial domains within a tissue boundary.
cx <- task_x[2]
tissue_shape(cx, 0.322, sx = 0.54, sy = 0.52, fill = "#FAFBFC", col = "#BFC9D6")
cluster_xy <- expand.grid(x = cx + c(-0.025, -0.008, 0.009, 0.026), y = 0.322 + c(-0.022, 0, 0.022))
for (i in seq_len(nrow(cluster_xy))) {
    cc <- if (cluster_xy$x[i] < cx - 0.004) pal[["blue"]] else if (cluster_xy$y[i] > 0.320) pal[["coral"]] else pal[["green"]]
    add(pt(cluster_xy$x[i], cluster_xy$y[i], cc, size = 0.90, col = NA))
}
add(tx("Clustering", cx, 0.240, size = 5.2, col = pal[["coral"]], face = "plain"))

# Alignment: offset slices are registered and linked by spatial correspondences.
ax <- task_x[3]
tissue_shape(ax - 0.018, 0.337, sx = 0.35, sy = 0.35, fill = pal[["blue_bg"]], col = pal[["blue"]])
tissue_shape(ax + 0.018, 0.302, sx = 0.35, sy = 0.35, fill = pal[["teal_bg"]], col = pal[["teal"]])
for (k in -1:1) {
    add(ln(ax - 0.027 + k * 0.010, 0.337 + k * 0.004,
           ax + 0.009 + k * 0.010, 0.302 + k * 0.004,
           col = "#9AA7B7", lwd = 0.45, lty = 2))
}
add(tx("Alignment", ax, 0.240, size = 5.2, col = pal[["green"]], face = "plain"))

# The two evidence tracks are synthesized before guidance.
add(ln(eval[["x"]] + eval[["w"]] + 0.006, 0.500,
       guide[["x"]] - 0.008, 0.500, col = pal[["slate"]], lwd = 1.15, arrow_end = TRUE))

# Stage 4: an actual decision framework, not a legend.
diamond_x <- guide[["x"]] + guide[["w"]] / 2
add(poly(
    c(diamond_x, diamond_x + 0.052, diamond_x, diamond_x - 0.052),
    c(0.745, 0.690, 0.635, 0.690),
    fill = pal[["gold_bg"]], col = pal[["gold"]], lwd = 0.72
))
add(tx("Analysis goal", diamond_x, 0.690, size = 5.9, col = pal[["gold"]], face = "bold"))

decision_cards <- list(
    list(x = guide[["x"]] + 0.050, fill = pal[["coral_bg"]], col = pal[["coral"]], label = "Spatial\ndiscovery", icon = "hotspot"),
    list(x = guide[["x"]] + 0.130, fill = pal[["indigo_bg"]], col = pal[["indigo"]], label = "Integrated\nanalysis", icon = "embed")
)
for (card in decision_cards) {
    add(rr(card$x, 0.475, 0.070, 0.145,
           fill = card$fill, col = card$col, lwd = 0.68, r = 1.9))
    icon_x <- card$x
    if (card$icon == "hotspot") {
        for (i in seq_len(5)) {
            add(pt(icon_x + c(-0.010, 0, 0.010, -0.006, 0.007)[i],
                   0.502 + c(-0.010, -0.009, -0.004, 0.008, 0.012)[i],
                   if (i >= 4) card$col else "#F0B7AD", size = 1.0, col = NA))
        }
    } else {
        add(pt(icon_x - 0.010, 0.512, pal[["blue"]], size = 1.0, col = NA))
        add(pt(icon_x, 0.490, pal[["coral"]], size = 1.0, col = NA))
        add(pt(icon_x + 0.011, 0.507, pal[["green"]], size = 1.0, col = NA))
        add(ln(icon_x - 0.009, 0.510, icon_x - 0.001, 0.493, col = card$col, lwd = 0.55))
        add(ln(icon_x + 0.002, 0.493, icon_x + 0.010, 0.505, col = card$col, lwd = 0.55))
    }
    add(tx(card$label, card$x, 0.443,
           size = 4.8, col = card$col, face = "plain"))
}
add(ln(diamond_x, 0.635, diamond_x, 0.585, col = pal[["gold"]], lwd = 0.65))
add(ln(guide[["x"]] + 0.050, 0.585, guide[["x"]] + 0.130, 0.585,
       col = pal[["gold"]], lwd = 0.65))
add(ln(guide[["x"]] + 0.050, 0.585, guide[["x"]] + 0.050, 0.553,
       col = pal[["gold"]], lwd = 0.65, arrow_end = TRUE))
add(ln(guide[["x"]] + 0.130, 0.585, guide[["x"]] + 0.130, 0.553,
       col = pal[["gold"]], lwd = 0.65, arrow_end = TRUE))

add(rr(diamond_x, 0.215, guide[["w"]] - 0.030, 0.120,
       fill = pal[["gold_bg"]], col = pal[["gold"]], lwd = 0.76, r = 2.0))
add(tx("Objective-aware\nmethod choice", diamond_x, 0.215,
       size = 5.7, col = pal[["gold"]], face = "bold"))
add(ln(guide[["x"]] + 0.050, 0.402, guide[["x"]] + 0.050, 0.345,
       col = pal[["gold"]], lwd = 0.65))
add(ln(guide[["x"]] + 0.130, 0.402, guide[["x"]] + 0.130, 0.345,
       col = pal[["gold"]], lwd = 0.65))
add(ln(guide[["x"]] + 0.050, 0.345, guide[["x"]] + 0.130, 0.345,
       col = pal[["gold"]], lwd = 0.65))
add(ln(diamond_x, 0.345, diamond_x, 0.280,
       col = pal[["gold"]], lwd = 0.70, arrow_end = TRUE))

fig_grob <- do.call(grobTree, g)
fig <- as.ggplot(fig_grob) +
    theme_void() +
    theme(plot.margin = margin(1.2, 1.2, 1.2, 1.2, "mm"))

stem <- file.path(output_dir, "Figure1_workflow_review_v3")
width_mm <- 183
height_mm <- 94

ggsave(paste0(stem, ".svg"), fig, width = width_mm, height = height_mm, units = "mm",
       device = svglite::svglite, bg = "white")
svg <- readLines(paste0(stem, ".svg"), warn = FALSE)
svg <- gsub(" textLength='[^']+' lengthAdjust='spacingAndGlyphs'", "", svg)
writeLines(svg, paste0(stem, ".svg"), useBytes = TRUE)

ggsave(paste0(stem, ".pdf"), fig, width = width_mm, height = height_mm, units = "mm",
       device = grDevices::cairo_pdf, family = "Arial", bg = "white")
ggsave(paste0(stem, ".tiff"), fig, width = width_mm, height = height_mm, units = "mm",
       dpi = 600, device = ragg::agg_tiff, compression = "lzw", bg = "white")
ggsave(paste0(stem, ".png"), fig, width = width_mm, height = height_mm, units = "mm",
       dpi = 300, device = ragg::agg_png, bg = "white")

message("Figure 1 workflow review v3 written to ", output_dir)
