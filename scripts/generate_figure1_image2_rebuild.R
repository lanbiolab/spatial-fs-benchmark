#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ggplot2)
    library(ggplotify)
    library(grid)
})

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) >= 1) args[[1]] else file.path("results", "figure1_image2_rebuild")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

canvas_w <- 1536
canvas_h <- 1024
X <- function(x) x / canvas_w
Y <- function(y) 1 - y / canvas_h
W <- function(w) w / canvas_w
H <- function(h) h / canvas_h

pal <- c(
    blue = "#356FAE",
    blue2 = "#6F9ED0",
    teal = "#1F766F",
    coral = "#C45140",
    coral2 = "#EF8B73",
    indigo = "#5E63A8",
    green = "#3D8A60",
    gold = "#B98606",
    ink = "#202A35",
    mid = "#667487",
    line = "#B7C2CE",
    blue_bg = "#EEF4FA",
    teal_bg = "#EDF7F5",
    coral_bg = "#FFF5F2",
    indigo_bg = "#F3F3FB",
    gold_bg = "#FFF9E9",
    neutral_bg = "#F6F8FA"
)

items <- list()
add <- function(x) items[[length(items) + 1L]] <<- x

rr <- function(cx, cy, w, h, fill = "white", col = pal[["line"]], lwd = 0.55, r = 2.0) {
    roundrectGrob(
        x = unit(X(cx), "npc"), y = unit(Y(cy), "npc"),
        width = unit(W(w), "npc"), height = unit(H(h), "npc"),
        r = unit(r, "mm"), gp = gpar(fill = fill, col = col, lwd = lwd)
    )
}

rect_px <- function(cx, cy, w, h, fill = "white", col = NA, lwd = 0.3) {
    rectGrob(
        x = unit(X(cx), "npc"), y = unit(Y(cy), "npc"),
        width = unit(W(w), "npc"), height = unit(H(h), "npc"),
        gp = gpar(fill = fill, col = col, lwd = lwd)
    )
}

txt <- function(text, x, y, size = 6.0, col = pal[["ink"]], face = "plain",
                just = "centre", rot = 0, lineheight = 0.92) {
    textGrob(
        text, x = unit(X(x), "npc"), y = unit(Y(y), "npc"),
        just = just, rot = rot,
        gp = gpar(
            fontfamily = "Arial", fontsize = size, col = col,
            fontface = face, lineheight = lineheight
        )
    )
}

seg <- function(x0, y0, x1, y1, col = pal[["mid"]], lwd = 0.55,
                arrow_end = FALSE, lty = 1) {
    segmentsGrob(
        x0 = unit(X(x0), "npc"), y0 = unit(Y(y0), "npc"),
        x1 = unit(X(x1), "npc"), y1 = unit(Y(y1), "npc"),
        arrow = if (arrow_end) arrow(type = "closed", length = unit(1.12, "mm")) else NULL,
        gp = gpar(col = col, lwd = lwd, lty = lty, lineend = "round")
    )
}

point_px <- function(x, y, fill, size = 0.85, col = "white", lwd = 0.2, pch = 21) {
    pointsGrob(
        x = unit(X(x), "npc"), y = unit(Y(y), "npc"), pch = pch,
        size = unit(size, "mm"), gp = gpar(fill = fill, col = col, lwd = lwd)
    )
}

poly_px <- function(x, y, fill = "white", col = pal[["line"]], lwd = 0.45, lty = 1) {
    polygonGrob(
        x = unit(X(x), "npc"), y = unit(Y(y), "npc"),
        gp = gpar(fill = fill, col = col, lwd = lwd, lty = lty, linejoin = "round")
    )
}

stage_header <- function(number, title = NULL, x, colour, title_size = 7.4) {
    add(point_px(x, 77, colour, size = 4.25, col = colour, lwd = 0.25))
    add(txt(number, x, 77, size = 6.8, col = "white", face = "bold"))
    if (!is.null(title) && nzchar(title)) {
        add(txt(title, x + 34, 78, size = title_size, col = pal[["ink"]], face = "bold", just = "left"))
    }
}

heatmap_px <- function(x0, y0, ncol, nrow, cell_w, cell_h, colours,
                       selected_cols = integer(), outline = pal[["gold"]], border = NA) {
    for (j in seq_len(ncol)) {
        for (i in seq_len(nrow)) {
            idx <- 1L + ((i * 2L + j * 3L + i * j) %% length(colours))
            add(rect_px(
                x0 + (j - 0.5) * cell_w, y0 + (i - 0.5) * cell_h,
                cell_w * 0.94, cell_h * 0.94, fill = colours[idx], col = "white", lwd = 0.15
            ))
        }
        if (j %in% selected_cols) {
            add(rr(
                x0 + (j - 0.5) * cell_w, y0 + nrow * cell_h / 2,
                cell_w * 1.08, nrow * cell_h * 1.02,
                fill = NA, col = outline, lwd = 0.7, r = 0.5
            ))
        }
    }
    if (!is.na(border)) {
        add(rect_px(x0 + ncol * cell_w / 2, y0 + nrow * cell_h / 2,
                    ncol * cell_w, nrow * cell_h, fill = NA, col = border, lwd = 0.42))
    }
}

cell_px <- function(cx, cy, rx, ry, fill, col = "#FFFFFF", lwd = 0.16,
                    phase = 0, alpha = 1) {
    theta <- seq(0, 2 * pi, length.out = 7)[-7] + phase
    wobble <- c(1.00, 0.88, 1.08, 0.94, 1.04, 0.90)
    add(poly_px(cx + rx * wobble * cos(theta), cy + ry * wobble * sin(theta),
                fill = grDevices::adjustcolor(fill, alpha.f = alpha),
                col = col, lwd = lwd))
}

ellipse_spots <- function(cx, cy, rx, ry, nx, ny, colour_fun, size = 0.75,
                          tissue_fill = "#F7F5F4", tissue_col = "#BAC2C8",
                          cell_border = "#FFFFFF", cell_alpha = 0.96,
                          underlay = TRUE, glyph = c("cell", "spot"), variant = 0) {
    glyph <- match.arg(glyph)
    theta <- seq(0, 2 * pi, length.out = 25)[-25]
    wobble <- 1 + 0.038 * sin(theta * (3 + (variant %% 3)) + variant * 0.71) +
        0.024 * cos(theta * (5 + (variant %% 2)) - variant * 0.43)
    outline_x <- cx + rx * wobble * cos(theta)
    outline_y <- cy + ry * wobble * sin(theta)
    add(poly_px(outline_x, outline_y, fill = tissue_fill, col = tissue_col, lwd = 0.42))

    if (underlay) {
        # Dense histology-like cells remain visible below the analysis overlay.
        ux <- seq(cx - rx * 0.82, cx + rx * 0.82, length.out = 13)
        uy <- seq(cy - ry * 0.76, cy + ry * 0.76, length.out = 10)
        q <- 0L
        for (iy in seq_along(uy)) {
            for (ix in seq_along(ux)) {
                xx <- ux[ix] + if (iy %% 2) rx * 0.055 else 0
                yy <- uy[iy] + sin(ix * 1.7 + iy * 0.8) * ry * 0.018
                if (((xx - cx) / rx)^2 + ((yy - cy) / ry)^2 < 0.82) {
                    q <- q + 1L
                    cell_px(xx, yy, rx * 0.055, ry * 0.060,
                            c("#E9DCE5", "#E1D4E1", "#E1E7EB", "#EEE5EA")[1 + (q %% 4)],
                            col = "#D7CCD4", lwd = 0.08,
                            phase = (q %% 6) * pi / 18, alpha = 0.72)
                    if (q %% 4 == 0) add(point_px(xx, yy, "#7E7082", size = 0.14, col = NA))
                }
            }
        }
    }

    xs <- seq(cx - rx * 0.79, cx + rx * 0.79, length.out = nx)
    ys <- seq(cy - ry * 0.73, cy + ry * 0.73, length.out = ny)
    dx <- if (nx > 1) diff(range(xs)) / (nx - 1) else rx
    dy <- if (ny > 1) diff(range(ys)) / (ny - 1) else ry
    k <- 0L
    for (iy in seq_along(ys)) {
        for (ix in seq_along(xs)) {
            xx <- xs[ix] + if (iy %% 2) dx * 0.20 else -dx * 0.20
            yy <- ys[iy] + sin((ix * 1.7 + iy) * 0.9) * dy * 0.05
            if (((xx - cx) / rx)^2 + ((yy - cy) / ry)^2 < 0.76) {
                k <- k + 1L
                if (glyph == "cell") {
                    cell_px(xx, yy, min(dx * 0.48, size * 6.6), min(dy * 0.48, size * 6.0),
                            colour_fun(xx, yy, k), col = cell_border, lwd = 0.14,
                            phase = (k %% 6) * pi / 18, alpha = cell_alpha)
                    if (size >= 0.95 && k %% 3 == 0) {
                        add(point_px(xx, yy, "#6E6682", size = size * 0.24, col = NA))
                    }
                } else {
                    add(point_px(xx, yy, colour_fun(xx, yy, k), size = size,
                                 col = "#F8FAFC", lwd = 0.18))
                }
            }
        }
    }
    add(poly_px(outline_x, outline_y, fill = NA, col = tissue_col, lwd = 0.42))
}

curve_seg <- function(x0, y0, x1, y1, curvature = 0.25, col = pal[["mid"]],
                      lwd = 0.55, arrow_end = TRUE, lty = 1) {
    curveGrob(
        x1 = unit(X(x0), "npc"), y1 = unit(Y(y0), "npc"),
        x2 = unit(X(x1), "npc"), y2 = unit(Y(y1), "npc"),
        curvature = curvature,
        arrow = if (arrow_end) arrow(type = "closed", length = unit(1.15, "mm")) else NULL,
        gp = gpar(col = col, lwd = lwd, lty = lty, lineend = "round")
    )
}

brace_right <- function(x, y0, y1, col = pal[["teal"]], lwd = 0.6) {
    ym <- (y0 + y1) / 2
    xs <- c(x - 14, x - 4, x - 4, x, x - 4, x - 4, x - 14)
    ys <- c(y0, y0 + 14, ym - 18, ym, ym + 18, y1 - 14, y1)
    xsplineGrob(
        x = unit(X(xs), "npc"), y = unit(Y(ys), "npc"),
        shape = c(0, 1, 1, 0, 1, 1, 0), open = TRUE,
        gp = gpar(col = col, lwd = lwd, fill = NA, lineend = "round")
    )
}

cluster_cloud <- function(cx, cy, rx, ry, fills, nx = 7, ny = 9,
                          cell_rx = 4.2, cell_ry = 4.0) {
    xs <- seq(cx - rx * 0.78, cx + rx * 0.78, length.out = nx)
    ys <- seq(cy - ry * 0.78, cy + ry * 0.78, length.out = ny)
    k <- 0L
    for (iy in seq_along(ys)) {
        for (ix in seq_along(xs)) {
            xx <- xs[ix] + if (iy %% 2) cell_rx * 0.65 else 0
            yy <- ys[iy] + sin(ix * 1.3 + iy * 0.8) * 1.1
            if (((xx - cx) / rx)^2 + ((yy - cy) / ry)^2 < 0.82) {
                k <- k + 1L
                fill <- fills[1 + ((k + ix + 2 * iy) %% length(fills))]
                cell_px(xx, yy, cell_rx, cell_ry, fill, col = "white", lwd = 0.12,
                        phase = (k %% 6) * pi / 18, alpha = 0.98)
                if (k %% 5 == 0) add(point_px(xx, yy, "#665F78", size = 0.18, col = NA))
            }
        }
    }
}

selected_strip <- function(cx, cy, colour) {
    heatmap_px(cx - 10, cy - 45, 2, 5, 10, 18,
               c("#DCEAE8", colour, "#94C5BE", "#E8F2F0"), border = pal[["teal"]])
    add(txt("…", cx + 17, cy + 48, size = 5.5, col = pal[["mid"]], face = "bold"))
}

add(rectGrob(gp = gpar(fill = "white", col = NA)))

# Stage headers -------------------------------------------------------------
stage_header("1", "Spatial transcriptomics", 50, pal[["blue"]], 7.2)
stage_header("2", "Feature selection", 397, pal[["teal"]], 7.2)
stage_header("3", NULL, 795, pal[["coral"]])
stage_header("4", NULL, 1249, pal[["gold"]])

# Stage 1: layered spatial inputs ------------------------------------------
slice_centres <- c(185, 240, 295, 350)
slice_cols <- list(
    c("#C7B0DD", "#8DB4DF", "#D8C4E7"),
    c("#AFA0D3", "#5793CA", "#D6C0E8"),
    c("#6E9FD0", "#7C76BA", "#CFC4E5"),
    c("#6A9DCB", "#4D87C0", "#C8BFE2")
)
for (s in seq_along(slice_centres)) {
    cy <- slice_centres[s]
    x <- c(58, 318, 282, 52)
    y <- c(cy - 38, cy - 38, cy + 38, cy + 38)
    add(poly_px(x, y, fill = grDevices::adjustcolor("#FBFAFC", alpha.f = 0.92),
                col = "#8F9AA5", lwd = 0.42))
    theta <- seq(0, 2 * pi, length.out = 25)[-25]
    rx <- 94
    ry <- 28
    wobble <- 1 + rep(c(0.02, -0.04, 0.05, -0.01, 0.03, -0.03), 4)
    ox <- 184 + rx * wobble * cos(theta)
    oy <- cy + ry * wobble * sin(theta)
    add(poly_px(ox, oy, fill = "#F1E8F2", col = "#C8B9CB", lwd = 0.24))
    gx <- seq(101, 267, length.out = 16)
    gy <- seq(cy - 20, cy + 20, length.out = 7)
    k <- 0L
    for (iy in seq_along(gy)) {
        for (ix in seq_along(gx)) {
            xx <- gx[ix] + if (iy %% 2) 5 else 0
            yy <- gy[iy] + sin(ix * 1.4 + s) * 1.4
            if (((xx - 184) / 92)^2 + ((yy - cy) / 27)^2 < 0.91) {
                k <- k + 1L
                base_col <- c("#E0CDE5", "#D6BDDF", "#CBB1DA", "#E8D8EA")[1 + ((k + s) %% 4)]
                cell_px(xx, yy, 4.9, 3.8, base_col, col = "#F8F5FA", lwd = 0.08,
                        phase = (k %% 6) * pi / 18, alpha = 0.94)
                if (k %% 4 == 0) add(point_px(xx, yy, "#7D668A", size = 0.23, col = NA))
                if (s >= 2 && ((xx - 172) / 76)^2 + ((yy - cy) / 22)^2 < 0.74 && k %% 3 != 0) {
                    add(point_px(xx, yy, slice_cols[[s]][1 + (k %% 3)],
                                 size = 0.70, col = "#F8FAFC", lwd = 0.11))
                }
            }
        }
    }
    add(poly_px(ox, oy, fill = NA, col = "#C4B6C7", lwd = 0.25))
}
add(txt("⋮", 185, 407, size = 9, col = pal[["ink"]], face = "bold"))

# One slice expanded to coordinates and matrix.
ellipse_spots(
    190, 505, 63, 58, 6, 5,
    function(x, y, k) c(pal[["blue"]], pal[["blue2"]], "#7D73B8", "#9B79C6")[1 + (k %% 4)],
    size = 1.28, tissue_fill = "#F4EDF5", tissue_col = "#9FA8B2", glyph = "spot", variant = 1
)
add(seg(98, 558, 98, 470, col = pal[["ink"]], lwd = 0.45, arrow_end = TRUE))
add(seg(90, 566, 280, 566, col = pal[["ink"]], lwd = 0.45, arrow_end = TRUE))
add(txt("y", 84, 510, size = 4.7, col = pal[["ink"]]))
add(txt("x", 187, 580, size = 4.7, col = pal[["ink"]]))
add(txt("Spots (spatial locations)", 107, 646, size = 5.2, col = pal[["ink"]], face = "bold", just = "left"))

spot_x <- seq(88, 292, length.out = 11)
spot_cols <- c("#356FAE", "#4D86C0", "#6A9ED0", "#8BB4DA", "#B2D0E8",
               "#D8E4EF", "#C9C8E5", "#B6A8D8", "#9C83C8", "#825EB4", "#6C4DA7")
for (i in seq_along(spot_x)) add(point_px(spot_x[i], 686, spot_cols[i], size = 1.55, col = "white", lwd = 0.25))
add(seg(148, 558, 88, 676, col = "#8C98A5", lwd = 0.45, lty = 2))
add(seg(229, 558, 292, 676, col = "#8C98A5", lwd = 0.45, lty = 2))

heatmap_px(75, 704, 10, 8, 23, 23,
           c("#E9F0F7", "#C8D9EB", "#93B6D8", "#5A91C6", "#356FAE"),
           border = pal[["blue"]])
add(txt("Genes\n(features)", 47, 796, size = 5.1, col = pal[["ink"]], face = "bold", rot = 90))
add(txt("…", 316, 802, size = 7, col = pal[["ink"]], face = "bold"))
add(txt("…", 316, 886, size = 7, col = pal[["ink"]], face = "bold"))

# Stage 2: common matrix and three selectors -------------------------------
add(txt("Common expression matrix\nor ranked gene candidates", 403, 153,
        size = 5.3, col = pal[["ink"]], face = "bold", just = "left"))
heatmap_px(408, 196, 8, 6, 21, 22,
           c("#E6F2F0", "#B9DCD7", "#6EB2AA", "#2C8B82", "#1F766F"),
           border = pal[["teal"]])
add(txt("…", 585, 250, size = 7, col = pal[["ink"]], face = "bold"))
add(txt("…", 585, 316, size = 7, col = pal[["ink"]], face = "bold"))

# Shared selector trunk.
add(seg(304, 800, 360, 800, col = pal[["ink"]], lwd = 0.65))
add(seg(360, 800, 360, 331, col = pal[["ink"]], lwd = 0.65))
add(seg(360, 331, 398, 331, col = pal[["ink"]], lwd = 0.65, arrow_end = TRUE))

# Expression-driven selector.
add(txt("Expression-driven", 410, 401, size = 6.1, col = pal[["teal"]], face = "bold", just = "left"))
add(seg(412, 506, 412, 432, col = pal[["ink"]], lwd = 0.45))
add(seg(412, 506, 535, 506, col = pal[["ink"]], lwd = 0.45, arrow_end = TRUE))
bar_h <- c(43, 65, 50, 36, 25)
for (i in seq_along(bar_h)) {
    add(rect_px(429 + (i - 1) * 20, 506 - bar_h[i] / 2, 13, bar_h[i],
                fill = c("#6FB6AD", "#4E9F96", "#73B5AC", "#8AC3BC", "#A4D0CB")[i],
                col = pal[["teal"]], lwd = 0.35))
}
add(seg(430, 445, 520, 486, col = pal[["teal"]], lwd = 0.45, lty = 2))
add(seg(360, 464, 397, 464, col = pal[["teal"]], lwd = 0.65, arrow_end = TRUE))
add(seg(545, 464, 598, 464, col = pal[["ink"]], lwd = 0.65, arrow_end = TRUE))
selected_strip(620, 458, "#54A99F")
add(brace_right(654, 413, 503, col = pal[["teal"]], lwd = 0.55))

# Spatially informed selector.
add(txt("Spatially informed", 410, 584, size = 6.1, col = pal[["teal"]], face = "bold", just = "left"))
ellipse_spots(
    475, 663, 54, 48, 9, 8,
    function(x, y, k) {
        d <- sqrt(((x - 462) / 54)^2 + ((y - 656) / 48)^2)
        if (d < 0.42) pal[["coral2"]] else if (d < 0.67) "#F4B48F" else "#B8CFE1"
    },
    size = 0.72, tissue_fill = "#F6F2F2", tissue_col = "#ABB4BE", variant = 2
)
add(seg(360, 655, 397, 655, col = pal[["teal"]], lwd = 0.65, arrow_end = TRUE))
add(seg(540, 655, 598, 655, col = pal[["ink"]], lwd = 0.65, arrow_end = TRUE))
selected_strip(620, 655, "#5FAFA6")
add(brace_right(654, 610, 700, col = pal[["teal"]], lwd = 0.55))

# Controls.
add(txt("Controls", 442, 770, size = 6.1, col = pal[["teal"]], face = "bold", just = "left"))
for (i in seq_len(5)) {
    for (j in seq_len(4)) {
        fill <- if (j %in% c(2, 3) && i %in% c(2, 3, 4)) "#85BDB6" else "white"
        add(rect_px(433 + (i - 1) * 21, 802 + (j - 1) * 21, 20, 20,
                    fill = fill, col = pal[["ink"]], lwd = 0.5))
    }
}
add(seg(360, 815, 397, 815, col = pal[["teal"]], lwd = 0.65, arrow_end = TRUE))
add(seg(540, 835, 598, 835, col = pal[["ink"]], lwd = 0.65, arrow_end = TRUE))
selected_strip(620, 835, "#5FAFA6")
add(brace_right(654, 790, 880, col = pal[["teal"]], lwd = 0.55))

# Selector outputs feed two evidence lanes independently.
add(seg(642, 458, 660, 458, col = pal[["teal"]], lwd = 0.65))
add(seg(642, 655, 660, 655, col = pal[["teal"]], lwd = 0.65))
add(seg(642, 835, 660, 835, col = pal[["teal"]], lwd = 0.65))
add(seg(660, 330, 660, 835, col = pal[["teal"]], lwd = 0.65))
add(seg(660, 330, 692, 330, col = pal[["teal"]], lwd = 0.75, arrow_end = TRUE))
add(seg(660, 740, 692, 740, col = pal[["teal"]], lwd = 0.75, arrow_end = TRUE))

# Stage 3a: spatial signal recovery ----------------------------------------
add(rr(907, 314, 414, 305, fill = pal[["coral_bg"]], col = pal[["coral2"]], lwd = 0.62, r = 2.5))
add(rect_px(907, 162, 218, 24, fill = "white", col = NA))
add(txt("Spatial signal recovery", 907, 157, size = 6.6, col = pal[["coral"]], face = "bold"))
add(txt("Simulation truth", 791, 216, size = 5.4, col = pal[["ink"]], face = "bold"))
add(txt("Held-out validation", 998, 216, size = 5.4, col = pal[["ink"]], face = "bold"))

ellipse_spots(
    791, 334, 77, 69, 10, 9,
    function(x, y, k) {
        d <- sqrt(((x - 770) / 77)^2 + ((y - 335) / 69)^2)
        if (d < 0.30) pal[["coral"]] else if (d < 0.56) "#F2A469" else "#78A8D2"
    }, size = 0.70, tissue_fill = "#F7F4F2", tissue_col = "#B8B6B4", variant = 3
)
ellipse_spots(
    999, 334, 77, 69, 10, 9,
    function(x, y, k) {
        d <- sqrt(((x - 980) / 77)^2 + ((y - 334) / 69)^2)
        if (d < 0.24) "#F4C58D" else if (d < 0.47) "#A7BEDB" else "#7BA5CF"
    }, size = 0.70, tissue_fill = "#F7F4F2", tissue_col = "#B8B6B4", variant = 4
)
add(point_px(1114, 331, pal[["coral2"]], size = 2.5, col = pal[["coral2"]]))

# Stage 3b: downstream utility ---------------------------------------------
add(rr(907, 764, 414, 357, fill = pal[["indigo_bg"]], col = pal[["indigo"]], lwd = 0.62, r = 2.5))
add(rect_px(907, 586, 188, 24, fill = "white", col = NA))
add(txt("Downstream utility", 907, 580, size = 6.6, col = pal[["indigo"]], face = "bold"))

add(txt("Integration", 761, 638, size = 5.2, col = pal[["ink"]], face = "bold"))
cluster_cloud(738, 723, 28, 45, c("#2F70B5", "#4E89C7", "#79A8D3"),
              nx = 6, ny = 9, cell_rx = 4.1, cell_ry = 3.8)
cluster_cloud(805, 723, 28, 45, c("#8B5BB5", "#A875C7", "#C29AD8"),
              nx = 6, ny = 9, cell_rx = 4.1, cell_ry = 3.8)
add(seg(771, 774, 771, 807, col = pal[["ink"]], lwd = 0.50, arrow_end = TRUE))
cluster_cloud(771, 854, 42, 38,
              c("#2F70B5", "#5D91C8", "#8B5BB5", "#B487CD"),
              nx = 9, ny = 8, cell_rx = 3.9, cell_ry = 3.6)

add(txt("Clustering", 910, 638, size = 5.2, col = pal[["ink"]], face = "bold"))
ellipse_spots(
    910, 765, 66, 93, 9, 11,
    function(x, y, k) {
        if (x < 900 && y < 760) "#70A8C8" else if (x > 922 && y < 765) "#7CB783" else if (y > 790) "#D889A8" else "#E4B153"
    }, size = 0.66, tissue_fill = "#F5F2EF", tissue_col = "#AEB6BC",
    underlay = FALSE, variant = 5
)

add(txt("Alignment", 1041, 638, size = 5.2, col = pal[["ink"]], face = "bold"))
ellipse_spots(
    1023, 714, 45, 47, 7, 7,
    function(x, y, k) c("#BFD2E8", "#9DBCE0", "#739FD0", "#D2DCEE")[1 + (k %% 4)],
    size = 0.62, tissue_fill = "#E8EFF7", tissue_col = pal[["blue"]],
    cell_border = "#F7FAFD", underlay = TRUE, variant = 6
)
ellipse_spots(
    1064, 846, 45, 47, 7, 7,
    function(x, y, k) c("#EBC8D8", "#DDA8C2", "#C982A5", "#F2DCE6")[1 + (k %% 4)],
    size = 0.62, tissue_fill = "#F6E8EE", tissue_col = "#B86B91",
    cell_border = "#FFF9FB", underlay = TRUE, variant = 7
)
for (k in -2:2) {
    add(curve_seg(1008 + k * 10, 708 + k * 4, 1047 + k * 10, 839 + k * 4,
                  curvature = -0.08, col = "#586473", lwd = 0.36,
                  arrow_end = FALSE, lty = 2))
}
add(point_px(1114, 741, pal[["indigo"]], size = 2.5, col = pal[["indigo"]]))

# Both evidence lanes merge before the objective decision.
add(seg(1120, 331, 1192, 331, col = pal[["ink"]], lwd = 0.65))
add(seg(1192, 331, 1192, 523, col = pal[["ink"]], lwd = 0.65))
add(seg(1120, 741, 1192, 741, col = pal[["ink"]], lwd = 0.65))
add(seg(1192, 741, 1192, 523, col = pal[["ink"]], lwd = 0.65))
add(point_px(1192, 523, pal[["ink"]], size = 1.8, col = "white", lwd = 0.25))
add(seg(1202, 523, 1222, 523, col = pal[["ink"]], lwd = 0.75, arrow_end = TRUE))

# Stage 4: objective-aware guidance ----------------------------------------
add(poly_px(c(1305, 1386, 1305, 1224), c(438, 523, 608, 523),
            fill = pal[["gold_bg"]], col = pal[["gold"]], lwd = 0.62))
add(txt("Analysis\ngoal", 1305, 523, size = 6.2, col = pal[["ink"]], face = "bold"))

add(txt("Spatial discovery", 1264, 174, size = 5.8, col = pal[["ink"]], face = "bold"))
add(rr(1305, 273, 154, 154, fill = "white", col = pal[["gold"]], lwd = 0.52, r = 0.8))
ellipse_spots(
    1305, 273, 64, 61, 8, 8,
    function(x, y, k) {
        d <- sqrt(((x - 1283) / 64)^2 + ((y - 276) / 61)^2)
        if (d < 0.35) pal[["coral"]] else if (d < 0.62) "#F5A15B" else "#D8BFC2"
    }, size = 0.62, tissue_fill = "#F8F4F2", tissue_col = "#C0B7B5", variant = 8
)
add(seg(1305, 438, 1305, 358, col = pal[["ink"]], lwd = 0.6, arrow_end = TRUE))

add(rr(1305, 795, 170, 160, fill = "white", col = pal[["gold"]], lwd = 0.52, r = 0.8))
ellipse_spots(
    1275, 797, 34, 46, 6, 7,
    function(x, y, k) c("#B7CCE5", "#8EAFD6", "#6F9ED0", "#D5E0EE")[1 + (k %% 4)],
    size = 0.52, tissue_fill = "#E9F0F8", tissue_col = pal[["blue"]],
    cell_border = "#F8FAFD", underlay = FALSE, variant = 9
)
ellipse_spots(
    1338, 797, 34, 46, 6, 7,
    function(x, y, k) c("#E7BCD0", "#D896B5", "#C9779E", "#F0D8E3")[1 + (k %% 4)],
    size = 0.52, tissue_fill = "#F6E8EF", tissue_col = "#B7698D",
    cell_border = "#FFF9FC", underlay = FALSE, variant = 10
)
for (k in -1:1) {
    add(curve_seg(1268 + k * 8, 785 + k * 12, 1330 + k * 8, 785 + k * 12,
                  curvature = if (k == 0) -0.10 else 0.06,
                  col = "#586473", lwd = 0.34, arrow_end = FALSE, lty = 2))
}
add(curve_seg(1276, 735, 1337, 735, curvature = 0.46,
              col = pal[["ink"]], lwd = 0.45, arrow_end = TRUE))
add(curve_seg(1337, 859, 1276, 859, curvature = 0.46,
              col = pal[["ink"]], lwd = 0.45, arrow_end = TRUE))
add(seg(1305, 608, 1305, 708, col = pal[["ink"]], lwd = 0.6, arrow_end = TRUE))
add(rect_px(1305, 684, 170, 22, fill = "white", col = NA))
add(txt("Integrated analysis", 1305, 684, size = 5.8, col = pal[["ink"]], face = "bold"))

# Target icon and overall choice.
add(seg(1386, 523, 1414, 523, col = pal[["ink"]], lwd = 0.65, arrow_end = TRUE))
add(pointsGrob(x = unit(X(1460), "npc"), y = unit(Y(523), "npc"), pch = 1,
               size = unit(8.0, "mm"), gp = gpar(col = pal[["gold"]], lwd = 1.0)))
add(pointsGrob(x = unit(X(1460), "npc"), y = unit(Y(523), "npc"), pch = 1,
               size = unit(4.8, "mm"), gp = gpar(col = pal[["gold"]], lwd = 0.8)))
add(point_px(1460, 523, pal[["gold"]], size = 1.6, col = pal[["gold"]]))
add(seg(1460, 523, 1490, 493, col = pal[["gold"]], lwd = 0.75, arrow_end = TRUE))
add(txt("Objective-aware\nmethod choice", 1460, 590, size = 5.3, col = pal[["ink"]], face = "bold"))

figure_grob <- do.call(grobTree, items)
figure <- as.ggplot(figure_grob) +
    theme_void() +
    theme(plot.margin = margin(1.0, 1.0, 1.0, 1.0, "mm"))

stem <- file.path(output_dir, "Figure1")
width_mm <- 183
height_mm <- 122

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
ggsave(file.path(output_dir, "Figure1_print_qa.png"), figure,
       width = width_mm, height = height_mm, units = "mm",
       dpi = 150, device = ragg::agg_png, bg = "white")

message("Figure 1 SVG rebuild written to ", output_dir)
