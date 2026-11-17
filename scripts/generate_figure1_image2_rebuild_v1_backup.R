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

rr <- function(cx, cy, w, h, fill = "white", col = pal[["line"]], lwd = 0.65, r = 2.0) {
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

seg <- function(x0, y0, x1, y1, col = pal[["mid"]], lwd = 0.65,
                arrow_end = FALSE, lty = 1) {
    segmentsGrob(
        x0 = unit(X(x0), "npc"), y0 = unit(Y(y0), "npc"),
        x1 = unit(X(x1), "npc"), y1 = unit(Y(y1), "npc"),
        arrow = if (arrow_end) arrow(type = "closed", length = unit(1.35, "mm")) else NULL,
        gp = gpar(col = col, lwd = lwd, lty = lty, lineend = "round")
    )
}

point_px <- function(x, y, fill, size = 0.85, col = "white", lwd = 0.2, pch = 21) {
    pointsGrob(
        x = unit(X(x), "npc"), y = unit(Y(y), "npc"), pch = pch,
        size = unit(size, "mm"), gp = gpar(fill = fill, col = col, lwd = lwd)
    )
}

poly_px <- function(x, y, fill = "white", col = pal[["line"]], lwd = 0.6, lty = 1) {
    polygonGrob(
        x = unit(X(x), "npc"), y = unit(Y(y), "npc"),
        gp = gpar(fill = fill, col = col, lwd = lwd, lty = lty, linejoin = "round")
    )
}

stage_header <- function(number, title, x, colour, title_size = 7.4) {
    add(point_px(x, 77, colour, size = 4.5, col = colour, lwd = 0.3))
    add(txt(number, x, 77, size = 7.1, col = "white", face = "bold"))
    add(txt(title, x + 34, 78, size = title_size, col = pal[["ink"]], face = "bold", just = "left"))
}

heatmap_px <- function(x0, y0, ncol, nrow, cell_w, cell_h, colours,
                       selected_cols = integer(), outline = pal[["gold"]]) {
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
}

ellipse_spots <- function(cx, cy, rx, ry, nx, ny, colour_fun, size = 0.75,
                          tissue_fill = "#F7F5F4", tissue_col = "#BAC2C8") {
    theta <- seq(0, 2 * pi, length.out = 17)[-17]
    wobble <- 1 + c(0.00, 0.06, -0.02, 0.05, -0.04, 0.02, 0.05, -0.03,
                    0.02, -0.04, 0.05, -0.01, 0.04, -0.05, 0.03, -0.01)
    add(poly_px(cx + rx * wobble * cos(theta), cy + ry * wobble * sin(theta),
                fill = tissue_fill, col = tissue_col, lwd = 0.55))
    xs <- seq(cx - rx * 0.78, cx + rx * 0.78, length.out = nx)
    ys <- seq(cy - ry * 0.72, cy + ry * 0.72, length.out = ny)
    k <- 0L
    for (yy in ys) {
        for (xx in xs) {
            if (((xx - cx) / rx)^2 + ((yy - cy) / ry)^2 < 0.76) {
                k <- k + 1L
                add(point_px(xx, yy, colour_fun(xx, yy, k), size = size, col = NA))
            }
        }
    }
}

selected_strip <- function(cx, cy, colour) {
    heatmap_px(cx - 10, cy - 45, 2, 5, 10, 18,
               c("#DCEAE8", colour, "#94C5BE", "#E8F2F0"))
    add(txt("…", cx + 17, cy + 48, size = 5.5, col = pal[["mid"]], face = "bold"))
}

add(rectGrob(gp = gpar(fill = "white", col = NA)))

# Stage headers -------------------------------------------------------------
stage_header("1", "Spatial transcriptomics", 50, pal[["blue"]], 7.2)
stage_header("2", "Feature selection", 397, pal[["teal"]], 7.2)
stage_header("3", "Two-track evaluation", 795, pal[["coral"]], 7.0)
stage_header("4", "Objective-aware guidance", 1249, pal[["gold"]], 6.6)

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
    add(poly_px(x, y, fill = "#F8F6FA", col = "#9EA8B3", lwd = 0.55))
    # Irregular tissue field within each perspective slice.
    theta <- seq(0, 2 * pi, length.out = 19)[-19]
    rx <- 93
    ry <- 28
    wobble <- 1 + rep(c(0.02, -0.04, 0.05, -0.01, 0.03, -0.03), 3)
    add(poly_px(184 + rx * wobble * cos(theta), cy + ry * wobble * sin(theta),
                fill = "#F0E8F4", col = "#D3C7DD", lwd = 0.35))
    for (j in seq_len(35)) {
        xx <- 100 + ((j * 31 + s * 17) %% 170)
        yy <- cy - 23 + ((j * 19 + s * 11) %% 46)
        if (((xx - 184) / 92)^2 + ((yy - cy) / 28)^2 < 0.95) {
            add(point_px(xx, yy, slice_cols[[s]][1 + (j %% 3)], size = 0.72, col = NA))
        }
    }
}
add(txt("⋮", 185, 407, size = 9, col = pal[["ink"]], face = "bold"))

# One slice expanded to coordinates and matrix.
ellipse_spots(
    190, 505, 63, 58, 6, 5,
    function(x, y, k) c(pal[["blue"]], pal[["blue2"]], "#7D73B8", "#9B79C6")[1 + (k %% 4)],
    size = 1.05, tissue_fill = "#F4EDF5", tissue_col = "#9FA8B2"
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
           c("#E9F0F7", "#C8D9EB", "#93B6D8", "#5A91C6", "#356FAE"))
add(txt("Genes\n(features)", 47, 796, size = 5.1, col = pal[["ink"]], face = "bold", rot = 90))
add(txt("…", 316, 802, size = 7, col = pal[["ink"]], face = "bold"))
add(txt("…", 316, 886, size = 7, col = pal[["ink"]], face = "bold"))

# Stage 2: common matrix and three selectors -------------------------------
add(txt("Common expression matrix\nor ranked gene candidates", 403, 153,
        size = 5.3, col = pal[["ink"]], face = "bold", just = "left"))
heatmap_px(408, 196, 8, 6, 21, 22,
           c("#E6F2F0", "#B9DCD7", "#6EB2AA", "#2C8B82", "#1F766F"))
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

# Spatially informed selector.
add(txt("Spatially informed", 410, 584, size = 6.1, col = pal[["teal"]], face = "bold", just = "left"))
ellipse_spots(
    475, 663, 54, 48, 7, 6,
    function(x, y, k) {
        d <- sqrt(((x - 462) / 54)^2 + ((y - 656) / 48)^2)
        if (d < 0.42) pal[["coral2"]] else if (d < 0.67) "#F4B48F" else "#B8CFE1"
    },
    size = 0.68, tissue_fill = "#F6F2F2", tissue_col = "#ABB4BE"
)
add(seg(360, 655, 397, 655, col = pal[["teal"]], lwd = 0.65, arrow_end = TRUE))
add(seg(540, 655, 598, 655, col = pal[["ink"]], lwd = 0.65, arrow_end = TRUE))
selected_strip(620, 655, "#5FAFA6")

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

# Selector outputs feed two evidence lanes independently.
add(seg(642, 458, 660, 458, col = pal[["teal"]], lwd = 0.65))
add(seg(642, 655, 660, 655, col = pal[["teal"]], lwd = 0.65))
add(seg(642, 835, 660, 835, col = pal[["teal"]], lwd = 0.65))
add(seg(660, 330, 660, 835, col = pal[["teal"]], lwd = 0.65))
add(seg(660, 330, 692, 330, col = pal[["teal"]], lwd = 0.75, arrow_end = TRUE))
add(seg(660, 740, 692, 740, col = pal[["teal"]], lwd = 0.75, arrow_end = TRUE))

# Stage 3a: spatial signal recovery ----------------------------------------
add(rr(907, 314, 414, 305, fill = pal[["coral_bg"]], col = pal[["coral2"]], lwd = 0.8, r = 2.5))
add(rect_px(907, 162, 218, 24, fill = "white", col = NA))
add(txt("Spatial signal recovery", 907, 157, size = 6.6, col = pal[["coral"]], face = "bold"))
add(txt("Simulation truth", 791, 216, size = 5.4, col = pal[["ink"]], face = "bold"))
add(txt("Held-out validation", 998, 216, size = 5.4, col = pal[["ink"]], face = "bold"))

ellipse_spots(
    791, 334, 77, 69, 10, 9,
    function(x, y, k) {
        d <- sqrt(((x - 770) / 77)^2 + ((y - 335) / 69)^2)
        if (d < 0.30) pal[["coral"]] else if (d < 0.56) "#F2A469" else "#78A8D2"
    }, size = 0.70, tissue_fill = "#F7F4F2", tissue_col = "#B8B6B4"
)
ellipse_spots(
    999, 334, 77, 69, 10, 9,
    function(x, y, k) {
        d <- sqrt(((x - 980) / 77)^2 + ((y - 334) / 69)^2)
        if (d < 0.24) "#F4C58D" else if (d < 0.47) "#A7BEDB" else "#7BA5CF"
    }, size = 0.70, tissue_fill = "#F7F4F2", tissue_col = "#B8B6B4"
)
add(point_px(1114, 331, pal[["coral2"]], size = 2.5, col = pal[["coral2"]]))

# Stage 3b: downstream utility ---------------------------------------------
add(rr(907, 764, 414, 357, fill = pal[["indigo_bg"]], col = pal[["indigo"]], lwd = 0.8, r = 2.5))
add(rect_px(907, 586, 188, 24, fill = "white", col = NA))
add(txt("Downstream utility", 907, 580, size = 6.6, col = pal[["indigo"]], face = "bold"))

add(txt("Integration", 761, 638, size = 5.2, col = pal[["ink"]], face = "bold"))
for (i in seq_len(42)) {
    dx <- ((i * 19) %% 47) - 23
    dy <- ((i * 31) %% 73) - 36
    if ((dx / 25)^2 + (dy / 39)^2 < 0.92) {
        add(point_px(738 + dx, 724 + dy, pal[["blue"]], size = 0.72, col = NA))
        add(point_px(804 + dx, 724 + dy, "#9B6CC1", size = 0.72, col = NA))
    }
}
add(seg(771, 772, 771, 807, col = pal[["ink"]], lwd = 0.55, arrow_end = TRUE))
for (i in seq_len(30)) {
    a <- 2 * pi * i / 30
    add(point_px(771 + 33 * cos(a) * (0.5 + (i %% 7) / 12),
                 851 + 37 * sin(a) * (0.5 + (i %% 6) / 12),
                 if (i %% 2) pal[["blue"]] else "#9B6CC1", size = 0.68, col = NA))
}

add(txt("Clustering", 910, 638, size = 5.2, col = pal[["ink"]], face = "bold"))
ellipse_spots(
    910, 765, 66, 93, 9, 11,
    function(x, y, k) {
        if (x < 900 && y < 760) "#70A8C8" else if (x > 922 && y < 765) "#7CB783" else if (y > 790) "#D889A8" else "#E4B153"
    }, size = 0.58, tissue_fill = "#F5F2EF", tissue_col = "#AEB6BC"
)

add(txt("Alignment", 1041, 638, size = 5.2, col = pal[["ink"]], face = "bold"))
theta <- seq(0, 2 * pi, length.out = 13)[-13]
add(poly_px(1024 + 52 * cos(theta), 730 + 56 * sin(theta), fill = "#DDE8F5", col = pal[["blue"]], lwd = 0.55))
add(poly_px(1065 + 52 * cos(theta), 822 + 56 * sin(theta), fill = "#F3DCE6", col = "#B86B91", lwd = 0.55))
for (k in -2:2) {
    add(seg(1004 + k * 12, 720 + k * 5, 1045 + k * 12, 812 + k * 5,
            col = pal[["ink"]], lwd = 0.42, lty = 2))
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
            fill = pal[["gold_bg"]], col = pal[["gold"]], lwd = 0.8))
add(txt("Analysis\ngoal", 1305, 523, size = 6.2, col = pal[["ink"]], face = "bold"))

add(txt("Spatial discovery", 1264, 174, size = 5.8, col = pal[["ink"]], face = "bold"))
add(rr(1305, 273, 154, 154, fill = "white", col = pal[["gold"]], lwd = 0.65, r = 0.8))
ellipse_spots(
    1305, 273, 64, 61, 8, 8,
    function(x, y, k) {
        d <- sqrt(((x - 1283) / 64)^2 + ((y - 276) / 61)^2)
        if (d < 0.35) pal[["coral"]] else if (d < 0.62) "#F5A15B" else "#D8BFC2"
    }, size = 0.62, tissue_fill = "#F8F4F2", tissue_col = "#C0B7B5"
)
add(seg(1305, 438, 1305, 358, col = pal[["ink"]], lwd = 0.6, arrow_end = TRUE))

add(rr(1305, 795, 170, 160, fill = "white", col = pal[["gold"]], lwd = 0.65, r = 0.8))
add(poly_px(c(1260, 1312, 1300, 1250), c(770, 750, 818, 835),
            fill = "#DCE8F5", col = pal[["blue"]], lwd = 0.55))
add(poly_px(c(1320, 1370, 1362, 1310), c(758, 780, 842, 820),
            fill = "#F2DDE5", col = "#B7698D", lwd = 0.55))
for (k in -1:1) {
    add(seg(1280 + k * 13, 782 + k * 4, 1338 + k * 13, 798 + k * 4,
            col = pal[["ink"]], lwd = 0.42, lty = 2))
}
add(seg(1305, 608, 1305, 708, col = pal[["ink"]], lwd = 0.6, arrow_end = TRUE))
add(rect_px(1255, 684, 172, 22, fill = "white", col = NA))
add(txt("Integrated analysis", 1255, 684, size = 5.8, col = pal[["ink"]], face = "bold"))

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

message("Figure 1 SVG rebuild written to ", output_dir)
