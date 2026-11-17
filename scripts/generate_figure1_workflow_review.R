#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ggplot2)
    library(ggplotify)
    library(grid)
})

source(file.path("external", "atlas-feature-selection-benchmark", "analysis", "R", "plotting.R"))

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) >= 1) args[[1]] else file.path("results", "figure1_workflow_review")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

pal <- c(
    blue = "#2F6DB0",
    teal = "#248F83",
    coral = "#DF624E",
    indigo = "#6B5AA8",
    gold = "#D9A329",
    green = "#4F9B67",
    ink = "#263241",
    line = "#C9D2DE",
    blue_bg = "#EEF4FB",
    teal_bg = "#ECF7F5",
    coral_bg = "#FDF0EC",
    indigo_bg = "#F2EFF9",
    gold_bg = "#FCF6E6",
    neutral_bg = "#F6F8FA"
)

items <- list()
add <- function(g) items[[length(items) + 1]] <<- g

box <- function(x, y, w, h, fill = "white", col = pal[["line"]], lwd = 0.8, r = 2.0) {
    roundrectGrob(
        x = unit(x, "npc"), y = unit(y, "npc"),
        width = unit(w, "npc"), height = unit(h, "npc"),
        r = unit(r, "mm"),
        gp = gpar(fill = fill, col = col, lwd = lwd)
    )
}

label <- function(text, x, y, size = 7, col = pal[["ink"]], face = "plain",
                  just = "centre", lineheight = 0.92) {
    textGrob(
        text, x = unit(x, "npc"), y = unit(y, "npc"), just = just,
        gp = gpar(
            fontfamily = "Helvetica", fontsize = size, col = col,
            fontface = face, lineheight = lineheight
        )
    )
}

arrow_line <- function(x0, y0, x1, y1, col = pal[["ink"]], lwd = 1.0, arrow_end = TRUE) {
    segmentsGrob(
        x0 = unit(x0, "npc"), y0 = unit(y0, "npc"),
        x1 = unit(x1, "npc"), y1 = unit(y1, "npc"),
        arrow = if (arrow_end) arrow(type = "closed", length = unit(1.5, "mm")) else NULL,
        gp = gpar(col = col, lwd = lwd, lineend = "round")
    )
}

dot <- function(x, y, fill, size = 1.25, col = "white") {
    pointsGrob(
        x = unit(x, "npc"), y = unit(y, "npc"), pch = 21,
        size = unit(size, "mm"), gp = gpar(fill = fill, col = col, lwd = 0.35)
    )
}

shape <- function(x, y, fill, col, lwd = 0.8) {
    polygonGrob(
        x = unit(x, "npc"), y = unit(y, "npc"),
        gp = gpar(fill = fill, col = col, lwd = lwd, linejoin = "round")
    )
}

# White canvas.
add(rectGrob(gp = gpar(fill = "white", col = NA)))

# Main geometry.
input <- c(x = 0.025, y = 0.10, w = 0.145, h = 0.80)
select <- c(x = 0.215, y = 0.10, w = 0.205, h = 0.80)
eval <- c(x = 0.465, y = 0.10, w = 0.335, h = 0.80)
guide <- c(x = 0.845, y = 0.10, w = 0.130, h = 0.80)

add(box(input["x"] + input["w"] / 2, 0.50, input["w"], input["h"],
        fill = "white", col = pal[["blue"]], lwd = 1.2, r = 2.6))
add(box(select["x"] + select["w"] / 2, 0.50, select["w"], select["h"],
        fill = "white", col = pal[["teal"]], lwd = 1.2, r = 2.6))
add(box(eval["x"] + eval["w"] / 2, 0.50, eval["w"], eval["h"],
        fill = "white", col = "#556273", lwd = 1.2, r = 2.6))
add(box(guide["x"] + guide["w"] / 2, 0.50, guide["w"], guide["h"],
        fill = "white", col = pal[["gold"]], lwd = 1.2, r = 2.6))

add(label("Spatial\ntranscriptomics", input["x"] + input["w"] / 2, 0.825,
          size = 8.4, col = pal[["blue"]], face = "bold"))
add(label("Feature\nselection", select["x"] + select["w"] / 2, 0.825,
          size = 8.4, col = pal[["teal"]], face = "bold"))
add(label("Evaluation", eval["x"] + eval["w"] / 2, 0.825,
          size = 8.4, col = pal[["ink"]], face = "bold"))
add(label("Guidance", guide["x"] + guide["w"] / 2, 0.825,
          size = 8.4, col = pal[["gold"]], face = "bold"))

# Input icon: two clean spatial slices.
slice_x <- c(-0.056, -0.028, 0.012, 0.050, 0.058, 0.030, -0.014, -0.054)
slice_y <- c(0.000, 0.035, 0.042, 0.025, -0.010, -0.038, -0.042, -0.025)
for (i in seq_len(2)) {
    cy <- c(0.635, 0.450)[i]
    add(shape(input["x"] + input["w"] / 2 + slice_x, cy + slice_y,
              fill = pal[["blue_bg"]], col = "#AFC3DA", lwd = 0.75))
    xs <- seq(input["x"] + 0.042, input["x"] + 0.103, length.out = 5)
    ys <- seq(cy - 0.020, cy + 0.020, length.out = 3)
    pos <- expand.grid(x = xs, y = ys)
    fills <- if (i == 1) {
        rep(c(pal[["blue"]], "#6FA8DC", pal[["teal"]]), length.out = nrow(pos))
    } else {
        rep(c(pal[["green"]], pal[["blue"]], "#80B86E"), length.out = nrow(pos))
    }
    for (j in seq_len(nrow(pos))) add(dot(pos$x[j], pos$y[j], fills[j], size = 1.10, col = NA))
}
# Feature-selection cards, with one label and one icon each.
selection_cards <- list(
    list(y = 0.670, fill = pal[["blue_bg"]], col = pal[["blue"]], text = "Expression-driven", icon = "bars"),
    list(y = 0.485, fill = pal[["teal_bg"]], col = pal[["teal"]], text = "Spatially informed", icon = "spots"),
    list(y = 0.300, fill = pal[["neutral_bg"]], col = "#697586", text = "Controls", icon = "cross")
)
for (card in selection_cards) {
    add(box(select["x"] + select["w"] / 2, card$y, select["w"] - 0.035, 0.145,
            fill = card$fill, col = card$col, lwd = 0.85, r = 2.0))
    ix <- select["x"] + 0.042
    if (card$icon == "bars") {
        heights <- c(0.030, 0.052, 0.040, 0.064)
        for (k in seq_along(heights)) {
            add(rectGrob(
                x = unit(ix - 0.020 + 0.012 * k, "npc"),
                y = unit(card$y - 0.027 + heights[k] / 2, "npc"),
                width = unit(0.008, "npc"), height = unit(heights[k], "npc"),
                gp = gpar(fill = card$col, col = NA)
            ))
        }
    } else if (card$icon == "spots") {
        coords <- expand.grid(x = ix + c(-0.020, 0, 0.020), y = card$y + c(-0.026, 0, 0.026))
        fills <- c("#A8D8D2", card$col, "#A8D8D2", card$col, card$col,
                   "#A8D8D2", "#A8D8D2", card$col, "#A8D8D2")
        for (k in seq_len(nrow(coords))) add(dot(coords$x[k], coords$y[k], fills[k], size = 1.35))
    } else {
        add(arrow_line(ix - 0.024, card$y - 0.028, ix + 0.024, card$y + 0.028,
                       col = card$col, lwd = 0.75, arrow_end = FALSE))
        add(arrow_line(ix - 0.024, card$y + 0.028, ix + 0.024, card$y - 0.028,
                       col = card$col, lwd = 0.75, arrow_end = FALSE))
        add(dot(ix, card$y, "white", size = 2.1, col = card$col))
    }
    add(label(card$text, select["x"] + 0.083, card$y, size = 6.4,
              col = card$col, face = "bold", just = "left"))
}

# Evaluation: two equal, independent tracks within one stage.
spatial_y <- 0.505
downstream_y <- 0.195
track_h <- 0.255

add(box(eval["x"] + eval["w"] / 2, spatial_y + track_h / 2, eval["w"], track_h,
        fill = pal[["coral_bg"]], col = pal[["coral"]], lwd = 1.0, r = 2.4))
add(box(eval["x"] + eval["w"] / 2, downstream_y + track_h / 2, eval["w"], track_h,
        fill = pal[["indigo_bg"]], col = pal[["indigo"]], lwd = 1.0, r = 2.4))

add(label("Spatial signal recovery", eval["x"] + 0.018, 0.715,
          size = 7.2, col = pal[["coral"]], face = "bold", just = "left"))
add(label("Downstream utility", eval["x"] + 0.018, 0.405,
          size = 7.2, col = pal[["indigo"]], face = "bold", just = "left"))

# Spatial evidence cards.
spatial_cards <- list(
    list(x = eval["x"] + 0.105, text = "Simulation truth", icon = "truth"),
    list(x = eval["x"] + 0.245, text = "Held-out validation", icon = "holdout")
)
for (card in spatial_cards) {
    add(box(card$x, 0.600, 0.120, 0.125, fill = "white", col = "#E8B8AE", lwd = 0.75, r = 1.8))
    if (card$icon == "truth") {
        coords <- expand.grid(x = card$x + c(-0.026, 0, 0.026), y = 0.620 + c(-0.016, 0, 0.016))
        fills <- ifelse(coords$x >= card$x & coords$y >= 0.615, pal[["coral"]], "#F4C7BE")
        for (k in seq_len(nrow(coords))) add(dot(coords$x[k], coords$y[k], fills[k], size = 1.25))
    } else {
        add(box(card$x - 0.018, 0.625, 0.046, 0.048, fill = "white", col = pal[["coral"]], lwd = 0.75, r = 1.0))
        add(box(card$x + 0.020, 0.610, 0.046, 0.048, fill = "#F7D8D1", col = pal[["coral"]], lwd = 0.75, r = 1.0))
        add(arrow_line(card$x - 0.003, 0.625, card$x + 0.006, 0.617,
                       col = pal[["coral"]], lwd = 0.75))
    }
    add(label(card$text, card$x, 0.560, size = 5.7, col = pal[["coral"]], face = "bold"))
}

# Downstream task cards.
downstream_cards <- list(
    list(x = eval["x"] + 0.070, text = "Integration", col = pal[["coral"]]),
    list(x = eval["x"] + 0.168, text = "Clustering", col = pal[["blue"]]),
    list(x = eval["x"] + 0.266, text = "Alignment", col = pal[["green"]])
)
for (card in downstream_cards) {
    add(box(card$x, 0.295, 0.086, 0.105, fill = "white", col = "#CEC7E4", lwd = 0.75, r = 1.7))
    add(dot(card$x - 0.020, 0.318, card$col, size = 1.30))
    add(dot(card$x, 0.332, card$col, size = 1.30))
    add(dot(card$x + 0.020, 0.313, card$col, size = 1.30))
    add(label(card$text, card$x, 0.263, size = 5.8, col = card$col, face = "bold"))
}

# Flow arrows: input -> selection -> split evaluation.
add(arrow_line(input["x"] + input["w"] + 0.007, 0.500, select["x"] - 0.010, 0.500,
               col = pal[["blue"]], lwd = 1.25))
add(arrow_line(select["x"] + select["w"] + 0.007, 0.500, eval["x"] - 0.026, 0.500,
               col = pal[["ink"]], lwd = 0.95, arrow_end = FALSE))
add(arrow_line(eval["x"] - 0.026, 0.500, eval["x"] - 0.026, 0.630,
               col = pal[["coral"]], lwd = 0.95, arrow_end = FALSE))
add(arrow_line(eval["x"] - 0.026, 0.500, eval["x"] - 0.026, 0.325,
               col = pal[["indigo"]], lwd = 0.95, arrow_end = FALSE))
add(arrow_line(eval["x"] - 0.026, 0.630, eval["x"] - 0.008, 0.630,
               col = pal[["coral"]], lwd = 0.95))
add(arrow_line(eval["x"] - 0.026, 0.325, eval["x"] - 0.008, 0.325,
               col = pal[["indigo"]], lwd = 0.95))

# Guidance is one endpoint, not a third legend or another evidence hierarchy.
add(box(guide["x"] + guide["w"] / 2, 0.500, guide["w"] - 0.024, 0.245,
        fill = pal[["gold_bg"]], col = pal[["gold"]], lwd = 0.95, r = 2.2))
add(label("Task-aware\nmethod choice", guide["x"] + guide["w"] / 2, 0.500,
          size = 7.0, col = pal[["gold"]], face = "bold"))

add(arrow_line(eval["x"] + eval["w"] + 0.007, 0.630, guide["x"] - 0.010, 0.555,
               col = pal[["coral"]], lwd = 1.0))
add(arrow_line(eval["x"] + eval["w"] + 0.007, 0.325, guide["x"] - 0.010, 0.445,
               col = pal[["indigo"]], lwd = 1.0))

figure_grob <- do.call(grobTree, items)
figure <- as.ggplot(figure_grob) +
    theme_void() +
    theme(plot.margin = margin(1.5, 1.5, 1.5, 1.5, "mm"))

stem <- file.path(output_dir, "Figure1_workflow_review")
width_mm <- 183
height_mm <- 82

ggsave(
    paste0(stem, ".svg"), figure,
    width = width_mm, height = height_mm, units = "mm",
    device = svglite::svglite, bg = "white"
)
svg_text <- readLines(paste0(stem, ".svg"), warn = FALSE)
svg_text <- gsub(" textLength='[^']+' lengthAdjust='spacingAndGlyphs'", "", svg_text)
svg_text <- gsub(
    'font-family: "Liberation Sans";',
    "font-family: Helvetica, Arial, sans-serif;",
    svg_text, fixed = TRUE
)
writeLines(svg_text, paste0(stem, ".svg"), useBytes = TRUE)

ggsave(
    paste0(stem, ".pdf"), figure,
    width = width_mm, height = height_mm, units = "mm",
    device = grDevices::cairo_pdf, family = "Arial", bg = "white"
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

message("Figure 1 review draft written to ", output_dir)
