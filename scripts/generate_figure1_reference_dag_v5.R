#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ggplot2)
    library(ggplotify)
    library(grid)
})

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) >= 1) args[[1]] else file.path("results", "figure1_reference_dag_v5")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

canvas_w <- 2300
canvas_h <- 720
X <- function(x) x / canvas_w
Y <- function(y) 1 - y / canvas_h
W <- function(w) w / canvas_w
H <- function(h) h / canvas_h

pal <- c(
    ink = "#26333F",
    arrow = "#44546A",
    node = "#62666A",
    node_border = "#454A4E",
    teal = "#86CDBF",
    teal_dark = "#177C73",
    blue = "#90BBDD",
    blue_dark = "#2F6FB1",
    green = "#8FD0AC",
    green_dark = "#2E8058",
    gold = "#F1BE72",
    gold_dark = "#A97200",
    coral = "#E78B78",
    indigo = "#8776BA",
    pale_gold = "#FFF7E7",
    white = "#FFFFFF"
)

grobs <- list()
add <- function(g) grobs[[length(grobs) + 1L]] <<- g

rr <- function(cx, cy, w, h, fill = "white", col = pal[["node_border"]],
               lwd = 0.65, radius = 2.0) {
    roundrectGrob(
        x = unit(X(cx), "npc"), y = unit(Y(cy), "npc"),
        width = unit(W(w), "npc"), height = unit(H(h), "npc"),
        r = unit(radius, "mm"),
        gp = gpar(fill = fill, col = col, lwd = lwd)
    )
}

txt <- function(label, x, y, size = 6.0, col = pal[["ink"]], face = "plain",
                just = "centre", lineheight = 0.92) {
    textGrob(
        label,
        x = unit(X(x), "npc"), y = unit(Y(y), "npc"), just = just,
        gp = gpar(fontfamily = "Helvetica", fontsize = size, col = col,
                  fontface = face, lineheight = lineheight)
    )
}

seg <- function(x0, y0, x1, y1, col = pal[["arrow"]], lwd = 0.9,
                arrow_end = FALSE, lty = 1) {
    segmentsGrob(
        x0 = unit(X(x0), "npc"), y0 = unit(Y(y0), "npc"),
        x1 = unit(X(x1), "npc"), y1 = unit(Y(y1), "npc"),
        arrow = if (arrow_end) arrow(type = "closed", length = unit(1.5, "mm")) else NULL,
        gp = gpar(col = col, lwd = lwd, lty = lty, lineend = "round")
    )
}

pt <- function(x, y, fill = pal[["arrow"]], size = 1.6, col = "white", lwd = 0.25) {
    pointsGrob(
        x = unit(X(x), "npc"), y = unit(Y(y), "npc"), pch = 21,
        size = unit(size, "mm"), gp = gpar(fill = fill, col = col, lwd = lwd)
    )
}

group_box <- function(cx, cy, w, h, fill, title, title_size = 7.0) {
    add(rr(cx, cy, w, h, fill = fill, col = fill, lwd = 0.7, radius = 3.0))
    add(txt(title, cx, cy - h / 2 + 34, size = title_size, col = "white", face = "bold"))
}

node <- function(cx, cy, w, h, label, size = 5.9, fill = pal[["node"]],
                 col = pal[["node_border"]], text_col = "white", lwd = 0.65) {
    add(rr(cx, cy, w, h, fill = fill, col = col, lwd = lwd, radius = 2.0))
    add(txt(label, cx, cy, size = size, col = text_col, face = "plain"))
}

add(rectGrob(gp = gpar(fill = "white", col = NA)))

# Input ---------------------------------------------------------------------
node(120, 355, 230, 72, "Spatial transcriptomics\ndatasets", size = 4.2,
     fill = "white", col = "#7B838A", text_col = pal[["ink"]], lwd = 0.75)

# Feature-selection group ---------------------------------------------------
group_box(455, 355, 350, 540, pal[["teal"]], "Select features", title_size = 6.2)
fs_y <- c(210, 355, 500)
node(455, fs_y[1], 270, 76, "Expression-driven / HVG", size = 4.2)
node(455, fs_y[2], 270, 76, "Spatially informed / SVG", size = 4.2)
node(455, fs_y[3], 270, 76, "Controls", size = 5.0)

# Input fans out to all feature-selection routes.
add(seg(235, 355, 255, 355, lwd = 1.15))
add(seg(255, fs_y[1], 255, fs_y[3], lwd = 1.15))
for (yy in fs_y) add(seg(255, yy, 316, yy, lwd = 1.15, arrow_end = TRUE))

# Workflow group ------------------------------------------------------------
group_box(825, 355, 300, 540, pal[["blue"]], "Run workflows", title_size = 6.2)
node(765, 255, 105, 66, "scVI", size = 5.3)
node(885, 255, 120, 66, "CellCharter", size = 4.2)
node(825, 355, 150, 66, "GraphST", size = 5.0)
node(825, 475, 230, 84, "Feature budgets\nand repeated seeds", size = 4.9)
add(txt("100 to 10,000 genes", 825, 535, size = 5.0, col = "white", face = "bold"))

# Feature routes merge before all workflow settings.
merge_fs_x <- 650
for (yy in fs_y) add(seg(590, yy, merge_fs_x, yy, lwd = 1.0))
add(seg(merge_fs_x, fs_y[1], merge_fs_x, fs_y[3], lwd = 1.0))
add(pt(merge_fs_x, 355, fill = pal[["arrow"]], size = 1.6))
add(seg(merge_fs_x, 355, 674, 355, lwd = 1.2, arrow_end = TRUE))

# Internal workflow relationship.
for (xx in c(765, 885)) add(seg(xx, 289, xx, 315, col = "#EAF2F8", lwd = 0.7))
add(seg(765, 315, 885, 315, col = "#EAF2F8", lwd = 0.7))
add(seg(825, 315, 825, 322, col = "#EAF2F8", lwd = 0.7, arrow_end = TRUE))
add(seg(825, 388, 825, 432, col = "#EAF2F8", lwd = 0.7, arrow_end = TRUE))

# Evaluation group ----------------------------------------------------------
group_box(1315, 355, 610, 620, pal[["green"]], "Evaluate", title_size = 6.2)
add(txt("Spatial signal evidence", 1315, 125, size = 6.2, col = pal[["green_dark"]], face = "bold"))
add(txt("Downstream representation utility", 1315, 372, size = 6.2,
        col = pal[["green_dark"]], face = "bold"))

eval_x <- c(1135, 1315, 1495)
top_y <- 235
bot_y <- 485
node(eval_x[1], top_y, 150, 78, "Simulation\ntruth", size = 4.6)
node(eval_x[2], top_y, 150, 78, "Held-out\nslices", size = 4.6)
node(eval_x[3], top_y, 150, 78, "Cross-statistic", size = 4.3)
node(eval_x[1], bot_y, 150, 78, "Integration", size = 4.8)
node(eval_x[2], bot_y, 150, 78, "Clustering", size = 4.8)
node(eval_x[3], bot_y, 150, 78, "Alignment", size = 4.8)
add(txt("truth and reproducibility", 1315, 305, size = 4.7, col = "white", face = "bold"))
add(txt("task-specific metrics", 1315, 610, size = 4.7, col = "white", face = "bold"))

# Workflow outputs split into the two evidence tracks.
run_out_x <- 995
split_x <- 1028
add(seg(run_out_x, 355, split_x, 355, lwd = 1.15))
add(seg(split_x, top_y, split_x, bot_y, lwd = 1.05))
add(pt(split_x, 355, fill = pal[["arrow"]], size = 1.6))
top_bus_y <- 158
bot_bus_y <- 567
add(seg(split_x, top_y, split_x, top_bus_y, col = pal[["coral"]], lwd = 1.05))
add(seg(split_x, top_bus_y, eval_x[3], top_bus_y, col = pal[["coral"]], lwd = 1.05))
for (xx in eval_x) add(seg(xx, top_bus_y, xx, top_y - 40,
                           col = pal[["coral"]], lwd = 0.9, arrow_end = TRUE))
add(seg(split_x, bot_y, split_x, bot_bus_y, col = pal[["indigo"]], lwd = 1.05))
add(seg(split_x, bot_bus_y, eval_x[3], bot_bus_y, col = pal[["indigo"]], lwd = 1.05))
for (xx in eval_x) add(seg(xx, bot_bus_y, xx, bot_y + 40,
                           col = pal[["indigo"]], lwd = 0.9, arrow_end = TRUE))

# Benchmark-output group ----------------------------------------------------
group_box(1825, 355, 340, 540, pal[["gold"]], "Summarize benchmark", title_size = 5.4)
node(1748, 240, 145, 76, "Method\nranking", size = 4.7)
node(1902, 240, 145, 76, "Feature-number\neffects", size = 4.4)
node(1748, 405, 145, 76, "Metric\nrelationships", size = 4.5)
node(1902, 405, 145, 76, "Integration-method\ncomparison", size = 4.3)
add(txt("scores, ranks and\nsensitivity analyses", 1825, 520,
        size = 4.1, col = "white", face = "bold"))

# The evaluation module produces one summarized benchmark stream.
add(pt(1620, 355, fill = pal[["arrow"]], size = 1.6))
add(seg(1620, 355, 1654, 355, lwd = 1.2, arrow_end = TRUE))

# Final guidance ------------------------------------------------------------
node(2180, 355, 220, 96, "Objective-aware\nmethod guidance", size = 4.8,
     fill = pal[["pale_gold"]], col = pal[["gold_dark"]],
     text_col = pal[["gold_dark"]], lwd = 0.9)

merge_sum_x <- 2028
for (yy in c(240, 405)) {
    add(seg(1975, yy, merge_sum_x, yy, lwd = 0.9))
}
add(seg(merge_sum_x, 240, merge_sum_x, 405, lwd = 0.9))
add(pt(merge_sum_x, 355, fill = pal[["arrow"]], size = 1.6))
add(seg(merge_sum_x, 355, 2069, 355, lwd = 1.2, arrow_end = TRUE))

figure_grob <- do.call(grobTree, grobs)
figure <- as.ggplot(figure_grob) +
    theme_void() +
    theme(plot.margin = margin(1, 1, 1, 1, "mm"))

stem <- file.path(output_dir, "Figure1_reference_dag_v5")
width_mm <- 183
height_mm <- 70

ggsave(paste0(stem, ".svg"), figure, width = width_mm, height = height_mm, units = "mm",
       device = svglite::svglite, bg = "white")
svg <- readLines(paste0(stem, ".svg"), warn = FALSE)
svg <- gsub(" textLength='[^']+' lengthAdjust='spacingAndGlyphs'", "", svg)
writeLines(svg, paste0(stem, ".svg"), useBytes = TRUE)

ggsave(paste0(stem, ".pdf"), figure, width = width_mm, height = height_mm, units = "mm",
       device = grDevices::pdf, family = "Helvetica", useDingbats = FALSE, bg = "white")
ggsave(paste0(stem, ".tiff"), figure, width = width_mm, height = height_mm, units = "mm",
       dpi = 600, device = ragg::agg_tiff, compression = "lzw", bg = "white")
ggsave(paste0(stem, ".png"), figure, width = width_mm, height = height_mm, units = "mm",
       dpi = 300, device = ragg::agg_png, bg = "white")
ggsave(file.path(output_dir, "Figure1_reference_dag_v5_print_qa.png"), figure,
       width = width_mm, height = height_mm, units = "mm",
       dpi = 150, device = ragg::agg_png, bg = "white")
saveRDS(figure, paste0(stem, ".rds"))

message("Reference-DAG Figure 1 written to ", output_dir)
