#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(tidyr)
})

source("analysis/manuscript_rebuild_v2/R/theme.R")
source("analysis/manuscript_rebuild_v2/R/data.R")
ensure_figure_dirs()

arrow_spec <- grid::arrow(length = unit(1.8, "mm"), type = "closed")

box_layer <- function(xmin, xmax, ymin, ymax, fill, colour = BENCH_COLOURS$light,
                      linewidth = 0.35) {
  annotate(
    "rect", xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
    fill = fill, colour = colour, linewidth = linewidth
  )
}

# Panel a --------------------------------------------------------------------
tissue <- expand.grid(col = 0:5, row = 0:5) %>%
  mutate(
    x = 8.1 + col * 1.15 + if_else(row %% 2 == 0, 0, 0.56),
    y = 19.2 + row * 0.82,
    keep = ((col - 2.5)^2 / 9 + (row - 2.5)^2 / 8) < 1,
    zone = if_else(x + 1.2 * y > 36, "domain 2", "domain 1")
  ) %>%
  filter(keep)

matrix_cells <- tidyr::expand_grid(route = 1:3, col = 1:7, row = 1:4) %>%
  mutate(
    x = 48.0 + col * 1.0,
    y0 = recode(route, `1` = 29.2, `2` = 20.6, `3` = 12.0),
    y = y0 + row * 0.64,
    signal = case_when(
      route == 1 ~ ((col + 2 * row) %% 7) / 6,
      route == 2 ~ ((2 * col + row + 2) %% 7) / 6,
      TRUE ~ ((3 * col + row + 1) %% 7) / 6
    )
  )

route_df <- tibble(
  route = c("Expression-based", "Spatially informed", "Controls / references"),
  detail = c("24 selectors", "5 SVG selectors", "4 controls + 1 oracle"),
  y = c(32, 23.4, 14.8),
  fill = c("#E8EEF3", "#F8E9E4", "#EFEFEF"),
  colour = c(BENCH_COLOURS$expression, BENCH_COLOURS$spatial, BENCH_COLOURS$mid)
)

task_df <- tibble(
  task = c("Integration", "Clustering", "Alignment"),
  detail = c("mixing + conservation", "domains + label agreement", "paired slices only"),
  y = c(30.4, 22.2, 14.0),
  fill = c("#E5F0F7", "#F8E8E7", "#E8F2EC"),
  colour = c(BENCH_COLOURS$integration, BENCH_COLOURS$clustering, BENCH_COLOURS$alignment)
)

p_a <- ggplot() +
  box_layer(2.2, 17.1, 10.0, 37.2, BENCH_COLOURS$panel) +
  annotate("text", x = 3.3, y = 35.6, label = "Shared spatial input", hjust = 0,
           family = BENCH_FONT, fontface = "bold", size = 3.0) +
  annotate("text", x = 3.3, y = 33.5, label = "counts  +  coordinates  +  slice", hjust = 0,
           family = BENCH_FONT, colour = BENCH_COLOURS$mid, size = 2.45) +
  geom_point(
    data = tissue,
    aes(x, y, fill = I(if_else(zone == "domain 1", "#88A8BE", "#D89579"))),
    shape = 21, size = 2.25,
    stroke = 0.18, colour = "white", show.legend = FALSE
  ) +
  annotate("rect", xmin = 5.2, xmax = 13.8, ymin = 12.3, ymax = 17.3,
           fill = "white", colour = "#BDBDBD", linewidth = 0.3) +
  geom_tile(
    data = expand.grid(x = seq(5.7, 13.3, length.out = 10), y = seq(12.8, 16.8, length.out = 6)) %>%
      mutate(v = ((row_number() * 7) %% 13) / 12),
    aes(x, y, fill = v), width = 0.65, height = 0.52, colour = NA,
    inherit.aes = FALSE, show.legend = FALSE
  ) +
  annotate("text", x = 9.5, y = 11.25, label = "spot × gene matrix", family = BENCH_FONT,
           size = 2.25, colour = BENCH_COLOURS$mid) +
  annotate("segment", x = 17.6, xend = 19.5, y = 23.5, yend = 23.5,
           colour = BENCH_COLOURS$ink, linewidth = 0.42) +
  annotate("segment", x = 19.5, xend = 19.5, y = 14.8, yend = 32,
           colour = BENCH_COLOURS$ink, linewidth = 0.42) +
  annotate("segment", x = 19.5, xend = 21.4, y = 32, yend = 32,
           colour = BENCH_COLOURS$ink, linewidth = 0.42, arrow = arrow_spec) +
  annotate("segment", x = 19.5, xend = 21.4, y = 23.4, yend = 23.4,
           colour = BENCH_COLOURS$ink, linewidth = 0.42, arrow = arrow_spec) +
  annotate("segment", x = 19.5, xend = 21.4, y = 14.8, yend = 14.8,
           colour = BENCH_COLOURS$ink, linewidth = 0.42, arrow = arrow_spec) +
  annotate("text", x = 31.8, y = 39.0, label = "Feature selection", hjust = 0.5,
           family = BENCH_FONT, fontface = "bold", size = 2.75) +
  annotate("text", x = 31.8, y = 36.9, label = "feature budgets: 100 to 10,000", hjust = 0.5,
           family = BENCH_FONT, colour = BENCH_COLOURS$mid, size = 2.15) +
  geom_rect(
    data = route_df,
    aes(xmin = 22.0, xmax = 41.5, ymin = y - 3.0, ymax = y + 3.0,
        fill = I(fill), colour = I(colour)), linewidth = 0.45
  ) +
  geom_text(
    data = route_df, aes(x = 23.4, y = y + 0.75, label = route, colour = I(colour)),
    hjust = 0, family = BENCH_FONT, fontface = "bold", size = 2.7
  ) +
  geom_text(
    data = route_df, aes(x = 23.4, y = y - 1.25, label = detail),
    hjust = 0, family = BENCH_FONT, colour = BENCH_COLOURS$mid, size = 2.3
  ) +
  annotate("segment", x = 42.0, xend = 45.2, y = 32, yend = 32,
           colour = BENCH_COLOURS$expression, linewidth = 0.5, arrow = arrow_spec) +
  annotate("segment", x = 42.0, xend = 45.2, y = 23.4, yend = 23.4,
           colour = BENCH_COLOURS$spatial, linewidth = 0.5, arrow = arrow_spec) +
  annotate("segment", x = 42.0, xend = 45.2, y = 14.8, yend = 14.8,
           colour = BENCH_COLOURS$mid, linewidth = 0.5, arrow = arrow_spec) +
  annotate("rect", xmin = 45.7, xmax = 58.0, ymin = 9.8, ymax = 37.2,
           fill = "#FAFAF8", colour = "#CFCFCF", linewidth = 0.35) +
  annotate("text", x = 51.85, y = 39.0,
           label = "Selected\nrepresentation", family = BENCH_FONT,
           fontface = "bold", lineheight = 0.88, size = 2.65) +
  geom_tile(
    data = matrix_cells,
    aes(x, y, fill = signal), width = 0.82, height = 0.49,
    colour = "white", linewidth = 0.08, show.legend = FALSE
  ) +
  scale_fill_gradient(low = "#F4F2EC", high = "#3F647A", guide = "none") +
  geom_path(
    data = tibble(
      x = rep(seq(47.2, 56.8, length.out = 12), 3),
      route = rep(1:3, each = 12),
      y = c(
        27.8 + c(0.1, 0.2, 0.8, 1.7, 1.1, 0.4, 0.6, 1.3, 1.8, 1.1, 0.3, 0.2),
        19.2 + c(1.5, 1.3, 0.7, 0.2, 0.4, 1.2, 1.8, 1.4, 0.8, 0.5, 0.9, 1.4),
        10.6 + c(0.5, 0.8, 0.6, 0.9, 0.7, 1.0, 0.8, 0.6, 0.9, 0.7, 0.8, 0.7)
      )
    ),
    aes(x, y, group = route, colour = factor(route)), linewidth = 0.7,
    show.legend = FALSE
  ) +
  scale_colour_manual(
    values = c(`1` = BENCH_COLOURS$expression, `2` = BENCH_COLOURS$spatial,
               `3` = BENCH_COLOURS$mid), guide = "none"
  ) +
  annotate("text", x = 51.85, y = 7.4,
           label = "same input, different selected representations",
           family = BENCH_FONT, fontface = "italic", size = 2.15,
           colour = BENCH_COLOURS$ink) +
  annotate("segment", x = 58.5, xend = 60.0, y = 23.5, yend = 23.5,
           colour = BENCH_COLOURS$ink, linewidth = 0.42) +
  annotate("segment", x = 60.0, xend = 60.0, y = 18.9, yend = 28.1,
           colour = BENCH_COLOURS$ink, linewidth = 0.42) +
  annotate("segment", x = 60.0, xend = 62.0, y = 28.1, yend = 28.1,
           colour = BENCH_COLOURS$ink, linewidth = 0.42, arrow = arrow_spec) +
  annotate("segment", x = 60.0, xend = 62.0, y = 18.9, yend = 18.9,
           colour = BENCH_COLOURS$ink, linewidth = 0.42, arrow = arrow_spec) +
  annotate("text", x = 69.1, y = 39.0, label = "Integration", hjust = 0.5,
           family = BENCH_FONT, fontface = "bold", size = 2.75) +
  annotate("rect", xmin = 62.6, xmax = 75.5, ymin = 25.0, ymax = 31.2,
           fill = "#E7EDF3", colour = BENCH_COLOURS$scvi, linewidth = 0.45) +
  annotate("text", x = 69.05, y = 28.1, label = "scVI", family = BENCH_FONT,
           fontface = "bold", size = 2.9, colour = BENCH_COLOURS$scvi) +
  annotate("rect", xmin = 62.6, xmax = 75.5, ymin = 15.8, ymax = 22.0,
           fill = "#F8EEDF", colour = BENCH_COLOURS$cellcharter, linewidth = 0.45) +
  annotate("text", x = 69.05, y = 18.9, label = "CellCharter", family = BENCH_FONT,
           fontface = "bold", size = 2.8, colour = BENCH_COLOURS$cellcharter) +
  annotate("text", x = 69.05, y = 12.8, label = "fixed parameters\nmatched seeds",
           family = BENCH_FONT, lineheight = 0.9, size = 2.0, colour = BENCH_COLOURS$mid) +
  annotate("segment", x = 76.0, xend = 79.0, y = 23.5, yend = 23.5,
           colour = BENCH_COLOURS$ink, linewidth = 0.45, arrow = arrow_spec) +
  annotate("text", x = 89.2, y = 39.0, label = "Task evidence", hjust = 0.5,
           family = BENCH_FONT, fontface = "bold", size = 2.75) +
  geom_rect(
    data = task_df,
    aes(xmin = 79.6, xmax = 98.8, ymin = y - 2.9, ymax = y + 2.9,
        fill = I(fill), colour = I(colour)), linewidth = 0.45
  ) +
  geom_text(
    data = task_df, aes(x = 81.0, y = y + 0.7, label = task, colour = I(colour)),
    hjust = 0, family = BENCH_FONT, fontface = "bold", size = 2.7
  ) +
  geom_text(
    data = task_df, aes(x = 81.0, y = y - 1.15, label = detail),
    hjust = 0, family = BENCH_FONT, colour = BENCH_COLOURS$mid, size = 2.2
  ) +
  annotate("text", x = 89.2, y = 7.4,
           label = "frozen scaling; missing values remain NA",
           family = BENCH_FONT, size = 2.05, colour = BENCH_COLOURS$mid) +
  coord_cartesian(xlim = c(0, 101), ylim = c(5.5, 41.3), clip = "off") +
  theme_void(base_family = BENCH_FONT) +
  theme(plot.margin = margin(1, 2, 1, 2, "mm"))

# Panel b --------------------------------------------------------------------
metadata <- read_tsv("results/current_rank/data/datasets-metadata.tsv", show_col_types = FALSE) %>%
  filter(Dataset != "STOmicsVisium5Samples") %>%
  transmute(
    dataset = Dataset,
    Dataset = display_dataset(Dataset),
    Platform,
    Slices = NSlices,
    Spots
  )

coverage <- metadata %>%
  select(dataset, Dataset) %>%
  crossing(Task = c("Integration", "Clustering", "Alignment")) %>%
  mutate(
    Available = Task != "Alignment" |
      dataset %in% c("DLPFC", "MouseBrainSerialSections"),
    Task = factor(Task, levels = c("Integration", "Clustering", "Alignment")),
    Dataset = factor(Dataset, levels = rev(metadata$Dataset))
  )

p_b <- ggplot(coverage, aes(Task, Dataset)) +
  geom_tile(fill = "#F4F4F1", colour = "white", linewidth = 0.6) +
  geom_point(
    aes(fill = Available), shape = 21, size = 3.1, stroke = 0.35,
    colour = ifelse(coverage$Available, BENCH_COLOURS$ink, "#C8C8C8")
  ) +
  scale_fill_manual(values = c(`TRUE` = BENCH_COLOURS$ink, `FALSE` = "white"), guide = "none") +
  scale_x_discrete(position = "top") +
  labs(
    x = NULL, y = NULL,
    title = "Seven datasets with integer counts",
    subtitle = "Alignment is evaluated only where explicit slice pairing is defined"
  ) +
  theme_benchmark(6.4) +
  theme(
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.text.x = element_text(face = "bold", margin = margin(b = 2)),
    axis.text.y = element_text(face = "bold"),
    plot.title = element_text(size = 7.6),
    plot.subtitle = element_text(size = 6.0),
    plot.margin = margin(3, 3, 2, 2, "mm")
  )

# Panel c --------------------------------------------------------------------
p_c <- ggplot() +
  annotate("text", x = 0.5, y = 9.8, label = "Score construction", hjust = 0,
           family = BENCH_FONT, fontface = "bold", size = 3.0) +
  annotate("text", x = 0.5, y = 8.7,
           label = "dataset × metric frozen ranges; higher is always better",
           hjust = 0, family = BENCH_FONT, size = 2.25, colour = BENCH_COLOURS$mid) +
  annotate("rect", xmin = 0.5, xmax = 4.1, ymin = 5.3, ymax = 7.5,
           fill = "#E5F0F7", colour = BENCH_COLOURS$integration, linewidth = 0.4) +
  annotate("text", x = 2.3, y = 6.7, label = "Integration", family = BENCH_FONT,
           fontface = "bold", size = 2.7, colour = BENCH_COLOURS$integration) +
  annotate("text", x = 2.3, y = 5.85, label = "mixing + conservation", family = BENCH_FONT,
           size = 2.05, colour = BENCH_COLOURS$mid) +
  annotate("rect", xmin = 0.5, xmax = 4.1, ymin = 2.4, ymax = 4.6,
           fill = "#F8E8E7", colour = BENCH_COLOURS$clustering, linewidth = 0.4) +
  annotate("text", x = 2.3, y = 3.8, label = "Clustering", family = BENCH_FONT,
           fontface = "bold", size = 2.7, colour = BENCH_COLOURS$clustering) +
  annotate("text", x = 2.3, y = 2.95, label = "structure + agreement", family = BENCH_FONT,
           size = 2.05, colour = BENCH_COLOURS$mid) +
  annotate("segment", x = 4.35, xend = 5.55, y = 5.0, yend = 5.0,
           colour = BENCH_COLOURS$ink, linewidth = 0.4, arrow = arrow_spec) +
  annotate("rect", xmin = 5.8, xmax = 10.0, ymin = 3.55, ymax = 6.45,
           fill = "#EEEEEC", colour = BENCH_COLOURS$core, linewidth = 0.5) +
  annotate("text", x = 7.9, y = 5.5, label = "CoreOverall", family = BENCH_FONT,
           fontface = "bold", size = 2.9) +
  annotate("text", x = 7.9, y = 4.55,
           label = "(Integration + Clustering) / 2", family = BENCH_FONT, size = 2.0) +
  annotate("rect", xmin = 5.8, xmax = 10.0, ymin = 0.45, ymax = 2.45,
           fill = "#E8F2EC", colour = BENCH_COLOURS$alignment, linewidth = 0.4) +
  annotate("text", x = 7.9, y = 1.65, label = "Alignment reported separately",
           family = BENCH_FONT, fontface = "bold", size = 2.35,
           colour = BENCH_COLOURS$alignment) +
  coord_cartesian(xlim = c(0.2, 10.3), ylim = c(0.1, 10.2), clip = "off") +
  theme_void(base_family = BENCH_FONT) +
  theme(plot.margin = margin(3, 2, 2, 3, "mm"))

figure <- p_a + (p_b + p_c + plot_layout(widths = c(0.96, 1.04))) +
  plot_layout(heights = c(1.55, 1.0)) +
  plot_annotation(tag_levels = "a") &
  panel_label_theme

write_source_data(metadata, "figure1_dataset_metadata.tsv")
write_source_data(coverage %>% mutate(Available = as.logical(Available)), "figure1_task_coverage.tsv")
write_source_data(
  tibble(
    FeatureSelectionGroup = c("Expression-based", "Spatially informed", "Control/reference", "Label-informed oracle"),
    NMethods = c(24L, 5L, 4L, 1L),
    IncludedInCompetitiveRank = c(TRUE, TRUE, FALSE, FALSE)
  ),
  "figure1_method_scope.tsv"
)

save_figure_set(
  figure,
  file.path(FIGURE_ROOT, "figures", "figure1_benchmark_framework"),
  width_mm = 183,
  height_mm = 128
)

message("Figure 1 written to ", file.path(FIGURE_ROOT, "figures"))
