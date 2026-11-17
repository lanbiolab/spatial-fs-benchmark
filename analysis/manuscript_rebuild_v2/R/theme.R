suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

BENCH_FONT <- "Liberation Sans"

BENCH_COLOURS <- list(
  expression = "#466B8A",
  spatial = "#D46A4C",
  control = "#A6A6A6",
  oracle = "#7A5AA6",
  integration = "#2C7FB8",
  clustering = "#D95F59",
  alignment = "#4F9A6D",
  core = "#252525",
  scvi = "#294C73",
  cellcharter = "#D78745",
  ink = "#252525",
  mid = "#717171",
  light = "#D9D9D9",
  paper = "#FFFFFF",
  panel = "#F5F5F3"
)

theme_benchmark <- function(base_size = 7, base_family = BENCH_FONT) {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      text = element_text(colour = BENCH_COLOURS$ink),
      axis.title = element_text(size = base_size + 0.5, face = "plain"),
      axis.text = element_text(size = base_size, colour = BENCH_COLOURS$ink),
      axis.line = element_line(linewidth = 0.32, colour = BENCH_COLOURS$ink),
      axis.ticks = element_line(linewidth = 0.28, colour = BENCH_COLOURS$ink),
      axis.ticks.length = unit(1.4, "mm"),
      panel.grid = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size + 0.5, face = "bold"),
      legend.title = element_text(size = base_size, face = "bold"),
      legend.text = element_text(size = base_size),
      legend.key.height = unit(3.4, "mm"),
      legend.key.width = unit(4.2, "mm"),
      plot.title = element_text(size = base_size + 1.5, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size, colour = BENCH_COLOURS$mid),
      plot.caption = element_text(size = base_size - 0.5, colour = BENCH_COLOURS$mid),
      plot.margin = margin(3, 3, 3, 3, unit = "mm")
    )
}

panel_label_theme <- theme(
  plot.tag = element_text(
    family = BENCH_FONT,
    size = 9,
    face = "bold",
    colour = BENCH_COLOURS$ink
  ),
  plot.tag.position = c(0, 1)
)

method_group_colours <- c(
  expression_driven = BENCH_COLOURS$expression,
  spatially_informed = BENCH_COLOURS$spatial,
  control = BENCH_COLOURS$control,
  label_informed = BENCH_COLOURS$oracle
)

task_colours <- c(
  CoreOverall = BENCH_COLOURS$core,
  Integration = BENCH_COLOURS$integration,
  Clustering = BENCH_COLOURS$clustering,
  Alignment = BENCH_COLOURS$alignment
)

integration_colours <- c(
  scvi = BENCH_COLOURS$scvi,
  cellcharter = BENCH_COLOURS$cellcharter
)

save_figure_set <- function(plot, stem, width_mm, height_mm, dpi = 600) {
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)

  if (requireNamespace("svglite", quietly = TRUE)) {
    ggplot2::ggsave(
      paste0(stem, ".svg"), plot = plot,
      width = width_mm, height = height_mm, units = "mm",
      device = svglite::svglite, bg = "white"
    )
  } else {
    ggplot2::ggsave(
      paste0(stem, ".svg"), plot = plot,
      width = width_mm, height = height_mm, units = "mm",
      device = grDevices::svg, bg = "white"
    )
  }

  svg_path <- paste0(stem, ".svg")
  svg_text <- readLines(svg_path, warn = FALSE)
  svg_text <- gsub(
    " textLength='[^']+' lengthAdjust='spacingAndGlyphs'",
    "",
    svg_text
  )
  writeLines(svg_text, svg_path, useBytes = TRUE)

  ggplot2::ggsave(
    paste0(stem, ".pdf"), plot = plot,
    width = width_mm, height = height_mm, units = "mm",
    device = grDevices::cairo_pdf, bg = "white"
  )
  ggplot2::ggsave(
    paste0(stem, ".tiff"), plot = plot,
    width = width_mm, height = height_mm, units = "mm", dpi = dpi,
    device = ragg::agg_tiff, compression = "lzw", bg = "white"
  )
  ggplot2::ggsave(
    paste0(stem, ".png"), plot = plot,
    width = width_mm, height = height_mm, units = "mm", dpi = 300,
    device = ragg::agg_png, bg = "white"
  )
}
