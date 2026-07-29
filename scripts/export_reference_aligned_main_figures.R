#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) args[[1]] else file.path("results", "reference_aligned_v2")
output_dir <- if (length(args) >= 2) args[[2]] else file.path(root, "main_figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

figure_specs <- list(
    Figure1 = list(
        path = file.path(root, "figure1", "figure1_complete_reference_aligned.rds"),
        width_mm = 284.5, height_mm = 294.6
    ),
    Figure2 = list(
        path = file.path(root, "baselines", "figure_baselines.rds"),
        width_mm = 279.4, height_mm = 330.2
    ),
    Figure3 = list(
        path = file.path(root, "num_features", "figure_num_features_benchmark.rds"),
        width_mm = 198.1, height_mm = 172.7
    ),
    Figure4 = list(
        path = file.path(root, "figure4", "figure4_complete_reference_aligned.rds"),
        width_mm = 208.3, height_mm = 259.1
    ),
    Figure5 = list(
        path = file.path(root, "fig6", "figure_fig6_spatial_integration.rds"),
        width_mm = 208.3, height_mm = 144.8
    )
)

export_one <- function(name, spec) {
    plot <- readRDS(spec$path)
    stem <- file.path(output_dir, paste0(name, "_reference_aligned"))

    ggsave(
        paste0(stem, ".svg"), plot,
        width = spec$width_mm, height = spec$height_mm, units = "mm",
        device = svglite::svglite, bg = "white"
    )
    svg_path <- paste0(stem, ".svg")
    svg_text <- readLines(svg_path, warn = FALSE)
    svg_text <- gsub(
        " textLength='[^']+' lengthAdjust='spacingAndGlyphs'",
        "",
        svg_text
    )
    svg_text <- gsub(
        'font-family: "Liberation Sans";',
        'font-family: Arial, Helvetica, sans-serif;',
        svg_text,
        fixed = TRUE
    )
    writeLines(svg_text, svg_path, useBytes = TRUE)
    if (!"Arial" %in% names(grDevices::pdfFonts())) {
        grDevices::pdfFonts(Arial = grDevices::pdfFonts("Helvetica")[[1]])
    }
    ggsave(
        paste0(stem, ".pdf"), plot,
        width = spec$width_mm, height = spec$height_mm, units = "mm",
        device = grDevices::pdf, family = "Arial",
        useDingbats = FALSE, bg = "white"
    )
    ggsave(
        paste0(stem, ".tiff"), plot,
        width = spec$width_mm, height = spec$height_mm, units = "mm",
        dpi = 600, device = ragg::agg_tiff, compression = "lzw", bg = "white"
    )
    ggsave(
        paste0(stem, ".png"), plot,
        width = spec$width_mm, height = spec$height_mm, units = "mm",
        dpi = 300, device = ragg::agg_png, bg = "white"
    )

    data.frame(
        Figure = name,
        SourceRDS = spec$path,
        WidthMM = spec$width_mm,
        HeightMM = spec$height_mm,
        stringsAsFactors = FALSE
    )
}

manifest <- do.call(rbind, Map(export_one, names(figure_specs), figure_specs))
readr::write_tsv(manifest, file.path(output_dir, "figure_manifest.tsv"))
