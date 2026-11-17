suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(patchwork)
    library(readr)
    library(scales)
    library(tidyr)
})

source(file.path("external", "atlas-feature-selection-benchmark", "analysis", "R", "plotting.R"))

supplement_root <- file.path("results", "manuscript_rebuild_v2", "supplemental")
supplement_figure_dir <- file.path(supplement_root, "figures")
supplement_source_dir <- file.path(supplement_root, "source_data")
dir.create(supplement_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supplement_source_dir, recursive = TRUE, showWarnings = FALSE)

supplement_method_labels <- c(
    all_features = "All features",
    random = "Random",
    TFs = "Transcription factors",
    scsegindex = "scSEGIndex",
    Brennecke = "Brennecke",
    anticor = "Anticor",
    dubstepr = "DUBStepR",
    hotspot = "Hotspot",
    nbumi = "NBumi",
    osca = "OSCA",
    scPNMF = "scPNMF",
    scanpy_cell_ranger = "Scanpy Cell Ranger",
    scanpy_cell_ranger_batch = "Scanpy Cell Ranger, batch",
    scanpy_pearson = "Scanpy Pearson",
    scanpy_pearson_batch = "Scanpy Pearson, batch",
    scanpy_seurat = "Scanpy Seurat",
    scanpy_seurat_batch = "Scanpy Seurat, batch",
    scanpy_seurat_v3 = "Scanpy Seurat v3",
    scanpy_seurat_v3_batch = "Scanpy Seurat v3, batch",
    scry = "scry",
    seurat_disp = "Seurat Dispersion",
    seurat_mvp = "Seurat MVP",
    seurat_sct = "Seurat sctransform",
    seurat_vst = "Seurat VST",
    singleCellHaystack = "singleCellHaystack",
    statistic_mean = "Mean expression",
    statistic_variance = "Expression variance",
    triku = "triku",
    wilcoxon = "Wilcoxon (label-informed)",
    morans_i = "Moran's I",
    nnsvg = "nnSVG",
    somde = "SOMDE",
    sparkx = "SPARK-X",
    spatialde = "SpatialDE"
)

supplement_dataset_labels <- c(
    DLPFC = "DLPFC",
    MouseBrainSerialSections = "Mouse Brain",
    STOmics0212 = "STOmics-0212",
    STOmics0218 = "STOmics-0218",
    STOmics0224 = "STOmics-0224",
    E8p5Embryo = "E8.5 Embryo",
    E9p5Embryo = "E9.5 Embryo"
)

supplement_dataset_palette <- c(
    "DLPFC" = "#8DD3C7",
    "Mouse Brain" = "#FB8072",
    "STOmics-0212" = "#80B1D3",
    "STOmics-0218" = "#FDB462",
    "STOmics-0224" = "#B3DE69",
    "E8.5 Embryo" = "#FFFFB3",
    "E9.5 Embryo" = "#BEBADA"
)

supplement_task_palette <- c(
    "Core overall" = "#F781BF",
    "Integration" = "#E41A1C",
    "Clustering" = "#377EB8",
    "Alignment" = "#4DAF4A"
)

supplement_integrator_labels <- c(
    scvi = "scVI",
    cellcharter = "CellCharter"
)

display_supplement_method <- function(x) {
    labels <- unname(supplement_method_labels[x])
    labels[is.na(labels)] <- x[is.na(labels)]
    labels
}

display_supplement_dataset <- function(x) {
    labels <- unname(supplement_dataset_labels[x])
    labels[is.na(labels)] <- x[is.na(labels)]
    labels
}

reference_supplement_theme <- function() {
    theme_features_pub() +
        theme(
            panel.grid.minor = element_blank(),
            strip.background = element_rect(fill = "black", colour = "black", linewidth = 0.35),
            strip.text = element_text(colour = "white", face = "bold", size = 6.2),
            panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.35),
            panel.spacing = unit(0.06, "cm"),
            plot.margin = margin(0.05, 0.05, 0.05, 0.05, "cm"),
            plot.tag = element_text(size = 8, face = "bold", colour = "black"),
            plot.tag.position = "topleft"
        )
}

save_reference_supplement <- function(figure, stem, width_mm = 208.3, height_mm = 180) {
    ggsave(
        paste0(stem, ".svg"), figure,
        width = width_mm, height = height_mm, units = "mm",
        device = svglite::svglite, bg = "white"
    )
    svg_path <- paste0(stem, ".svg")
    svg_text <- readLines(svg_path, warn = FALSE)
    svg_text <- gsub(" textLength='[^']+' lengthAdjust='spacingAndGlyphs'", "", svg_text)
    svg_text <- gsub(
        'font-family: "Liberation Sans";',
        "font-family: Arial, Helvetica, sans-serif;",
        svg_text,
        fixed = TRUE
    )
    writeLines(svg_text, svg_path, useBytes = TRUE)

    if (!"Arial" %in% names(grDevices::pdfFonts())) {
        grDevices::pdfFonts(Arial = grDevices::pdfFonts("Helvetica")[[1]])
    }
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
}
