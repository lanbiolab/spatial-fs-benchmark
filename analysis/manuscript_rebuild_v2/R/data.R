suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
})

SCORE_ROOT <- "results/spatial_svg_rebuild_v1/frozen_scores"
FIGURE_ROOT <- "results/manuscript_rebuild_v2"

ensure_figure_dirs <- function() {
  invisible(lapply(
    file.path(FIGURE_ROOT, c("figures", "source_data", "previews")),
    dir.create,
    recursive = TRUE,
    showWarnings = FALSE
  ))
}

read_score_table <- function(name) {
  readr::read_tsv(file.path(SCORE_ROOT, name), show_col_types = FALSE)
}

method_labels <- c(
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

dataset_labels <- c(
  DLPFC = "DLPFC",
  MouseBrainSerialSections = "Mouse brain sections",
  STOmics0212 = "STOmics 0212",
  STOmics0218 = "STOmics 0218",
  STOmics0224 = "STOmics 0224",
  E8p5Embryo = "E8.5 embryo",
  E9p5Embryo = "E9.5 embryo",
  STOmicsVisium5Samples = "STOmics Visium (5 samples)"
)

display_method <- function(x) {
  out <- unname(method_labels[x])
  out[is.na(out)] <- stringr::str_to_title(stringr::str_replace_all(x[is.na(out)], "_", " "))
  out
}

display_dataset <- function(x) {
  out <- unname(dataset_labels[x])
  out[is.na(out)] <- x[is.na(out)]
  out
}

display_integrator <- function(x) {
  dplyr::recode(x, scvi = "scVI", cellcharter = "CellCharter", .default = x)
}

group_label <- function(x) {
  dplyr::recode(
    x,
    expression_driven = "Expression-based",
    spatially_informed = "Spatially informed",
    control = "Control/reference",
    label_informed = "Label-informed oracle",
    .default = x
  )
}

write_source_data <- function(x, filename) {
  ensure_figure_dirs()
  readr::write_tsv(x, file.path(FIGURE_ROOT, "source_data", filename), na = "NA")
}

