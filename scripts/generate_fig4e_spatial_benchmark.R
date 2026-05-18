#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[[1]] else file.path("results", "current_rank", "data")
output_dir <- if (length(args) >= 2) args[[2]] else file.path("results", "fig4e_spatial_benchmark", "figures")
base_dir <- normalizePath(file.path(data_dir, ".."), mustWork = FALSE)
ranges_path <- file.path(base_dir, "output", "baseline-ranges.tsv")

source(file.path("external", "atlas-feature-selection-benchmark", "analysis", "R", "plotting.R"))
source(file.path("analysis", "atlas_style", "R", "spatial_summarisation.R"))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

metrics <- read_tsv(file.path(data_dir, "benchmark.tsv"), show_col_types = FALSE) |>
    filter(.data$Type %in% c("Integration", "Clustering", "Alignment"))

metric_ranges <- read_tsv(ranges_path, show_col_types = FALSE) |>
    filter(.data$Type %in% c("Integration", "Clustering", "Alignment"))

datasets_meta <- read_tsv(file.path(data_dir, "datasets-metadata.tsv"), show_col_types = FALSE)
dataset_names <- datasets_meta$Name
names(dataset_names) <- datasets_meta$Dataset

dataset_keep <- c(
    "DLPFC",
    "MouseBrainSerialSections",
    "STOmics0212",
    "STOmics0218",
    "STOmics0224",
    "STOmicsVisium5Samples"
)

metrics_summary <- summarise_spatial_metrics(
    metrics,
    metric_ranges,
    type_weights = c("Integration" = 1/3, "Clustering" = 1/3, "Alignment" = 1/3),
    require_types_for_overall = c("Integration", "Clustering", "Alignment")
) |>
    filter(.data$IntegrationLabel == "scVI", .data$Dataset %in% dataset_keep)

scanpy_pairs <- tribble(
    ~MethodBase, ~Standard, ~BatchAware, ~Label,
    "scanpy_seurat", "scanpy_seurat-N2000", "scanpy_seurat_batch-N2000", "scanpy-Seurat",
    "scanpy_pearson", "scanpy_pearson-N2000", "scanpy_pearson_batch-N2000", "scanpy-Pearson",
    "scanpy_cell_ranger", "scanpy_cell_ranger-N2000", "scanpy_cell_ranger_batch-N2000", "scanpy-Cell Ranger",
    "scanpy_seurat_v3", "scanpy_seurat_v3-N2000", "scanpy_seurat_v3_batch-N2000", "scanpy-SeuratV3"
)

plotting <- metrics_summary |>
    select(.data$Dataset, .data$Method, .data$Overall, .data$Integration, .data$Clustering, .data$Alignment) |>
    pivot_longer(
        cols = c("Overall", "Integration", "Clustering", "Alignment"),
        names_to = "Type",
        values_to = "Value"
    ) |>
    inner_join(
        scanpy_pairs |>
            pivot_longer(
                cols = c("Standard", "BatchAware"),
                names_to = "Version",
                values_to = "Method"
            ),
        by = "Method"
    ) |>
    select(.data$Dataset, .data$Type, .data$MethodBase, .data$Label, .data$Version, .data$Value) |>
    pivot_wider(
        names_from = "Version",
        values_from = "Value"
    ) |>
    mutate(Difference = .data$BatchAware - .data$Standard)

method_order <- c("scanpy-Seurat", "scanpy-Pearson", "scanpy-Cell Ranger", "scanpy-SeuratV3")
dataset_order <- dataset_keep
type_order <- c("Overall", "Integration", "Clustering", "Alignment")
y_positions <- c(
    "scanpy-Seurat" = 4.00,
    "scanpy-Pearson" = 3.35,
    "scanpy-Cell Ranger" = 2.70,
    "scanpy-SeuratV3" = 2.05
)

plotting <- plotting |>
    mutate(
        Label = factor(.data$Label, levels = rev(method_order)),
        Y = unname(y_positions[as.character(.data$Label)]),
        DatasetLabel = factor(.data$Dataset, levels = dataset_order, labels = dataset_names[dataset_order]),
        Type = factor(.data$Type, levels = type_order)
    )

p <- ggplot(plotting, aes(x = .data$DatasetLabel, y = .data$Y, fill = .data$Difference)) +
    geom_tile(width = 0.97, height = 0.54) +
    colorspace::scale_fill_continuous_diverging(
        palette = "Purple-Green",
        name = "Batch aware -\nstandard",
        na.value = "#D9D9D9",
        limits = c(-0.35, 0.65),
        oob = scales::squish
    ) +
    scale_y_continuous(
        breaks = unname(y_positions[rev(method_order)]),
        labels = rev(method_order),
        expand = expansion(mult = c(0.03, 0.03))
    ) +
    facet_grid(. ~ .data$Type) +
    labs(x = NULL, y = NULL) +
    theme_features_pub() +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "bold"),
        axis.text.y = element_text(face = "bold"),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "black", colour = "black"),
        strip.text = element_text(colour = "white", face = "bold"),
        legend.position = "right",
        legend.title.position = "top",
        panel.spacing.x = unit(0.06, "cm"),
        plot.margin = margin(0.08, 0.10, 0.08, 0.08, "cm")
    )

save_figure_files(p, file.path(output_dir, "figure_fig4e_panel_e"), width = 7.6, height = 2.35)
write_tsv(plotting, file.path(output_dir, "fig4e_batch_aware_differences.tsv"))
