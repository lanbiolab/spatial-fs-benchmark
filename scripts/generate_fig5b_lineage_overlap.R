#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(jsonlite)
})

source(file.path("external", "atlas-feature-selection-benchmark", "analysis", "R", "plotting.R"))

args <- commandArgs(trailingOnly = TRUE)
methods_path <- if (length(args) >= 1) args[[1]] else file.path("results", "fig5a_spatial_lineages", "figures", "fig5a_lineage_methods.tsv")
output_dir <- if (length(args) >= 2) args[[2]] else file.path("results", "fig5b_spatial_lineages", "figures")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

methods_tbl <- read_tsv(methods_path, show_col_types = FALSE)

method_levels <- rev(methods_tbl$Name)
combo_levels <- c("Full - Immune", "Full - Epithelial", "Immune - Epithelial")

read_feature_set <- function(dataset_key, method_base, n_features) {
    path <- switch(
        dataset_key,
        "Full" = file.path(
            "results", "spatial_main_native_seed0_fix3", "stomics_0212_wilcoxon_spatial_main",
            "stomics0212", "feature_selection", method_base, paste0("n", n_features), "seed0", "selected_features.json"
        ),
        "Immune" = file.path(
            "results", "lineage_subsets", "stomics_0212_immune_subset_scvi",
            "stomics0212immune", "feature_selection", method_base, paste0("n", n_features), "seed0", "selected_features.json"
        ),
        "Epithelial" = file.path(
            "results", "lineage_subsets", "stomics_0212_epithelial_subset_scvi",
            "stomics0212epithelial", "feature_selection", method_base, paste0("n", n_features), "seed0", "selected_features.json"
        )
    )

    obj <- fromJSON(path)
    unique(obj$feature_names)
}

jaccard <- function(x, y) {
    both <- length(intersect(x, y))
    either <- length(union(x, y))
    if (either == 0) {
        return(NA_real_)
    }
    both / either
}

overlap_df <- bind_rows(lapply(seq_len(nrow(methods_tbl)), function(i) {
    method_base <- methods_tbl$MethodBase[[i]]
    method_name <- methods_tbl$Name[[i]]
    method_id <- methods_tbl$Method[[i]]
    n_features <- sub(".*-N", "", method_id)

    full_set <- read_feature_set("Full", method_base, n_features)
    immune_set <- read_feature_set("Immune", method_base, n_features)
    epithelial_set <- read_feature_set("Epithelial", method_base, n_features)

    tibble(
        Method = method_name,
        Combination = combo_levels,
        Jaccard = c(
            jaccard(full_set, immune_set),
            jaccard(full_set, epithelial_set),
            jaccard(immune_set, epithelial_set)
        )
    )
})) |>
    mutate(
        Method = factor(.data$Method, levels = method_levels),
        Combination = factor(.data$Combination, levels = combo_levels)
    )

write_tsv(overlap_df, file.path(output_dir, "fig5b_lineage_overlap.tsv"))

fig <- ggplot(overlap_df, aes(x = .data$Combination, y = .data$Method, fill = .data$Jaccard)) +
    geom_tile(colour = "white", linewidth = 0.35) +
    scale_fill_viridis_c(option = "C", limits = c(0, 1), name = "Jaccard\nindex") +
    labs(x = NULL, y = NULL) +
    theme_features_pub() +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "bold", colour = "black"),
        axis.text.y = element_text(face = "bold", colour = "black", size = 6),
        panel.grid = element_blank(),
        plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"),
        legend.position = "bottom"
    )

save_figure_files(fig, file.path(output_dir, "figure_fig5b_spatial_lineages"), width = 3.2, height = 7.8)
