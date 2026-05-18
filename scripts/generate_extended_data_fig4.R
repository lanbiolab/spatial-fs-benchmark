#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(patchwork)
    library(cowplot)
    library(jsonlite)
    library(purrr)
    library(stringr)
    library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[[1]] else file.path("results", "current_rank", "data")
results_root <- if (length(args) >= 2) args[[2]] else file.path("results", "spatial_main_native_seed0_fix3")
output_dir <- if (length(args) >= 3) args[[3]] else file.path("results", "extended_data_fig4", "figures")

source(file.path("external", "atlas-feature-selection-benchmark", "analysis", "R", "plotting.R"))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

dataset_dir_map <- c(
    "DLPFC" = file.path(results_root, "dlpfc_spatial_main", "dlpfc"),
    "E8p5Embryo" = file.path(results_root, "e8p5_embryo_region_spatial_main", "e8p5embryo"),
    "E9p5Embryo" = file.path(results_root, "e9p5_embryo_region_spatial_main", "e9p5embryo"),
    "MouseBrainSerialSections" = file.path(results_root, "mouse_brain_serial_sections_spatial_main", "mousebrainserialsections"),
    "STOmics0212" = file.path(results_root, "stomics_0212_wilcoxon_spatial_main", "stomics0212"),
    "STOmics0218" = file.path(results_root, "stomics_0218_wilcoxon_spatial_main", "stomics0218"),
    "STOmics0224" = file.path(results_root, "stomics_0224_wilcoxon_spatial_main", "stomics0224"),
    "STOmicsVisium5Samples" = file.path(results_root, "stomics_visium_5samples_spatial_main", "stomicsvisium5samples")
)

benchmark <- read_tsv(file.path(data_dir, "benchmark.tsv"), show_col_types = FALSE)
methods_meta <- read_tsv(file.path(data_dir, "methods-metadata.tsv"), show_col_types = FALSE)
datasets_meta <- read_tsv(file.path(data_dir, "datasets-metadata.tsv"), show_col_types = FALSE)
selected_meta <- read_tsv(
    file.path("results", "fig4a_spatial_benchmark", "figures", "representative_settings_selected.tsv"),
    show_col_types = FALSE
)

dataset_names <- datasets_meta$Name
names(dataset_names) <- datasets_meta$Dataset

available_by_dataset <- benchmark |>
    distinct(.data$Dataset, .data$Method) |>
    group_by(.data$Dataset) |>
    summarise(Methods = list(.data$Method), .groups = "drop")

common_methods <- selected_meta$Method
for (i in seq_len(nrow(available_by_dataset))) {
    common_methods <- intersect(common_methods, available_by_dataset$Methods[[i]])
}
common_methods <- selected_meta$Method[selected_meta$Method %in% common_methods]

method_names <- methods_meta$Name
names(method_names) <- methods_meta$Method

method_labels <- method_names[common_methods]
method_labels[is.na(method_labels)] <- common_methods[is.na(method_labels)]
method_labels <- str_remove(method_labels, "\\s*\\(N=2000\\)")
method_labels <- str_replace(method_labels, "All features \\(N=all\\)", "All features")
method_labels <- str_replace_all(method_labels, "\\s*\\(\\)", "")

method_folder_info <- function(method) {
    if (method == "all_features-Nall") {
        return(list(base = "all_features", n_dir = "n1"))
    }
    m <- str_match(method, "^(.*)-N([^\\-]+)$")
    list(base = m[, 2], n_dir = paste0("n", m[, 3]))
}

load_feature_set <- function(dataset, method) {
    info <- method_folder_info(method)
    path <- file.path(
        dataset_dir_map[[dataset]],
        "feature_selection",
        info$base,
        info$n_dir,
        "seed0",
        "selected_features.json"
    )
    if (!file.exists(path)) return(NULL)
    obj <- jsonlite::fromJSON(path)
    unique(obj$feature_names)
}

selected_list <- map(set_names(names(dataset_dir_map)), function(dataset) {
    feats <- map(common_methods, ~ load_feature_set(dataset, .x))
    feats <- feats[!vapply(feats, is.null, logical(1))]
    kept_methods <- common_methods[!vapply(map(common_methods, ~ load_feature_set(dataset, .x)), is.null, logical(1))]
    names(feats) <- kept_methods
    feats
})
selected_list <- selected_list[lengths(selected_list) > 0]

overlaps <- map_dfr(names(selected_list), function(.dataset) {
    selected <- selected_list[[.dataset]]
    methods_here <- names(selected)
    expand_grid(Method1 = methods_here, Method2 = methods_here) |>
        filter(
            str_detect(.data$Method1, "^random", negate = TRUE),
            str_detect(.data$Method2, "^random", negate = TRUE)
        ) |>
        pmap_dfr(function(Method1, Method2) {
            selected1 <- selected[[Method1]]
            selected2 <- selected[[Method2]]
            both <- length(intersect(selected1, selected2))
            either <- length(union(selected1, selected2))
            tibble(
                Dataset = .dataset,
                Method1 = Method1,
                Method2 = Method2,
                Jaccard = ifelse(either > 0, both / either, NA_real_)
            )
        })
})

overlap_means <- overlaps |>
    group_by(.data$Method1, .data$Method2) |>
    summarise(
        MeanJaccard = mean(.data$Jaccard, na.rm = TRUE),
        SDJaccard = sd(.data$Jaccard, na.rm = TRUE),
        .groups = "drop"
    )

clust <- overlap_means |>
    select(.data$Method1, .data$Method2, .data$MeanJaccard) |>
    pivot_wider(names_from = .data$Method2, values_from = .data$MeanJaccard) |>
    column_to_rownames("Method1") |>
    as.matrix() |>
    dist() |>
    hclust()
method_order <- clust$labels[clust$order]

label_map <- method_labels[method_order]

overall_plot_df <- overlap_means |>
    mutate(
        Method1 = factor(.data$Method1, levels = method_order, labels = label_map),
        Method2 = factor(.data$Method2, levels = rev(method_order), labels = rev(label_map))
    )

dataset_plot_df <- overlaps |>
    mutate(
        Dataset = factor(.data$Dataset, levels = datasets_meta$Dataset, labels = dataset_names[datasets_meta$Dataset]),
        Method1 = factor(.data$Method1, levels = method_order, labels = label_map),
        Method2 = factor(.data$Method2, levels = rev(method_order), labels = rev(label_map))
    )

make_title_bar <- function(text) {
    ggplot() +
        annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1, fill = "black", colour = "black") +
        annotate("text", x = 0.5, y = 0.5, label = text, colour = "white", size = 2.4, fontface = "plain") +
        coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE, clip = "off") +
        theme_void() +
        theme(plot.margin = margin(0, 0, 0, 0))
}

mean_plot <- ggplot(overall_plot_df, aes(x = .data$Method1, y = .data$Method2, fill = .data$MeanJaccard)) +
    geom_tile() +
    scale_fill_viridis_c(option = "plasma", limits = c(0, 1)) +
    coord_equal() +
    labs(title = NULL, fill = "Mean\nJaccard\nIndex") +
    theme_features_pub() +
    theme(
        legend.position = "bottom",
        legend.title.position = "top",
        axis.text.x = element_text(size = 5.2, angle = 90, hjust = 0, vjust = 0.5, face = "bold", colour = "black"),
        axis.text.y = element_text(size = 5.2, face = "bold", colour = "black"),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        plot.margin = margin(0, 0, 0, 0),
        legend.key.width = unit(0.55, "cm"),
        legend.key.height = unit(0.18, "cm"),
        legend.text = element_text(size = 4.8),
        legend.title = element_text(size = 5.0)
    )

sd_plot <- ggplot(overall_plot_df, aes(x = .data$Method1, y = .data$Method2, fill = .data$SDJaccard)) +
    geom_tile() +
    scale_fill_viridis_c(option = "cividis", limits = c(0, max(overall_plot_df$SDJaccard, na.rm = TRUE))) +
    coord_equal() +
    labs(title = NULL, fill = "Jaccard\nIndex\nSD") +
    theme_features_pub() +
    theme(
        legend.position = "bottom",
        legend.title.position = "top",
        axis.text.x = element_text(size = 5.2, angle = 90, hjust = 0, vjust = 0.5, face = "bold", colour = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        plot.margin = margin(0, 0, 0, 0),
        legend.key.width = unit(0.55, "cm"),
        legend.key.height = unit(0.18, "cm"),
        legend.text = element_text(size = 4.8),
        legend.title = element_text(size = 5.0)
    )

dataset_heatmaps <- map(levels(dataset_plot_df$Dataset), function(ds) {
    base_plot <- ggplot(filter(dataset_plot_df, .data$Dataset == ds), aes(x = .data$Method1, y = .data$Method2, fill = .data$Jaccard)) +
        geom_tile() +
        scale_fill_viridis_c(option = "plasma", limits = c(0, 1), guide = "none") +
        coord_equal() +
        theme_features_pub() +
        theme(
            panel.background = element_rect(fill = "black", colour = NA),
            panel.grid = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            plot.margin = margin(0, 0, 0, 0)
        )
    wrap_plots(
        make_title_bar(ds),
        base_plot,
        ncol = 1,
        heights = c(0.13, 1)
    )
})

top_row <- wrap_plots(dataset_heatmaps[1:4], nrow = 1)
bottom_row <- wrap_plots(dataset_heatmaps[5:8], nrow = 1)

legend_row <- wrap_plots(
    cowplot::get_legend(
        mean_plot +
            theme(
                legend.position = "bottom",
                legend.margin = margin(0, 0, 0, 0),
                legend.box.margin = margin(0, 0, 0, 0)
            )
    ),
    cowplot::get_legend(
        sd_plot +
            theme(
                legend.position = "bottom",
                legend.margin = margin(0, 0, 0, 0),
                legend.box.margin = margin(0, 0, 0, 0)
            )
    ),
    nrow = 1,
    widths = c(1, 1)
)

panel_a <- wrap_plots(
    mean_plot + theme(legend.position = "none"),
    sd_plot + theme(legend.position = "none"),
    nrow = 1,
    widths = c(1, 1)
)
panel_b <- wrap_plots(
    top_row,
    bottom_row,
    ncol = 1,
    heights = c(1, 1)
) &
    theme(plot.margin = margin(0.05, 0, 0.05, 0, "cm"))

figure <- wrap_plots(
    panel_a,
    plot_spacer(),
    panel_b,
    legend_row,
    ncol = 1,
    heights = c(1.02, 0.10, 1.60, 0.16)
) &
    theme(plot.margin = margin(0, 0, 0, 0))

write_tsv(overlap_means, file.path(output_dir, "extended_fig4_overlap_means.tsv"))
write_tsv(overlaps, file.path(output_dir, "extended_fig4_overlap_by_dataset.tsv"))

save_figure_files(
    figure,
    file.path(output_dir, "extended_data_fig4"),
    width = 8.3,
    height = 8.9
)
