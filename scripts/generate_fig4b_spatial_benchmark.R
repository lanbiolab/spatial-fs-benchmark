#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(patchwork)
    library(jsonlite)
    library(purrr)
    library(stringr)
    library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[[1]] else file.path("results", "current_rank", "data")
results_root <- if (length(args) >= 2) args[[2]] else file.path("results", "spatial_main_native_seed0_fix3")
output_dir <- if (length(args) >= 3) args[[3]] else file.path("results", "fig4b_spatial_benchmark", "figures")

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
    if (!file.exists(path)) {
        return(NULL)
    }
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
    expand_grid(
        method1 = methods_here,
        method2 = methods_here
    ) |>
        filter(
            str_detect(.data$method1, "^random", negate = TRUE),
            str_detect(.data$method2, "^random", negate = TRUE)
        ) |>
        pmap_dfr(function(method1, method2) {
            selected1 <- selected[[method1]]
            selected2 <- selected[[method2]]
            both <- length(intersect(selected1, selected2))
            either <- length(union(selected1, selected2))
            subset <- length(intersect(selected1, selected2))
            tibble(
                Dataset = .dataset,
                Method1 = method1,
                Method2 = method2,
                Both = both,
                Either = either,
                Jaccard = ifelse(either > 0, both / either, NA_real_),
                Prop = ifelse(length(selected1) > 0, subset / length(selected1), NA_real_)
            )
        })
})

overlap_means <- overlaps |>
    group_by(.data$Method1, .data$Method2) |>
    summarise(
        JaccardSD = sd(.data$Jaccard, na.rm = TRUE),
        Jaccard = mean(.data$Jaccard, na.rm = TRUE),
        PropSD = sd(.data$Prop, na.rm = TRUE),
        Prop = mean(.data$Prop, na.rm = TRUE),
        .groups = "drop"
    )

overlaps_clust <- overlap_means |>
    select(.data$Method1, .data$Method2, .data$Jaccard) |>
    pivot_wider(names_from = .data$Method2, values_from = .data$Jaccard) |>
    column_to_rownames("Method1") |>
    as.matrix() |>
    dist() |>
    hclust()
overlaps_clust_order <- overlaps_clust$labels[overlaps_clust$order]

overlaps_plotting <- overlap_means |>
    mutate(
        Method1 = factor(.data$Method1, levels = overlaps_clust_order, labels = method_labels[overlaps_clust_order]),
        Method2 = factor(.data$Method2, levels = overlaps_clust_order, labels = method_labels[overlaps_clust_order])
    )

mean_overlaps_plot <- ggplot(
    overlaps_plotting,
    aes(
        x = .data$Method1, y = .data$Method2,
        fill = .data$Jaccard, colour = .data$Jaccard > 0.5, size = .data$JaccardSD
    )
) +
    geom_point(shape = 22, stroke = 0.25) +
    scale_fill_viridis_c(
        option = "plasma",
        begin = 0.08,
        end = 0.98,
        limits = c(0, 0.4),
        oob = scales::squish
    ) +
    scale_colour_manual(
        values = c("FALSE" = "black", "TRUE" = "white"),
        breaks = c("FALSE", "TRUE"),
        labels = c("FALSE", "TRUE")
    ) +
    scale_size_continuous(trans = "reverse", range = c(0.1, 1.8), breaks = scales::breaks_extended(n = 3)) +
    coord_equal() +
    labs(
        title = NULL,
        fill = "Mean\nJaccard\nIndex",
        colour = "JI > 0.5",
        size = "Standard\ndeviation"
    ) +
    guides(
        fill = guide_colourbar(
            order = 1,
            barheight = unit(0.22, "cm"),
            barwidth = unit(1.7, "cm"),
            title.position = "top"
        ),
        colour = guide_legend(
            order = 2,
            theme = theme(legend.text.position = "bottom"),
            override.aes = list(
                shape = 22,
                fill = c("black", "white"),
                colour = "grey25",
                size = 3.2,
                stroke = 0.35
            )
        ),
        size = guide_legend(
            order = 3,
            theme = theme(legend.text.position = "bottom"),
            override.aes = list(
                shape = 22,
                fill = "white",
                colour = "grey25",
                stroke = 0.35
            )
        )
    ) +
    theme_features_pub() +
    theme(
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.title.position = "top",
        legend.key.width = unit(0.32, "cm"),
        legend.key.height = unit(0.26, "cm"),
        legend.spacing.x = unit(0.16, "cm"),
        legend.box.spacing = unit(0.25, "cm"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(face = "bold"),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = "#1A1024", colour = NA),
        panel.grid = element_blank()
    )

num_selected <- map_dfr(names(selected_list), function(.dataset) {
    selected <- selected_list[[.dataset]]
    tibble(
        Dataset = .dataset,
        Method = names(selected),
        Selected = lengths(selected)
    )
})

num_selected_means <- num_selected |>
    group_by(.data$Method) |>
    summarise(
        MeanSelected = mean(.data$Selected),
        SDSelected = sd(.data$Selected),
        .groups = "drop"
    )

fig4d_methods <- c(
    "all_features-Nall",
    "TFs-N2000",
    "seurat_mvp-N2000",
    "dubstepr-N5000",
    "scsegindex-N10000",
    "osca-N10000",
    "wilcoxon-N10000",
    "triku-N2000"
)

fig4d_method_labels <- c(
    "all_features-Nall" = "All\nfeatures",
    "TFs-N2000" = "Transcription\nfactors",
    "seurat_mvp-N2000" = "Seurat-MVP",
    "dubstepr-N5000" = "DUBStepR",
    "scsegindex-N10000" = "scSEGIndex",
    "osca-N10000" = "OSCA",
    "wilcoxon-N10000" = "Wilcoxon",
    "triku-N2000" = "triku"
)

fig4d_display_order <- c(
    "All\nfeatures",
    "Transcription\nfactors",
    "Seurat-MVP",
    "DUBStepR",
    "scSEGIndex",
    "OSCA",
    "Wilcoxon",
    "triku"
)

num_selected_d <- map_dfr(names(dataset_dir_map), function(.dataset) {
    map_dfr(fig4d_methods, function(.method) {
        feats <- load_feature_set(.dataset, .method)
        if (is.null(feats)) {
            return(tibble())
        }
        tibble(
            Dataset = .dataset,
            Method = .method,
            Selected = length(feats)
        )
    })
})

num_selected_means_d <- num_selected_d |>
    group_by(.data$Method) |>
    summarise(
        MeanSelected = mean(.data$Selected),
        SDSelected = sd(.data$Selected),
        .groups = "drop"
    )

num_methods <- map_dfr(names(selected_list), function(.dataset) {
    selected <- selected_list[[.dataset]]
    counts <- table(unlist(selected))
    count_table <- table(counts)
    tibble(
        Dataset = .dataset,
        NumMethods = as.numeric(names(count_table)),
        NumFeatures = as.numeric(count_table)
    )
})

cum_num_methods <- num_methods |>
    arrange(.data$Dataset, desc(.data$NumMethods)) |>
    group_by(.data$Dataset) |>
    mutate(CumFeatures = cumsum(.data$NumFeatures)) |>
    ungroup()

cum_num_methods_plotting <- cum_num_methods |>
    filter(.data$NumMethods %in% c(5, 10, 15, 20)) |>
    mutate(
        Dataset = factor(.data$Dataset, levels = names(dataset_names), labels = dataset_names),
        NumMethodsLegend = factor(
            .data$NumMethods,
            levels = c(5, 10, 15),
            labels = c("5", "10", "15")
        )
    )

num_methods_plot <- ggplot(
    cum_num_methods_plotting,
    aes(y = .data$Dataset, x = .data$CumFeatures)
) +
    geom_point(
        aes(fill = .data$NumMethodsLegend),
        data = dplyr::filter(cum_num_methods_plotting, .data$NumMethods %in% c(5, 10, 15)),
        size = 2.1, shape = 21, colour = "grey20", stroke = 0.25
    ) +
    geom_point(
        aes(fill = factor(.data$NumMethods)),
        data = dplyr::filter(cum_num_methods_plotting, .data$NumMethods == 20),
        size = 2.1, shape = 21, colour = "grey20", stroke = 0.25,
        show.legend = FALSE
    ) +
    scale_fill_viridis_d(option = "plasma", breaks = c("5", "10", "15"), drop = FALSE) +
    scale_x_log10(breaks = c(10, 100, 500, 1000, 5000, 10000), limits = c(10, NA)) +
    labs(
        x = "Number of selected features (log scale)",
        fill = "Number of\nmethods"
    ) +
    theme_features_pub() +
    theme(
        axis.text.y = element_text(face = "bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.y = element_blank(),
        panel.grid.minor = element_blank()
    )

num_selected_means_plotting <- num_selected_means_d |>
    mutate(
        MethodLabel = factor(
            .data$Method,
            levels = names(fig4d_method_labels),
            labels = unname(fig4d_method_labels)
        )
    )

num_selected_means_plotting$MethodLabel <- factor(
    num_selected_means_plotting$MethodLabel,
    levels = rev(fig4d_display_order)
)

num_selected_plotting <- num_selected_d |>
    filter(.data$Method %in% num_selected_means_plotting$Method) |>
    mutate(
        MethodLabel = factor(
            .data$Method,
            levels = num_selected_means_plotting$Method,
            labels = fig4d_method_labels[num_selected_means_plotting$Method]
        ),
        Dataset = factor(.data$Dataset, levels = names(dataset_names), labels = dataset_names)
    )

num_selected_plotting$MethodLabel <- factor(
    num_selected_plotting$MethodLabel,
    levels = rev(fig4d_display_order)
)

num_selected_plot <- ggplot(
    num_selected_plotting,
    aes(x = .data$Selected, y = .data$MethodLabel)
) +
    geom_vline(xintercept = 2000, colour = "red") +
    geom_point(
        aes(fill = .data$Dataset),
        position = position_jitter(width = 0, height = 0.2, seed = 1),
        size = 2.0, shape = 21, alpha = 0.9, colour = "grey20", stroke = 0.25
    ) +
    geom_point(
        data = num_selected_means_plotting,
        aes(x = .data$MeanSelected),
        shape = "|", size = 4, colour = "blue"
    ) +
    scale_x_log10(breaks = c(100, 200, 500, 1000, 2000, 5000, 10000, 20000)) +
    scale_fill_brewer(palette = "Set3") +
    labs(x = "Number of selected\nfeatures (log scale)") +
    theme_features_pub() +
    theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.minor.x = element_blank()
    )

figure <- wrap_plots(
    mean_overlaps_plot,
    num_methods_plot,
    num_selected_plot,
    nrow = 1,
    widths = c(1.2, 0.75, 0.95),
    guides = "collect"
) &
    theme(
        legend.position = "bottom",
        legend.title.position = "top"
    )

write_tsv(overlap_means, file.path(output_dir, "fig4b_overlap_means.tsv"))
write_tsv(num_selected, file.path(output_dir, "fig4b_num_selected.tsv"))
write_tsv(cum_num_methods, file.path(output_dir, "fig4b_cum_num_methods.tsv"))
write_tsv(tibble(Method = common_methods, Label = method_labels[common_methods]), file.path(output_dir, "fig4b_methods_used.tsv"))

save_figure_files(
    mean_overlaps_plot +
        theme(
            legend.position = "bottom",
            legend.box = "horizontal",
            legend.title.position = "top"
        ),
    file.path(output_dir, "figure_fig4b_panel_b"),
    width = 4.3,
    height = 4.2
)

save_figure_files(
    num_methods_plot +
        theme(
            legend.position = "right",
            legend.box = "vertical",
            legend.title.position = "top"
        ),
    file.path(output_dir, "figure_fig4c_panel_c"),
    width = 3.1,
    height = 3.6
)

save_figure_files(
    num_selected_plot +
        theme(
            legend.position = "right",
            legend.box = "vertical",
            legend.title.position = "top"
        ),
    file.path(output_dir, "figure_fig4d_panel_d"),
    width = 3.2,
    height = 3.6
)

save_figure_files(figure, file.path(output_dir, "figure_fig4b_spatial_benchmark"), width = 8.2, height = 4.7)
