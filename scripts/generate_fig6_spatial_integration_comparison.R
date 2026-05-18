#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(colorspace)
    library(patchwork)
    library(cowplot)
})

source(file.path("external", "atlas-feature-selection-benchmark", "analysis", "R", "plotting.R"))

args <- commandArgs(trailingOnly = TRUE)
results_path <- if (length(args) >= 1) args[[1]] else file.path("results", "spatial_main_native_seed0_fix3", "results_with_all_wilcoxon_datasets.csv")
output_dir <- if (length(args) >= 2) args[[2]] else file.path("results", "fig6_spatial_integration", "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

type_names <- c(
    "Overall" = "Overall",
    "Integration" = "Integration",
    "Clustering" = "Clustering",
    "Alignment" = "Alignment"
)

integration_labels <- c(
    "scvi" = "scVI",
    "cellcharter" = "Cell\nCharter",
    "gpsa" = "GPSA",
    "staligner" = "STAligner"
)

safe_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
safe_sd <- function(x) if (sum(!is.na(x)) <= 1) NA_real_ else sd(x, na.rm = TRUE)

score_direction <- function(metric_name) {
    lower_better <- c("chaos", "pas", "ratio")
    ifelse(metric_name %in% lower_better, -1, 1)
}

metrics <- read_csv(results_path, show_col_types = FALSE) |>
    transmute(
        Dataset = .data$dataset,
        MethodBase = .data$fs_method,
        Method = case_when(
            .data$fs_method == "all_features" ~ "All features",
            .data$fs_method == "random" ~ "Random",
            .data$fs_method == "TFs" ~ "Transcription factors",
            .data$fs_method == "scPNMF" ~ "scPNMF",
            .data$fs_method == "scsegindex" ~ "scSEGIndex",
            .data$fs_method == "seurat_disp" ~ "Seurat-Dispersion",
            .data$fs_method == "seurat_mvp" ~ "Seurat-MVP",
            .data$fs_method == "seurat_sct" ~ "Seurat-scTransform",
            .data$fs_method == "seurat_vst" ~ "Seurat-VST",
            .data$fs_method == "scanpy_cell_ranger" ~ "scanpy-CellRanger",
            .data$fs_method == "scanpy_cell_ranger_batch" ~ "scanpy-CellRanger (batch)",
            .data$fs_method == "scanpy_pearson" ~ "scanpy-Pearson",
            .data$fs_method == "scanpy_pearson_batch" ~ "scanpy-Pearson (batch)",
            .data$fs_method == "scanpy_seurat" ~ "scanpy-Seurat",
            .data$fs_method == "scanpy_seurat_batch" ~ "scanpy-Seurat (batch)",
            .data$fs_method == "scanpy_seurat_v3" ~ "scanpy-SeuratV3",
            .data$fs_method == "scanpy_seurat_v3_batch" ~ "scanpy-SeuratV3 (batch)",
            .data$fs_method == "singleCellHaystack" ~ "singleCellHaystack",
            .data$fs_method == "statistic_mean" ~ "Statistic mean",
            .data$fs_method == "statistic_variance" ~ "Statistic variance",
            .data$fs_method == "dubstepr" ~ "DUBStepR",
            .data$fs_method == "hotspot" ~ "Hotspot",
            .data$fs_method == "nbumi" ~ "NBumi",
            .data$fs_method == "osca" ~ "OSCA",
            .data$fs_method == "triku" ~ "triku",
            .data$fs_method == "wilcoxon" ~ "Wilcoxon",
            .data$fs_method == "anticor" ~ "Anticor",
            TRUE ~ .data$fs_method
        ),
        IntegrationMethod = .data$integration_method,
        N = as.character(.data$n_features),
        Task = case_when(
            .data$task == "integration_eval" ~ "Integration",
            .data$task == "clustering_eval" ~ "Clustering",
            .data$task == "alignment_eval" ~ "Alignment",
            TRUE ~ NA_character_
        ),
        Metric = tolower(.data$metric_name),
        ValueRaw = .data$metric_value
    ) |>
    filter(!is.na(.data$Task), .data$IntegrationMethod %in% names(integration_labels), !is.na(.data$ValueRaw)) |>
    mutate(
        Keep = case_when(
            .data$MethodBase == "all_features" ~ .data$N == "all",
            .data$MethodBase == "random" ~ .data$N == "500",
            TRUE ~ .data$N == "2000"
        )
    ) |>
    filter(.data$Keep) |>
    mutate(
        Value = .data$ValueRaw * score_direction(.data$Metric)
    )

common_methods <- metrics |>
    distinct(.data$MethodBase, .data$Method, .data$IntegrationMethod) |>
    count(.data$MethodBase, .data$Method, name = "n_integrations") |>
    filter(.data$n_integrations == 4) |>
    arrange(.data$MethodBase) |>
    pull(.data$Method)

metrics <- metrics |>
    filter(.data$Method %in% common_methods)

scaled <- metrics |>
    group_by(.data$Dataset, .data$Metric, .data$Task) |>
    filter(any(!is.na(.data$Value))) |>
    mutate(
        Lower = min(.data$Value, na.rm = TRUE),
        Upper = max(.data$Value, na.rm = TRUE),
        Range = ifelse(.data$Upper > .data$Lower, .data$Upper - .data$Lower, 1),
        Scaled = (.data$Value - .data$Lower) / .data$Range
    ) |>
    ungroup() |>
    filter(is.finite(.data$Scaled))

summary_scores <- scaled |>
    group_by(.data$Dataset, .data$Method, .data$IntegrationMethod, .data$Task) |>
    summarise(Value = mean(.data$Scaled, na.rm = TRUE), .groups = "drop") |>
    pivot_wider(names_from = "Task", values_from = "Value")

for (nm in c("Integration", "Clustering", "Alignment")) {
    if (!nm %in% names(summary_scores)) {
        summary_scores[[nm]] <- NA_real_
    }
}

summary_scores <- summary_scores |>
    mutate(
        Overall = if_else(
            !is.na(.data$Integration) & !is.na(.data$Clustering) & !is.na(.data$Alignment),
            (.data$Integration + .data$Clustering + .data$Alignment) / 3,
            NA_real_
        )
    ) |>
    pivot_longer(cols = c("Overall", "Integration", "Clustering", "Alignment"), names_to = "Type", values_to = "Value")

methods_order <- summary_scores |>
    filter(.data$IntegrationMethod == "scvi") |>
    group_by(.data$Method, .data$Type) |>
    summarise(Mean = mean(.data$Value, na.rm = TRUE), .groups = "drop") |>
    filter(.data$Type == "Overall") |>
    arrange(desc(.data$Mean)) |>
    pull(.data$Method)

summary_scores_means <- summary_scores |>
    group_by(.data$Method, .data$IntegrationMethod, .data$Type) |>
    summarise(
        Mean = safe_mean(.data$Value),
        SD = safe_sd(.data$Value),
        .groups = "drop"
    ) |>
    mutate(
        Method = factor(.data$Method, levels = rev(methods_order)),
        Integration = factor(.data$IntegrationMethod, levels = names(integration_labels), labels = integration_labels),
        Type = factor(.data$Type, levels = names(type_names), labels = type_names)
    )

score_diffs <- summary_scores |>
    select(.data$Dataset, .data$Method, .data$IntegrationMethod, .data$Type, .data$Value) |>
    left_join(
        summary_scores |>
            filter(.data$IntegrationMethod == "scvi") |>
            select(.data$Dataset, .data$Method, .data$Type, ScviValue = .data$Value),
        by = c("Dataset", "Method", "Type")
    ) |>
    filter(.data$IntegrationMethod %in% c("cellcharter", "gpsa", "staligner")) |>
    mutate(Diff = .data$Value - .data$ScviValue) |>
    group_by(.data$Method, .data$IntegrationMethod, .data$Type) |>
    summarise(
        Mean = safe_mean(.data$Diff),
        SD = safe_sd(.data$Diff),
        .groups = "drop"
    ) |>
    mutate(
        Method = factor(.data$Method, levels = rev(methods_order)),
        Integration = factor(.data$IntegrationMethod, levels = c("cellcharter", "gpsa", "staligner"), labels = integration_labels[c("cellcharter", "gpsa", "staligner")]),
        Type = factor(.data$Type, levels = names(type_names), labels = type_names)
    ) |>
    select(-.data$IntegrationMethod)

rank_data <- summary_scores |>
    group_by(.data$Dataset, .data$IntegrationMethod, .data$Type) |>
    mutate(Rank = rank(-.data$Value, ties.method = "average")) |>
    ungroup()

rank_means <- rank_data |>
    group_by(.data$Method, .data$IntegrationMethod, .data$Type) |>
    summarise(
        MeanRank = safe_mean(.data$Rank),
        SDRank = safe_sd(.data$Rank),
        .groups = "drop"
    ) |>
    mutate(
        Method = factor(.data$Method, levels = rev(methods_order)),
        Integration = factor(.data$IntegrationMethod, levels = names(integration_labels), labels = integration_labels),
        Type = factor(.data$Type, levels = names(type_names), labels = type_names)
    )

rank_diffs <- rank_data |>
    select(.data$Dataset, .data$Method, .data$IntegrationMethod, .data$Type, .data$Rank) |>
    left_join(
        rank_data |>
            filter(.data$IntegrationMethod == "scvi") |>
            select(.data$Dataset, .data$Method, .data$Type, ScviRank = .data$Rank),
        by = c("Dataset", "Method", "Type")
    ) |>
    filter(.data$IntegrationMethod %in% c("cellcharter", "gpsa", "staligner")) |>
    mutate(Diff = .data$Rank - .data$ScviRank) |>
    group_by(.data$Method, .data$IntegrationMethod, .data$Type) |>
    summarise(
        Mean = safe_mean(.data$Diff),
        SD = safe_sd(.data$Diff),
        .groups = "drop"
    ) |>
    mutate(
        Method = factor(.data$Method, levels = rev(methods_order)),
        Integration = factor(.data$IntegrationMethod, levels = c("cellcharter", "gpsa", "staligner"), labels = integration_labels[c("cellcharter", "gpsa", "staligner")]),
        Type = factor(.data$Type, levels = names(type_names), labels = type_names)
    ) |>
    select(-.data$IntegrationMethod)

theme_integration <- theme_features_pub() +
    theme(
        axis.title = element_blank(),
        axis.text.x = element_text(size = 5.2, angle = 90, hjust = 1, vjust = 0.5, face = "bold", colour = "black"),
        axis.text.y = element_text(size = 5.4, face = "bold", colour = "black"),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "black", colour = "black"),
        strip.text = element_text(colour = "white", face = "bold")
    )

plot_scores <- ggplot(summary_scores_means, aes(x = .data$Type, y = .data$Method, colour = .data$Mean, size = .data$SD)) +
    geom_point(shape = 15) +
    scale_colour_viridis_c(option = "magma", limits = c(min(0, min(summary_scores_means$Mean, na.rm = TRUE)), max(summary_scores_means$Mean, na.rm = TRUE)), na.value = "grey90") +
    scale_size_continuous(trans = "reverse", range = c(0.5, 3.2)) +
    facet_wrap(~ .data$Integration, nrow = 1) +
    labs(colour = "Mean score", size = "s.d.") +
    theme_integration

plot_score_diffs <- ggplot(score_diffs, aes(x = .data$Type, y = .data$Method, colour = .data$Mean, size = .data$SD)) +
    geom_point(shape = 15) +
    colorspace::scale_colour_continuous_diverging(palette = "Purple-Green", limits = c(-max(abs(score_diffs$Mean), na.rm = TRUE), max(abs(score_diffs$Mean), na.rm = TRUE)), na.value = "grey90") +
    scale_size_continuous(trans = "reverse", range = c(0.5, 3.2)) +
    facet_wrap(~ .data$Integration, nrow = 1) +
    labs(colour = "Mean difference", size = "s.d.") +
    theme_integration +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

plot_ranks <- ggplot(rank_means, aes(x = .data$Type, y = .data$Method, colour = .data$Type, alpha = .data$MeanRank, size = .data$SDRank)) +
    geom_point(shape = 15) +
    scale_colour_manual(values = c("Overall" = "#f781bf", "Integration" = "#e41a1c", "Clustering" = "#377eb8", "Alignment" = "#4daf4a"), guide = "none") +
    scale_alpha_continuous(
        limits = c(max(rank_means$MeanRank, na.rm = TRUE), 1),
        range = c(0.18, 1.0),
        trans = "reverse",
        breaks = c(5, 10, 15, 20, 25)
    ) +
    scale_size_continuous(trans = "reverse", range = c(0.5, 3.2)) +
    facet_wrap(~ .data$Integration, nrow = 1) +
    labs(alpha = "Mean rank", size = "Rank s.d.") +
    theme_integration +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

plot_rank_diffs <- ggplot(rank_diffs, aes(x = .data$Type, y = .data$Method, colour = .data$Mean, size = .data$SD)) +
    geom_point(shape = 15) +
    colorspace::scale_colour_continuous_diverging(palette = "Tropic", rev = TRUE, limits = c(-max(abs(rank_diffs$Mean), na.rm = TRUE), max(abs(rank_diffs$Mean), na.rm = TRUE)), na.value = "grey90") +
    scale_size_continuous(trans = "reverse", range = c(0.5, 3.2)) +
    facet_wrap(~ .data$Integration, nrow = 1) +
    labs(colour = "Mean difference", size = "s.d.") +
    theme_integration +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

main_fig <- wrap_plots(
    plot_scores + theme(legend.position = "none"),
    plot_score_diffs + theme(legend.position = "none"),
    plot_ranks + theme(legend.position = "none"),
    plot_rank_diffs + theme(legend.position = "none"),
    nrow = 1,
    widths = c(3, 2, 3, 2)
)

legend_common <- theme(
    legend.position = "bottom",
    legend.title.position = "top",
    legend.box = "horizontal",
    legend.box.just = "left",
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7, face = "bold")
)

legend_score_colour <- cowplot::get_legend(
    plot_scores +
        guides(
            colour = guide_colourbar(title = "Mean score", barwidth = unit(1.5, "cm"), barheight = unit(0.22, "cm")),
            size = "none"
        ) +
        legend_common
)

legend_score_sd <- cowplot::get_legend(
    plot_scores +
        guides(
            colour = "none",
            size = guide_legend(title = "s.d.", nrow = 1)
        ) +
        legend_common
)

legend_diff_colour <- cowplot::get_legend(
    plot_score_diffs +
        guides(
            colour = guide_colourbar(title = "Mean difference to scVI", barwidth = unit(1.5, "cm"), barheight = unit(0.22, "cm")),
            size = "none"
        ) +
        legend_common
)

legend_diff_sd <- cowplot::get_legend(
    plot_score_diffs +
        guides(
            colour = "none",
            size = guide_legend(title = "s.d.", nrow = 1)
        ) +
        legend_common
)

legend_rank_alpha <- cowplot::get_legend(
    plot_ranks +
        guides(
            alpha = guide_legend(title = "Mean rank", nrow = 1),
            size = "none"
        ) +
        legend_common
)

legend_rank_sd <- cowplot::get_legend(
    plot_ranks +
        guides(
            alpha = "none",
            size = guide_legend(title = "Rank s.d.", nrow = 1)
        ) +
        legend_common
)

legend_rankdiff_colour <- cowplot::get_legend(
    plot_rank_diffs +
        guides(
            colour = guide_colourbar(title = "Mean difference to scVI", barwidth = unit(2.1, "cm"), barheight = unit(0.22, "cm")),
            size = "none"
        ) +
        legend_common
)

legend_rankdiff_sd <- cowplot::get_legend(
    plot_rank_diffs +
        guides(
            colour = "none",
            size = guide_legend(title = "s.d.", nrow = 1)
        ) +
        legend_common
)

legend_grid_row1 <- wrap_plots(
    patchwork::wrap_elements(legend_score_colour),
    patchwork::wrap_elements(legend_diff_colour),
    patchwork::wrap_elements(legend_rank_alpha),
    patchwork::wrap_elements(legend_rankdiff_colour),
    nrow = 1,
    widths = c(3, 2, 3, 2)
)

legend_grid_row2 <- wrap_plots(
    patchwork::wrap_elements(legend_score_sd),
    patchwork::wrap_elements(legend_diff_sd),
    patchwork::wrap_elements(legend_rank_sd),
    patchwork::wrap_elements(legend_rankdiff_sd),
    nrow = 1,
    widths = c(3, 2, 3, 2)
)

fig <- wrap_plots(
    main_fig,
    legend_grid_row1,
    legend_grid_row2,
    ncol = 1,
    heights = c(1, 0.085, 0.075)
) &
    theme(
        plot.margin = margin(2, 2, 2, 2)
    )

write_tsv(summary_scores_means, file.path(output_dir, "fig6_scores_means.tsv"))
write_tsv(score_diffs, file.path(output_dir, "fig6_score_diffs.tsv"))
write_tsv(rank_means, file.path(output_dir, "fig6_rank_means.tsv"))
write_tsv(rank_diffs, file.path(output_dir, "fig6_rank_diffs.tsv"))

save_figure_files(fig, file.path(output_dir, "figure_fig6_spatial_integration"), width = 8.2, height = 4.8)
