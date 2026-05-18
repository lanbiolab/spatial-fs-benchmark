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
output_dir <- if (length(args) >= 2) args[[2]] else file.path("results", "extended_data_fig6", "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

type_names <- c(
    "Overall" = "Overall",
    "Integration" = "Integration",
    "Clustering" = "Clustering",
    "Alignment" = "Alignment"
)

type_palette <- c(
    "Overall" = "#f781bf",
    "Integration" = "#e41a1c",
    "Clustering" = "#377eb8",
    "Alignment" = "#4daf4a"
)

integration_labels <- c(
    "scvi" = "scVI",
    "cellcharter" = "CellCharter",
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
    mutate(Value = .data$ValueRaw * score_direction(.data$Metric))

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
    if (!nm %in% names(summary_scores)) summary_scores[[nm]] <- NA_real_
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

scores_means <- summary_scores |>
    group_by(.data$Method, .data$IntegrationMethod, .data$Type) |>
    summarise(Mean = safe_mean(.data$Value), SD = safe_sd(.data$Value), .groups = "drop") |>
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
    summarise(Mean = safe_mean(.data$Diff), SD = safe_sd(.data$Diff), .groups = "drop") |>
    mutate(
        Method = factor(.data$Method, levels = rev(methods_order)),
        Integration = factor(.data$IntegrationMethod, levels = c("cellcharter", "gpsa", "staligner"), labels = integration_labels[c("cellcharter", "gpsa", "staligner")]),
        Type = factor(.data$Type, levels = names(type_names), labels = type_names)
    )

rank_data <- summary_scores |>
    group_by(.data$Dataset, .data$IntegrationMethod, .data$Type) |>
    mutate(Rank = rank(-.data$Value, ties.method = "average")) |>
    ungroup()

rank_means <- rank_data |>
    group_by(.data$Method, .data$IntegrationMethod, .data$Type) |>
    summarise(MeanRank = safe_mean(.data$Rank), SDRank = safe_sd(.data$Rank), .groups = "drop") |>
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
    summarise(Mean = safe_mean(.data$Diff), SD = safe_sd(.data$Diff), .groups = "drop") |>
    mutate(
        Method = factor(.data$Method, levels = rev(methods_order)),
        Integration = factor(.data$IntegrationMethod, levels = c("cellcharter", "gpsa", "staligner"), labels = integration_labels[c("cellcharter", "gpsa", "staligner")]),
        Type = factor(.data$Type, levels = names(type_names), labels = type_names)
    )

theme_ext6 <- theme_features_pub() +
    theme(
        axis.title = element_blank(),
        axis.text.x = element_text(size = 5.3, angle = 90, hjust = 1, vjust = 0.5, face = "bold", colour = "black"),
        axis.text.y = element_text(size = 5.0, face = "bold", colour = "black"),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "black", colour = "black"),
        strip.text = element_text(colour = "white", face = "bold", size = 5.8),
        plot.margin = margin(0, -2, 0, -2)
    )

score_lim <- c(min(0, min(scores_means$Mean, na.rm = TRUE)), max(scores_means$Mean, na.rm = TRUE))
score_diff_lim <- max(abs(score_diffs$Mean), na.rm = TRUE)
rank_diff_lim <- max(abs(rank_diffs$Mean), na.rm = TRUE)
sd_score_lim <- c(0, max(c(scores_means$SD, score_diffs$SD), na.rm = TRUE))
sd_rank_lim <- c(0, max(c(rank_means$SDRank, rank_diffs$SD), na.rm = TRUE))

plot_scores_mean <- ggplot(scores_means, aes(x = .data$Type, y = .data$Method, colour = .data$Mean)) +
    geom_point(shape = 15, size = 1.55) +
    scale_colour_viridis_c(option = "magma", limits = score_lim, na.value = "grey90") +
    facet_wrap(~ .data$Integration, nrow = 1) +
    labs(title = "Scores", colour = "Mean score") +
    theme_ext6

plot_score_diff_mean <- ggplot(score_diffs, aes(x = .data$Type, y = .data$Method, colour = .data$Mean)) +
    geom_point(shape = 15, size = 1.55) +
    colorspace::scale_colour_continuous_diverging(palette = "Purple-Green", limits = c(-score_diff_lim, score_diff_lim), na.value = "grey90") +
    facet_wrap(~ .data$Integration, nrow = 1) +
    labs(title = "Difference in scores", colour = "Mean difference") +
    theme_ext6 +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

plot_ranks_mean <- ggplot(rank_means, aes(x = .data$Type, y = .data$Method, colour = .data$Type, alpha = .data$MeanRank)) +
    geom_point(shape = 15, size = 1.55) +
    scale_colour_manual(values = type_palette, guide = "none") +
    scale_alpha_continuous(
        limits = c(max(rank_means$MeanRank, na.rm = TRUE), 1),
        range = c(0.12, 1.0),
        trans = "reverse"
    ) +
    facet_wrap(~ .data$Integration, nrow = 1) +
    labs(title = "Ranks", alpha = "Mean rank") +
    theme_ext6 +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

plot_rank_diff_mean <- ggplot(rank_diffs, aes(x = .data$Type, y = .data$Method, colour = .data$Mean)) +
    geom_point(shape = 15, size = 1.55) +
    colorspace::scale_colour_continuous_diverging(palette = "Tropic", rev = TRUE, limits = c(-rank_diff_lim, rank_diff_lim), na.value = "grey90") +
    facet_wrap(~ .data$Integration, nrow = 1) +
    labs(title = "Difference in ranks", colour = "Mean difference") +
    theme_ext6 +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

plot_scores_sd <- ggplot(scores_means, aes(x = .data$Type, y = .data$Method, colour = .data$SD)) +
    geom_point(shape = 15, size = 1.55) +
    scale_colour_viridis_c(option = "cividis", limits = sd_score_lim, na.value = "grey90") +
    facet_wrap(~ .data$Integration, nrow = 1) +
    labs(colour = "Standard deviation") +
    theme_ext6

plot_score_diff_sd <- ggplot(score_diffs, aes(x = .data$Type, y = .data$Method, colour = .data$SD)) +
    geom_point(shape = 15, size = 1.55) +
    scale_colour_viridis_c(option = "cividis", limits = sd_score_lim, na.value = "grey90") +
    facet_wrap(~ .data$Integration, nrow = 1) +
    labs(colour = "Difference SD") +
    theme_ext6 +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

plot_ranks_sd <- ggplot(rank_means, aes(x = .data$Type, y = .data$Method, colour = .data$SDRank)) +
    geom_point(shape = 15, size = 1.55) +
    scale_colour_viridis_c(option = "cividis", limits = sd_rank_lim, na.value = "grey90") +
    facet_wrap(~ .data$Integration, nrow = 1) +
    labs(colour = "Rank SD") +
    theme_ext6 +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

plot_rank_diff_sd <- ggplot(rank_diffs, aes(x = .data$Type, y = .data$Method, colour = .data$SD)) +
    geom_point(shape = 15, size = 1.55) +
    scale_colour_viridis_c(option = "cividis", limits = sd_rank_lim, na.value = "grey90") +
    facet_wrap(~ .data$Integration, nrow = 1) +
    labs(colour = "Difference SD") +
    theme_ext6 +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

row_a <- wrap_plots(
    plot_scores_mean + theme(legend.position = "none"),
    plot_score_diff_mean + theme(legend.position = "none"),
    plot_ranks_mean + theme(legend.position = "none"),
    plot_rank_diff_mean + theme(legend.position = "none"),
    nrow = 1,
    widths = c(3.2, 2.2, 3.2, 2.2)
)

row_b <- wrap_plots(
    plot_scores_sd + theme(legend.position = "none"),
    plot_score_diff_sd + theme(legend.position = "none"),
    plot_ranks_sd + theme(legend.position = "none"),
    plot_rank_diff_sd + theme(legend.position = "none"),
    nrow = 1,
    widths = c(3.2, 2.2, 3.2, 2.2)
)

legend_row <- wrap_plots(
    cowplot::get_legend(plot_scores_mean + theme(legend.position = "bottom")),
    cowplot::get_legend(plot_score_diff_mean + theme(legend.position = "bottom")),
    cowplot::get_legend(plot_ranks_mean + theme(legend.position = "bottom")),
    cowplot::get_legend(plot_rank_diff_mean + theme(legend.position = "bottom")),
    cowplot::get_legend(plot_scores_sd + theme(legend.position = "bottom")),
    cowplot::get_legend(plot_score_diff_sd + theme(legend.position = "bottom")),
    cowplot::get_legend(plot_ranks_sd + theme(legend.position = "bottom")),
    cowplot::get_legend(plot_rank_diff_sd + theme(legend.position = "bottom")),
    nrow = 2,
    widths = c(3.2, 2.2, 3.2, 2.2)
)

figure <- wrap_plots(
    row_a,
    row_b,
    legend_row,
    ncol = 1,
    heights = c(1, 1, 0.28)
) &
    theme(plot.margin = margin(0, 0, 0, 0))

write_tsv(scores_means, file.path(output_dir, "extended_fig6_scores_means.tsv"))
write_tsv(score_diffs, file.path(output_dir, "extended_fig6_score_diffs.tsv"))
write_tsv(rank_means, file.path(output_dir, "extended_fig6_rank_means.tsv"))
write_tsv(rank_diffs, file.path(output_dir, "extended_fig6_rank_diffs.tsv"))

save_figure_files(
    figure,
    file.path(output_dir, "extended_data_fig6"),
    width = 8.3,
    height = 8.3
)
