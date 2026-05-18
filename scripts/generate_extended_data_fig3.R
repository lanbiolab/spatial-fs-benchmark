#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(patchwork)
    library(stringr)
    library(colorspace)
})

args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[[1]] else file.path("results", "current_rank", "data")
output_dir <- if (length(args) >= 2) args[[2]] else file.path("results", "extended_data_fig3", "figures")
base_dir <- normalizePath(file.path(data_dir, ".."), mustWork = FALSE)
ranges_path <- file.path(base_dir, "output", "baseline-ranges.tsv")

source(file.path("external", "atlas-feature-selection-benchmark", "analysis", "R", "plotting.R"))
source(file.path("analysis", "atlas_style", "R", "spatial_summarisation.R"))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

types_palette <- c(
    "Overall" = "#f781bf",
    "Integration" = "#e41a1c",
    "Clustering" = "#377eb8",
    "Alignment" = "#4daf4a"
)

metrics <- read_tsv(file.path(data_dir, "benchmark.tsv"), show_col_types = FALSE) |>
    filter(.data$Type %in% c("Integration", "Clustering", "Alignment"))

metric_ranges <- read_tsv(ranges_path, show_col_types = FALSE) |>
    filter(.data$Type %in% c("Integration", "Clustering", "Alignment"))

methods_meta_all <- read_tsv(file.path(data_dir, "methods-metadata.tsv"), show_col_types = FALSE)
datasets_meta <- read_tsv(file.path(data_dir, "datasets-metadata.tsv"), show_col_types = FALSE)

dataset_names <- datasets_meta$Name
names(dataset_names) <- datasets_meta$Dataset

metrics_summary_all <- summarise_spatial_metrics(
    metrics,
    metric_ranges,
    type_weights = c("Integration" = 1/3, "Clustering" = 1/3, "Alignment" = 1/3),
    require_types_for_overall = c("Integration", "Clustering", "Alignment")
)

pick_representative_methods <- function(metrics_summary_all) {
    available <- metrics_summary_all |>
        distinct(.data$MethodBase, .data$Method)

    family_pick <- function(family, methods) {
        methods <- sort(unique(methods))
        if (family == "all_features") {
            return("all_features-Nall")
        }
        if (family == "random") {
            target <- paste0(family, "-N500")
            return(if (target %in% methods) target else methods[[1]])
        }
        target <- paste0(family, "-N2000")
        if (target %in% methods) return(target)
        preferred <- c(
            paste0(family, "-N1000"),
            paste0(family, "-N500"),
            paste0(family, "-N5000"),
            paste0(family, "-N10000"),
            paste0(family, "-N200"),
            paste0(family, "-N100")
        )
        hit <- preferred[preferred %in% methods]
        if (length(hit) > 0) hit[[1]] else methods[[1]]
    }

    available |>
        group_by(.data$MethodBase) |>
        summarise(
            Method = family_pick(.data$MethodBase[[1]], .data$Method),
            .groups = "drop"
        ) |>
        mutate(
            IsBaseline = .data$Method %in% c(
                "all_features-Nall",
                "random-N500",
                "scanpy_cell_ranger_batch-N2000",
                "scsegindex-N2000"
            )
        )
}

methods_meta <- pick_representative_methods(metrics_summary_all) |>
    left_join(methods_meta_all |> select(.data$Method, .data$Name), by = "Method") |>
    mutate(
        Name = case_when(
            !is.na(.data$Name) ~ .data$Name,
            .data$Method == "all_features-Nall" ~ "All features",
            TRUE ~ .data$Method
        )
    )

metrics_summary <- metrics_summary_all |>
    filter(.data$Method %in% methods_meta$Method, .data$IntegrationLabel == "scVI")

method_names <- methods_meta$Name
names(method_names) <- methods_meta$Method

plotting <- metrics_summary |>
    pivot_longer(
        cols = c("Overall", "Integration", "Clustering", "Alignment"),
        names_to = "Type",
        values_to = "Value"
    ) |>
    filter(!is.na(.data$Value)) |>
    mutate(
        Type = factor(.data$Type, levels = c("Overall", "Integration", "Clustering", "Alignment")),
        Dataset = factor(.data$Dataset, levels = datasets_meta$Dataset, labels = dataset_names[datasets_meta$Dataset])
    )

method_means <- plotting |>
    group_by(.data$Method, .data$Type) |>
    summarise(Value = mean(.data$Value), .groups = "drop")

methods_order <- method_means |>
    filter(.data$Type == "Overall") |>
    arrange(desc(.data$Value)) |>
    pull(.data$Method)

format_method <- function(x) {
    str_remove(x, " \\(N=2000\\)") |>
        str_remove("N=2000, ") |>
        str_remove(" \\(N=500\\)") |>
        str_remove(" \\(N=1000\\)") |>
        str_remove(" \\(N=5000\\)") |>
        str_remove(" \\(N=10000\\)") |>
        str_remove(" \\(N=all\\)") |>
        str_wrap(width = 16)
}

plotting <- plotting |>
    mutate(
        Method = factor(.data$Method, levels = rev(methods_order), labels = rev(format_method(method_names[methods_order])))
    )

baseline_methods <- methods_meta |>
    filter(.data$Method %in% methods_order) |>
    mutate(
        Method = factor(.data$Method, levels = rev(methods_order), labels = rev(format_method(method_names[methods_order])))
    )

score_limits <- max(abs(plotting$Value), na.rm = TRUE)

category_scores_plot <- ggplot(
    plotting,
    aes(x = .data$Dataset, y = .data$Method)
) +
    geom_linerange(
        data = baseline_methods |> filter(.data$IsBaseline),
        aes(x = NULL, y = .data$Method),
        xmin = -Inf, xmax = Inf,
        linewidth = 2.4, colour = "grey55", alpha = 0.28
    ) +
    geom_point(aes(colour = .data$Value), shape = 15, size = 1.18) +
    scale_y_discrete(drop = FALSE) +
    colorspace::scale_colour_continuous_divergingx(
        palette = "Zissou 1",
        limits = c(-score_limits, score_limits),
        name = "Metric category score"
    ) +
    facet_grid(. ~ .data$Type) +
    labs(title = "Metric category scores") +
    theme_features_pub() +
    theme(
        legend.position = "bottom",
        legend.title.position = "top",
        axis.text.x = element_text(size = 4.1, angle = 90, hjust = 1, vjust = 0.5, face = "bold", colour = "black"),
        axis.text.y = element_text(size = 4.2, face = "bold", colour = "black"),
        axis.title = element_blank(),
        plot.title = element_text(size = 6.6, margin = margin(b = 2)),
        panel.spacing = unit(0.05, "cm")
    )

metric_ranks <- metrics_summary |>
    pivot_longer(
        cols = c("Overall", "Integration", "Clustering", "Alignment"),
        names_to = "Type",
        values_to = "Value"
    ) |>
    filter(!is.na(.data$Value)) |>
    group_by(.data$Dataset, .data$Type) |>
    mutate(Rank = rank(-.data$Value, ties.method = "average")) |>
    ungroup() |>
    mutate(
        Type = factor(.data$Type, levels = c("Overall", "Integration", "Clustering", "Alignment")),
        Dataset = factor(.data$Dataset, levels = datasets_meta$Dataset, labels = dataset_names[datasets_meta$Dataset]),
        Method = factor(.data$Method, levels = rev(methods_order), labels = rev(format_method(method_names[methods_order])))
    )

n_methods <- length(unique(metric_ranks$Method))

category_ranks_plot <- ggplot(
    metric_ranks,
    aes(x = .data$Dataset, y = .data$Method, colour = .data$Type)
) +
    geom_linerange(
        data = baseline_methods |> filter(.data$IsBaseline),
        aes(x = NULL, y = .data$Method),
        xmin = -Inf, xmax = Inf,
        linewidth = 2.4, colour = "grey55", alpha = 0.28
    ) +
    geom_point(
        aes(alpha = n_methods - .data$Rank),
        shape = 15,
        size = 1.18
    ) +
    scale_y_discrete(drop = FALSE) +
    scale_colour_manual(values = types_palette, guide = "none") +
    scale_alpha_continuous(
        limits = c(n_methods, 1),
        range = c(0.10, 1.0),
        trans = "reverse",
        breaks = c(5, 10, 15, 20, 25, 30),
        name = "Mean category rank"
    ) +
    guides(
        alpha = guide_legend(theme = theme(legend.text.position = "bottom"), nrow = 1)
    ) +
    facet_grid(. ~ .data$Type) +
    labs(title = "Metric category ranks") +
    theme_features_pub() +
    theme(
        legend.position = "bottom",
        legend.title.position = "top",
        axis.text.x = element_text(size = 4.1, angle = 90, hjust = 1, vjust = 0.5, face = "bold", colour = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 6.6, margin = margin(b = 2)),
        panel.spacing = unit(0.05, "cm")
    )

figure <- wrap_plots(
    category_scores_plot,
    category_ranks_plot,
    nrow = 1,
    widths = c(1, 1),
    guides = "collect"
) &
    theme(
        legend.position = "bottom",
        legend.title.position = "top"
    )

write_tsv(plotting |> mutate(Method = as.character(.data$Method), Dataset = as.character(.data$Dataset)), file.path(output_dir, "extended_fig3_scores.tsv"))
write_tsv(metric_ranks |> mutate(Method = as.character(.data$Method), Dataset = as.character(.data$Dataset)), file.path(output_dir, "extended_fig3_ranks.tsv"))

save_figure_files(
    figure,
    file.path(output_dir, "extended_data_fig3"),
    width = 8.3,
    height = 6.3
)
