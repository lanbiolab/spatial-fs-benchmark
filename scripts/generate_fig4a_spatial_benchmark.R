#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(patchwork)
    library(ggforce)
})

args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[[1]] else file.path("results", "current_rank", "data")
output_dir <- if (length(args) >= 2) args[[2]] else file.path("results", "fig4a_spatial_benchmark", "figures")
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

methods_meta_all <- read_tsv(file.path(data_dir, "methods-metadata.tsv"), show_col_types = FALSE)

metric_ranges <- read_tsv(ranges_path, show_col_types = FALSE) |>
    filter(.data$Type %in% c("Integration", "Clustering", "Alignment"))

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
        if (target %in% methods) {
            return(target)
        }
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

    selected <- available |>
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

    selected
}

methods_meta <- pick_representative_methods(metrics_summary_all) |>
    left_join(methods_meta_all |> select(.data$Method, .data$Name), by = "Method") |>
    mutate(
        Name = case_when(
            !is.na(.data$Name) ~ .data$Name,
            .data$Method == "all_features-Nall" ~ "All",
            TRUE ~ .data$Method
        )
    )

metrics_summary <- metrics_summary_all |>
    filter(.data$Method %in% methods_meta$Method, .data$IntegrationLabel == "scVI")

method_names <- methods_meta$Name
names(method_names) <- methods_meta$Method

metrics_summary_plotting <- metrics_summary |>
    pivot_longer(
        cols = c("Overall", "Integration", "Clustering", "Alignment"),
        names_to = "Type",
        values_to = "Value"
    ) |>
    filter(!is.na(.data$Value))

metric_means <- metrics_summary_plotting |>
    group_by(.data$Method, .data$Type) |>
    summarise(Value = mean(.data$Value), .groups = "drop")

methods_order <- metric_means |>
    filter(.data$Type == "Overall") |>
    arrange(desc(.data$Value)) |>
    pull(.data$Method)

metrics_summary_plotting <- metrics_summary_plotting |>
    mutate(
        Type = factor(.data$Type, levels = c("Overall", "Integration", "Clustering", "Alignment")),
        Method = factor(.data$Method, levels = rev(methods_order), labels = rev(method_names[methods_order]))
    )

metric_means <- metric_means |>
    mutate(
        Type = factor(.data$Type, levels = c("Overall", "Integration", "Clustering", "Alignment")),
        Method = factor(.data$Method, levels = rev(methods_order), labels = rev(method_names[methods_order]))
    )

baseline_methods <- methods_meta |>
    filter(.data$Method %in% methods_order) |>
    mutate(
        Method = factor(.data$Method, levels = rev(methods_order), labels = rev(method_names[methods_order]))
    )

overall_values_figure <- ggplot(
    metrics_summary_plotting,
    aes(x = .data$Value, y = .data$Method, fill = .data$Type)
) +
    annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.08) +
    annotate("rect", xmin = 1, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "blue", alpha = 0.08) +
    geom_vline(xintercept = 0, colour = "red", linewidth = 0.3) +
    geom_vline(xintercept = 1, colour = "blue", linewidth = 0.3) +
    geom_linerange(
        data = baseline_methods,
        aes(x = NULL, y = .data$Method, fill = NULL, alpha = .data$IsBaseline),
        xmin = -Inf, xmax = Inf,
        linewidth = 2.5, colour = "grey55"
    ) +
    ggforce::geom_sina(
        size = 1.0, alpha = 0.55, shape = 21, stroke = 0, colour = "white", maxwidth = 0.55
    ) +
    geom_point(
        data = metric_means,
        aes(fill = .data$Type),
        shape = 23, size = 2.15, colour = "white", stroke = 0.38
    ) +
    scale_fill_manual(values = types_palette) +
    scale_alpha_manual(values = c(`FALSE` = 0, `TRUE` = 0.25), guide = "none") +
    facet_grid(. ~ .data$Type, scales = "free_x") +
    labs(x = NULL, y = NULL) +
    theme_features_pub() +
    theme(
        legend.position = "bottom",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "black", colour = "black"),
        strip.text = element_text(colour = "white", face = "bold"),
        axis.text.y = element_text(face = "bold", colour = "black", size = 5.1),
        plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm")
    )

metric_ranks_means <- metrics_summary_plotting |>
    mutate(MethodRaw = as.character(.data$Method)) |>
    group_by(.data$Dataset, .data$Type) |>
    mutate(Rank = rank(-.data$Value, ties.method = "average")) |>
    ungroup() |>
    group_by(.data$MethodRaw, .data$Type) |>
    summarise(
        MeanRank = mean(.data$Rank),
        SDRank = sd(.data$Rank),
        .groups = "drop"
    ) |>
    mutate(
        Type = factor(.data$Type, levels = c("Overall", "Integration", "Clustering", "Alignment")),
        Method = factor(.data$MethodRaw, levels = rev(method_names[methods_order]))
    )

n_methods <- length(methods_order)

overall_ranks_figure <- ggplot(
    metric_ranks_means,
    aes(x = .data$Type, y = .data$Method, colour = .data$Type)
) +
    geom_linerange(
        data = baseline_methods |> filter(.data$IsBaseline),
        aes(x = NULL, y = .data$Method),
        xmin = -Inf, xmax = Inf,
        linewidth = 2.5, colour = "grey55", alpha = 0.25
    ) +
    geom_point(
        aes(alpha = .data$MeanRank, size = .data$SDRank),
        shape = 15
    ) +
    scale_y_discrete(drop = FALSE) +
    scale_colour_manual(values = types_palette, guide = "none") +
    scale_size_continuous(
        trans = "reverse",
        limits = c(max(metric_ranks_means$SDRank, na.rm = TRUE), 0),
        range = c(0.7, 3.0)
    ) +
    scale_alpha_continuous(
        limits = c(n_methods, 1),
        range = c(0.18, 1.0),
        trans = "reverse",
        breaks = c(10, 20, n_methods),
        labels = c("10", "20", "30")
    ) +
    guides(
        alpha = guide_legend(
            override.aes = list(shape = 15, size = 3.8, colour = "black"),
            order = 2
        ),
        size = guide_legend(
            override.aes = list(shape = 15, colour = "black", alpha = 1),
            order = 3
        ),
        fill = guide_legend(order = 1)
    ) +
    labs(
        title = "Mean ranks",
        alpha = "Mean rank",
        size = "Rank s.d.",
        x = NULL,
        y = NULL
    ) +
    theme_features_pub() +
    theme(
        plot.title = element_text(hjust = 0.5, vjust = -0.2, margin = margin(t = 3, b = -1)),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = margin(0.2, 0.2, 0.2, 0.05, "cm")
    )

figure <- wrap_plots(
    overall_values_figure,
    overall_ranks_figure,
    nrow = 1,
    widths = c(5, 0.88),
    guides = "collect"
) &
    theme(
        legend.position = "bottom",
        legend.title.position = "top",
        legend.box = "horizontal",
        legend.direction = "horizontal",
        legend.justification = "left",
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.text = element_text(size = 5.2),
        legend.title = element_text(size = 5.8),
        legend.key.size = unit(0.34, "cm"),
        legend.spacing.x = unit(0.18, "cm"),
        legend.box.spacing = unit(0.24, "cm")
    )

readr::write_tsv(metric_means, file.path(output_dir, "fig4a_metric_means.tsv"))
readr::write_tsv(metric_ranks_means, file.path(output_dir, "fig4a_metric_ranks.tsv"))

save_figure_files(figure, file.path(output_dir, "figure_fig4a_spatial_benchmark"), width = 8.3, height = 5.95)
