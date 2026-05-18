#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(purrr)
})

source(file.path("external", "atlas-feature-selection-benchmark", "analysis", "R", "plotting.R"))
source(file.path("analysis", "atlas_style", "R", "spatial_summarisation.R"))

args <- commandArgs(trailingOnly = TRUE)
benchmark_path <- if (length(args) >= 1) args[[1]] else file.path("results", "current_rank", "data", "benchmark.tsv")
epi_path <- if (length(args) >= 2) args[[2]] else file.path("results", "lineage_subsets", "stomics_0212_epithelial_subset_scvi", "results.csv")
imm_path <- if (length(args) >= 3) args[[3]] else file.path("results", "lineage_subsets", "stomics_0212_immune_subset_scvi", "results.csv")
output_dir <- if (length(args) >= 4) args[[4]] else file.path("results", "fig5a_spatial_lineages", "figures")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

types_palette <- c(
    "Overall" = "#f781bf",
    "Integration" = "#e41a1c",
    "Clustering" = "#377eb8",
    "Alignment" = "#4daf4a"
)

task_to_type <- c(
    "integration_eval" = "Integration",
    "clustering_eval" = "Clustering",
    "alignment_eval" = "Alignment"
)

pretty_dataset_names <- c(
    "STOmics0212" = "Full",
    "STOmics0212Immune" = "Immune",
    "STOmics0212Epithelial" = "Epithelial"
)

benchmark <- read_tsv(benchmark_path, show_col_types = FALSE)
methods_meta_all <- read_tsv(file.path(dirname(benchmark_path), "methods-metadata.tsv"), show_col_types = FALSE)

metric_meta <- benchmark |>
    distinct(Metric, MetricName, HigherBetter)

full_metrics <- benchmark |>
    filter(
        .data$Dataset == "STOmics0212",
        .data$IntegrationLabel == "scVI",
        .data$Type %in% c("Integration", "Clustering", "Alignment")
    ) |>
    select(
        Dataset, Platform, NSlices, Method, MethodBase, SelFeatures,
        Integration, IntegrationLabel, Task, TaskLabel, Type,
        Metric, MetricName, Value, ValueRaw, HigherBetter,
        Seed, Runtime, EffectiveFeatures, Notes
    )

subset_results_to_metrics <- function(path, dataset_label) {
    raw <- read_csv(path, show_col_types = FALSE)

    raw |>
        mutate(
            Dataset = dataset_label,
            Platform = .data$platform,
            NSlices = .data$n_slices,
            Method = if_else(
                .data$fs_method == "all_features",
                "all_features-Nall",
                paste0(.data$fs_method, "-N", .data$n_features)
            ),
            MethodBase = .data$fs_method,
            SelFeatures = as.character(.data$n_features),
            Integration = .data$integration_method,
            IntegrationLabel = case_when(
                .data$integration_method == "scvi" ~ "scVI",
                .default = .data$integration_method
            ),
            Task = .data$task,
            TaskLabel = recode(
                .data$task,
                "integration_eval" = "Integration",
                "clustering_eval" = "Clustering",
                "alignment_eval" = "Alignment"
            ),
            Type = unname(task_to_type[.data$task]),
            Metric = .data$metric_name,
            Seed = .data$random_seed,
            Runtime = .data$runtime,
            EffectiveFeatures = .data$effective_n_features,
            Notes = .data$notes
        ) |>
        left_join(metric_meta, by = "Metric") |>
        mutate(
            ValueRaw = .data$metric_value,
            Value = if_else(.data$HigherBetter, .data$ValueRaw, -.data$ValueRaw)
        ) |>
        select(
            Dataset, Platform, NSlices, Method, MethodBase, SelFeatures,
            Integration, IntegrationLabel, Task, TaskLabel, Type,
            Metric, MetricName, Value, ValueRaw, HigherBetter,
            Seed, Runtime, EffectiveFeatures, Notes
        ) |>
        filter(.data$Type %in% c("Integration", "Clustering", "Alignment"))
}

epi_metrics <- subset_results_to_metrics(epi_path, "STOmics0212Epithelial")
imm_metrics <- subset_results_to_metrics(imm_path, "STOmics0212Immune")

all_metrics <- bind_rows(full_metrics, imm_metrics, epi_metrics)

pick_canonical_methods <- function(metrics) {
    available <- metrics |>
        distinct(.data$Dataset, .data$MethodBase, .data$Method)

    common_families <- available |>
        group_by(.data$MethodBase) |>
        summarise(NDatasets = n_distinct(.data$Dataset), .groups = "drop") |>
        filter(.data$NDatasets == 3) |>
        pull(.data$MethodBase)

    family_pick <- function(family, methods) {
        methods <- sort(unique(methods))
        if (family == "random") {
            target <- paste0(family, "-N500")
            return(if (target %in% methods) target else methods[[1]])
        }
        preferred <- c(
            paste0(family, "-N2000"),
            paste0(family, "-N1000"),
            paste0(family, "-N500")
        )
        hit <- preferred[preferred %in% methods]
        if (length(hit) > 0) hit[[1]] else methods[[1]]
    }

    available_common <- available |>
        filter(.data$MethodBase %in% common_families)

    family_dfs <- split(available_common, available_common$MethodBase)
    rows <- lapply(names(family_dfs), function(family) {
        if (family == "Brennecke") {
            return(NULL)
        }
        df <- family_dfs[[family]]
        tibble(
            MethodBase = family,
            Method = family_pick(family, df$Method)
        )
    })

    bind_rows(rows)
}

canonical_methods <- pick_canonical_methods(all_metrics)

baseline_methods <- c(
    "random-N500",
    "scanpy_cell_ranger_batch-N2000",
    "scsegindex-N2000"
)

metrics_filtered <- all_metrics |>
    filter(.data$Method %in% canonical_methods$Method)

baseline_ranges <- metrics_filtered |>
    filter(.data$Method %in% baseline_methods) |>
    group_by(.data$Dataset, .data$Metric, .data$Type) |>
    summarise(
        Lower = min(.data$Value, na.rm = TRUE),
        Upper = max(.data$Value, na.rm = TRUE),
        Range = Upper - Lower,
        .groups = "drop"
    ) |>
    mutate(
        Range = if_else(.data$Range <= 0 | is.na(.data$Range), 1, .data$Range)
    )

metrics_summary <- summarise_spatial_metrics(
    metrics_filtered,
    baseline_ranges,
    type_weights = c("Integration" = 1/3, "Clustering" = 1/3, "Alignment" = 1/3),
    require_types_for_overall = c("Integration", "Clustering", "Alignment")
)

methods_meta <- canonical_methods |>
    left_join(methods_meta_all |> select(.data$Method, .data$Name), by = "Method") |>
    mutate(Name = coalesce(.data$Name, .data$Method))

metrics_ranks <- metrics_summary |>
    pivot_longer(
        cols = c("Integration", "Clustering", "Alignment", "Overall"),
        names_to = "Type",
        values_to = "Value"
    ) |>
    filter(!is.na(.data$Value)) |>
    group_by(.data$Dataset, .data$Type) |>
    mutate(Rank = rank(-.data$Value, ties.method = "average")) |>
    ungroup()

method_order <- metrics_ranks |>
    filter(.data$Type == "Overall") |>
    group_by(.data$Method) |>
    summarise(MeanRank = mean(.data$Rank, na.rm = TRUE), .groups = "drop") |>
    arrange(.data$MeanRank) |>
    pull(.data$Method)

plot_df <- metrics_ranks |>
    mutate(
        Dataset = factor(.data$Dataset, levels = c("STOmics0212", "STOmics0212Immune", "STOmics0212Epithelial"),
                         labels = c("Full", "Immune", "Epithelial")),
        MethodLabel = methods_meta$Name[match(.data$Method, methods_meta$Method)],
        MethodLabel = factor(.data$MethodLabel, levels = rev(methods_meta$Name[match(method_order, methods_meta$Method)])),
        Type = factor(.data$Type, levels = c("Overall", "Integration", "Clustering", "Alignment"))
    )

write_tsv(
    plot_df |>
        select(.data$Dataset, .data$Method, .data$MethodLabel, .data$Type, .data$Rank, .data$Value),
    file.path(output_dir, "fig5a_lineage_ranks.tsv")
)
write_tsv(methods_meta, file.path(output_dir, "fig5a_lineage_methods.tsv"))

fig <- ggplot(plot_df, aes(x = .data$Dataset, y = .data$MethodLabel, fill = .data$Rank)) +
    geom_tile(aes(fill = .data$Type, alpha = .data$Rank), colour = "white", linewidth = 0.35) +
    facet_wrap(~ .data$Type, nrow = 1) +
    scale_fill_manual(values = types_palette, guide = "none") +
    scale_alpha_continuous(
        limits = c(max(plot_df$Rank, na.rm = TRUE), 1),
        range = c(0.22, 1.0),
        trans = "reverse",
        breaks = c(5, 10, 15, 20, 25),
        name = "Rank"
    ) +
    labs(x = NULL, y = NULL) +
    theme_features_pub() +
    theme(
        strip.background = element_rect(fill = "black", colour = "black"),
        strip.text = element_text(colour = "white", face = "bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "bold", colour = "black"),
        axis.text.y = element_text(face = "bold", colour = "black", size = 6),
        panel.grid = element_blank(),
        legend.position = "bottom",
        plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm")
    )

save_figure_files(fig, file.path(output_dir, "figure_fig5a_spatial_lineages"), width = 7.8, height = 7.8)
