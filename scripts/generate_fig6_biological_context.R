#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(tidyr)
    library(stringr)
    library(ggplot2)
    library(patchwork)
    library(jsonlite)
    library(colorspace)
    library(grid)
})

source(file.path("external", "atlas-feature-selection-benchmark", "analysis", "R", "plotting.R"))
source(file.path("analysis", "manuscript_rebuild_v2", "R", "data.R"))

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) >= 1) args[[1]] else file.path("results", "fig6_biological_context")
source_dir <- file.path(output_dir, "source_data")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

full_scores_path <- file.path(
    "results", "spatial_svg_rebuild_v1", "frozen_scores",
    "dataset_representative_ranks.tsv"
)
subset_scores_path <- file.path(
    "results", "biological_subsets_rebuild_v1", "frozen_scores",
    "dataset_representative_ranks.tsv"
)
representatives_path <- file.path(
    "results", "biological_subsets_rebuild_v1", "frozen_scores",
    "representative_settings.tsv"
)
overlap_path <- file.path(
    "results", "biological_subsets_rebuild_v1", "feature_overlap_summary.tsv"
)
marker_overlap_path <- file.path(
    "results", "biological_subsets_rebuild_v1", "marker_overlap_summary.tsv"
)
marker_genes_path <- file.path(
    "results", "biological_subsets_rebuild_v1", "marker_gene_sets_top100.tsv"
)

context_labels <- c(
    STOmics0212 = "Full",
    STOmics0212Immune = "Immune",
    STOmics0212Epithelial = "Epithelial"
)
context_levels <- unname(context_labels)

representatives <- read_tsv(representatives_path, show_col_types = FALSE) |>
    filter(.data$dataset == "STOmics0212Immune") |>
    distinct(.data$fs_method, .data$n_features, .data$MethodGroup)

format_n <- function(value) {
    if (as.character(value) == "all") return("all")
    format(as.integer(value), big.mark = ",", scientific = FALSE, trim = TRUE)
}

method_display <- setNames(
    vapply(seq_len(nrow(representatives)), function(i) {
        method <- representatives$fs_method[[i]]
        base <- display_method(method)
        if (method == "all_features") return(base)
        paste0(base, " (N=", format_n(representatives$n_features[[i]]), ")")
    }, character(1)),
    representatives$fs_method
)

full_ranks <- read_tsv(full_scores_path, show_col_types = FALSE) |>
    filter(.data$dataset == "STOmics0212", .data$integration_method == "scvi")
subset_ranks <- read_tsv(subset_scores_path, show_col_types = FALSE) |>
    filter(.data$integration_method == "scvi")

rank_wide <- bind_rows(full_ranks, subset_ranks) |>
    filter(.data$dataset %in% names(context_labels)) |>
    mutate(Context = unname(context_labels[.data$dataset]))

method_order <- rank_wide |>
    group_by(.data$fs_method) |>
    summarise(MeanOverallRank = mean(.data$CoreOverallRank), .groups = "drop") |>
    arrange(.data$MeanOverallRank, .data$fs_method) |>
    pull(.data$fs_method)

method_factor <- function(methods) {
    factor(
        methods,
        levels = rev(method_order),
        labels = rev(unname(method_display[method_order]))
    )
}

rank_df <- rank_wide |>
    select(
        .data$dataset, .data$Context, .data$fs_method, .data$MethodGroup,
        Overall = .data$CoreOverallRank,
        Integration = .data$IntegrationRank,
        Clustering = .data$ClusteringRank
    ) |>
    pivot_longer(
        cols = c("Overall", "Integration", "Clustering"),
        names_to = "ScoreType",
        values_to = "Rank"
    ) |>
    mutate(
        Context = factor(.data$Context, levels = context_levels),
        ScoreType = factor(.data$ScoreType, levels = c("Overall", "Integration", "Clustering")),
        Method = method_factor(.data$fs_method)
    )

overlap_df <- read_tsv(overlap_path, show_col_types = FALSE) |>
    group_by(.data$method, .data$group, .data$pair) |>
    summarise(
        Jaccard = mean(.data$jaccard),
        JaccardSD = sd(.data$jaccard),
        NSeeds = n(),
        .groups = "drop"
    ) |>
    mutate(
        Jaccard = if_else(.data$method == "random", NA_real_, .data$Jaccard),
        Combination = recode(
            .data$pair,
            "Full vs Immune" = "Full vs Imm.",
            "Full vs Epithelial" = "Full vs Epi.",
            "Immune vs Epithelial" = "Epi. vs Imm."
        ),
        Combination = factor(
            .data$Combination,
            levels = c("Full vs Imm.", "Full vs Epi.", "Epi. vs Imm.")
        ),
        Method = method_factor(.data$method)
    )

marker_lineage_df <- read_tsv(marker_overlap_path, show_col_types = FALSE) |>
    group_by(.data$method, .data$group, .data$context, .data$lineage) |>
    summarise(
        MeanProp = mean(.data$mean_marker_recovery),
        SDProp = mean(.data$sd_celltypes),
        NSeeds = n(),
        .groups = "drop"
    ) |>
    mutate(
        MeanProp = if_else(.data$method %in% c("all_features", "random"), NA_real_, .data$MeanProp),
        SDProp = if_else(.data$method %in% c("all_features", "random"), NA_real_, .data$SDProp),
        Context = factor(.data$context, levels = context_levels),
        Lineage = factor(.data$lineage, levels = c("Epithelial", "Immune")),
        Method = method_factor(.data$method)
    )

marker_genes <- read_tsv(marker_genes_path, show_col_types = FALSE) |>
    group_by(.data$label) |>
    summarise(Markers = list(unique(.data$gene)), .groups = "drop")
marker_lookup <- setNames(marker_genes$Markers, marker_genes$label)

feature_roots <- c(
    Full = file.path(
        "results", "spatial_svg_rebuild_v1", "stomics_0212", "stomics0212",
        "feature_selection"
    ),
    Immune = file.path(
        "results", "biological_subsets_rebuild_v1", "stomics_0212_immune_subset",
        "stomics0212immune", "feature_selection"
    ),
    Epithelial = file.path(
        "results", "biological_subsets_rebuild_v1", "stomics_0212_epithelial_subset",
        "stomics0212epithelial", "feature_selection"
    )
)
context_celltypes <- list(
    Full = c(
        "Langerhans Cell", "T Cell", "Keratinocyte", "Merkel Cell",
        "Stem Cell", "Progenitor Cell"
    ),
    Immune = c("Langerhans Cell", "T Cell"),
    Epithelial = c("Keratinocyte", "Merkel Cell", "Stem Cell", "Progenitor Cell")
)

read_selected_features <- function(path) {
    unique(jsonlite::fromJSON(path)$feature_names)
}

celltype_rows <- list()
row_index <- 1L
for (i in seq_len(nrow(representatives))) {
    method <- representatives$fs_method[[i]]
    n_value <- if (as.character(representatives$n_features[[i]]) == "all") {
        "1"
    } else {
        as.character(as.integer(representatives$n_features[[i]]))
    }
    for (context in names(feature_roots)) {
        files <- Sys.glob(file.path(
            feature_roots[[context]], method, paste0("n", n_value),
            "seed*", "selected_features.json"
        ))
        if (length(files) == 0) {
            stop("No selected feature files for ", context, "/", method)
        }
        for (path in files) {
            selected <- read_selected_features(path)
            seed <- as.integer(str_match(path, "seed([0-9]+)")[, 2])
            for (label in context_celltypes[[context]]) {
                markers <- marker_lookup[[label]]
                celltype_rows[[row_index]] <- tibble(
                    fs_method = method,
                    Context = context,
                    Label = label,
                    Seed = seed,
                    MarkerOverlap = length(intersect(selected, markers)) / length(markers)
                )
                row_index <- row_index + 1L
            }
        }
    }
}

celltype_df <- bind_rows(celltype_rows) |>
    group_by(.data$fs_method, .data$Context, .data$Label) |>
    summarise(
        MarkerOverlap = mean(.data$MarkerOverlap),
        MarkerOverlapSD = sd(.data$MarkerOverlap),
        NSeeds = n(),
        .groups = "drop"
    )

full_celltype <- celltype_df |>
    filter(.data$Context == "Full") |>
    select(.data$fs_method, .data$Label, FullMarkerOverlap = .data$MarkerOverlap)

celltype_diff_df <- celltype_df |>
    filter(.data$Context %in% c("Immune", "Epithelial")) |>
    left_join(full_celltype, by = c("fs_method", "Label")) |>
    mutate(Difference = .data$MarkerOverlap - .data$FullMarkerOverlap)

label_levels <- c(
    "Langerhans Cell", "T Cell", "Keratinocyte", "Merkel Cell",
    "Stem Cell", "Progenitor Cell"
)

celltype_plot_df <- celltype_df |>
    mutate(
        MarkerOverlap = if_else(
            .data$fs_method %in% c("all_features", "random"),
            NA_real_, .data$MarkerOverlap
        ),
        Method = method_factor(.data$fs_method),
        Context = factor(.data$Context, levels = context_levels),
        Label = factor(.data$Label, levels = label_levels)
    )

celltype_diff_plot_df <- celltype_diff_df |>
    mutate(
        Difference = if_else(
            .data$fs_method %in% c("all_features", "random"),
            NA_real_, .data$Difference
        ),
        Method = method_factor(.data$fs_method),
        Context = factor(.data$Context, levels = c("Immune", "Epithelial")),
        Label = factor(.data$Label, levels = label_levels)
    )

write_tsv(rank_df, file.path(source_dir, "fig6a_context_ranks.tsv"), na = "NA")
write_tsv(overlap_df, file.path(source_dir, "fig6b_feature_overlap.tsv"), na = "NA")
write_tsv(marker_lineage_df, file.path(source_dir, "fig6c_lineage_marker_overlap.tsv"), na = "NA")
write_tsv(celltype_df, file.path(source_dir, "fig6d_celltype_marker_overlap.tsv"), na = "NA")
write_tsv(celltype_diff_df, file.path(source_dir, "fig6d_celltype_marker_difference.tsv"), na = "NA")
mean_rank_lookup <- rank_wide |>
    group_by(.data$fs_method) |>
    summarise(Value = mean(.data$CoreOverallRank), .groups = "drop")
write_tsv(
    tibble(
        Method = method_order,
        DisplayName = unname(method_display[method_order]),
        MeanOverallRank = mean_rank_lookup$Value[
            match(method_order, mean_rank_lookup$fs_method)
        ]
    ),
    file.path(source_dir, "fig6_method_order.tsv")
)

task_palette <- c(
    Overall = "#F781BF",
    Integration = "#E64B4B",
    Clustering = "#4E91C4"
)

panel_theme <- theme_features_pub() +
    theme(
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(
            angle = 90, hjust = 1, vjust = 0.5,
            colour = "black", size = 5.0
        ),
        axis.text.y = element_text(colour = "black", size = 4.7),
        strip.background = element_rect(fill = "black", colour = "black", linewidth = 0.3),
        strip.text = element_text(colour = "white", size = 5.6, face = "plain"),
        panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.3),
        panel.spacing = unit(0.045, "cm"),
        plot.margin = margin(0.02, 0.04, 0.04, 0.04, "cm"),
        plot.tag = element_text(size = 8, face = "bold", colour = "black"),
        plot.tag.position = "topleft"
    )

max_rank <- max(rank_df$Rank, na.rm = TRUE)
rank_breaks <- c(5, 10, 15, 20, 25, 30)
rank_breaks <- rank_breaks[rank_breaks <= max_rank]

plot_a <- ggplot(rank_df, aes(x = .data$Context, y = .data$Method)) +
    geom_tile(aes(fill = .data$ScoreType, alpha = .data$Rank)) +
    facet_wrap(~ .data$ScoreType, nrow = 1) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0), drop = FALSE) +
    scale_fill_manual(values = task_palette, guide = "none") +
    scale_alpha_continuous(
        limits = c(max_rank, 1),
        breaks = rank_breaks,
        range = c(0.14, 1),
        trans = "reverse",
        name = "Rank",
        guide = guide_legend(
            order = 1,
            nrow = 1,
            byrow = TRUE,
            title.position = "top",
            override.aes = list(fill = "#555555", colour = NA)
        )
    ) +
    labs(tag = "a") +
    panel_theme +
    theme(axis.text.y = element_text(colour = "black", size = 4.7))

plot_b <- ggplot(overlap_df, aes(x = .data$Combination, y = .data$Method, fill = .data$Jaccard)) +
    geom_tile() +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0), drop = FALSE) +
    scale_fill_viridis_c(
        option = "C", limits = c(0, 1), breaks = c(0, 0.5, 1), na.value = "white",
        name = "Jaccard index",
        guide = guide_colourbar(order = 2, title.position = "top", barwidth = unit(0.95, "cm"))
    ) +
    labs(tag = "b") +
    panel_theme +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

sd_limit <- max(marker_lineage_df$SDProp, na.rm = TRUE)
plot_c <- ggplot(marker_lineage_df, aes(x = .data$Context, y = .data$Method)) +
    geom_point(
        aes(colour = .data$MeanProp, size = .data$SDProp),
        shape = 15, stroke = 0, na.rm = TRUE
    ) +
    facet_wrap(~ .data$Lineage, nrow = 1) +
    scale_x_discrete(expand = expansion(mult = c(0.06, 0.06))) +
    scale_y_discrete(expand = expansion(add = 0.55), drop = FALSE) +
    scale_colour_viridis_c(
        option = "D", limits = c(0, 1), breaks = c(0, 0.5, 1),
        name = "Mean proportion",
        guide = guide_colourbar(order = 3, title.position = "top", barwidth = unit(1.15, "cm"))
    ) +
    scale_size_continuous(
        trans = "reverse",
        limits = c(sd_limit, 0),
        range = c(0.55, 2.20),
        breaks = scales::breaks_extended(n = 3),
        name = "s.d. proportion",
        guide = guide_legend(order = 4, nrow = 1, byrow = TRUE, title.position = "top")
    ) +
    labs(tag = "c") +
    panel_theme +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

top_panel <- (plot_a | plot_b | plot_c) +
    plot_layout(widths = c(1.08, 0.34, 0.62), guides = "collect") &
    theme(
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.box.just = "left",
        legend.title.position = "top",
        legend.margin = margin(0, 0, 0, 0),
        legend.spacing.x = unit(0.28, "cm"),
        legend.key.height = unit(0.22, "cm"),
        legend.key.width = unit(0.36, "cm"),
        legend.title = element_text(size = 5.8),
        legend.text = element_text(size = 5.0)
    )

plot_d_left <- ggplot(
    celltype_plot_df,
    aes(x = .data$Label, y = .data$Method, fill = .data$MarkerOverlap)
) +
    geom_tile() +
    facet_grid(. ~ .data$Context, scales = "free_x", space = "free_x") +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0), drop = FALSE) +
    scale_fill_viridis_c(
        option = "magma", direction = -1,
        limits = c(0, 1), breaks = c(0, 0.5, 1), na.value = "white",
        name = "Marker overlap",
        guide = guide_colourbar(order = 1, title.position = "top", barwidth = unit(1.25, "cm"))
    ) +
    labs(tag = "d") +
    panel_theme +
    theme(
        axis.text.x = element_text(
            angle = 90, hjust = 1, vjust = 0.5,
            lineheight = 0.82, colour = "black", size = 4.6
        ),
        axis.text.y = element_text(colour = "black", size = 4.7)
    )

diff_limit <- max(abs(celltype_diff_plot_df$Difference), na.rm = TRUE)
diff_limit <- min(1, max(0.3, ceiling(diff_limit * 10) / 10))
plot_d_right <- ggplot(
    celltype_diff_plot_df,
    aes(x = .data$Label, y = .data$Method, fill = .data$Difference)
) +
    geom_tile() +
    facet_grid(. ~ .data$Context, scales = "free_x", space = "free_x") +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0), drop = FALSE) +
    colorspace::scale_fill_continuous_diverging(
        palette = "Purple-Green",
        limits = c(-diff_limit, diff_limit),
        breaks = c(-diff_limit, 0, diff_limit),
        na.value = "white",
        name = "Difference from Full",
        guide = guide_colourbar(order = 2, title.position = "top", barwidth = unit(1.25, "cm"))
    ) +
    labs(tag = " ") +
    panel_theme +
    theme(
        axis.text.x = element_text(
            angle = 90, hjust = 1, vjust = 0.5,
            lineheight = 0.82, colour = "black", size = 4.6
        ),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
    )

bottom_panel <- (plot_d_left | plot_d_right) +
    plot_layout(widths = c(1, 0.48), guides = "collect") &
    theme(
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.box.just = "left",
        legend.title.position = "top",
        legend.margin = margin(0, 0, 0, 0),
        legend.spacing.x = unit(0.20, "cm"),
        legend.key.height = unit(0.22, "cm"),
        legend.key.width = unit(0.40, "cm"),
        legend.title = element_text(size = 5.8),
        legend.text = element_text(size = 5.0)
    )

figure <- top_panel / bottom_panel +
    plot_layout(heights = c(1.04, 0.96)) &
    theme(
        text = element_text(family = "Arial", face = "plain", colour = "black"),
        plot.background = element_rect(fill = "white", colour = NA)
    )

base_path <- file.path(output_dir, "figure_fig6_biological_context")
width <- 8.2
height <- 7.75

ragg::agg_png(
    paste0(base_path, ".png"),
    width = width, height = height, units = "in", res = 450,
    background = "white"
)
print(figure)
dev.off()

ragg::agg_tiff(
    paste0(base_path, ".tiff"),
    width = width, height = height, units = "in", res = 600,
    background = "white", compression = "lzw"
)
print(figure)
dev.off()

register_arial_pdf_font()
grDevices::pdf(
    paste0(base_path, ".pdf"),
    width = width, height = height,
    family = "Arial", useDingbats = FALSE,
    bg = "white"
)
print(figure)
dev.off()

svglite::svglite(
    paste0(base_path, ".svg"),
    width = width, height = height,
    bg = "white"
)
print(figure)
dev.off()

svg_path <- paste0(base_path, ".svg")
svg_text <- readLines(svg_path, warn = FALSE)
svg_text <- gsub(
    " textLength='[^']+' lengthAdjust='spacingAndGlyphs'",
    "",
    svg_text
)
svg_text <- gsub(
    'font-family: "Liberation Sans";',
    "font-family: Arial, Helvetica, sans-serif;",
    svg_text,
    fixed = TRUE
)
writeLines(svg_text, svg_path, useBytes = TRUE)

saveRDS(figure, paste0(base_path, ".rds"))

qa <- tibble(
    Check = c(
        "methods in all panels",
        "rank cells",
        "feature-overlap cells",
        "lineage-marker observations",
        "cell-type marker observations",
        "cell-type difference observations",
        "rank values outside expected range",
        "Jaccard values outside [0,1]",
        "marker values outside [0,1]"
    ),
    Observed = c(
        length(method_order),
        nrow(rank_df),
        nrow(overlap_df),
        sum(!is.na(marker_lineage_df$MeanProp)),
        sum(!is.na(celltype_plot_df$MarkerOverlap)),
        sum(!is.na(celltype_diff_plot_df$Difference)),
        sum(rank_df$Rank < 1 | rank_df$Rank > length(method_order), na.rm = TRUE),
        sum(overlap_df$Jaccard < 0 | overlap_df$Jaccard > 1, na.rm = TRUE),
        sum(celltype_df$MarkerOverlap < 0 | celltype_df$MarkerOverlap > 1, na.rm = TRUE)
    ),
    Expected = c(34, 34 * 3 * 3, 34 * 3, NA, NA, NA, 0, 0, 0)
)
write_tsv(qa, file.path(output_dir, "qa_data_integrity.tsv"), na = "NA")

writeLines(
    c(
        "Core conclusion: lineage restriction reshapes selected features and method ranks without guaranteeing improved marker recovery.",
        "Archetype: quantitative grid, aligned to the reference lineage-subset figure.",
        "Panel a: method ranks for Full, Immune, and Epithelial contexts.",
        "Panel b: same-method feature-set Jaccard overlap across contexts.",
        "Panel c: lineage-level marker recovery; square size encodes variation across cell types.",
        "Panel d: cell-type marker recovery and subset-minus-Full differences.",
        "Scaling note: rank values are context-specific; absolute scaled scores are not compared across contexts.",
        "Caveat: Immune contains 1,533 Langerhans cells and 23 T cells; it is a stress-test context."
    ),
    file.path(output_dir, "figure_contract.txt")
)
