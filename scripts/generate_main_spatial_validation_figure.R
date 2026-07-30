#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(patchwork)
    library(readr)
    library(tidyr)
})

source(file.path("external", "atlas-feature-selection-benchmark", "analysis", "R", "plotting.R"))

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) >= 1) {
    args[[1]]
} else {
    file.path("results", "manuscript_rebuild_v2", "main_figure2_validation")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

method_labels <- c(
    "all_features" = "All features",
    "random" = "Random",
    "scsegindex" = "scSEGIndex",
    "Brennecke" = "Brennecke",
    "anticor" = "Anticor",
    "dubstepr" = "DUBStepR",
    "hotspot" = "Hotspot",
    "nbumi" = "NBumi",
    "osca" = "OSCA",
    "scPNMF" = "scPNMF",
    "scanpy_cell_ranger" = "Scanpy Cell Ranger",
    "scanpy_cell_ranger_batch" = "Scanpy Cell Ranger, batch",
    "scanpy_pearson" = "Scanpy Pearson",
    "scanpy_pearson_batch" = "Scanpy Pearson, batch",
    "scanpy_seurat" = "Scanpy Seurat",
    "scanpy_seurat_batch" = "Scanpy Seurat, batch",
    "scanpy_seurat_v3" = "Scanpy Seurat v3",
    "scanpy_seurat_v3_batch" = "Scanpy Seurat v3, batch",
    "scry" = "scry",
    "seurat_disp" = "Seurat Dispersion",
    "seurat_mvp" = "Seurat MVP",
    "seurat_sct" = "Seurat sctransform",
    "seurat_vst" = "Seurat VST",
    "singleCellHaystack" = "singleCellHaystack",
    "statistic_mean" = "Mean expression",
    "statistic_variance" = "Expression variance",
    "triku" = "triku",
    "wilcoxon" = "Wilcoxon, label-informed",
    "morans_i" = "Moran's I",
    "nnsvg" = "nnSVG",
    "somde" = "SOMDE",
    "sparkx" = "SPARK-X",
    "spatialde" = "SpatialDE"
)

method_label <- function(x) {
    recode(x, !!!method_labels, .default = x)
}

base_theme <- theme_features_pub() +
    theme(
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.32),
        strip.background = element_rect(fill = "black", colour = "black", linewidth = 0.32),
        strip.text = element_text(colour = "white", face = "bold", size = 6.0),
        panel.spacing = unit(0.05, "cm"),
        plot.margin = margin(0.8, 1.2, 0.8, 0.8, "mm"),
        plot.tag = element_text(size = 9, face = "bold", colour = "black"),
        plot.tag.position = "topleft",
        axis.text = element_text(colour = "black")
    )

semi <- read_tsv(
    file.path("results", "validation_v1", "semi_synthetic_metrics.tsv"),
    show_col_types = FALSE
)
heldout <- read_tsv(
    file.path("results", "validation_v1", "heldout_slice_metrics.tsv"),
    show_col_types = FALSE
)
rank_frame <- read_tsv(
    file.path(
        "results", "submission_validations_v1", "figures",
        "Supplemental_Figure_S7c_rank_concordance.tsv"
    ),
    show_col_types = FALSE
)
concordance <- read_tsv(
    file.path("results", "submission_validations_v1", "heldout_cross_reference_concordance.tsv"),
    show_col_types = FALSE
)
association <- read_tsv(
    file.path(
        "results", "submission_validations_v1", "figures",
        "Supplemental_Figure_S7d_ari_association.tsv"
    ),
    show_col_types = FALSE
)

# a, Semi-synthetic truth recovery.
auprc_raw <- semi |>
    filter(.data$metric == "AUPRC")
semi_order <- auprc_raw |>
    group_by(.data$method) |>
    summarise(MeanAUPRC = mean(.data$value), .groups = "drop") |>
    arrange(desc(.data$MeanAUPRC), .data$method) |>
    pull(.data$method)
semi_levels <- rev(semi_order)

prevalence_palette <- c("5% SVGs" = "#377eb8", "15% SVGs" = "#e41a1c")
auprc_plot <- auprc_raw |>
    mutate(
        Prevalence = factor(
            .data$prevalence, levels = c(0.05, 0.15),
            labels = c("5% SVGs", "15% SVGs")
        ),
        Method = factor(
            method_label(.data$method),
            levels = method_label(semi_levels)
        )
    )
auprc_summary <- auprc_plot |>
    group_by(.data$method, .data$Method, .data$Prevalence) |>
    summarise(
        Mean = mean(.data$value),
        SD = sd(.data$value),
        Lower = pmax(0, .data$Mean - .data$SD),
        Upper = pmin(1, .data$Mean + .data$SD),
        .groups = "drop"
    )
control_bands_a <- auprc_plot |>
    filter(.data$method %in% c("all_features", "random", "scsegindex")) |>
    distinct(.data$Method) |>
    mutate(y = as.numeric(.data$Method))

panel_a <- ggplot(auprc_plot, aes(x = .data$value, y = .data$Method)) +
    geom_rect(
        data = control_bands_a,
        aes(xmin = -Inf, xmax = Inf, ymin = .data$y - 0.47, ymax = .data$y + 0.47),
        inherit.aes = FALSE, fill = "grey70", alpha = 0.28, colour = NA
    ) +
    geom_point(
        aes(colour = .data$Prevalence), size = 0.72, alpha = 0.38,
        position = position_jitter(height = 0.09, width = 0, seed = 17)
    ) +
    geom_linerange(
        data = auprc_summary,
        aes(xmin = .data$Lower, xmax = .data$Upper, y = .data$Method),
        inherit.aes = FALSE, linewidth = 0.38, colour = "grey25"
    ) +
    geom_point(
        data = auprc_summary,
        aes(x = .data$Mean, y = .data$Method, fill = .data$Prevalence),
        inherit.aes = FALSE,
        shape = 23, size = 1.75, stroke = 0.30, colour = "white"
    ) +
    facet_grid(. ~ .data$Prevalence) +
    scale_colour_manual(values = prevalence_palette, guide = "none") +
    scale_fill_manual(values = prevalence_palette, guide = "none") +
    scale_x_continuous(
        limits = c(0, 1), breaks = c(0, 0.5, 1),
        expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_y_discrete(drop = FALSE) +
    labs(x = "AUPRC", y = NULL, tag = "a") +
    base_theme +
    theme(
        axis.text.y = element_text(size = 4.65),
        axis.text.x = element_text(size = 4.8),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.spacing.x = unit(0.25, "cm")
    )

# b, Recovery by spatial pattern and effect size.
recall_plot <- semi |>
    filter(
        .data$metric == "Recall",
        .data$summary %in% c("pattern", "effect")
    ) |>
    group_by(.data$method, .data$summary, .data$stratum) |>
    summarise(MeanRecall = mean(.data$value), .groups = "drop") |>
    complete(
        method = semi_order,
        nesting(summary, stratum)
    ) |>
    mutate(
        Method = factor(method_label(.data$method), levels = method_label(semi_levels)),
        Group = factor(.data$summary, levels = c("pattern", "effect"), labels = c("Pattern", "Effect")),
        Stratum = recode(
            .data$stratum,
            "domain" = "Domain", "gradient" = "Gradient", "focal" = "Focal",
            "periodic" = "Periodic", "weak" = "Weak", "moderate" = "Moderate",
            "strong" = "Strong"
        ),
        Stratum = factor(
            .data$Stratum,
            levels = c("Domain", "Gradient", "Focal", "Periodic", "Weak", "Moderate", "Strong")
        )
    )

panel_b <- ggplot(recall_plot, aes(x = 1, y = .data$Method, fill = .data$MeanRecall)) +
    geom_tile(colour = "white", linewidth = 0.30) +
    facet_grid(. ~ .data$Stratum, scales = "free_x", space = "free_x") +
    scale_fill_viridis_c(
        option = "C", limits = c(0, 1), breaks = c(0, 0.5, 1),
        na.value = "grey85", name = "Mean recall",
        guide = guide_colourbar(
            title.position = "top", barwidth = unit(1.45, "cm"),
            barheight = unit(0.20, "cm")
        )
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_discrete(drop = FALSE, expand = c(0, 0)) +
    labs(x = NULL, y = NULL, tag = "b") +
    base_theme +
    theme(
        axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
        strip.text.x = element_text(size = 4.05, angle = 0, hjust = 0.5),
        legend.position = "bottom",
        legend.title.position = "top",
        legend.title = element_text(face = "bold")
    )

# c, Generalization of training-selected features to unseen slices.
heldout_long <- heldout |>
    filter(
        .data$method != "all_features",
        .data$metric %in% c("heldout_moran_percentile", "heldout_topn_jaccard")
    ) |>
    mutate(
        Metric = recode(
            .data$metric,
            "heldout_moran_percentile" = "Moran percentile",
            "heldout_topn_jaccard" = "Top-N Jaccard"
        )
    )
heldout_order <- heldout_long |>
    filter(.data$Metric == "Moran percentile") |>
    group_by(.data$method) |>
    summarise(Mean = mean(.data$value), .groups = "drop") |>
    arrange(desc(.data$Mean), .data$method) |>
    pull(.data$method)
heldout_levels <- rev(heldout_order)
heldout_long <- heldout_long |>
    mutate(
        Metric = factor(.data$Metric, levels = c("Moran percentile", "Top-N Jaccard")),
        Method = factor(method_label(.data$method), levels = method_label(heldout_levels))
    )
heldout_summary <- heldout_long |>
    group_by(.data$method, .data$Method, .data$Metric) |>
    summarise(
        Mean = mean(.data$value), SD = sd(.data$value),
        Lower = pmax(0, .data$Mean - .data$SD),
        Upper = pmin(1, .data$Mean + .data$SD),
        .groups = "drop"
    )
control_bands_c <- heldout_long |>
    filter(.data$method %in% c("random", "scsegindex")) |>
    distinct(.data$Method) |>
    mutate(y = as.numeric(.data$Method))
spatial_palette <- c("Moran percentile" = "#4daf4a", "Top-N Jaccard" = "#984ea3")

panel_c <- ggplot(heldout_long, aes(x = .data$value, y = .data$Method)) +
    geom_rect(
        data = control_bands_c,
        aes(xmin = -Inf, xmax = Inf, ymin = .data$y - 0.47, ymax = .data$y + 0.47),
        inherit.aes = FALSE, fill = "grey70", alpha = 0.28, colour = NA
    ) +
    geom_point(
        aes(colour = .data$Metric), size = 0.76, alpha = 0.40,
        position = position_jitter(height = 0.09, width = 0, seed = 23)
    ) +
    geom_linerange(
        data = heldout_summary,
        aes(xmin = .data$Lower, xmax = .data$Upper, y = .data$Method),
        inherit.aes = FALSE, linewidth = 0.38, colour = "grey25"
    ) +
    geom_point(
        data = heldout_summary,
        aes(x = .data$Mean, y = .data$Method, fill = .data$Metric),
        inherit.aes = FALSE,
        shape = 23, size = 1.80, stroke = 0.30, colour = "white"
    ) +
    facet_grid(. ~ .data$Metric) +
    scale_colour_manual(values = spatial_palette, guide = "none") +
    scale_fill_manual(values = spatial_palette, guide = "none") +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    scale_y_discrete(drop = FALSE) +
    labs(x = "Held-out spatial reproducibility", y = NULL, tag = "c") +
    base_theme +
    theme(
        axis.text.y = element_text(size = 4.8),
        axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        panel.spacing.x = unit(0.22, "cm")
    )

# d, Agreement when the held-out statistic is changed.
metric_labels <- c(
    "Jaccard" = "Top-N Jaccard",
    "MeanReferencePercentile" = "Percentile",
    "RankSpearman" = "Rank agreement"
)
rho_labels <- concordance |>
    filter(.data$Metric %in% names(metric_labels)) |>
    transmute(
        Metric = .data$Metric,
        Label = sprintf("rho = %.2f", .data$SpearmanRho),
        x = 1.5, y = 12.4
    )
panel_d <- ggplot(rank_frame, aes(x = .data[["Moran's I"]], y = .data$nnSVG)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.32, colour = "grey65") +
    geom_point(shape = 21, fill = "#f781bf", colour = "white", stroke = 0.28, size = 1.65) +
    geom_text(
        data = rho_labels,
        aes(x = .data$x, y = .data$y, label = .data$Label),
        inherit.aes = FALSE, hjust = 0, size = 1.85
    ) +
    facet_wrap(
        ~ .data$Metric, nrow = 1,
        labeller = as_labeller(metric_labels)
    ) +
    scale_x_continuous(limits = c(1, 13), breaks = c(1, 7, 13)) +
    scale_y_continuous(limits = c(1, 13), breaks = c(1, 7, 13)) +
    labs(
        x = "Rank using Moran's I reference",
        y = "Rank using nnSVG reference",
        tag = "d"
    ) +
    base_theme +
    theme(
        axis.title = element_text(size = 5.5),
        axis.text = element_text(size = 4.7),
        strip.text = element_text(size = 4.35, colour = "white", face = "bold"),
        panel.grid.major = element_line(linewidth = 0.20, colour = "grey90")
    )

# e, Association of spatial reproducibility with held-out clustering.
association <- association |>
    mutate(
        Measure = factor(
            recode(
                .data$Measure,
                "Full-ranking agreement" = "Rank agreement",
                "Mean reference percentile" = "Mean percentile"
            ),
            levels = c("Rank agreement", "Mean percentile", "Top-N Jaccard")
        ),
        Panel = "Association with ARI"
    )
panel_e <- ggplot(
    association,
    aes(x = .data$SpearmanRho, y = .data$Measure, colour = .data$Reference, shape = .data$Reference)
) +
    geom_vline(xintercept = 0, linewidth = 0.28, colour = "grey70") +
    geom_errorbar(
        aes(xmin = .data$Lower95, xmax = .data$Upper95),
        orientation = "y", width = 0.15, linewidth = 0.45,
        position = position_dodge(width = 0.36)
    ) +
    geom_point(size = 1.75, position = position_dodge(width = 0.36)) +
    facet_wrap(~ .data$Panel, nrow = 1) +
    scale_colour_manual(
        values = c("Moran's I" = "#377eb8", "nnSVG" = "#4daf4a"),
        name = "Held-out reference"
    ) +
    scale_shape_manual(
        values = c("Moran's I" = 16, "nnSVG" = 17),
        name = "Held-out reference"
    ) +
    scale_x_continuous(limits = c(-0.25, 0.72), breaks = c(-0.2, 0.2, 0.6)) +
    labs(x = "Spearman correlation with ARI", y = NULL, tag = "e") +
    base_theme +
    theme(
        axis.text.y = element_text(size = 4.55),
        strip.text = element_text(size = 5.0, colour = "white", face = "bold"),
        legend.position = "bottom",
        legend.title.position = "top",
        legend.title = element_text(face = "bold"),
        legend.key.size = unit(0.28, "cm")
    )

top_row <- panel_a + panel_b + plot_layout(widths = c(0.62, 0.38))
bottom_row <- panel_c + panel_d + panel_e + plot_layout(widths = c(0.47, 0.31, 0.22))
figure <- top_row / bottom_row +
    plot_layout(heights = c(1.36, 0.64))

save_figure <- function(plot, stem, width_mm = 183, height_mm = 190) {
    ggsave(
        paste0(stem, ".svg"), plot,
        width = width_mm, height = height_mm, units = "mm",
        device = svglite::svglite, bg = "white"
    )
    svg_path <- paste0(stem, ".svg")
    svg_text <- readLines(svg_path, warn = FALSE)
    svg_text <- gsub(" textLength='[^']+' lengthAdjust='spacingAndGlyphs'", "", svg_text)
    svg_text <- gsub(
        'font-family: "Liberation Sans";',
        "font-family: Arial, Helvetica, sans-serif;",
        svg_text, fixed = TRUE
    )
    writeLines(svg_text, svg_path, useBytes = TRUE)

    register_arial_pdf_font()
    ggsave(
        paste0(stem, ".pdf"), plot,
        width = width_mm, height = height_mm, units = "mm",
        device = grDevices::pdf, family = "Arial",
        useDingbats = FALSE, bg = "white"
    )
    ggsave(
        paste0(stem, ".tiff"), plot,
        width = width_mm, height = height_mm, units = "mm", dpi = 600,
        device = ragg::agg_tiff, compression = "lzw", bg = "white"
    )
    ggsave(
        paste0(stem, ".png"), plot,
        width = width_mm, height = height_mm, units = "mm", dpi = 300,
        device = ragg::agg_png, bg = "white"
    )
    saveRDS(plot, paste0(stem, ".rds"))
}

stem <- file.path(output_dir, "Figure2_spatial_validation")
save_figure(figure, stem)

write_tsv(auprc_plot, file.path(output_dir, "Figure2a_semi_synthetic_auprc.tsv"), na = "NA")
write_tsv(recall_plot, file.path(output_dir, "Figure2b_pattern_effect_recall.tsv"), na = "NA")
write_tsv(heldout_long, file.path(output_dir, "Figure2c_heldout_reproducibility.tsv"), na = "NA")
write_tsv(rank_frame, file.path(output_dir, "Figure2d_cross_statistic_ranks.tsv"), na = "NA")
write_tsv(association, file.path(output_dir, "Figure2e_ari_association.tsv"), na = "NA")

qa <- tibble(
    Check = c(
        "semi-synthetic methods", "semi-synthetic seeds", "held-out methods",
        "held-out folds", "cross-statistic methods", "association rows"
    ),
    Observed = c(
        n_distinct(auprc_plot$method), n_distinct(auprc_plot$seed),
        n_distinct(heldout_long$method),
        n_distinct(paste(heldout_long$dataset, heldout_long$heldout_slice)),
        n_distinct(rank_frame$Method), nrow(association)
    ),
    Expected = c(33, 3, 13, 6, 13, 6)
) |>
    mutate(Pass = .data$Observed == .data$Expected)
write_tsv(qa, file.path(output_dir, "Figure2_QA.tsv"))

if (!all(qa$Pass)) {
    stop("Figure 2 source-data QA failed")
}

message("Main Figure 2 written to ", output_dir)
