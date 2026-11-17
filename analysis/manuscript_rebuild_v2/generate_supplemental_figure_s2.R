#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(patchwork)
    library(readr)
    library(scales)
    library(tidyr)
})

source(file.path("external", "atlas-feature-selection-benchmark", "analysis", "R", "plotting.R"))

output_root <- file.path("results", "manuscript_rebuild_v2", "supplemental")
figure_dir <- file.path(output_root, "figures")
source_dir <- file.path(output_root, "source_data")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

method_label <- function(x, include_n = TRUE) {
    labels <- dplyr::recode(
        x,
        "all_features" = "All features",
        "random" = "Random",
        "scsegindex" = "scSEGIndex",
        "scanpy_cell_ranger_batch" = "scanpy-CellRanger (batch)",
        "scanpy_seurat_v3_batch" = "scanpy-SeuratV3 (batch)",
        "seurat_vst" = "Seurat-VST",
        "singleCellHaystack" = "singleCellHaystack",
        "triku" = "triku",
        "wilcoxon" = "Wilcoxon (label-informed)",
        "morans_i" = "Moran's I",
        "nnsvg" = "nnSVG",
        "somde" = "SOMDE",
        "sparkx" = "SPARK-X",
        "spatialde" = "SpatialDE",
        .default = x
    )
    if (include_n) {
        labels <- ifelse(x == "all_features", labels, paste0(labels, " (N=2,000)"))
    }
    labels
}

save_reference_files <- function(figure, stem, width_mm, height_mm) {
    ggsave(
        paste0(stem, ".svg"), figure,
        width = width_mm, height = height_mm, units = "mm",
        device = svglite::svglite, bg = "white"
    )
    svg_path <- paste0(stem, ".svg")
    svg_text <- readLines(svg_path, warn = FALSE)
    svg_text <- gsub(
        " textLength='[^']+' lengthAdjust='spacingAndGlyphs'",
        "", svg_text
    )
    svg_text <- gsub(
        'font-family: "Liberation Sans";',
        "font-family: Arial, Helvetica, sans-serif;",
        svg_text, fixed = TRUE
    )
    writeLines(svg_text, svg_path, useBytes = TRUE)

    if (!"Arial" %in% names(grDevices::pdfFonts())) {
        grDevices::pdfFonts(Arial = grDevices::pdfFonts("Helvetica")[[1]])
    }
    ggsave(
        paste0(stem, ".pdf"), figure,
        width = width_mm, height = height_mm, units = "mm",
        device = grDevices::pdf, family = "Arial",
        useDingbats = FALSE, bg = "white"
    )
    ggsave(
        paste0(stem, ".tiff"), figure,
        width = width_mm, height = height_mm, units = "mm", dpi = 600,
        device = ragg::agg_tiff, compression = "lzw", bg = "white"
    )
    ggsave(
        paste0(stem, ".png"), figure,
        width = width_mm, height = height_mm, units = "mm", dpi = 300,
        device = ragg::agg_png, bg = "white"
    )
    saveRDS(figure, paste0(stem, ".rds"))
}

metrics <- read_tsv(
    file.path("results", "validation_v1", "heldout_slice_metrics.tsv"),
    show_col_types = FALSE
)

wide <- metrics |>
    select("dataset", "heldout_slice", "method", "metric", "value") |>
    pivot_wider(names_from = "metric", values_from = "value")

method_order <- wide |>
    filter(.data$method != "all_features") |>
    group_by(.data$method) |>
    summarise(MeanSpatial = mean(.data$heldout_moran_percentile), .groups = "drop") |>
    arrange(desc(.data$MeanSpatial), .data$method) |>
    pull(.data$method)
method_order <- c(method_order, "all_features")
method_levels <- rev(method_order)

fold_labels <- c(
    "151507" = "151507",
    "151673" = "151673",
    "section1" = "Section 1",
    "section2" = "Section 2",
    "GSM5907080_AD_4_LS" = "AD 4 LS",
    "GSM5907093_HE_3" = "HE 3"
)
dataset_labels <- c(
    "DLPFC" = "DLPFC",
    "MouseBrainSerialSections" = "Mouse Brain",
    "STOmics0212" = "STOmics-0212"
)
dataset_palette <- c(
    "DLPFC" = "#8DD3C7",
    "Mouse Brain" = "#FB8072",
    "STOmics-0212" = "#80B1D3"
)

common_theme <- theme_features_pub() +
    theme(
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "black", colour = "black", linewidth = 0.35),
        strip.text = element_text(colour = "white", face = "bold", size = 6.2),
        panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.35),
        panel.spacing = unit(0.05, "cm"),
        plot.margin = margin(0.05, 0.05, 0.05, 0.05, "cm"),
        plot.tag = element_text(size = 8, face = "bold", colour = "black"),
        plot.tag.position = "topleft"
    )

spatial_palette <- c(
    "Moran percentile" = "#4DAF4A",
    "Top-N Jaccard" = "#984EA3"
)

spatial_long <- metrics |>
    filter(.data$metric %in% c("heldout_moran_percentile", "heldout_topn_jaccard")) |>
    mutate(
        Metric = recode(
            .data$metric,
            "heldout_moran_percentile" = "Moran percentile",
            "heldout_topn_jaccard" = "Top-N Jaccard"
        ),
        Metric = factor(.data$Metric, levels = names(spatial_palette)),
        Method = factor(
            .data$method,
            levels = method_levels,
            labels = method_label(method_levels)
        )
    )

spatial_summary <- spatial_long |>
    filter(!is.na(.data$value)) |>
    group_by(.data$method, .data$Method, .data$Metric) |>
    summarise(
        Mean = mean(.data$value), SD = sd(.data$value),
        Lower = pmax(0, .data$Mean - .data$SD),
        Upper = pmin(1, .data$Mean + .data$SD),
        .groups = "drop"
    )

control_bands <- spatial_long |>
    filter(.data$method %in% c("all_features", "random", "scsegindex")) |>
    distinct(.data$Method) |>
    mutate(y = as.numeric(.data$Method))

plot_a <- ggplot(
    spatial_long |> filter(!is.na(.data$value)),
    aes(x = .data$value, y = .data$Method)
) +
    geom_rect(
        data = control_bands,
        aes(xmin = -Inf, xmax = Inf, ymin = .data$y - 0.47, ymax = .data$y + 0.47),
        inherit.aes = FALSE, fill = "grey70", alpha = 0.28, colour = NA
    ) +
    geom_point(
        aes(colour = .data$Metric), size = 0.85, alpha = 0.42,
        position = position_jitter(height = 0.10, width = 0, seed = 23)
    ) +
    geom_linerange(
        data = spatial_summary,
        aes(xmin = .data$Lower, xmax = .data$Upper, y = .data$Method),
        inherit.aes = FALSE, linewidth = 0.42, colour = "grey25"
    ) +
    geom_point(
        data = spatial_summary,
        aes(x = .data$Mean, y = .data$Method, fill = .data$Metric),
        inherit.aes = FALSE,
        shape = 23, size = 2.05, stroke = 0.32, colour = "white"
    ) +
    facet_grid(. ~ .data$Metric) +
    scale_colour_manual(values = spatial_palette, guide = "none") +
    scale_fill_manual(values = spatial_palette, guide = "none") +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    scale_y_discrete(drop = FALSE) +
    labs(x = "Held-out spatial reproducibility", y = NULL, tag = "a") +
    common_theme +
    theme(
        axis.text.y = element_text(face = "plain", colour = "black", size = 4.9),
        axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.spacing.x = unit(0.25, "cm")
    )

jaccard_heat <- wide |>
    mutate(
        Dataset = factor(
            unname(dataset_labels[.data$dataset]),
            levels = unname(dataset_labels)
        ),
        Fold = factor(
            unname(fold_labels[.data$heldout_slice]),
            levels = unname(fold_labels)
        ),
        Method = factor(
            .data$method,
            levels = method_levels,
            labels = method_label(method_levels)
        )
    )

plot_b <- ggplot(jaccard_heat, aes(x = .data$Fold, y = .data$Method, fill = .data$heldout_topn_jaccard)) +
    geom_tile(colour = "white", linewidth = 0.35) +
    facet_grid(. ~ .data$Dataset, scales = "free_x", space = "free_x") +
    scale_fill_viridis_c(
        option = "C", limits = c(0, 1), breaks = c(0, 0.5, 1),
        na.value = "grey85", name = "Jaccard index",
        guide = guide_colourbar(
            title.position = "top", barwidth = unit(1.35, "cm"),
            barheight = unit(0.20, "cm"), order = 1
        )
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(drop = FALSE, expand = c(0, 0)) +
    labs(x = NULL, y = NULL, tag = "b") +
    common_theme +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 4.8),
        axis.text.y = element_text(face = "plain", colour = "black", size = 4.7),
        axis.ticks = element_blank(), axis.line = element_blank(), panel.grid = element_blank()
    )

cluster_palette <- c("ARI" = "#F781BF", "NMI" = "#377EB8")
cluster_long <- metrics |>
    filter(.data$metric %in% c("heldout_ari", "heldout_nmi")) |>
    mutate(
        Metric = recode(.data$metric, "heldout_ari" = "ARI", "heldout_nmi" = "NMI"),
        Metric = factor(.data$Metric, levels = names(cluster_palette)),
        Method = factor(
            .data$method,
            levels = method_levels,
            labels = method_label(method_levels)
        )
    )

cluster_summary <- cluster_long |>
    group_by(.data$method, .data$Method, .data$Metric) |>
    summarise(
        Mean = mean(.data$value), SD = sd(.data$value),
        Lower = pmax(0, .data$Mean - .data$SD),
        Upper = pmin(1, .data$Mean + .data$SD),
        .groups = "drop"
    )

plot_c <- ggplot(cluster_long, aes(x = .data$value, y = .data$Method)) +
    geom_rect(
        data = control_bands,
        aes(xmin = -Inf, xmax = Inf, ymin = .data$y - 0.47, ymax = .data$y + 0.47),
        inherit.aes = FALSE, fill = "grey70", alpha = 0.28, colour = NA
    ) +
    geom_point(
        aes(colour = .data$Metric), size = 0.85, alpha = 0.42,
        position = position_jitter(height = 0.10, width = 0, seed = 29)
    ) +
    geom_linerange(
        data = cluster_summary,
        aes(xmin = .data$Lower, xmax = .data$Upper, y = .data$Method),
        inherit.aes = FALSE, linewidth = 0.42, colour = "grey25"
    ) +
    geom_point(
        data = cluster_summary,
        aes(x = .data$Mean, y = .data$Method, fill = .data$Metric),
        inherit.aes = FALSE,
        shape = 23, size = 2.05, stroke = 0.32, colour = "white"
    ) +
    facet_grid(. ~ .data$Metric) +
    scale_colour_manual(values = cluster_palette, guide = "none") +
    scale_fill_manual(values = cluster_palette, guide = "none") +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    scale_y_discrete(drop = FALSE) +
    labs(x = "Held-out clustering agreement", y = NULL, tag = "c") +
    common_theme +
    theme(
        axis.text.y = element_text(face = "plain", colour = "black", size = 4.9),
        axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.spacing.x = unit(0.25, "cm")
    )

association <- wide |>
    filter(.data$method != "all_features") |>
    group_by(.data$dataset, .data$heldout_slice) |>
    mutate(
        JaccardPercentile = percent_rank(.data$heldout_topn_jaccard),
        ARIPercentile = percent_rank(.data$heldout_ari)
    ) |>
    ungroup() |>
    mutate(
        Dataset = factor(
            unname(dataset_labels[.data$dataset]),
            levels = unname(dataset_labels)
        )
    )

rho <- cor(
    association$JaccardPercentile,
    association$ARIPercentile,
    method = "spearman"
)

plot_d <- ggplot(association, aes(x = .data$JaccardPercentile, y = .data$ARIPercentile)) +
    geom_smooth(
        method = "lm", formula = y ~ x, se = TRUE,
        colour = "grey45", fill = "grey85", linewidth = 0.45
    ) +
    geom_point(
        aes(fill = .data$Dataset),
        shape = 21, size = 1.8, stroke = 0.28, colour = "grey35", alpha = 0.86
    ) +
    facet_wrap(~ factor("Within-fold association")) +
    scale_fill_manual(values = dataset_palette, name = "Dataset") +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    annotate(
        "text", x = 0.97, y = 0.06,
        label = sprintf("Spearman rho = %.2f", rho),
        hjust = 1, family = "Arial", size = 2.15, fontface = "bold"
    ) +
    labs(
        x = "Jaccard percentile", y = "ARI percentile", tag = "d"
    ) +
    common_theme +
    theme(
        legend.position = "bottom",
        legend.title.position = "top",
        legend.box.just = "left",
        legend.key.size = unit(0.28, "cm"),
        legend.text = element_text(size = 4.8),
        legend.title = element_text(size = 5.5, face = "bold")
    )

top_row <- plot_a + plot_b + plot_layout(widths = c(0.57, 0.43), guides = "collect")
bottom_row <- plot_c + plot_d + plot_layout(widths = c(0.57, 0.43), guides = "collect")
figure <- top_row / bottom_row +
    plot_layout(heights = c(1, 1), guides = "collect") &
    theme(
        legend.position = "bottom",
        legend.title.position = "top",
        legend.box = "horizontal",
        legend.box.just = "left",
        legend.margin = margin(0, 0, 0, 0),
        legend.spacing.x = unit(0.12, "cm"),
        legend.text = element_text(size = 5.0),
        legend.title = element_text(size = 5.7, face = "bold")
    )

write_tsv(spatial_long, file.path(source_dir, "Supplemental_Fig_S2a_spatial_reproducibility.tsv"), na = "NA")
write_tsv(jaccard_heat, file.path(source_dir, "Supplemental_Fig_S2b_fold_jaccard.tsv"), na = "NA")
write_tsv(cluster_long, file.path(source_dir, "Supplemental_Fig_S2c_clustering.tsv"), na = "NA")
write_tsv(association, file.path(source_dir, "Supplemental_Fig_S2d_association.tsv"), na = "NA")

qa <- tibble(
    Check = c(
        "methods", "held-out folds", "spatial rows", "clustering rows",
        "association rows", "plotted values outside [0,1]", "expected missing values"
    ),
    Observed = c(
        n_distinct(metrics$method),
        n_distinct(paste(metrics$dataset, metrics$heldout_slice)),
        nrow(spatial_long), nrow(cluster_long), nrow(association),
        sum(c(spatial_long$value, cluster_long$value) < 0 |
                c(spatial_long$value, cluster_long$value) > 1,
            na.rm = TRUE),
        sum(is.na(metrics$value))
    ),
    Expected = c(14, 6, 14 * 6 * 2, 14 * 6 * 2, 13 * 6, 0, 18)
) |>
    mutate(Pass = .data$Observed == .data$Expected)
write_tsv(qa, file.path(output_root, "Supplemental_Fig_S2_QA.tsv"))

writeLines(
    paste(
        "Supplemental Fig. S2 | Generalization of training-selected features to unseen spatial slices.",
        "a, Held-out Moran percentile and top-N Jaccard index across six folds.",
        "Small points show folds; diamonds and horizontal bars show mean and s.d.",
        "b, Fold-specific Jaccard indices between training-selected features and the top 2,000 held-out Moran-ranked genes.",
        "c, ARI and NMI after clustering held-out slices using features selected only from training slices.",
        "d, Descriptive association between within-fold Jaccard and ARI percentiles (Spearman rho = 0.37).",
        "No inferential P value is reported because method-fold observations are not independent."
    ),
    file.path(output_root, "Supplemental_Fig_S2_caption.txt")
)

save_reference_files(
    figure,
    file.path(figure_dir, "Supplemental_Fig_S2"),
    width_mm = 208.3,
    height_mm = 182
)

message("Reference-aligned Supplemental Figure S2 written to ", figure_dir)
