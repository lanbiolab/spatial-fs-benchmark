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

method_label <- function(x) {
    dplyr::recode(
        x,
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
        "scanpy_cell_ranger" = "scanpy-CellRanger",
        "scanpy_cell_ranger_batch" = "scanpy-CellRanger (batch)",
        "scanpy_pearson" = "scanpy-Pearson",
        "scanpy_pearson_batch" = "scanpy-Pearson (batch)",
        "scanpy_seurat" = "scanpy-Seurat",
        "scanpy_seurat_batch" = "scanpy-Seurat (batch)",
        "scanpy_seurat_v3" = "scanpy-SeuratV3",
        "scanpy_seurat_v3_batch" = "scanpy-SeuratV3 (batch)",
        "scry" = "scry",
        "seurat_disp" = "Seurat-Dispersion",
        "seurat_mvp" = "Seurat-MVP",
        "seurat_sct" = "Seurat-scTransform",
        "seurat_vst" = "Seurat-VST",
        "singleCellHaystack" = "singleCellHaystack",
        "statistic_mean" = "Statistic mean",
        "statistic_variance" = "Statistic variance",
        "triku" = "triku",
        "wilcoxon" = "Wilcoxon (label-informed)",
        "morans_i" = "Moran's I",
        "nnsvg" = "nnSVG",
        "somde" = "SOMDE",
        "sparkx" = "SPARK-X",
        "spatialde" = "SpatialDE",
        .default = x
    )
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
    file.path("results", "validation_v1", "semi_synthetic_metrics.tsv"),
    show_col_types = FALSE
)

auprc_raw <- metrics |>
    filter(.data$metric == "AUPRC")

method_order <- auprc_raw |>
    group_by(.data$method) |>
    summarise(MeanAUPRC = mean(.data$value), .groups = "drop") |>
    arrange(desc(.data$MeanAUPRC), .data$method) |>
    pull(.data$method)
method_levels <- rev(method_order)

prevalence_palette <- c("5% SVGs" = "#377EB8", "15% SVGs" = "#E41A1C")
control_methods <- c("all_features", "random", "scsegindex")

auprc_plot <- auprc_raw |>
    mutate(
        Prevalence = factor(
            .data$prevalence,
            levels = c(0.05, 0.15),
            labels = c("5% SVGs", "15% SVGs")
        ),
        Method = factor(
            .data$method,
            levels = method_levels,
            labels = method_label(method_levels)
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

control_bands <- auprc_plot |>
    filter(.data$method %in% control_methods) |>
    distinct(.data$Method) |>
    mutate(y = as.numeric(.data$Method))

common_theme <- theme_features_pub() +
    theme(
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "black", colour = "black", linewidth = 0.35),
        strip.text = element_text(colour = "white", face = "bold", size = 6.2),
        panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.35),
        panel.spacing = unit(0.05, "cm"),
        plot.margin = margin(0.05, 0.05, 0.05, 0.05, "cm"),
        plot.tag = element_text(size = 8, face = "bold", colour = "black"),
        plot.tag.position = "topleft"
    )

plot_a <- ggplot(auprc_plot, aes(x = .data$value, y = .data$Method)) +
    geom_rect(
        data = control_bands,
        aes(xmin = -Inf, xmax = Inf, ymin = .data$y - 0.47, ymax = .data$y + 0.47),
        inherit.aes = FALSE, fill = "grey70", alpha = 0.28, colour = NA
    ) +
    geom_point(
        aes(colour = .data$Prevalence),
        size = 0.85, alpha = 0.42,
        position = position_jitter(height = 0.10, width = 0, seed = 17)
    ) +
    geom_linerange(
        data = auprc_summary,
        aes(xmin = .data$Lower, xmax = .data$Upper, y = .data$Method),
        inherit.aes = FALSE, linewidth = 0.42, colour = "grey25"
    ) +
    geom_point(
        data = auprc_summary,
        aes(x = .data$Mean, y = .data$Method, fill = .data$Prevalence),
        inherit.aes = FALSE,
        shape = 23, size = 2.05, stroke = 0.32, colour = "white"
    ) +
    facet_grid(. ~ .data$Prevalence) +
    scale_colour_manual(values = prevalence_palette, guide = "none") +
    scale_fill_manual(values = prevalence_palette, guide = "none") +
    scale_x_continuous(
        limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1),
        expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_y_discrete(drop = FALSE) +
    labs(x = "AUPRC", y = NULL, tag = "a") +
    common_theme +
    theme(
        axis.text.y = element_text(face = "plain", colour = "black", size = 5.1),
        axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(colour = "black", size = 5.0),
        panel.spacing.x = unit(0.35, "cm")
    )

make_recall <- function(summary_name, strata, labels) {
    metrics |>
        filter(.data$summary == summary_name, .data$metric == "Recall") |>
        group_by(.data$method, .data$stratum) |>
        summarise(MeanRecall = mean(.data$value), SDRecall = sd(.data$value), .groups = "drop") |>
        complete(method = method_order, stratum = strata) |>
        mutate(
            Method = factor(
                .data$method,
                levels = method_levels,
                labels = method_label(method_levels)
            ),
            Stratum = factor(.data$stratum, levels = strata, labels = labels)
        )
}

pattern_recall <- make_recall(
    "pattern",
    c("domain", "gradient", "focal", "periodic"),
    c("Domain", "Gradient", "Focal", "Periodic")
)
effect_recall <- make_recall(
    "effect",
    c("weak", "moderate", "strong"),
    c("Weak", "Moderate", "Strong")
)

heat_scale <- scale_fill_viridis_c(
    option = "C", limits = c(0, 1), breaks = c(0, 0.5, 1),
    na.value = "grey85", name = "Mean recall",
    guide = guide_colourbar(
        title.position = "top", barwidth = unit(1.45, "cm"),
        barheight = unit(0.20, "cm"), order = 1
    )
)

plot_b <- ggplot(pattern_recall, aes(x = 1, y = .data$Method, fill = .data$MeanRecall)) +
    geom_tile(colour = "white", linewidth = 0.35) +
    facet_grid(. ~ .data$Stratum) +
    heat_scale +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_discrete(drop = FALSE, expand = c(0, 0)) +
    labs(x = NULL, y = NULL, tag = "b") +
    common_theme +
    theme(
        axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank()
    )

plot_c <- ggplot(effect_recall, aes(x = 1, y = .data$Method, fill = .data$MeanRecall)) +
    geom_tile(colour = "white", linewidth = 0.35) +
    facet_grid(. ~ .data$Stratum) +
    heat_scale +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_discrete(drop = FALSE, expand = c(0, 0)) +
    labs(x = NULL, y = NULL, tag = "c") +
    common_theme +
    theme(
        axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank()
    )

figure <- wrap_plots(
    plot_a, plot_b, plot_c,
    nrow = 1, widths = c(0.58, 0.24, 0.18), guides = "collect"
) &
    theme(
        legend.position = "bottom",
        legend.title.position = "top",
        legend.box = "horizontal",
        legend.box.just = "left",
        legend.margin = margin(0, 0, 0, 0),
        legend.key.height = unit(0.20, "cm"),
        legend.key.width = unit(0.42, "cm"),
        legend.text = element_text(size = 5.2),
        legend.title = element_text(size = 5.8, face = "bold")
    )

write_tsv(auprc_plot, file.path(source_dir, "Supplemental_Fig_S1a_AUPRC.tsv"), na = "NA")
write_tsv(pattern_recall, file.path(source_dir, "Supplemental_Fig_S1b_pattern_recall.tsv"), na = "NA")
write_tsv(effect_recall, file.path(source_dir, "Supplemental_Fig_S1c_effect_recall.tsv"), na = "NA")

qa <- tibble(
    Check = c(
        "methods", "prevalence levels", "pattern cells", "effect cells",
        "AUPRC outside [0,1]", "recall outside [0,1]"
    ),
    Observed = c(
        n_distinct(auprc_plot$method), n_distinct(auprc_plot$Prevalence),
        nrow(pattern_recall), nrow(effect_recall),
        sum(auprc_plot$value < 0 | auprc_plot$value > 1),
        sum(c(pattern_recall$MeanRecall, effect_recall$MeanRecall) < 0 |
                c(pattern_recall$MeanRecall, effect_recall$MeanRecall) > 1,
            na.rm = TRUE)
    ),
    Expected = c(33, 2, 33 * 4, 33 * 3, 0, 0)
) |>
    mutate(Pass = .data$Observed == .data$Expected)
write_tsv(qa, file.path(output_root, "Supplemental_Fig_S1_QA.tsv"))

writeLines(
    paste(
        "Supplemental Fig. S1 | Recovery of known spatially variable genes in semi-synthetic data.",
        "a, AUPRC for 33 feature-selection methods when 5% or 15% of genes were simulated as spatially variable.",
        "Small points show three independent simulation seeds; diamonds and horizontal bars show mean and s.d.",
        "b, Mean recall at N = 500 for domain, gradient, focal and periodic spatial patterns.",
        "c, Mean recall at N = 500 for weak, moderate and strong effects.",
        "Gray cells indicate that ranking-based recall is not defined for the all-features control."
    ),
    file.path(output_root, "Supplemental_Fig_S1_caption.txt")
)

save_reference_files(
    figure,
    file.path(figure_dir, "Supplemental_Fig_S1"),
    width_mm = 208.3,
    height_mm = 178
)

message("Reference-aligned Supplemental Figure S1 written to ", figure_dir)
