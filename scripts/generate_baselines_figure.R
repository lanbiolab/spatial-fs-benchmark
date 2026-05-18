#!/usr/bin/env Rscript
# Fig2-style figure: baseline visualization (panel a) + pipeline illustration (panel b)

suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(patchwork)

})

args       <- commandArgs(trailingOnly = TRUE)
data_dir   <- if (length(args) >= 1) args[[1]] else file.path("results", "current_rank", "data")
output_dir <- if (length(args) >= 2) args[[2]] else file.path("results", "metric_overview", "figures")

source(file.path("external", "atlas-feature-selection-benchmark", "analysis", "R", "plotting.R"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ── Constants ─────────────────────────────────────────────────────────────────
baseline_methods_fixed <- c(
    "scanpy_cell_ranger_batch-N2000",
    "all_features-Nall",
    "random-N500",
    "scsegindex-N500"
)

type_levels <- c("Integration", "Clustering", "Alignment")
type_colors <- c(
    "Integration" = "#6a3d9a",
    "Clustering"  = "#e84a8a",
    "Alignment"   = "#f0b429"
)

metric_order <- c(
    "bASW", "dASW", "iLISI", "dLISI", "ILL", "GC",
    "ARI", "NMI", "CHAOS", "PAS", "Silhouette",
    "Accuracy", "Ratio"
)

dataset_labels <- c(
    "DLPFC"                    = "DLPFC",
    "E8p5Embryo"               = "E8.5 Embryo",
    "E9p5Embryo"               = "E9.5 Embryo",
    "MouseBrainSerialSections"  = "Mouse Brain",
    "STOmics0212"              = "STOmics-0212",
    "STOmics0218"              = "STOmics-0218",
    "STOmics0224"              = "STOmics-0224",
    "STOmicsVisium5Samples"    = "STOmics-V5"
)

baseline_labels <- c(
    "scanpy_cell_ranger_batch-N2000" = "Scanpy CellRanger (N=2000, batch)",
    "all_features-Nall"              = "All features",
    "scsegindex-N500"                = "scSEGIndex (N=500)",
    "random-N500"                    = "Random (N=500)"
)
baseline_colors <- c(
    "Scanpy CellRanger (N=2000, batch)" = "#7570b3",
    "All features"                     = "#1b9e77",
    "scSEGIndex (N=500)"               = "#66a61e",
    "Random (N=500)"                   = "#d95f02"
)
baseline_shapes <- c(21, 23, 22, 24)

group_order <- c(
    '"Best" method',
    "All features", "Random (N=500)",
    "Scanpy CellRanger (N=2000, batch)", "scSEGIndex (N=500)",
    '"Worst" method'
)
group_colors <- c(
    '"Best" method'              = "#27ae60",
    "All features"               = "#1b9e77",
    "Random (N=500)"             = "#d95f02",
    "Scanpy CellRanger (N=2000, batch)" = "#7570b3",
    "scSEGIndex (N=500)"         = "#66a61e",
    '"Worst" method'             = "#c0392b"
)

EXAMPLE_DATASET <- "DLPFC"

# ── Load & compute baselines ──────────────────────────────────────────────────
metrics_raw <- readr::read_tsv(
    file.path(data_dir, "benchmark.tsv"), show_col_types = FALSE
) |>
    dplyr::filter(.data$Type %in% type_levels)

baselines <- metrics_raw |>
    dplyr::filter(.data$Method %in% baseline_methods_fixed, !is.na(.data$Value)) |>
    dplyr::group_by(.data$Dataset, .data$Method, .data$MetricName, .data$Type) |>
    dplyr::summarise(MeanValue = mean(.data$Value, na.rm = TRUE), .groups = "drop")

baseline_ranges <- baselines |>
    dplyr::group_by(.data$Dataset, .data$MetricName, .data$Type) |>
    dplyr::summarise(
        Lower = min(.data$MeanValue),
        Upper = max(.data$MeanValue),
        .groups = "drop"
    )

# ── Panel a: grid baseline dot plot (Dataset × Type) ─────────────────────────
# Order metrics by type (Integration → Clustering → Alignment) for grouped facet_wrap
metric_order_typed <- c(
    "bASW", "dASW", "iLISI", "dLISI", "ILL", "GC",         # Integration
    "ARI", "NMI", "CHAOS", "PAS", "Silhouette",              # Clustering
    "Accuracy", "Ratio"                                        # Alignment
)
type_for_metric <- c(
    bASW = "Integration", dASW = "Integration", iLISI = "Integration",
    dLISI = "Integration", ILL = "Integration", GC = "Integration",
    ARI = "Clustering", NMI = "Clustering", CHAOS = "Clustering",
    PAS = "Clustering", Silhouette = "Clustering",
    Accuracy = "Alignment", Ratio = "Alignment"
)
strip_fills <- type_colors[type_for_metric[metric_order_typed]]

baselines_plt <- baselines |>
    dplyr::mutate(
        Dataset    = dplyr::recode(.data$Dataset, !!!dataset_labels),
        Type       = factor(.data$Type,       levels = type_levels),
        MetricName = factor(.data$MetricName, levels = metric_order_typed),
        Method = factor(.data$Method, levels = names(baseline_labels), labels = baseline_labels)
    )

ranges_plt <- baseline_ranges |>
    dplyr::mutate(
        Dataset    = dplyr::recode(.data$Dataset, !!!dataset_labels),
        Type       = factor(.data$Type,       levels = type_levels),
        MetricName = factor(.data$MetricName, levels = metric_order_typed)
    )

panel_a <- ggplot2::ggplot(
    ranges_plt,
    ggplot2::aes(y = .data$Dataset)
) +
    ggplot2::geom_linerange(
        ggplot2::aes(xmin = .data$Lower, xmax = .data$Upper, colour = .data$Type),
        linewidth = 3.0, alpha = 0.35
    ) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.25, linetype = "dashed") +
    ggplot2::geom_point(
        data = baselines_plt,
        ggplot2::aes(x = .data$MeanValue, shape = .data$Method, fill = .data$Method,
                     group = .data$Method),
        size = 2.0, colour = "white",
        position = ggplot2::position_dodge(width = 0.85)
    ) +
    ggplot2::scale_x_continuous(
        breaks = scales::pretty_breaks(n = 4),
        expand = ggplot2::expansion(mult = c(0.08, 0.08))
    ) +
    ggplot2::scale_colour_manual(values = type_colors, name = "Metric type") +
    ggplot2::scale_shape_manual(
        values = baseline_shapes,
        name = "Baseline method",
        guide = ggplot2::guide_legend(
            title.position = "top",
            order = 1,
            nrow = 2,
            byrow = TRUE,
            override.aes = list(
                fill = unname(baseline_colors),
                colour = "white",
                size = 2.2
            )
        )
    ) +
    ggplot2::scale_fill_manual(values = baseline_colors, name = "Baseline method", guide = "none") +
    ggplot2::facet_wrap(~ .data$MetricName, nrow = 3, scales = "free_x") +
    ggplot2::labs(x = "Metric value (higher is better)", y = NULL, colour = "Metric type", fill = "Baseline method") +
    ggplot2::guides(
        colour = ggplot2::guide_legend(
            title.position = "top",
            order = 2,
            nrow = 1,
            override.aes = list(linewidth = 3, alpha = 0.8, shape = NA)
        )
    ) +
    theme_features_pub() +
    ggplot2::theme(
        legend.position    = c(0.61, 0.22),
        legend.box         = "vertical",
        legend.justification = c(0, 0.5),
        legend.box.just    = "left",
        legend.spacing.x   = ggplot2::unit(0.15, "cm"),
        legend.spacing.y   = ggplot2::unit(0.03, "cm"),
        legend.margin      = ggplot2::margin(0, 0, 0, 0),
        legend.box.margin  = ggplot2::margin(0, 0, 0, 0),
        legend.key.width   = ggplot2::unit(0.45, "cm"),
        axis.title.y       = ggplot2::element_blank(),
        axis.text.y        = ggplot2::element_text(size = 5.5, face = "bold", colour = "black"),
        axis.text.x        = ggplot2::element_text(size = 5),
        strip.text         = ggplot2::element_text(face = "bold", size = 6.5, colour = "white",
                                                    margin = ggplot2::margin(1, 0, 1, 0)),
        strip.background   = ggplot2::element_rect(fill = "#2d2d2d", colour = "#2d2d2d"),
        panel.border       = ggplot2::element_rect(fill = NA, colour = "grey70", linewidth = 0.3),
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.minor   = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_line(colour = "grey88", linewidth = 0.25),
        panel.spacing      = ggplot2::unit(0.02, "cm"),
        plot.margin        = ggplot2::margin(1, 2, 0.5, 2)
    )

# ── Panel b: pipeline illustration ───────────────────────────────────────────
# Identify best & worst non-baseline methods on the example dataset
ranges_ex <- baseline_ranges |> dplyr::filter(.data$Dataset == EXAMPLE_DATASET)

ranked <- metrics_raw |>
    dplyr::filter(
        .data$Dataset == EXAMPLE_DATASET,
        !(.data$Method %in% baseline_methods_fixed),
        !is.na(.data$Value)
    ) |>
    dplyr::group_by(.data$Method, .data$MetricName, .data$Type) |>
    dplyr::summarise(Value = mean(.data$Value, na.rm = TRUE), .groups = "drop") |>
    dplyr::left_join(ranges_ex, by = c("MetricName", "Type")) |>
    dplyr::mutate(Scaled = dplyr::if_else(
        is.finite(.data$Lower) & is.finite(.data$Upper) & (.data$Upper - .data$Lower) > 0,
        (.data$Value - .data$Lower) / (.data$Upper - .data$Lower),
        NA_real_
    )) |>
    dplyr::group_by(.data$Method) |>
    dplyr::summarise(Score = mean(.data$Scaled, na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(.data$Score))

# Pick methods with scores in a moderate range for clean display
n <- nrow(ranked)
reasonable <- ranked |> dplyr::filter(.data$Score >= 0.3, .data$Score <= 1.2)
best_method  <- if (nrow(reasonable) >= 2) reasonable$Method[1] else ranked$Method[max(1, round(n * 0.1))]
worst_method <- if (nrow(reasonable) >= 2) reasonable$Method[nrow(reasonable)] else ranked$Method[min(n, round(n * 0.85))]

make_method_rows <- function(method_str, group_label) {
    metrics_raw |>
        dplyr::filter(.data$Dataset == EXAMPLE_DATASET, .data$Method == method_str, !is.na(.data$Value)) |>
        dplyr::group_by(.data$MetricName, .data$Type) |>
        dplyr::summarise(Value = mean(.data$Value, na.rm = TRUE), .groups = "drop") |>
        dplyr::mutate(Group = group_label)
}

baselines_ex <- baselines |>
    dplyr::filter(.data$Dataset == EXAMPLE_DATASET) |>
    dplyr::rename(Value = .data$MeanValue) |>
    dplyr::mutate(Group = dplyr::recode(.data$Method, !!!baseline_labels))

example_values <- dplyr::bind_rows(
    make_method_rows(best_method,  '"Best" method'),
    baselines_ex |> dplyr::select("MetricName", "Type", "Value", "Group"),
    make_method_rows(worst_method, '"Worst" method')
) |>
    dplyr::mutate(
        MetricName = factor(.data$MetricName, levels = metric_order),
        Type       = factor(.data$Type,       levels = type_levels),
        Group      = factor(.data$Group,      levels = group_order)
    )

# Step 2: scaled
example_scaled <- example_values |>
    dplyr::left_join(ranges_ex, by = c("MetricName", "Type")) |>
    dplyr::mutate(ScaledValue = dplyr::if_else(
        is.finite(.data$Lower) & is.finite(.data$Upper) & (.data$Upper - .data$Lower) > 0,
        (.data$Value - .data$Lower) / (.data$Upper - .data$Lower),
        NA_real_
    ))

# Step 3: type means
example_means <- example_scaled |>
    dplyr::group_by(.data$Group, .data$Type) |>
    dplyr::summarise(TypeMean = mean(.data$ScaledValue, na.rm = TRUE), .groups = "drop") |>
    dplyr::mutate(Group = factor(.data$Group, levels = group_order))

# Step 4: overall score
example_overall <- example_means |>
    dplyr::group_by(.data$Group) |>
    dplyr::summarise(Overall = mean(.data$TypeMean, na.rm = TRUE), .groups = "drop") |>
    dplyr::mutate(Group = factor(.data$Group, levels = group_order))

# Shared theme for pipeline panels
pipe_theme <- theme_features_pub() +
    ggplot2::theme(
        panel.border        = ggplot2::element_blank(),
        panel.grid.major.y  = ggplot2::element_blank(),
        panel.grid.minor.x  = ggplot2::element_blank(),
        panel.grid.major.x  = ggplot2::element_line(colour = "grey88", linewidth = 0.25),
        axis.title.y        = ggplot2::element_blank(),
        strip.text          = ggplot2::element_blank(),
        strip.background    = ggplot2::element_blank(),
        legend.position     = "none",
        plot.title          = ggplot2::element_text(face = "bold", size = 7, hjust = 0,
                                                    margin = ggplot2::margin(b = 3))
    )

step1 <- ggplot2::ggplot(
    example_values,
    ggplot2::aes(y = .data$Group, x = .data$Value, fill = .data$Type, group = .data$MetricName)
) +
    ggplot2::geom_col(position = "dodge", width = 0.75) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
    ggplot2::scale_x_continuous(
        breaks = scales::pretty_breaks(n = 3),
        expand = ggplot2::expansion(mult = c(0.04, 0.04))
    ) +
    ggplot2::scale_fill_manual(values = type_colors) +
    ggplot2::facet_grid(.data$Group ~ .data$Type, scales = "free", space = "free_y") +
    ggplot2::labs(title = "1  Measure metrics", x = "Value") +
    pipe_theme +
    ggplot2::theme(
        axis.text.y   = ggplot2::element_text(size = 5.5, face = "bold"),
        strip.text.x  = ggplot2::element_text(size = 5.5, face = "bold", colour = "white",
                                               margin = ggplot2::margin(1.5, 0, 1.5, 0)),
        strip.text.y  = ggplot2::element_blank(),
        strip.background.x = ggplot2::element_rect(fill = "#2d2d2d", colour = "#2d2d2d"),
        strip.background.y = ggplot2::element_blank()
    )

step2 <- ggplot2::ggplot(
    example_scaled,
    ggplot2::aes(y = .data$Group, x = .data$ScaledValue, fill = .data$Type, group = .data$MetricName)
) +
    ggplot2::geom_col(position = "dodge", width = 0.75) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
    ggplot2::scale_x_continuous(
        breaks = c(-0.25, 0, 0.25, 0.5, 0.75, 1.0, 1.25),
        labels = c("", "0", ".25", ".50", ".75", "1", ""),
        expand = ggplot2::expansion(mult = c(0.02, 0.04))
    ) +
    ggplot2::coord_cartesian(xlim = c(-0.4, 1.35)) +
    ggplot2::scale_fill_manual(values = type_colors) +
    ggplot2::facet_grid(.data$Group ~ ., scales = "free_y", space = "free_y") +
    ggplot2::labs(title = "2  Scale using baselines", x = "Scaled value") +
    pipe_theme +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank())

step3 <- ggplot2::ggplot(
    example_means,
    ggplot2::aes(y = .data$Group, x = .data$TypeMean, fill = .data$Type)
) +
    ggplot2::geom_col(position = "dodge", width = 0.7) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
    ggplot2::scale_x_continuous(
        breaks = c(-0.25, 0, 0.25, 0.5, 0.75, 1.0, 1.25),
        labels = c("", "0", ".25", ".50", ".75", "1", ""),
        expand = ggplot2::expansion(mult = c(0.02, 0.04))
    ) +
    ggplot2::coord_cartesian(xlim = c(-0.4, 1.35)) +
    ggplot2::scale_fill_manual(values = type_colors) +
    ggplot2::facet_grid(.data$Group ~ ., scales = "free_y", space = "free_y") +
    ggplot2::labs(title = "3  Average by metric type", x = "Type mean value") +
    pipe_theme +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank())

step4 <- ggplot2::ggplot(
    example_overall,
    ggplot2::aes(y = .data$Group, x = .data$Overall, fill = .data$Group)
) +
    ggplot2::geom_col(width = 0.6) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
    ggplot2::scale_x_continuous(
        breaks = c(0, 0.25, 0.5, 0.75, 1.0),
        labels = c("0", ".25", ".50", ".75", "1"),
        expand = ggplot2::expansion(mult = c(0.02, 0.04))
    ) +
    ggplot2::coord_cartesian(xlim = c(-0.1, 1.25)) +
    ggplot2::scale_fill_manual(values = group_colors) +
    ggplot2::facet_grid(.data$Group ~ ., scales = "free_y", space = "free_y") +
    ggplot2::labs(title = "4  Overall score", x = "Overall score") +
    pipe_theme +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank())

panel_b <- patchwork::wrap_plots(
    step1, step2, step3, step4,
    nrow = 1, widths = c(2.2, 1.5, 1.5, 1.3)
)

# ── Overall score annotation for panel b ─────────────────────────────────────
formula_note <- paste0(
    "Overall = \u00bd \u00d7 Integration + \u00bd \u00d7 Clustering + \u00bd \u00d7 Alignment   ",
    "  (example dataset: ", EXAMPLE_DATASET,
    " | best: ", best_method, " | worst: ", worst_method, ")"
)

# ── Combine ───────────────────────────────────────────────────────────────────
figure <- (
    (panel_a   + patchwork::plot_annotation(tag_levels = list(c("a")))) /
    (panel_b   + patchwork::plot_annotation(tag_levels = list(c("b"))))
) +
    patchwork::plot_layout(heights = c(1.6, 1)) +
    patchwork::plot_annotation(caption = formula_note) &
    ggplot2::theme(
        plot.tag     = ggplot2::element_text(face = "bold", size = 10),
        plot.caption = ggplot2::element_text(size = 5.5, colour = "grey40",
                                             hjust = 0, margin = ggplot2::margin(t = 4))
    )

# ── Save ──────────────────────────────────────────────────────────────────────
panel_a_figure <- panel_a +
    patchwork::plot_annotation(tag_levels = list(c("a"))) &
    ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 10))

panel_b_figure <- panel_b +
    patchwork::plot_annotation(tag_levels = list(c("b"))) &
    ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 10))

upper_png_path <- file.path(output_dir, "figure_baselines_upper.png")
upper_pdf_path <- file.path(output_dir, "figure_baselines_upper.pdf")
lower_png_path <- file.path(output_dir, "figure_baselines_lower.png")
lower_pdf_path <- file.path(output_dir, "figure_baselines_lower.pdf")
png_path <- file.path(output_dir, "figure_baselines.png")
pdf_path <- file.path(output_dir, "figure_baselines.pdf")
ggplot2::ggsave(upper_png_path, panel_a_figure, width = 8.5, height = 6.6, dpi = 600)
ggplot2::ggsave(upper_pdf_path, panel_a_figure, width = 8.5, height = 6.6, device = grDevices::cairo_pdf)
ggplot2::ggsave(lower_png_path, panel_b_figure, width = 11, height = 4.8, dpi = 600)
ggplot2::ggsave(lower_pdf_path, panel_b_figure, width = 11, height = 4.8, device = grDevices::cairo_pdf)
ggplot2::ggsave(png_path, figure, width = 11, height = 13, dpi = 600)
ggplot2::ggsave(pdf_path, figure, width = 11, height = 13, device = grDevices::cairo_pdf)

cat("Saved:", upper_png_path, "\n")
cat("Saved:", lower_png_path, "\n")
cat("Saved:", png_path, "\n")
cat("Best method :", best_method,  "(score =", round(ranked$Score[1], 3), ")\n")
cat("Worst method:", worst_method, "(score =", round(ranked$Score[nrow(ranked)], 3), ")\n")
