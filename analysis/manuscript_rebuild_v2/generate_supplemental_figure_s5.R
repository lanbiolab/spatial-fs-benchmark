#!/usr/bin/env Rscript

source(file.path("analysis", "manuscript_rebuild_v2", "R", "reference_supplement.R"))
suppressPackageStartupMessages(library(ggrepel))

score_root <- file.path("results", "spatial_svg_rebuild_v1", "frozen_scores")
representative_scores <- read_tsv(file.path(score_root, "representative_task_scores.tsv"), show_col_types = FALSE)
dataset_ranks <- read_tsv(file.path(score_root, "dataset_representative_ranks.tsv"), show_col_types = FALSE)
global_ranks <- read_tsv(file.path(score_root, "global_method_ranks.tsv"), show_col_types = FALSE)

competitive_groups <- c("expression_driven", "spatially_informed")
dataset_levels <- unname(supplement_dataset_labels)
task_levels <- c("Core overall", "Integration", "Clustering", "Alignment")

paired_scores <- representative_scores |>
    filter(.data$MethodGroup %in% competitive_groups) |>
    select(
        "dataset", "fs_method", "n_features", "MethodGroup", "integration_method",
        "CoreOverallMean", "IntegrationMean", "ClusteringMean", "AlignmentMean"
    ) |>
    pivot_wider(
        names_from = "integration_method",
        values_from = c("CoreOverallMean", "IntegrationMean", "ClusteringMean", "AlignmentMean")
    ) |>
    mutate(
        Dataset = factor(display_supplement_dataset(.data$dataset), levels = dataset_levels),
        Method = display_supplement_method(.data$fs_method),
        `Core overall` = .data$CoreOverallMean_cellcharter - .data$CoreOverallMean_scvi,
        Integration = .data$IntegrationMean_cellcharter - .data$IntegrationMean_scvi,
        Clustering = .data$ClusteringMean_cellcharter - .data$ClusteringMean_scvi,
        Alignment = .data$AlignmentMean_cellcharter - .data$AlignmentMean_scvi
    )

paired_ranks <- dataset_ranks |>
    filter(.data$MethodGroup %in% competitive_groups) |>
    select(
        "dataset", "fs_method", "n_features", "MethodGroup", "integration_method",
        "CoreOverallRank", "IntegrationRank", "ClusteringRank", "AlignmentRank"
    ) |>
    pivot_wider(
        names_from = "integration_method",
        values_from = c("CoreOverallRank", "IntegrationRank", "ClusteringRank", "AlignmentRank")
    ) |>
    mutate(
        Dataset = factor(display_supplement_dataset(.data$dataset), levels = dataset_levels),
        Method = display_supplement_method(.data$fs_method),
        `Core overall` = .data$CoreOverallRank_cellcharter - .data$CoreOverallRank_scvi,
        Integration = .data$IntegrationRank_cellcharter - .data$IntegrationRank_scvi,
        Clustering = .data$ClusteringRank_cellcharter - .data$ClusteringRank_scvi,
        Alignment = .data$AlignmentRank_cellcharter - .data$AlignmentRank_scvi
    )

method_order <- dataset_ranks |>
    filter(.data$MethodGroup %in% competitive_groups) |>
    group_by(.data$fs_method) |>
    summarise(MeanRank = mean(.data$CoreOverallRank, na.rm = TRUE), .groups = "drop") |>
    mutate(Method = display_supplement_method(.data$fs_method)) |>
    arrange(desc(.data$MeanRank), .data$Method) |>
    pull(.data$Method)

score_long <- paired_scores |>
    select("dataset", "fs_method", "n_features", "MethodGroup", "Dataset", "Method", all_of(task_levels)) |>
    pivot_longer(cols = all_of(task_levels), names_to = "Task", values_to = "Difference") |>
    mutate(
        Task = factor(.data$Task, levels = task_levels),
        Method = factor(.data$Method, levels = method_order),
        Dataset = factor(.data$Dataset, levels = dataset_levels)
    )

rank_long <- paired_ranks |>
    select("dataset", "fs_method", "n_features", "MethodGroup", "Dataset", "Method", all_of(task_levels)) |>
    pivot_longer(cols = all_of(task_levels), names_to = "Task", values_to = "Difference") |>
    mutate(
        Task = factor(.data$Task, levels = task_levels),
        Method = factor(.data$Method, levels = method_order),
        Dataset = factor(.data$Dataset, levels = dataset_levels)
    )

score_limit <- max(abs(score_long$Difference), na.rm = TRUE)
rank_limit <- max(abs(rank_long$Difference), na.rm = TRUE)

score_heatmap <- score_long |>
    filter(is.finite(.data$Difference))
rank_heatmap <- rank_long |>
    filter(is.finite(.data$Difference))

plot_a <- ggplot(score_heatmap, aes(x = .data$Dataset, y = .data$Method, fill = .data$Difference)) +
    geom_tile(colour = "white", linewidth = 0.34) +
    facet_grid(. ~ .data$Task, scales = "free_x", space = "free_x") +
    scale_fill_gradient2(
        low = "#762A83", mid = "#F7F7F7", high = "#1B7837",
        midpoint = 0, limits = c(-score_limit, score_limit),
        breaks = c(-score_limit, 0, score_limit), labels = label_number(accuracy = 0.01),
        oob = squish, na.value = "grey85",
        name = "Score difference\nto scVI",
        guide = guide_colourbar(
            title.position = "top", barwidth = unit(1.8, "cm"),
            barheight = unit(0.20, "cm"), order = 1
        )
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0), drop = FALSE) +
    labs(x = NULL, y = NULL, tag = "a") +
    reference_supplement_theme() +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 4.2, face = "bold"),
        axis.text.y = element_text(face = "bold", colour = "black", size = 4.5),
        axis.ticks = element_blank(), axis.line = element_blank(), panel.grid = element_blank(),
        legend.position = "bottom", legend.title.position = "top", legend.box.just = "left"
    )

plot_b <- ggplot(rank_heatmap, aes(x = .data$Dataset, y = .data$Method, fill = .data$Difference)) +
    geom_tile(colour = "white", linewidth = 0.34) +
    facet_grid(. ~ .data$Task, scales = "free_x", space = "free_x") +
    scale_fill_gradient2(
        low = "#1B7837", mid = "#F7F7F7", high = "#762A83",
        midpoint = 0, limits = c(-rank_limit, rank_limit),
        breaks = c(-rank_limit, 0, rank_limit), labels = label_number(accuracy = 1),
        oob = squish, na.value = "grey85",
        name = "Rank difference\nto scVI",
        guide = guide_colourbar(
            title.position = "top", barwidth = unit(1.8, "cm"),
            barheight = unit(0.20, "cm"), order = 1
        )
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0), drop = FALSE) +
    labs(x = NULL, y = NULL, tag = "b") +
    reference_supplement_theme() +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 4.2, face = "bold"),
        axis.text.y = element_text(face = "bold", colour = "black", size = 4.5),
        axis.ticks = element_blank(), axis.line = element_blank(), panel.grid = element_blank(),
        legend.position = "bottom", legend.title.position = "top", legend.box.just = "left"
    )

global_pair <- global_ranks |>
    filter(.data$MethodGroup %in% competitive_groups) |>
    select("fs_method", "MethodGroup", "integration_method", "MeanCoreOverallRank") |>
    pivot_wider(names_from = "integration_method", values_from = "MeanCoreOverallRank") |>
    mutate(
        Method = display_supplement_method(.data$fs_method),
        Group = recode(
            .data$MethodGroup,
            expression_driven = "Expression-based",
            spatially_informed = "Spatially informed"
        ),
        Label = if_else(.data$MethodGroup == "spatially_informed", .data$Method, NA_character_)
    )

rank_rho <- cor(global_pair$scvi, global_pair$cellcharter, method = "spearman")
rank_max <- ceiling(max(c(global_pair$scvi, global_pair$cellcharter), na.rm = TRUE))
group_palette <- c("Expression-based" = "#466B8A", "Spatially informed" = "#D46A4C")

plot_c <- ggplot(global_pair, aes(x = .data$scvi, y = .data$cellcharter)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey55", linewidth = 0.42, linetype = "dashed") +
    geom_point(
        aes(fill = .data$Group), shape = 21, size = 2.0,
        stroke = 0.3, colour = "grey25", alpha = 0.9
    ) +
    geom_text_repel(
        aes(label = .data$Label), family = "Arial", size = 1.95,
        min.segment.length = 0, segment.size = 0.25,
        box.padding = 0.25, point.padding = 0.18,
        max.overlaps = Inf, seed = 71, na.rm = TRUE
    ) +
    annotate(
        "text", x = 1.2, y = rank_max - 0.5,
        label = sprintf("Spearman rho = %.2f", rank_rho),
        hjust = 0, family = "Arial", fontface = "bold", size = 2.05
    ) +
    facet_wrap(~ factor("Global mean rank")) +
    scale_fill_manual(values = group_palette, name = "Method type") +
    scale_x_continuous(limits = c(1, rank_max), breaks = pretty_breaks(5)) +
    scale_y_continuous(limits = c(1, rank_max), breaks = pretty_breaks(5)) +
    coord_equal() +
    labs(x = "Mean rank under scVI", y = "Mean rank under CellCharter", tag = "c") +
    reference_supplement_theme() +
    theme(
        legend.position = "bottom", legend.title.position = "top", legend.box.just = "left",
        panel.grid.major = element_line(colour = "grey92", linewidth = 0.25)
    )

score_distribution <- score_long |>
    filter(is.finite(.data$Difference))

plot_d <- ggplot(score_distribution, aes(x = .data$Task, y = .data$Difference, fill = .data$Task)) +
    geom_hline(yintercept = 0, colour = "grey45", linewidth = 0.38, linetype = "dashed") +
    geom_violin(width = 0.82, scale = "width", colour = NA, alpha = 0.32, trim = TRUE) +
    geom_boxplot(
        width = 0.24, outlier.shape = NA, fill = "white",
        colour = "grey25", linewidth = 0.35
    ) +
    geom_point(
        aes(colour = .data$Task), size = 0.62, alpha = 0.32,
        position = position_jitter(width = 0.14, height = 0, seed = 79)
    ) +
    facet_wrap(~ factor("CellCharter minus scVI")) +
    scale_fill_manual(values = supplement_task_palette, guide = "none") +
    scale_colour_manual(values = supplement_task_palette, guide = "none") +
    scale_y_continuous(breaks = pretty_breaks(6)) +
    labs(x = NULL, y = "Score difference", tag = "d") +
    reference_supplement_theme() +
    theme(
        axis.text.x = element_text(face = "bold", size = 5.0),
        panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.25),
        legend.position = "none"
    )

figure <- plot_a / plot_b / (plot_c + plot_d + plot_layout(widths = c(0.46, 0.54), guides = "collect")) +
    plot_layout(heights = c(1.45, 1.45, 0.85), guides = "keep") &
    theme(
        legend.position = "bottom", legend.title.position = "top",
        legend.box = "horizontal", legend.box.just = "left",
        legend.margin = margin(0, 0, 0, 0),
        legend.spacing.x = unit(0.12, "cm"),
        legend.text = element_text(size = 5.0),
        legend.title = element_text(size = 5.7, face = "bold")
    )

write_tsv(score_long, file.path(supplement_source_dir, "Supplemental_Fig_S5a_score_differences.tsv"), na = "NA")
write_tsv(rank_long, file.path(supplement_source_dir, "Supplemental_Fig_S5b_rank_differences.tsv"), na = "NA")
write_tsv(global_pair, file.path(supplement_source_dir, "Supplemental_Fig_S5c_global_rank_agreement.tsv"), na = "NA")
write_tsv(score_distribution, file.path(supplement_source_dir, "Supplemental_Fig_S5d_difference_distributions.tsv"), na = "NA")

task_summary <- score_distribution |>
    group_by(.data$Task) |>
    summarise(MeanDifference = mean(.data$Difference), N = n(), .groups = "drop")

qa <- tibble(
    Check = c(
        "paired competitive settings", "competitive methods", "score heatmap rows",
        "rank heatmap rows", "plotted score cells", "plotted rank cells",
        "finite alignment differences", "global rank pairs",
        "global rank Spearman correlation", "score color limit", "rank color limit"
    ),
    Observed = c(
        nrow(paired_scores), n_distinct(paired_scores$fs_method), nrow(score_long),
        nrow(rank_long), nrow(score_heatmap), nrow(rank_heatmap),
        sum(is.finite(score_long$Difference[score_long$Task == "Alignment"])),
        nrow(global_pair), rank_rho, score_limit, rank_limit
    ),
    Expected = c(203, 29, 812, 812, 667, 667, 58, 29, NA, NA, NA),
    Pass = c(
        nrow(paired_scores) == 203,
        n_distinct(paired_scores$fs_method) == 29,
        nrow(score_long) == 812,
        nrow(rank_long) == 812,
        nrow(score_heatmap) == 667,
        nrow(rank_heatmap) == 667,
        sum(is.finite(score_long$Difference[score_long$Task == "Alignment"])) == 58,
        nrow(global_pair) == 29,
        TRUE, TRUE, TRUE
    )
)
write_tsv(qa, file.path(supplement_root, "Supplemental_Fig_S5_QA.tsv"), na = "NA")

summary_text <- paste0(
    as.character(task_summary$Task), " ", sprintf("%+.3f", task_summary$MeanDifference),
    collapse = "; "
)
writeLines(
    paste(
        "Supplemental Fig. S5 | Sensitivity of feature-selection performance to the integration strategy.",
        "a, Dataset- and task-specific score differences for CellCharter relative to scVI using representative feature-selection settings.",
        "b, Corresponding rank differences, defined as rank(CellCharter) minus rank(scVI); negative values therefore indicate improved ranks under CellCharter.",
        "c, Agreement of mean Core overall ranks between integration strategies. Spatially informed selectors are labelled and the dashed line indicates identity.",
        "d, Distribution of CellCharter-minus-scVI score differences across method-dataset combinations.",
        "Green indicates improvement under CellCharter in both heatmaps, whereas purple indicates deterioration.",
        paste0("Mean score differences: ", summary_text, ".")
    ),
    file.path(supplement_root, "Supplemental_Fig_S5_caption.txt")
)

writeLines(
    c(
        "Core conclusion: integration strategy changes absolute task performance and some method rankings; CellCharter tends to improve clustering and alignment while reducing integration mixing relative to scVI.",
        "Archetype: quantitative grid with paired difference heatmaps and two summary panels.",
        "Panels a-b: dataset-method heterogeneity; panel c: global rank agreement; panel d: direction and magnitude of task shifts.",
        "Review risk: score and rank differences use opposite numeric directions, so the rank color scale is reversed to keep green consistently meaning better under CellCharter."
    ),
    file.path(supplement_root, "Supplemental_Fig_S5_contract.txt")
)

save_reference_supplement(
    figure,
    file.path(supplement_figure_dir, "Supplemental_Fig_S5"),
    width_mm = 208.3,
    height_mm = 265
)

message("Supplemental Figure S5 written to ", supplement_figure_dir)
