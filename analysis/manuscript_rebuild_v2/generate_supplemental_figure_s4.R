#!/usr/bin/env Rscript

source(file.path("analysis", "manuscript_rebuild_v2", "R", "reference_supplement.R"))

score_root <- file.path("results", "spatial_svg_rebuild_v1", "frozen_scores")
setting_scores <- read_tsv(file.path(score_root, "setting_task_scores.tsv"), show_col_types = FALSE)

n_levels <- c(100, 200, 500, 1000, 2000, 5000, 10000)
dataset_levels <- unname(supplement_dataset_labels)
integrator_levels <- unname(supplement_integrator_labels)

competitive <- setting_scores |>
    filter(.data$MethodGroup %in% c("expression_driven", "spatially_informed")) |>
    mutate(
        NFeatures = as.integer(.data$n_features),
        Dataset = factor(display_supplement_dataset(.data$dataset), levels = dataset_levels),
        Integrator = factor(unname(supplement_integrator_labels[.data$integration_method]), levels = integrator_levels),
        Method = display_supplement_method(.data$fs_method)
    ) |>
    filter(.data$NFeatures %in% n_levels, is.finite(.data$CoreOverallMean))

global_summary <- competitive |>
    group_by(.data$Integrator, .data$NFeatures) |>
    summarise(
        Mean = mean(.data$CoreOverallMean),
        SD = sd(.data$CoreOverallMean),
        N = n(),
        Lower = .data$Mean - .data$SD,
        Upper = .data$Mean + .data$SD,
        .groups = "drop"
    )

plot_a <- ggplot(competitive, aes(x = .data$NFeatures, y = .data$CoreOverallMean)) +
    geom_point(
        colour = "grey55", size = 0.52, alpha = 0.22,
        position = position_jitter(width = 0.035, height = 0, seed = 101)
    ) +
    geom_linerange(
        data = global_summary,
        aes(x = .data$NFeatures, ymin = .data$Lower, ymax = .data$Upper),
        inherit.aes = FALSE, linewidth = 0.42, colour = "grey25"
    ) +
    geom_line(
        data = global_summary,
        aes(x = .data$NFeatures, y = .data$Mean, group = 1),
        inherit.aes = FALSE, linewidth = 0.52, colour = "grey20"
    ) +
    geom_point(
        data = global_summary,
        aes(x = .data$NFeatures, y = .data$Mean),
        inherit.aes = FALSE,
        shape = 23, size = 2.25, stroke = 0.35,
        fill = "#D95F59", colour = "white"
    ) +
    facet_grid(. ~ .data$Integrator) +
    scale_x_log10(breaks = n_levels, labels = label_number(big.mark = ",")) +
    scale_y_continuous(breaks = pretty_breaks(5)) +
    labs(x = "Number of selected features", y = "Core overall score", tag = "a") +
    reference_supplement_theme() +
    theme(
        axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 5.0),
        panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.25),
        panel.grid.major.x = element_blank()
    )

baseline <- competitive |>
    filter(.data$NFeatures == 2000) |>
    select(
        "dataset", "fs_method", "integration_method",
        Baseline = "CoreOverallMean"
    )

relative <- competitive |>
    left_join(baseline, by = c("dataset", "fs_method", "integration_method")) |>
    mutate(Delta = .data$CoreOverallMean - .data$Baseline)

dataset_delta <- relative |>
    group_by(.data$Dataset, .data$Integrator, .data$NFeatures) |>
    summarise(
        MeanDelta = mean(.data$Delta),
        SDDelta = sd(.data$Delta),
        NMethods = n(),
        .groups = "drop"
    ) |>
    mutate(
        Dataset = factor(.data$Dataset, levels = rev(dataset_levels)),
        NLabel = factor(.data$NFeatures, levels = n_levels, labels = label_number(big.mark = ",")(n_levels))
    )

method_order <- competitive |>
    filter(.data$NFeatures == 2000) |>
    group_by(.data$fs_method, .data$Method) |>
    summarise(ReferenceScore = mean(.data$CoreOverallMean), .groups = "drop") |>
    arrange(.data$ReferenceScore, .data$Method) |>
    pull(.data$Method)

method_delta <- relative |>
    group_by(.data$fs_method, .data$Method, .data$Integrator, .data$NFeatures) |>
    summarise(
        MeanDelta = mean(.data$Delta),
        SDDelta = sd(.data$Delta),
        NDatasets = n(),
        .groups = "drop"
    ) |>
    mutate(
        Method = factor(.data$Method, levels = method_order),
        NLabel = factor(.data$NFeatures, levels = n_levels, labels = label_number(big.mark = ",")(n_levels))
    )

delta_limit <- max(abs(c(dataset_delta$MeanDelta, method_delta$MeanDelta)), na.rm = TRUE)
delta_breaks <- c(-delta_limit, 0, delta_limit)

delta_scale <- scale_fill_gradient2(
    low = "#762A83", mid = "#F7F7F7", high = "#1B7837",
    midpoint = 0, limits = c(-delta_limit, delta_limit),
    breaks = delta_breaks, labels = label_number(accuracy = 0.01),
    oob = squish, na.value = "grey85",
    name = "Score difference\nfrom N=2,000",
    guide = guide_colourbar(
        title.position = "top", barwidth = unit(2.0, "cm"),
        barheight = unit(0.20, "cm"), order = 1
    )
)

plot_b <- ggplot(dataset_delta, aes(x = .data$NLabel, y = .data$Dataset, fill = .data$MeanDelta)) +
    geom_tile(colour = "white", linewidth = 0.45) +
    facet_grid(. ~ .data$Integrator) +
    delta_scale +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(x = "Number of selected features", y = NULL, tag = "b") +
    reference_supplement_theme() +
    theme(
        axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 4.9),
        axis.text.y = element_text(face = "bold", size = 5.0),
        axis.ticks = element_blank(), axis.line = element_blank(), panel.grid = element_blank(),
        legend.position = "bottom", legend.title.position = "top", legend.box.just = "left"
    )

plot_c <- ggplot(method_delta, aes(x = .data$NLabel, y = .data$Method, fill = .data$MeanDelta)) +
    geom_tile(colour = "white", linewidth = 0.38) +
    facet_grid(. ~ .data$Integrator) +
    delta_scale +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0), drop = FALSE) +
    labs(x = "Number of selected features", y = NULL, tag = "c") +
    reference_supplement_theme() +
    theme(
        axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 4.9),
        axis.text.y = element_text(face = "bold", colour = "black", size = 4.7),
        axis.ticks = element_blank(), axis.line = element_blank(), panel.grid = element_blank(),
        legend.position = "bottom", legend.title.position = "top", legend.box.just = "left"
    )

figure <- plot_a / plot_b / plot_c +
    plot_layout(heights = c(0.80, 0.75, 2.25), guides = "collect") &
    theme(
        legend.position = "bottom", legend.title.position = "top",
        legend.box.just = "left", legend.margin = margin(0, 0, 0, 0),
        legend.text = element_text(size = 5.0),
        legend.title = element_text(size = 5.7, face = "bold")
    )

write_tsv(competitive, file.path(supplement_source_dir, "Supplemental_Fig_S4a_feature_number_scores.tsv"), na = "NA")
write_tsv(dataset_delta, file.path(supplement_source_dir, "Supplemental_Fig_S4b_dataset_differences.tsv"), na = "NA")
write_tsv(method_delta, file.path(supplement_source_dir, "Supplemental_Fig_S4c_method_differences.tsv"), na = "NA")

best_n <- global_summary |>
    group_by(.data$Integrator) |>
    slice_max(.data$Mean, n = 1, with_ties = FALSE) |>
    ungroup()

qa <- tibble(
    Check = c(
        "competitive setting rows", "competitive methods", "feature-number settings",
        "datasets", "integrators", "dataset heatmap cells", "method heatmap cells",
        "maximum absolute plotted difference"
    ),
    Observed = c(
        nrow(competitive), n_distinct(competitive$fs_method), n_distinct(competitive$NFeatures),
        n_distinct(competitive$dataset), n_distinct(competitive$integration_method),
        nrow(dataset_delta), nrow(method_delta), delta_limit
    ),
    Expected = c(2842, 29, 7, 7, 2, 98, 406, NA),
    Pass = c(
        nrow(competitive) == 2842,
        n_distinct(competitive$fs_method) == 29,
        setequal(unique(competitive$NFeatures), n_levels),
        n_distinct(competitive$dataset) == 7,
        n_distinct(competitive$integration_method) == 2,
        nrow(dataset_delta) == 98,
        nrow(method_delta) == 406,
        TRUE
    )
)
write_tsv(qa, file.path(supplement_root, "Supplemental_Fig_S4_QA.tsv"), na = "NA")

best_text <- paste0(
    as.character(best_n$Integrator), " peaked at N=", label_number(big.mark = ",")(best_n$NFeatures),
    " (mean Core overall=", sprintf("%.3f", best_n$Mean), ")",
    collapse = "; "
)
writeLines(
    paste(
        "Supplemental Fig. S4 | Sensitivity of benchmark performance to the number of selected features.",
        "a, Core overall scores across 29 competitive expression-based and spatially informed selectors. Small points denote method-dataset combinations; diamonds and vertical bars show the mean and s.d.",
        "b, Dataset-level mean score difference relative to N=2,000, averaged over methods.",
        "c, Method-level mean score difference relative to N=2,000, averaged over datasets.",
        "Purple and green respectively indicate lower and higher scores than the N=2,000 reference.",
        best_text
    ),
    file.path(supplement_root, "Supplemental_Fig_S4_caption.txt")
)

writeLines(
    c(
        "Core conclusion: performance generally forms a broad plateau between 500 and 2,000 features, whereas very small and very large feature sets can reduce performance in a dataset- and method-dependent manner.",
        "Archetype: quantitative grid with a global trend and two heterogeneity heatmaps.",
        "Panel a: overall trend; panel b: dataset heterogeneity; panel c: method heterogeneity.",
        "Review risk: only competitive expression-driven and spatially informed methods are included; controls, all features and the label-informed oracle are excluded."
    ),
    file.path(supplement_root, "Supplemental_Fig_S4_contract.txt")
)

save_reference_supplement(
    figure,
    file.path(supplement_figure_dir, "Supplemental_Fig_S4"),
    width_mm = 208.3,
    height_mm = 235
)

message("Supplemental Figure S4 written to ", supplement_figure_dir)
