#!/usr/bin/env Rscript

source(file.path("analysis", "manuscript_rebuild_v2", "R", "reference_supplement.R"))

score_root <- file.path("results", "spatial_svg_rebuild_v1", "frozen_scores")

seed_scores <- read_tsv(file.path(score_root, "seed_task_scores.tsv"), show_col_types = FALSE)
representative <- read_tsv(file.path(score_root, "representative_settings.tsv"), show_col_types = FALSE)

dataset_levels <- unname(supplement_dataset_labels)
integrator_levels <- unname(supplement_integrator_labels)
task_levels <- names(supplement_task_palette)

representative_keys <- representative |>
    distinct(.data$dataset, .data$fs_method, .data$n_features)
seed_representative <- seed_scores |>
    inner_join(representative_keys, by = c("dataset", "fs_method", "n_features")) |>
    mutate(
        Dataset = factor(display_supplement_dataset(.data$dataset), levels = dataset_levels),
        Integrator = factor(unname(supplement_integrator_labels[.data$integration_method]), levels = integrator_levels)
    )

seed_sd <- seed_representative |>
    select(
        "dataset", "fs_method", "n_features", "integration_method", "Dataset", "Integrator",
        "random_seed", "Integration", "Clustering", "Alignment", "CoreOverall"
    ) |>
    pivot_longer(
        cols = c("Integration", "Clustering", "Alignment", "CoreOverall"),
        names_to = "Task", values_to = "Score"
    ) |>
    mutate(Task = recode(.data$Task, CoreOverall = "Core overall")) |>
    group_by(
        .data$dataset, .data$fs_method, .data$n_features,
        .data$integration_method, .data$Dataset, .data$Integrator, .data$Task
    ) |>
    summarise(SeedSD = sd(.data$Score, na.rm = TRUE), NSeeds = sum(is.finite(.data$Score)), .groups = "drop") |>
    filter(.data$NSeeds == 3, is.finite(.data$SeedSD)) |>
    mutate(Task = factor(.data$Task, levels = task_levels))

dataset_seed_sd <- seed_sd |>
    group_by(.data$Dataset, .data$Integrator, .data$Task) |>
    summarise(
        MedianSeedSD = median(.data$SeedSD),
        Q90SeedSD = quantile(.data$SeedSD, 0.9),
        NMethods = n(),
        .groups = "drop"
    ) |>
    complete(
        Dataset = factor(dataset_levels, levels = dataset_levels),
        Integrator = factor(integrator_levels, levels = integrator_levels),
        Task = factor(task_levels, levels = task_levels)
    ) |>
    mutate(
        Dataset = factor(.data$Dataset, levels = rev(dataset_levels)),
        Integrator = factor(.data$Integrator, levels = integrator_levels),
        Task = factor(.data$Task, levels = task_levels),
        ValueLabel = if_else(is.na(.data$MedianSeedSD), "NA", sprintf("%.3f", .data$MedianSeedSD))
    )

seed_sd_limit <- max(dataset_seed_sd$MedianSeedSD, na.rm = TRUE)

plot_a <- ggplot(dataset_seed_sd, aes(x = .data$Task, y = .data$Dataset, fill = .data$MedianSeedSD)) +
    geom_tile(colour = "white", linewidth = 0.55) +
    geom_text(
        aes(label = .data$ValueLabel),
        family = "Arial", fontface = "bold", size = 2.05,
        colour = ifelse(
            is.na(dataset_seed_sd$MedianSeedSD),
            "grey25",
            ifelse(dataset_seed_sd$MedianSeedSD > seed_sd_limit * 0.70, "grey20", "white")
        )
    ) +
    facet_grid(. ~ .data$Integrator) +
    scale_fill_viridis_c(
        option = "C", limits = c(0, seed_sd_limit),
        breaks = c(0, seed_sd_limit / 2, seed_sd_limit),
        labels = label_number(accuracy = 0.001), na.value = "grey85",
        name = "Median seed s.d.",
        guide = guide_colourbar(
            title.position = "top", barwidth = unit(1.4, "cm"),
            barheight = unit(0.20, "cm"), order = 1
        )
    ) +
    scale_x_discrete(
        expand = c(0, 0),
        labels = c("Core overall" = "Overall"),
        guide = guide_axis(angle = 45)
    ) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(x = NULL, y = NULL, tag = "a") +
    reference_supplement_theme() +
    theme(
        axis.text.x = element_text(face = "bold", size = 4.5, hjust = 1),
        axis.text.y = element_text(face = "bold", size = 5.0),
        axis.ticks = element_blank(), axis.line = element_blank(),
        panel.grid = element_blank(), legend.position = "bottom",
        legend.title.position = "top", legend.box.just = "left"
    )

core_wide <- seed_representative |>
    select("dataset", "fs_method", "n_features", "integration_method", "Dataset", "Integrator", "random_seed", "CoreOverall") |>
    pivot_wider(names_from = "random_seed", values_from = "CoreOverall", names_prefix = "Seed")

seed_pairs <- bind_rows(
    core_wide |> transmute(across(c("dataset", "fs_method", "n_features", "integration_method", "Dataset", "Integrator")),
                            Comparison = "Seed 0 vs 1", Seed0 = .data$Seed0, OtherSeed = .data$Seed1),
    core_wide |> transmute(across(c("dataset", "fs_method", "n_features", "integration_method", "Dataset", "Integrator")),
                            Comparison = "Seed 0 vs 2", Seed0 = .data$Seed0, OtherSeed = .data$Seed2)
) |>
    filter(is.finite(.data$Seed0), is.finite(.data$OtherSeed)) |>
    mutate(Comparison = factor(.data$Comparison, levels = c("Seed 0 vs 1", "Seed 0 vs 2")))

seed_cor <- seed_pairs |>
    group_by(.data$Integrator, .data$Comparison) |>
    summarise(
        rho = cor(.data$Seed0, .data$OtherSeed, method = "spearman"),
        x = min(.data$Seed0), y = max(.data$OtherSeed),
        label = sprintf("rho = %.2f", .data$rho),
        .groups = "drop"
    )

plot_b <- ggplot(seed_pairs, aes(x = .data$Seed0, y = .data$OtherSeed)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey55", linewidth = 0.4, linetype = "dashed") +
    geom_point(
        aes(fill = .data$Dataset), shape = 21, size = 1.35,
        stroke = 0.24, colour = "grey30", alpha = 0.72
    ) +
    geom_text(
        data = seed_cor,
        aes(x = .data$x, y = .data$y, label = .data$label),
        inherit.aes = FALSE, hjust = 0, vjust = 1,
        family = "Arial", fontface = "bold", size = 2.0
    ) +
    facet_grid(.data$Integrator ~ .data$Comparison) +
    scale_fill_manual(values = supplement_dataset_palette, name = "Dataset") +
    scale_x_continuous(breaks = pretty_breaks(4)) +
    scale_y_continuous(breaks = pretty_breaks(4)) +
    coord_equal() +
    labs(x = "Core overall, seed 0", y = "Core overall, comparison seed", tag = "b") +
    reference_supplement_theme() +
    theme(
        legend.position = "bottom", legend.title.position = "top",
        legend.box.just = "left", legend.key.size = unit(0.28, "cm"),
        panel.grid.major = element_line(colour = "grey92", linewidth = 0.25)
    )

plot_c <- ggplot(seed_sd, aes(x = .data$Task, y = .data$SeedSD, fill = .data$Task)) +
    geom_violin(width = 0.82, scale = "width", colour = NA, alpha = 0.32, trim = TRUE) +
    geom_boxplot(
        width = 0.24, outlier.shape = NA, fill = "white",
        colour = "grey25", linewidth = 0.35
    ) +
    geom_point(
        aes(colour = .data$Task), size = 0.62, alpha = 0.35,
        position = position_jitter(width = 0.14, height = 0, seed = 47)
    ) +
    facet_grid(. ~ .data$Integrator) +
    scale_fill_manual(values = supplement_task_palette, guide = "none") +
    scale_colour_manual(values = supplement_task_palette, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.06)), breaks = pretty_breaks(5)) +
    labs(x = NULL, y = "Score s.d. across three seeds", tag = "c") +
    reference_supplement_theme() +
    theme(
        axis.text.x = element_text(face = "bold", size = 5.2),
        panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.25),
        panel.grid.major.x = element_blank()
    )

figure <- (plot_a + plot_b + plot_layout(widths = c(0.42, 0.58), guides = "collect")) /
    plot_c +
    plot_layout(heights = c(1.15, 0.85), guides = "collect") &
    theme(
        legend.position = "bottom", legend.title.position = "top",
        legend.box = "horizontal", legend.box.just = "left",
        legend.margin = margin(0, 0, 0, 0),
        legend.spacing.x = unit(0.12, "cm"),
        legend.text = element_text(size = 5.0),
        legend.title = element_text(size = 5.7, face = "bold")
    )

write_tsv(dataset_seed_sd, file.path(supplement_source_dir, "Supplemental_Fig_S3a_dataset_seed_variability.tsv"), na = "NA")
write_tsv(seed_pairs, file.path(supplement_source_dir, "Supplemental_Fig_S3b_seed_consistency.tsv"), na = "NA")
write_tsv(seed_sd, file.path(supplement_source_dir, "Supplemental_Fig_S3c_seed_variability.tsv"), na = "NA")

qa <- tibble(
    Check = c(
        "representative seed rows", "seed comparisons", "dataset-task-integrator cells",
        "applicable dataset-task-integrator cells",
        "three-seed s.d. groups", "median CoreOverall seed s.d.", "maximum CoreOverall seed s.d."
    ),
    Observed = c(
        nrow(seed_representative), nrow(seed_pairs), nrow(dataset_seed_sd),
        sum(is.finite(dataset_seed_sd$MedianSeedSD)),
        nrow(seed_sd),
        median(seed_sd$SeedSD[seed_sd$Task == "Core overall"]),
        max(seed_sd$SeedSD[seed_sd$Task == "Core overall"])
    ),
    Expected = c(1416, 944, 56, 46, 1552, NA, NA),
    Pass = c(
        nrow(seed_representative) == 1416,
        nrow(seed_pairs) == 944,
        nrow(dataset_seed_sd) == 56,
        sum(is.finite(dataset_seed_sd$MedianSeedSD)) == 46,
        nrow(seed_sd) == 1552,
        TRUE, TRUE
    )
)
write_tsv(qa, file.path(supplement_root, "Supplemental_Fig_S3_QA.tsv"), na = "NA")

writeLines(
    paste(
        "Supplemental Fig. S3 | Random-seed stability of the rebuilt benchmark.",
        "a, Median score s.d. across three seeds for each dataset, task and integration strategy; gray cells denote tasks that were not applicable.",
        "b, Agreement of representative-setting Core overall scores from seed 0 with seeds 1 and 2. Points denote feature-selection method-dataset combinations and dashed lines indicate identity.",
        "c, Distribution of score s.d. across three seeds for each task. Violin envelopes show the distribution, boxes show the median and interquartile range, and points show method-dataset combinations.",
        "Alignment stability is shown only for DLPFC and Mouse Brain, where explicit slice-pair evaluation was available."
    ),
    file.path(supplement_root, "Supplemental_Fig_S3_caption.txt")
)

writeLines(
    c(
        "Core conclusion: canonical three-seed scores are generally stable, with the remaining variability concentrated in specific dataset-task combinations.",
        "Archetype: quantitative grid.",
        "Panel a: dataset- and task-specific variability; panel b: seed-to-seed agreement; panel c: global variability distribution.",
        "Review risk: inapplicable alignment cells are shown as NA rather than zero."
    ),
    file.path(supplement_root, "Supplemental_Fig_S3_contract.txt")
)

save_reference_supplement(
    figure,
    file.path(supplement_figure_dir, "Supplemental_Fig_S3"),
    width_mm = 208.3,
    height_mm = 175
)

message("Supplemental Figure S3 written to ", supplement_figure_dir)
