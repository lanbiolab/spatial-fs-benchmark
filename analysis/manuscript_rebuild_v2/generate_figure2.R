#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(readr)
  library(scales)
  library(tidyr)
})

source("analysis/manuscript_rebuild_v2/R/theme.R")
source("analysis/manuscript_rebuild_v2/R/data.R")
ensure_figure_dirs()

global <- read_score_table("global_competitive_method_ranks.tsv") %>%
  mutate(
    Method = display_method(fs_method),
    Integrator = display_integrator(integration_method),
    Family = group_label(MethodGroup)
  )

representative <- read_score_table("representative_task_scores.tsv") %>%
  mutate(
    Method = display_method(fs_method),
    Integrator = display_integrator(integration_method),
    Dataset = display_dataset(dataset),
    Family = group_label(MethodGroup)
  )

# Panel a: the global competitive result -------------------------------------
scvi_order <- global %>%
  filter(integration_method == "scvi") %>%
  arrange(GlobalCoreOverallRank) %>%
  pull(Method)

rank_long <- global %>%
  mutate(Method = factor(Method, levels = rev(scvi_order)))

rank_links <- rank_long %>%
  select(Method, Integrator, MeanCoreOverallRank, MethodGroup) %>%
  pivot_wider(names_from = Integrator, values_from = MeanCoreOverallRank)

spatial_bands <- rank_long %>%
  filter(MethodGroup == "spatially_informed") %>%
  distinct(Method) %>%
  mutate(y = as.numeric(Method))

p_a <- ggplot(rank_long, aes(MeanCoreOverallRank, Method)) +
  geom_rect(
    data = spatial_bands,
    aes(xmin = -Inf, xmax = Inf, ymin = y - 0.47, ymax = y + 0.47),
    inherit.aes = FALSE, fill = "#FBEEE9", colour = NA
  ) +
  geom_segment(
    data = rank_links,
    aes(x = scVI, xend = CellCharter, y = Method, yend = Method),
    inherit.aes = FALSE, colour = "#B8B8B8", linewidth = 0.38
  ) +
  geom_point(aes(fill = integration_method), shape = 21, size = 2.55,
             stroke = 0.35, colour = "white") +
  scale_fill_manual(
    values = integration_colours,
    breaks = c("scvi", "cellcharter"),
    labels = c("scVI", "CellCharter"),
    name = "Integration strategy"
  ) +
  scale_x_continuous(
    breaks = c(1, 5, 10, 15, 20, 25, 29),
    limits = c(0.4, 30.2), expand = expansion(mult = 0)
  ) +
  labs(
    x = "Mean CoreOverall rank across datasets",
    y = NULL,
    title = "Global competitive ranks",
    subtitle = "Methods are ordered by scVI; rank 1 is best"
  ) +
  theme_benchmark(6.4) +
  theme(
    axis.text.y = element_text(face = "bold", size = 5.9),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "top",
    legend.justification = "left",
    legend.box.margin = margin(0, 0, 1, 0),
    plot.title = element_text(size = 7.8),
    plot.subtitle = element_text(size = 6.2),
    plot.margin = margin(3, 3, 2, 2, "mm")
  )

# Panel b: task-specific rank structure --------------------------------------
combined_order <- global %>%
  group_by(fs_method, Method, MethodGroup) %>%
  summarise(CombinedRank = mean(GlobalCoreOverallRank), .groups = "drop") %>%
  arrange(CombinedRank)

heat_methods <- union(
  head(combined_order$fs_method, 13),
  combined_order %>% filter(MethodGroup == "spatially_informed") %>% pull(fs_method)
)

heat_order <- combined_order %>%
  filter(fs_method %in% heat_methods) %>%
  arrange(CombinedRank) %>%
  pull(Method)

task_ranks <- global %>%
  filter(fs_method %in% heat_methods) %>%
  select(fs_method, Method, Integrator, MethodGroup,
         MeanIntegrationRank, MeanClusteringRank, MeanAlignmentRank) %>%
  pivot_longer(
    starts_with("Mean"), names_to = "Task", values_to = "MeanRank"
  ) %>%
  mutate(
    Task = recode(
      Task,
      MeanIntegrationRank = "Integration",
      MeanClusteringRank = "Clustering",
      MeanAlignmentRank = "Alignment"
    ),
    Task = factor(Task, levels = c("Integration", "Clustering", "Alignment")),
    Integrator = factor(Integrator, levels = c("scVI", "CellCharter")),
    Method = factor(Method, levels = rev(heat_order))
  )

p_b <- ggplot(task_ranks, aes(Task, Method, fill = MeanRank)) +
  geom_tile(colour = "white", linewidth = 0.42) +
  geom_point(
    data = task_ranks %>% group_by(Integrator, Task) %>% slice_min(MeanRank, n = 1, with_ties = TRUE),
    shape = 21, size = 1.5, stroke = 0.32, fill = NA, colour = "white"
  ) +
  facet_grid(. ~ Integrator, scales = "free_x", space = "free_x") +
  scale_fill_gradient(
    low = "#244C68", high = "#F2EEE6", limits = c(1, 29),
    name = "Mean rank", breaks = c(1, 10, 20, 29)
  ) +
  scale_x_discrete(position = "top") +
  labs(
    x = NULL, y = NULL,
    title = "Task-specific performance",
    subtitle = "White rings mark the best mean rank within each task"
  ) +
  theme_benchmark(6.1) +
  theme(
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.text.x = element_text(angle = 42, hjust = 0, face = "bold", size = 5.8),
    axis.text.y = element_text(face = "bold", size = 5.6),
    strip.text = element_text(size = 6.5, face = "bold"),
    legend.position = "bottom",
    legend.key.width = unit(10, "mm"),
    plot.title = element_text(size = 7.5),
    plot.subtitle = element_text(size = 5.9),
    plot.margin = margin(3, 2, 2, 3, "mm")
  )

# Panel c: controls and label-informed oracle --------------------------------
competitive_dataset <- representative %>%
  filter(MethodGroup %in% c("expression_driven", "spatially_informed")) %>%
  group_by(dataset, Dataset, integration_method, Integrator) %>%
  summarise(CoreOverallMean = median(CoreOverallMean, na.rm = TRUE), .groups = "drop") %>%
  mutate(Reference = "Competitive median")

reference_scores <- representative %>%
  filter(fs_method %in% c("all_features", "random", "TFs", "scsegindex", "wilcoxon")) %>%
  mutate(
    Reference = recode(
      fs_method,
      all_features = "All features",
      random = "Random N=500",
      TFs = "TF set",
      scsegindex = "scSEGIndex",
      wilcoxon = "Wilcoxon oracle"
    )
  ) %>%
  select(dataset, Dataset, integration_method, Integrator, CoreOverallMean, Reference)

control_scores <- bind_rows(competitive_dataset, reference_scores) %>%
  mutate(
    Reference = factor(
      Reference,
      levels = c("Competitive median", "All features", "Random N=500", "TF set",
                 "scSEGIndex", "Wilcoxon oracle")
    ),
    Integrator = factor(Integrator, levels = c("scVI", "CellCharter"))
  )

p_c <- ggplot(control_scores, aes(CoreOverallMean, Reference)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, linewidth = 0.35,
               fill = "#E9E9E6", colour = BENCH_COLOURS$ink) +
  geom_point(
    aes(colour = Reference), position = position_jitter(width = 0, height = 0.10),
    size = 1.15, alpha = 0.78, show.legend = FALSE
  ) +
  facet_grid(Integrator ~ .) +
  scale_colour_manual(values = c(
    "Competitive median" = BENCH_COLOURS$core,
    "All features" = "#777777", "Random N=500" = "#B0B0B0",
    "TF set" = "#8F8F8F", "scSEGIndex" = "#696969",
    "Wilcoxon oracle" = BENCH_COLOURS$oracle
  )) +
  scale_x_continuous(limits = c(0.28, 0.86), breaks = c(0.3, 0.5, 0.7, 0.85)) +
  labs(
    x = "CoreOverall score", y = NULL,
    title = "Controls and label-informed oracle",
    subtitle = "Each point is one dataset"
  ) +
  theme_benchmark(5.8) +
  theme(
    axis.text.y = element_text(size = 4.9),
    strip.text = element_text(size = 6.2, face = "bold"),
    plot.title = element_text(size = 7.0),
    plot.subtitle = element_text(size = 5.6),
    plot.margin = margin(3, 2, 3, 3, "mm")
  )

# Panel d: rank association across integration strategies --------------------
rank_compare <- global %>%
  select(fs_method, Method, MethodGroup, Integrator, GlobalCoreOverallRank) %>%
  pivot_wider(names_from = Integrator, values_from = GlobalCoreOverallRank)

rho <- cor(rank_compare$scVI, rank_compare$CellCharter, method = "spearman")
rho_p <- cor.test(rank_compare$scVI, rank_compare$CellCharter, method = "spearman", exact = FALSE)$p.value

label_methods <- c("triku", "seurat_vst", "somde", "spatialde", "sparkx")

p_d <- ggplot(rank_compare, aes(scVI, CellCharter)) +
  geom_abline(slope = 1, intercept = 0, linetype = "22", linewidth = 0.4,
              colour = "#A5A5A5") +
  geom_point(aes(fill = MethodGroup), shape = 21, size = 2.25,
             stroke = 0.35, colour = "white") +
  ggrepel::geom_text_repel(
    data = rank_compare %>% filter(fs_method %in% label_methods),
    aes(label = Method), family = BENCH_FONT, size = 1.8,
    colour = BENCH_COLOURS$ink, max.overlaps = Inf,
    box.padding = 0.22, point.padding = 0.12, min.segment.length = 0,
    segment.size = 0.2, seed = 17, show.legend = FALSE
  ) +
  scale_fill_manual(
    values = method_group_colours,
    breaks = c("expression_driven", "spatially_informed"),
    labels = c("Expression-based", "Spatially informed"),
    name = "Method family"
  ) +
  scale_x_continuous(breaks = c(1, 10, 20, 29), limits = c(0.3, 30)) +
  scale_y_continuous(breaks = c(1, 10, 20, 29), limits = c(0.3, 30)) +
  coord_equal(clip = "off") +
  labs(
    x = "scVI rank", y = "CellCharter rank",
    title = "Rank association",
    subtitle = sprintf("rho = %.2f; P = %.3g\nspatial methods in coral", rho, rho_p)
  ) +
  theme_benchmark(5.8) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 7.0),
    plot.subtitle = element_text(size = 5.6),
    plot.margin = margin(3, 3, 3, 2, "mm")
  )

right_column <- p_b / (p_c + p_d + plot_layout(widths = c(1.18, 0.82))) +
  plot_layout(heights = c(1.18, 0.82))

figure <- p_a + right_column +
  plot_layout(widths = c(0.46, 0.54)) +
  plot_annotation(tag_levels = "a") &
  panel_label_theme

write_source_data(
  rank_long %>%
    transmute(fs_method, Method, MethodGroup, integration_method, Integrator,
              MeanCoreOverallScore, MeanCoreOverallRank, GlobalCoreOverallRank),
  "figure2_global_competitive_ranks.tsv"
)
write_source_data(task_ranks, "figure2_task_specific_ranks.tsv")
write_source_data(control_scores, "figure2_control_scores.tsv")
write_source_data(rank_compare, "figure2_cross_integrator_ranks.tsv")

save_figure_set(
  figure,
  file.path(FIGURE_ROOT, "figures", "figure2_global_benchmark"),
  width_mm = 183,
  height_mm = 178
)

message("Figure 2 written to ", file.path(FIGURE_ROOT, "figures"))
