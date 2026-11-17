#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(scales)
  library(tidyr)
})

source("analysis/manuscript_rebuild_v2/R/theme.R")
source("analysis/manuscript_rebuild_v2/R/data.R")
ensure_figure_dirs()

feature_levels <- c(100, 200, 500, 1000, 2000, 5000, 10000)

settings <- read_score_table("setting_task_scores.tsv") %>%
  filter(
    MethodGroup %in% c("expression_driven", "spatially_informed"),
    n_features != "all"
  ) %>%
  mutate(
    NFeatures = as.numeric(n_features),
    Integrator = display_integrator(integration_method),
    Dataset = display_dataset(dataset),
    Method = display_method(fs_method),
    Integrator = factor(Integrator, levels = c("scVI", "CellCharter"))
  )

# The dataset, rather than each method-setting row, is the replication unit.
dataset_budget <- settings %>%
  group_by(dataset, Dataset, integration_method, Integrator, NFeatures) %>%
  summarise(
    CoreOverall = mean(CoreOverallMean, na.rm = TRUE),
    Integration = mean(IntegrationMean, na.rm = TRUE),
    Clustering = mean(ClusteringMean, na.rm = TRUE),
    .groups = "drop"
  )

global_budget <- dataset_budget %>%
  group_by(integration_method, Integrator, NFeatures) %>%
  summarise(
    Mean = mean(CoreOverall),
    SD = sd(CoreOverall),
    NDatasets = n(),
    SE = SD / sqrt(NDatasets),
    .groups = "drop"
  )

peaks <- global_budget %>%
  group_by(Integrator) %>%
  slice_max(Mean, n = 1, with_ties = FALSE) %>%
  ungroup()

p_a <- ggplot(global_budget, aes(NFeatures, Mean, colour = integration_method,
                                  fill = integration_method)) +
  geom_ribbon(aes(ymin = Mean - SE, ymax = Mean + SE), alpha = 0.13,
              colour = NA, show.legend = FALSE) +
  geom_line(linewidth = 0.8) +
  geom_point(shape = 21, size = 2.7, stroke = 0.45, colour = "white") +
  geom_point(
    data = peaks, shape = 21, size = 4.0, stroke = 0.65,
    fill = "white", aes(colour = integration_method), show.legend = FALSE
  ) +
  geom_point(
    data = peaks, shape = 21, size = 2.4, stroke = 0,
    aes(fill = integration_method), colour = NA, show.legend = FALSE
  ) +
  geom_text(
    data = peaks,
    aes(label = paste0("peak: N=", comma(NFeatures), "\n", sprintf("%.3f", Mean))),
    family = BENCH_FONT, size = 2.1, fontface = "bold",
    nudge_y = 0.018, show.legend = FALSE
  ) +
  scale_colour_manual(values = integration_colours, breaks = c("scvi", "cellcharter"),
                      labels = c("scVI", "CellCharter"),
                      name = "Integration strategy") +
  scale_fill_manual(values = integration_colours, breaks = c("scvi", "cellcharter"),
                    labels = c("scVI", "CellCharter"),
                    name = "Integration strategy") +
  scale_x_log10(breaks = feature_levels, labels = scales::comma) +
  scale_y_continuous(expand = expansion(mult = c(0.07, 0.17))) +
  labs(
    x = "Number of selected features", y = "Mean CoreOverall score",
    title = "Intermediate feature budgets perform best on average",
    subtitle = "Mean ± s.e.m. across seven datasets; methods are averaged within each dataset"
  ) +
  theme_benchmark(6.5) +
  theme(
    legend.position = "top", legend.justification = "left",
    plot.title = element_text(size = 8.0),
    plot.subtitle = element_text(size = 6.1),
    plot.margin = margin(3, 4, 2, 3, "mm")
  )

# Task decomposition ---------------------------------------------------------
task_budget <- dataset_budget %>%
  select(dataset, Dataset, integration_method, Integrator, NFeatures,
         Integration, Clustering) %>%
  pivot_longer(c(Integration, Clustering), names_to = "Task", values_to = "Score") %>%
  group_by(integration_method, Integrator, NFeatures, Task) %>%
  summarise(
    Mean = mean(Score), SD = sd(Score), NDatasets = n(),
    SE = SD / sqrt(NDatasets), .groups = "drop"
  ) %>%
  mutate(Task = factor(Task, levels = c("Integration", "Clustering")))

p_b <- ggplot(task_budget, aes(NFeatures, Mean, colour = integration_method,
                                fill = integration_method)) +
  geom_ribbon(aes(ymin = Mean - SE, ymax = Mean + SE), alpha = 0.12,
              colour = NA, show.legend = FALSE) +
  geom_line(linewidth = 0.65) +
  geom_point(shape = 21, size = 1.8, stroke = 0.3, colour = "white") +
  facet_grid(. ~ Task) +
  scale_colour_manual(values = integration_colours, guide = "none") +
  scale_fill_manual(values = integration_colours, guide = "none") +
  scale_x_log10(breaks = c(100, 500, 2000, 10000), labels = scales::comma) +
  labs(
    x = "Number of selected features", y = "Mean task score",
    title = "Integration and clustering respond differently"
  ) +
  theme_benchmark(5.9) +
  theme(
    strip.text = element_text(size = 6.4, face = "bold"),
    plot.title = element_text(size = 7.1),
    plot.margin = margin(3, 2, 3, 3, "mm")
  )

# Dataset-specific deviations from each dataset's budget mean ----------------
dataset_heat <- dataset_budget %>%
  group_by(dataset, Dataset, integration_method, Integrator) %>%
  mutate(Deviation = CoreOverall - mean(CoreOverall)) %>%
  ungroup() %>%
  mutate(
    NLabel = factor(scales::comma(NFeatures), levels = scales::comma(feature_levels)),
    Dataset = factor(
      Dataset,
      levels = rev(c("DLPFC", "Mouse brain sections", "STOmics 0212", "STOmics 0218",
                     "STOmics 0224", "E8.5 embryo", "E9.5 embryo"))
    )
  )

heat_limit <- max(abs(dataset_heat$Deviation), na.rm = TRUE)

p_c <- ggplot(dataset_heat, aes(NLabel, Dataset, fill = Deviation)) +
  geom_tile(colour = "white", linewidth = 0.45) +
  facet_grid(. ~ Integrator) +
  scale_fill_gradient2(
    low = "#2A6F97", mid = "#F6F4EE", high = "#D87545",
    midpoint = 0, limits = c(-heat_limit, heat_limit),
    name = "Deviation\nfrom dataset mean",
    labels = label_number(accuracy = 0.01)
  ) +
  scale_x_discrete(position = "top") +
  labs(
    x = NULL, y = NULL,
    title = "The budget effect varies by dataset"
  ) +
  theme_benchmark(5.8) +
  theme(
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.text.x = element_text(angle = 48, hjust = 0, size = 5.2),
    axis.text.y = element_text(face = "bold", size = 5.4),
    strip.text = element_text(size = 6.2, face = "bold"),
    legend.position = "bottom", legend.key.width = unit(10, "mm"),
    plot.title = element_text(size = 7.1),
    plot.margin = margin(3, 2, 3, 3, "mm")
  )

# Distribution of method-specific optimal feature budgets --------------------
method_optima <- settings %>%
  group_by(fs_method, Method, MethodGroup, integration_method, Integrator, NFeatures) %>%
  summarise(MeanCoreOverall = mean(CoreOverallMean, na.rm = TRUE), .groups = "drop") %>%
  group_by(fs_method, Method, MethodGroup, integration_method, Integrator) %>%
  arrange(NFeatures) %>%
  slice_max(MeanCoreOverall, n = 1, with_ties = FALSE) %>%
  ungroup()

optimum_counts <- method_optima %>%
  count(Integrator, integration_method, MethodGroup, NFeatures, name = "NMethods") %>%
  group_by(Integrator, integration_method, MethodGroup) %>%
  complete(NFeatures = feature_levels, fill = list(NMethods = 0)) %>%
  ungroup() %>%
  mutate(
    NLabel = factor(scales::comma(NFeatures), levels = scales::comma(feature_levels)),
    Family = group_label(MethodGroup)
  )

p_d <- ggplot(optimum_counts, aes(NLabel, NMethods, fill = MethodGroup)) +
  geom_col(width = 0.74, colour = "white", linewidth = 0.2) +
  facet_grid(Integrator ~ .) +
  scale_fill_manual(
    values = method_group_colours,
    breaks = c("expression_driven", "spatially_informed"),
    labels = c("Expression-based", "Spatially informed"),
    name = "Method family"
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(4), expand = expansion(mult = c(0, 0.1))) +
  labs(
    x = "Method-specific optimal N", y = "Number of methods",
    title = "Method-specific optima",
    subtitle = "Spatial methods in coral"
  ) +
  theme_benchmark(5.8) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 5.0),
    strip.text = element_text(size = 6.2, face = "bold"),
    legend.position = "none",
    plot.title = element_text(size = 7.1),
    plot.subtitle = element_text(size = 5.6),
    plot.margin = margin(3, 3, 3, 2, "mm")
  )

figure <- p_a / (p_b + p_c + p_d + plot_layout(widths = c(0.32, 0.42, 0.26))) +
  plot_layout(heights = c(0.45, 0.55)) +
  plot_annotation(tag_levels = "a") &
  panel_label_theme

write_source_data(global_budget, "figure3_global_feature_budget.tsv")
write_source_data(task_budget, "figure3_task_feature_budget.tsv")
write_source_data(dataset_heat, "figure3_dataset_budget_deviation.tsv")
write_source_data(method_optima, "figure3_method_optimal_budget.tsv")

save_figure_set(
  figure,
  file.path(FIGURE_ROOT, "figures", "figure3_feature_budget"),
  width_mm = 183,
  height_mm = 148
)

message("Figure 3 written to ", file.path(FIGURE_ROOT, "figures"))
