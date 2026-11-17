#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(jsonlite)
  library(patchwork)
  library(readr)
  library(scales)
  library(tidyr)
})

source("analysis/manuscript_rebuild_v2/R/theme.R")
source("analysis/manuscript_rebuild_v2/R/data.R")
ensure_figure_dirs()

spatial_methods <- c("somde", "nnsvg", "spatialde", "morans_i", "sparkx")

dataset_ranks <- read_score_table("dataset_competitive_ranks.tsv") %>%
  mutate(
    Method = display_method(fs_method),
    Dataset = display_dataset(dataset),
    Integrator = display_integrator(integration_method)
  )

global <- read_score_table("global_competitive_method_ranks.tsv") %>%
  mutate(Method = display_method(fs_method), Integrator = display_integrator(integration_method))

# Panel a: SVG ranks by dataset ----------------------------------------------
svg_order <- global %>%
  filter(fs_method %in% spatial_methods) %>%
  group_by(fs_method, Method) %>%
  summarise(CombinedRank = mean(GlobalCoreOverallRank), .groups = "drop") %>%
  arrange(CombinedRank) %>%
  pull(Method)

dataset_order <- c(
  "DLPFC", "Mouse brain sections", "STOmics 0212", "STOmics 0218",
  "STOmics 0224", "E8.5 embryo", "E9.5 embryo"
)

svg_ranks <- dataset_ranks %>%
  filter(fs_method %in% spatial_methods) %>%
  mutate(
    Method = factor(Method, levels = rev(svg_order)),
    Dataset = factor(Dataset, levels = dataset_order),
    Integrator = factor(Integrator, levels = c("scVI", "CellCharter"))
  )

p_a <- ggplot(svg_ranks, aes(Dataset, Method, fill = CoreOverallRank)) +
  geom_tile(colour = "white", linewidth = 0.55) +
  geom_text(
    aes(label = number(CoreOverallRank, accuracy = 1)),
    family = BENCH_FONT, size = 2.0,
    colour = ifelse(svg_ranks$CoreOverallRank <= 14, "white", BENCH_COLOURS$ink)
  ) +
  facet_grid(Integrator ~ .) +
  scale_fill_gradient(
    low = "#B4482D", high = "#F4EEE7", limits = c(1, 29),
    breaks = c(1, 10, 20, 29), name = "CoreOverall rank"
  ) +
  scale_x_discrete(position = "top") +
  labs(
    x = NULL, y = NULL,
    title = "Spatial selectors are competitive in specific contexts",
    subtitle = "Dataset-level rank; rank 1 is best"
  ) +
  theme_benchmark(6.1) +
  theme(
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.text.x = element_text(angle = 43, hjust = 0, size = 5.4, face = "bold"),
    axis.text.y = element_text(size = 5.8, face = "bold"),
    strip.text = element_text(size = 6.5, face = "bold"),
    legend.position = "bottom", legend.key.width = unit(12, "mm"),
    plot.title = element_text(size = 7.5),
    plot.subtitle = element_text(size = 5.9),
    plot.margin = margin(3, 2, 3, 3, "mm")
  )

# Panel b: best spatial selector versus best expression selector -------------
best_family <- dataset_ranks %>%
  filter(MethodGroup %in% c("expression_driven", "spatially_informed")) %>%
  group_by(dataset, Dataset, integration_method, Integrator, MethodGroup) %>%
  slice_max(CoreOverallMean, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(dataset, Dataset, integration_method, Integrator, MethodGroup,
         fs_method, Method, CoreOverallMean) %>%
  pivot_wider(
    names_from = MethodGroup,
    values_from = c(fs_method, Method, CoreOverallMean)
  ) %>%
  mutate(
    Difference = CoreOverallMean_spatially_informed - CoreOverallMean_expression_driven,
    Direction = if_else(Difference >= 0, "Spatial higher", "Expression higher"),
    Dataset = factor(Dataset, levels = rev(dataset_order)),
    Integrator = factor(Integrator, levels = c("scVI", "CellCharter"))
  )

p_b <- ggplot(best_family, aes(Difference, Dataset)) +
  geom_vline(xintercept = 0, colour = "#9A9A9A", linewidth = 0.38) +
  geom_segment(aes(x = 0, xend = Difference, yend = Dataset, colour = Direction),
               linewidth = 1.0, lineend = "round") +
  geom_point(aes(fill = Direction), shape = 21, size = 2.3,
             stroke = 0.35, colour = "white") +
  facet_grid(. ~ Integrator) +
  scale_colour_manual(values = c(
    "Spatial higher" = BENCH_COLOURS$spatial,
    "Expression higher" = BENCH_COLOURS$expression
  ), guide = "none") +
  scale_fill_manual(values = c(
    "Spatial higher" = BENCH_COLOURS$spatial,
    "Expression higher" = BENCH_COLOURS$expression
  ), guide = "none") +
  scale_x_continuous(labels = label_number(accuracy = 0.01),
                     expand = expansion(mult = c(0.08, 0.08))) +
  labs(
    x = "Best spatial − best expression CoreOverall",
    y = NULL,
    title = "Spatial gains are dataset-dependent",
    subtitle = "Coral indicates a higher score for the best spatial selector"
  ) +
  theme_benchmark(5.9) +
  theme(
    axis.text.y = element_text(size = 5.4, face = "bold"),
    axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    strip.text = element_text(size = 6.2, face = "bold"),
    plot.title = element_text(size = 7.2),
    plot.subtitle = element_text(size = 5.6),
    plot.margin = margin(3, 2, 3, 3, "mm")
  )

# Feature-set overlap --------------------------------------------------------
representative <- read_score_table("representative_settings.tsv") %>%
  filter(MethodGroup %in% c("expression_driven", "spatially_informed")) %>%
  mutate(
    dataset_key = tolower(gsub("[^A-Za-z0-9]", "", dataset)),
    audit_n = as.numeric(n_features)
  )

feature_audit <- read_tsv(
  "results/spatial_svg_rebuild_v1/audit/selected_feature_audit.tsv",
  show_col_types = FALSE
) %>%
  filter(selection_seed == 0) %>%
  select(dataset_key, fs_method, audit_n = n_features, path)

feature_index <- representative %>%
  inner_join(feature_audit, by = c("dataset_key", "fs_method", "audit_n"))

if (nrow(feature_index) != nrow(representative)) {
  stop("Representative feature-set paths are incomplete: ", nrow(feature_index),
       " of ", nrow(representative), call. = FALSE)
}

feature_sets <- lapply(feature_index$path, function(path) {
  unique(as.character(jsonlite::fromJSON(path)$feature_names))
})
names(feature_sets) <- paste(feature_index$dataset, feature_index$fs_method, sep = "::")

pairwise_rows <- lapply(unique(feature_index$dataset), function(ds) {
  ds_index <- feature_index %>% filter(dataset == ds) %>% arrange(fs_method)
  method_pairs <- combn(ds_index$fs_method, 2, simplify = FALSE)
  bind_rows(lapply(method_pairs, function(pair) {
    set1 <- feature_sets[[paste(ds, pair[[1]], sep = "::")]]
    set2 <- feature_sets[[paste(ds, pair[[2]], sep = "::")]]
    group1 <- ds_index$MethodGroup[match(pair[[1]], ds_index$fs_method)]
    group2 <- ds_index$MethodGroup[match(pair[[2]], ds_index$fs_method)]
    pair_type <- case_when(
      group1 == "spatially_informed" && group2 == "spatially_informed" ~ "Spatial–spatial",
      group1 == "expression_driven" && group2 == "expression_driven" ~ "Expression–expression",
      TRUE ~ "Expression–spatial"
    )
    tibble(
      dataset = ds,
      Method1 = pair[[1]], Method2 = pair[[2]],
      Group1 = group1, Group2 = group2, PairType = pair_type,
      Jaccard = length(intersect(set1, set2)) / length(union(set1, set2))
    )
  }))
}) %>% bind_rows()

combined_order <- global %>%
  group_by(fs_method, Method, MethodGroup) %>%
  summarise(CombinedRank = mean(GlobalCoreOverallRank), .groups = "drop") %>%
  arrange(CombinedRank)

matrix_methods <- union(head(combined_order$fs_method, 10), spatial_methods)
matrix_order <- combined_order %>%
  filter(fs_method %in% matrix_methods) %>%
  arrange(CombinedRank) %>%
  pull(fs_method)

matrix_mean <- pairwise_rows %>%
  filter(Method1 %in% matrix_methods, Method2 %in% matrix_methods) %>%
  group_by(Method1, Method2) %>%
  summarise(MeanJaccard = mean(Jaccard), SDJaccard = sd(Jaccard), .groups = "drop")

matrix_full <- bind_rows(
  matrix_mean,
  matrix_mean %>% transmute(
    Method1 = Method2, Method2 = Method1,
    MeanJaccard = MeanJaccard, SDJaccard = SDJaccard
  ),
  tibble(Method1 = matrix_methods, Method2 = matrix_methods,
         MeanJaccard = 1, SDJaccard = 0)
) %>%
  distinct(Method1, Method2, .keep_all = TRUE) %>%
  mutate(
    IsDiagonal = Method1 == Method2,
    Label1 = paste0(display_method(Method1), if_else(Method1 %in% spatial_methods, " *", "")),
    Label2 = paste0(display_method(Method2), if_else(Method2 %in% spatial_methods, " *", "")),
    Label1 = factor(
      Label1,
      levels = paste0(display_method(matrix_order), ifelse(matrix_order %in% spatial_methods, " *", ""))
    ),
    Label2 = factor(
      Label2,
      levels = rev(paste0(display_method(matrix_order), ifelse(matrix_order %in% spatial_methods, " *", "")))
    )
  )

p_c <- ggplot(matrix_full, aes(Label1, Label2)) +
  geom_tile(
    data = matrix_full %>% filter(!IsDiagonal),
    aes(fill = MeanJaccard), colour = "white", linewidth = 0.42
  ) +
  geom_tile(
    data = matrix_full %>% filter(IsDiagonal),
    fill = "#E4B65B", colour = "white", linewidth = 0.42
  ) +
  scale_fill_gradient(
    low = "#F4F1E9", high = "#286878", limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1), name = "Mean Jaccard index"
  ) +
  scale_x_discrete(position = "top") +
  labs(
    x = NULL, y = NULL,
    title = "Representative feature sets are not interchangeable",
    subtitle = "Mean across seven datasets; spatial methods marked *"
  ) +
  theme_benchmark(5.7) +
  theme(
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.text.x = element_text(angle = 46, hjust = 0, size = 5.0, face = "bold"),
    axis.text.y = element_text(size = 5.0, face = "bold"),
    legend.position = "bottom", legend.key.width = unit(12, "mm"),
    plot.title = element_text(size = 7.3),
    plot.subtitle = element_text(size = 5.7),
    plot.margin = margin(3, 3, 3, 3, "mm")
  )

pair_summary <- pairwise_rows %>%
  group_by(PairType) %>%
  summarise(
    MedianJaccard = median(Jaccard),
    MeanJaccard = mean(Jaccard),
    ProportionAboveHalf = mean(Jaccard > 0.5),
    NPairs = n(), .groups = "drop"
  )

pairwise_plot <- pairwise_rows %>%
  mutate(PairType = factor(
    PairType,
    levels = c("Expression–expression", "Expression–spatial", "Spatial–spatial")
  ))

p_d <- ggplot(pairwise_plot, aes(PairType, Jaccard, fill = PairType)) +
  geom_violin(width = 0.86, trim = TRUE, scale = "width", alpha = 0.82,
              colour = "white", linewidth = 0.25) +
  geom_boxplot(width = 0.20, outlier.shape = NA, fill = "white",
               colour = BENCH_COLOURS$ink, linewidth = 0.32) +
  geom_text(
    data = pair_summary,
    aes(PairType, 0.98, label = paste0("median ", sprintf("%.2f", MedianJaccard))),
    inherit.aes = FALSE, family = BENCH_FONT, size = 2.0, fontface = "bold"
  ) +
  scale_fill_manual(values = c(
    "Expression–expression" = "#8DA8BA",
    "Expression–spatial" = "#C9A38F",
    "Spatial–spatial" = BENCH_COLOURS$spatial
  ), guide = "none") +
  scale_y_continuous(limits = c(0, 1.03), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(
    x = NULL, y = "Jaccard index",
    title = "Spatial selectors overlap more with each other",
    subtitle = "Descriptive distributions across dataset–method pairs"
  ) +
  theme_benchmark(5.8) +
  theme(
    axis.text.x = element_text(angle = 18, hjust = 0.5, size = 5.0),
    plot.title = element_text(size = 7.1),
    plot.subtitle = element_text(size = 5.5),
    plot.margin = margin(3, 3, 3, 3, "mm")
  )

figure <- (p_a / p_b + plot_layout(heights = c(0.57, 0.43))) |
  (p_c / p_d + plot_layout(heights = c(0.68, 0.32)))
figure <- figure +
  plot_layout(widths = c(0.49, 0.51)) +
  plot_annotation(tag_levels = "a") &
  panel_label_theme

write_source_data(svg_ranks, "figure4_svg_dataset_ranks.tsv")
write_source_data(best_family, "figure4_best_spatial_vs_expression.tsv")
write_source_data(pairwise_rows, "figure4_representative_pairwise_jaccard.tsv")
write_source_data(matrix_full, "figure4_selected_method_jaccard_matrix.tsv")
write_source_data(pair_summary, "figure4_jaccard_family_summary.tsv")

save_figure_set(
  figure,
  file.path(FIGURE_ROOT, "figures", "figure4_spatial_selector_context"),
  width_mm = 183,
  height_mm = 168
)

message("Figure 4 written to ", file.path(FIGURE_ROOT, "figures"))

