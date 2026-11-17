from pathlib import Path
import hashlib

ROOT = Path('/home/xd/y/Benchmark')
OUT = ROOT / 'results/revision_minor1_3_20261106'
source = ROOT / 'scripts/generate_num_features_benchmark_figure.R'
text = source.read_text()
anchor = 'overall_lineplot <- ggplot('
assert text.count(anchor) == 1
code = r'''
# R1 minor 3: add study-resampled uncertainty; retain the original plot objects,
# transformations, panels, palettes, typography and 7.8 x 6.8 inch dimensions.
ci <- readr::read_tsv(
    "results/revision_minor1_3_20261106/analysis/figure3_panelA_standardized_study_ci.tsv",
    show_col_types = FALSE
)
ci_check <- overall_means |>
    dplyr::mutate(N = as.integer(as.character(.data$SelFeatures))) |>
    dplyr::left_join(ci, by = c("Type", "N"))
stopifnot(nrow(ci_check) == 28, all(is.finite(ci_check$MeanStandardValue)),
          max(abs(ci_check$StandardValue - ci_check$MeanStandardValue)) < 1e-10)
ci <- ci |> dplyr::mutate(SelFeatures = factor(as.character(.data$N), levels = sel_levels))
readr::write_tsv(ci_check, file.path(output_dir, "Figure3_panelA_mean_reproduction_audit.tsv"))
'''
text = text.replace(anchor, code + '\n' + anchor)
anchor = '''    geom_line(
        data = overall_means,'''
assert text.count(anchor) == 1
errorbar = r'''    geom_errorbar(
        data = ci,
        aes(x = .data$SelFeatures, ymin = .data$Lower95, ymax = .data$Upper95, colour = .data$Type),
        inherit.aes = FALSE, width = 0.16, linewidth = 0.32
    ) +
'''
text = text.replace(anchor, errorbar + anchor)
extra = r'''

# Embed fonts in the new candidate; the source plotting theme is unchanged.
ggsave(file.path(output_dir, "figure_num_features_benchmark.pdf"), summary_plot,
       width = 7.8, height = 6.8, device = grDevices::cairo_pdf, family = "Arial", bg = "white")
ggsave(file.path(output_dir, "figure_num_features_benchmark.svg"), summary_plot,
       width = 7.8, height = 6.8, device = svglite::svglite, bg = "white")
ggsave(file.path(output_dir, "figure_num_features_benchmark.tiff"), summary_plot,
       width = 7.8, height = 6.8, dpi = 600, device = ragg::agg_tiff, compression = "lzw", bg = "white")

raw <- readr::read_tsv("results/revision_minor1_3_20261106/analysis/budget_curves_with_study_uncertainty.tsv", show_col_types = FALSE)
exact_count <- unique(raw$NMethods[raw$panel == "fixed_exact_budget_panel"])
raw <- raw |>
    mutate(BudgetPosition = match(.data$N, c(100, 200, 500, 1000, 2000, 5000, 10000)),
           Panel = if_else(.data$panel == "all_29_requested_budget", "All 29 selectors: requested budget",
                           paste0("Fixed ", exact_count, " selectors: exact budget")),
           Panel = factor(.data$Panel, levels = c("All 29 selectors: requested budget", paste0("Fixed ", exact_count, " selectors: exact budget"))),
           Integrator = factor(.data$integration_method, levels = c("scvi", "cellcharter"), labels = c("scVI", "CellCharter")))
palette <- setNames(c("#f781bf", "#404040"), levels(raw$Panel))
budget_plot <- ggplot(raw, aes(x = .data$BudgetPosition, y = .data$Mean, colour = .data$Panel,
                              fill = .data$Panel, group = .data$Panel)) +
    geom_ribbon(aes(ymin = .data$Lower95, ymax = .data$Upper95), alpha = 0.10, colour = NA) +
    geom_line(aes(linetype = .data$Panel), linewidth = 0.65) +
    geom_point(aes(shape = .data$Panel), size = 1.55, stroke = 0.5, fill = "white") +
    facet_grid(. ~ .data$Integrator) +
    scale_colour_manual(values = palette, name = NULL) +
    scale_fill_manual(values = palette, name = NULL) +
    scale_linetype_manual(values = c("solid", "dashed"), name = NULL) +
    scale_shape_manual(values = c(23, 21), name = NULL) +
    scale_x_continuous(breaks = 1:7, labels = c(100, 200, 500, 1000, 2000, 5000, 10000)) +
    labs(x = "Requested number of selected features", y = "Mean CoreOverall") +
    theme_features_pub() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.position = "bottom", legend.box = "vertical", panel.grid.minor = element_blank(),
          plot.margin = margin(4, 4, 4, 4, "mm"))
save_figure_files(budget_plot, file.path(output_dir, "budget_exact_panel_sensitivity"), width = 7.8, height = 3.6)
ggsave(file.path(output_dir, "budget_exact_panel_sensitivity.pdf"), budget_plot,
       width = 7.8, height = 3.6, device = grDevices::cairo_pdf, family = "Arial", bg = "white")
ggsave(file.path(output_dir, "budget_exact_panel_sensitivity.svg"), budget_plot,
       width = 7.8, height = 3.6, device = svglite::svglite, bg = "white")
ggsave(file.path(output_dir, "budget_exact_panel_sensitivity.tiff"), budget_plot,
       width = 7.8, height = 3.6, dpi = 600, device = ragg::agg_tiff, compression = "lzw", bg = "white")
write_tsv(raw, file.path(output_dir, "budget_exact_panel_plot_source.tsv"))
message("New Minor 3 candidates only; existing main figures were not overwritten.")
'''
(OUT / 'scripts/plot_minor3_original_style.R').write_text(text + extra)
(OUT / 'scripts/original_plot_source_sha256.txt').write_text(hashlib.sha256(source.read_bytes()).hexdigest() + '\n')
print('Prepared minimally extended original Figure 3 plotting script.')
