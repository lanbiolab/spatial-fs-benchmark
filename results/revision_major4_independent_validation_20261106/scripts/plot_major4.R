#!/usr/bin/env Rscript
# Reuse the submitted Figure 2 plotting setup and export function without evaluating its panels.
original_path <- "scripts/generate_main_spatial_validation_figure.R"
original <- readLines(original_path)
setup_end <- grep("^semi <-", original)[1] - 1L
eval(parse(text = paste(original[seq_len(setup_end)], collapse = "\n")))
export_start <- grep("^save_figure <-", original)[1]
export_end <- grep('^stem <- file.path', original)[1] - 1L
eval(parse(text = paste(original[export_start:export_end], collapse = "\n")))

root <- "results/revision_major4_independent_validation_20261106"
analysis <- file.path(root, "analysis")
baseline <- read_tsv(file.path(analysis, "original_protocol_seed0_fold_results.tsv"), show_col_types = FALSE)
valid <- read_tsv(file.path(analysis, "valid_annotation_three_seeds_fold_results.tsv"), show_col_types = FALSE)
method_order <- baseline |>
    group_by(.data$method) |>
    summarise(Mean = mean(.data$ARI), .groups = "drop") |>
    arrange(desc(.data$Mean), .data$method) |>
    pull(.data$method)
method_levels <- rev(method_label(method_order))
point_colours <- c("DLPFC" = "#377eb8", "MouseBrainSerialSections" = "#4daf4a", "STOmics0212" = "#e41a1c")

draw_ranked <- function(data, facet, xlab, palette, tag = NULL, reference_x = NULL) {
    data <- data |> mutate(Method = factor(method_label(.data$method), levels = method_levels))
    means <- data |> group_by(.data$method, .data$Method, .data$Panel) |>
        summarise(Mean = mean(.data$Value), SD = sd(.data$Value), .groups = "drop")
    bands <- data |> filter(.data$method %in% c("all_features", "random", "wilcoxon")) |>
        distinct(.data$Method, .data$Panel) |> mutate(y = as.numeric(.data$Method))
    p <- ggplot(data, aes(x = .data$Value, y = .data$Method)) +
        geom_rect(data = bands, aes(xmin = -Inf, xmax = Inf, ymin = .data$y - 0.47, ymax = .data$y + 0.47),
                  inherit.aes = FALSE, fill = "grey70", alpha = 0.28, colour = NA) +
        geom_point(aes(colour = .data$PointGroup), size = 1.05, alpha = 0.6,
                   position = position_jitter(height = 0.09, width = 0, seed = 17)) +
        geom_linerange(data = means, aes(xmin = .data$Mean - .data$SD, xmax = .data$Mean + .data$SD, y = .data$Method),
                      inherit.aes = FALSE, linewidth = 0.38, colour = "grey25") +
        geom_point(data = means, aes(x = .data$Mean, y = .data$Method), inherit.aes = FALSE,
                   shape = 23, size = 1.75, stroke = 0.3, fill = "grey20", colour = "white") +
        facet_wrap(~ .data$Panel, nrow = 1, scales = "free_x") +
        scale_colour_manual(values = palette, name = NULL) +
        scale_y_discrete(drop = FALSE) +
        labs(x = xlab, y = NULL, tag = tag) + base_theme +
        theme(axis.text.y = element_text(size = 5.5), axis.text.x = element_text(size = 5.5),
              axis.line.y = element_blank(), axis.ticks.y = element_blank(),
              strip.text = element_text(size = 6, colour = "white", face = "bold"),
              panel.spacing.x = unit(0.35, "cm"), legend.position = "bottom",
              legend.text = element_text(size = 6), legend.key.size = unit(0.28, "cm"))
    if (!is.null(reference_x)) p <- p + geom_vline(xintercept = reference_x, linewidth = 0.28, colour = "grey60", linetype = "dashed")
    p
}

ari_data <- bind_rows(
    baseline |> mutate(Panel = "Original held-out ARI (seed 0)", Value = .data$ARI),
    valid |> mutate(Panel = "Valid-annotation ARI (3-seed mean)", Value = .data$ARI)
) |> mutate(PointGroup = .data$dataset, Panel = factor(.data$Panel, levels = c("Original held-out ARI (seed 0)", "Valid-annotation ARI (3-seed mean)")))
ari_plot <- draw_ranked(ari_data, "Panel", "Adjusted Rand index", point_colours)
save_figure(ari_plot, file.path(output_dir, "Supplemental_Major4_heldout_ARI"), width_mm = 183, height_mm = 106)
# Embed fonts in these new supplemental PDFs so that repeated exports render all labels consistently.
ggsave(file.path(output_dir, "Supplemental_Major4_heldout_ARI.pdf"), ari_plot,
       width = 183, height = 106, units = "mm", device = grDevices::cairo_pdf, family = "Arial", bg = "white")
write_tsv(ari_data, file.path(output_dir, "Supplemental_Major4_heldout_ARI_source.tsv"))

markers <- read_tsv(file.path(analysis, "dlpfc_external_marker_recovery_by_fold.tsv"), show_col_types = FALSE) |>
    filter(.data$marker_set == "All human-supported layers") |>
    mutate(PointGroup = if_else(grepl("151507", .data$fold), "Held-out 151507", "Held-out 151673"))
marker_long <- bind_rows(
    markers |> mutate(Panel = "Marker recall", Value = .data$marker_recall),
    markers |> mutate(Panel = "Uniform-background enrichment", Value = .data$uniform_fold_enrichment),
    markers |> mutate(Panel = "Matched-background enrichment", Value = .data$matched_fold_enrichment)
) |> mutate(Panel = factor(.data$Panel, levels = c("Marker recall", "Uniform-background enrichment", "Matched-background enrichment")))
marker_plot <- draw_ranked(marker_long, "Panel", "Recall / fold enrichment", c("Held-out 151507" = "#377eb8", "Held-out 151673" = "#e41a1c"))
null_lines <- marker_long |> distinct(.data$Panel) |> filter(.data$Panel != "Marker recall") |> mutate(Null = 1)
marker_plot <- marker_plot + geom_vline(data = null_lines, aes(xintercept = .data$Null), linewidth = 0.28, colour = "grey60", linetype = "dashed")
save_figure(marker_plot, file.path(output_dir, "Supplemental_Major4_external_markers"), width_mm = 183, height_mm = 106)
ggsave(file.path(output_dir, "Supplemental_Major4_external_markers.pdf"), marker_plot,
       width = 183, height = 106, units = "mm", device = grDevices::cairo_pdf, family = "Arial", bg = "white")
write_tsv(marker_long, file.path(output_dir, "Supplemental_Major4_external_markers_source.tsv"))
write_tsv(tibble(Check = c("ARI methods", "ARI folds", "marker methods", "marker folds"),
                 Observed = c(n_distinct(ari_data$method), n_distinct(ari_data$fold), n_distinct(markers$method), n_distinct(markers$fold)),
                 Expected = c(14, 6, 14, 2)) |> mutate(Pass = .data$Observed == .data$Expected), file.path(output_dir, "plot_QA.tsv"))
message("Major 4 supplemental candidates exported; original Figure 2 was not modified.")
