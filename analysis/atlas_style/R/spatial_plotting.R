source(file.path("external", "atlas-feature-selection-benchmark", "analysis", "R", "plotting.R"))

type_names_spatial <- c(
    "Overall" = "Overall",
    "Integration" = "Integration",
    "Clustering" = "Clustering",
    "Alignment" = "Alignment",
    "SliceRepresentation" = "Slice Representation"
)

type_palette_spatial <- c(
    "Overall" = "#f781bf",
    "Integration" = "#e41a1c",
    "Clustering" = "#377eb8",
    "Alignment" = "#4daf4a",
    "SliceRepresentation" = "#984ea3"
)

plot_overview_heatmap <- function(metrics_summary, methods_meta) {
    method_order <- compute_method_order(metrics_summary, methods_meta)

    overview <- metrics_summary |>
        dplyr::group_by(.data$Method, .data$IntegrationLabel) |>
        dplyr::summarise(
            Integration = mean(.data$Integration, na.rm = TRUE),
            Clustering = mean(.data$Clustering, na.rm = TRUE),
            Alignment = mean(.data$Alignment, na.rm = TRUE),
            SliceRepresentation = mean(.data$SliceRepresentation, na.rm = TRUE),
            Overall = mean(.data$Overall, na.rm = TRUE),
            .groups = "drop"
        ) |>
        tidyr::pivot_longer(
            cols = c("Integration", "Clustering", "Alignment", "SliceRepresentation", "Overall"),
            names_to = "Type",
            values_to = "Value"
        ) |>
        dplyr::mutate(
            Value = dplyr::if_else(is.nan(.data$Value), NA_real_, .data$Value),
            Column = ifelse(.data$Type == "Overall", "Overall", paste(.data$IntegrationLabel, .data$Type, sep = "\n")),
            MethodLabel = methods_meta$Name[match(.data$Method, methods_meta$Method)]
        )

    overview |>
        dplyr::mutate(
            MethodLabel = factor(.data$MethodLabel, levels = rev(method_order)),
            Column = factor(.data$Column, levels = unique(.data$Column))
        ) |>
        ggplot2::ggplot(ggplot2::aes(x = .data$Column, y = .data$MethodLabel, fill = .data$Value)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_gradient2(
            low = "#762a83",
            mid = "white",
            high = "#1b7837",
            midpoint = 0.5,
            na.value = "grey92"
        ) +
        ggplot2::labs(x = NULL, y = NULL, fill = "Mean scaled\nvalue") +
        theme_features_pub() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
            panel.grid = ggplot2::element_blank()
        )
}

compute_method_order <- function(metrics_summary, methods_meta) {
    ranking <- metrics_summary |>
        dplyr::group_by(.data$Dataset, .data$IntegrationLabel) |>
        dplyr::mutate(
            RankOverall = rank(-.data$Overall)
        ) |>
        dplyr::ungroup() |>
        dplyr::group_by(.data$Method) |>
        dplyr::summarise(MeanRank = mean(.data$RankOverall, na.rm = TRUE), .groups = "drop") |>
        dplyr::mutate(MethodLabel = methods_meta$Name[match(.data$Method, methods_meta$Method)]) |>
        dplyr::arrange(.data$MeanRank)

    ranking$MethodLabel
}

plot_overview_dotpanels <- function(metrics_scaled, metrics_summary, methods_meta) {
    method_order <- compute_method_order(metrics_summary, methods_meta)

    dataset_scores <- metrics_scaled |>
        dplyr::group_by(.data$Dataset, .data$Method, .data$Type) |>
        dplyr::summarise(TaskScore = mean(.data$ScaledValue, na.rm = TRUE), .groups = "drop") |>
        dplyr::mutate(
            TaskScore = dplyr::if_else(is.nan(.data$TaskScore), NA_real_, .data$TaskScore),
            MethodLabel = methods_meta$Name[match(.data$Method, methods_meta$Method)],
            Type = factor(.data$Type, levels = names(type_names_spatial), labels = type_names_spatial),
            Dataset = factor(.data$Dataset, levels = sort(unique(.data$Dataset)))
        )

    summary_scores <- dataset_scores |>
        dplyr::group_by(.data$Method, .data$MethodLabel, .data$Type) |>
        dplyr::summarise(
            MeanScore = mean(.data$TaskScore, na.rm = TRUE),
            SDScore = sd(.data$TaskScore, na.rm = TRUE),
            .groups = "drop"
        ) |>
        dplyr::mutate(
            MeanScore = dplyr::if_else(is.nan(.data$MeanScore), NA_real_, .data$MeanScore),
            SDScore = dplyr::if_else(is.na(.data$SDScore), 0, .data$SDScore)
        )

    baseline_rows <- methods_meta |>
        dplyr::filter(.data$Kind == "selector-setting", .data$IsBaseline) |>
        dplyr::mutate(
            MethodLabel = .data$Name,
            y = match(.data$Name, method_order)
        ) |>
        dplyr::filter(!is.na(.data$y))

    baseline_bands <- tidyr::expand_grid(
        Type = factor(unname(type_names_spatial), levels = unname(type_names_spatial)),
        baseline_rows
    ) |>
        dplyr::transmute(
            Type = .data$Type,
            ymin = .data$y - 0.5,
            ymax = .data$y + 0.5,
            xmin = -Inf,
            xmax = Inf
        )

    dataset_scores |>
        dplyr::mutate(MethodLabel = factor(.data$MethodLabel, levels = rev(method_order))) |>
        ggplot2::ggplot(ggplot2::aes(x = .data$TaskScore, y = .data$MethodLabel)) +
        ggplot2::geom_rect(
            data = baseline_bands,
            ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax, ymin = .data$ymin, ymax = .data$ymax),
            inherit.aes = FALSE,
            fill = "grey96",
            colour = NA
        ) +
        ggplot2::geom_point(
            ggplot2::aes(colour = .data$Dataset),
            shape = 16,
            size = 1.2,
            alpha = 0.7,
            position = ggplot2::position_jitter(height = 0.14, width = 0)
        ) +
        ggplot2::geom_linerange(
            data = summary_scores |>
                dplyr::mutate(MethodLabel = factor(.data$MethodLabel, levels = rev(method_order))),
            ggplot2::aes(
                xmin = pmax(.data$MeanScore - .data$SDScore, 0),
                xmax = pmin(.data$MeanScore + .data$SDScore, 1),
                y = .data$MethodLabel
            ),
            inherit.aes = FALSE,
            linewidth = 0.35,
            colour = "black",
            alpha = 0.65
        ) +
        ggplot2::geom_point(
            data = summary_scores |>
                dplyr::mutate(MethodLabel = factor(.data$MethodLabel, levels = rev(method_order))),
            ggplot2::aes(x = .data$MeanScore, y = .data$MethodLabel),
            inherit.aes = FALSE,
            shape = 23,
            size = 2.1,
            stroke = 0.3,
            fill = "black",
            colour = "black"
        ) +
        ggplot2::facet_grid(. ~ .data$Type, scales = "free_x", space = "free_x") +
        ggplot2::scale_x_continuous(limits = c(0, 1), expand = ggplot2::expansion(mult = c(0.01, 0.02))) +
        ggplot2::labs(
            x = "Dataset-level aligned score",
            y = NULL,
            colour = "Dataset"
        ) +
        theme_features_pub() +
        ggplot2::theme(
            panel.grid.major.y = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            strip.text.x = ggplot2::element_text(face = "bold"),
            legend.position = "bottom"
        )
}

plot_setting_heatmap <- function(metrics_summary, methods_meta) {
    setting_meta <- methods_meta |>
        dplyr::filter(.data$Kind == "selector-setting")

    setting_summary <- metrics_summary |>
        dplyr::group_by(.data$Method, .data$IntegrationLabel) |>
        dplyr::summarise(
            Integration = mean(.data$Integration, na.rm = TRUE),
            Clustering = mean(.data$Clustering, na.rm = TRUE),
            Alignment = mean(.data$Alignment, na.rm = TRUE),
            SliceRepresentation = mean(.data$SliceRepresentation, na.rm = TRUE),
            Overall = mean(.data$Overall, na.rm = TRUE),
            .groups = "drop"
        ) |>
        tidyr::pivot_longer(
            cols = c("Integration", "Clustering", "Alignment", "SliceRepresentation", "Overall"),
            names_to = "Type",
            values_to = "Value"
        ) |>
        dplyr::mutate(
            Column = ifelse(.data$Type == "Overall", "Overall", paste(.data$IntegrationLabel, .data$Type, sep = "\n")),
            MethodLabel = setting_meta$Name[match(.data$Method, setting_meta$Method)]
        )

    method_order <- setting_summary |>
        dplyr::filter(.data$Type == "Overall") |>
        dplyr::group_by(.data$MethodLabel) |>
        dplyr::summarise(Value = mean(.data$Value, na.rm = TRUE), .groups = "drop") |>
        dplyr::arrange(.data$Value) |>
        dplyr::pull(.data$MethodLabel)

    setting_summary |>
        dplyr::mutate(
            MethodLabel = factor(.data$MethodLabel, levels = method_order),
            Column = factor(.data$Column, levels = unique(.data$Column))
        ) |>
        ggplot2::ggplot(ggplot2::aes(x = .data$Column, y = .data$MethodLabel, fill = .data$Value)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_gradient2(low = "#762a83", mid = "white", high = "#1b7837", midpoint = 0.5) +
        ggplot2::labs(x = NULL, y = NULL, fill = "Mean scaled\nvalue") +
        theme_features_pub() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
            panel.grid = ggplot2::element_blank()
        )
}

plot_overview_ranking <- function(metrics_summary, methods_meta) {
    ranking <- metrics_summary |>
        dplyr::group_by(.data$Dataset, .data$IntegrationLabel) |>
        dplyr::mutate(
            RankOverall = rank(-.data$Overall),
            RankIntegration = rank(-.data$Integration),
            RankClustering = rank(-.data$Clustering),
            RankAlignment = rank(-.data$Alignment),
            RankSliceRepresentation = rank(-.data$SliceRepresentation)
        ) |>
        dplyr::ungroup() |>
        dplyr::select(
            .data$Method,
            .data$RankOverall,
            .data$RankIntegration,
            .data$RankClustering,
            .data$RankAlignment,
            .data$RankSliceRepresentation
        ) |>
        tidyr::pivot_longer(
            cols = c("RankOverall", "RankIntegration", "RankClustering", "RankAlignment", "RankSliceRepresentation"),
            names_to = "Key",
            values_to = "Rank"
        ) |>
        dplyr::mutate(
            Type = dplyr::recode(
                .data$Key,
                "RankOverall" = "Overall",
                "RankIntegration" = "Integration",
                "RankClustering" = "Clustering",
                "RankAlignment" = "Alignment",
                "RankSliceRepresentation" = "SliceRepresentation"
            )
        ) |>
        dplyr::group_by(.data$Method, .data$Type) |>
        dplyr::summarise(
            MeanRank = mean(.data$Rank, na.rm = TRUE),
            SDRank = sd(.data$Rank, na.rm = TRUE),
            .groups = "drop"
        ) |>
        dplyr::mutate(
            MethodLabel = methods_meta$Name[match(.data$Method, methods_meta$Method)],
            Type = factor(.data$Type, levels = names(type_names_spatial), labels = type_names_spatial)
        )

    ranking |>
        dplyr::mutate(MethodLabel = factor(.data$MethodLabel, levels = rev(compute_method_order(metrics_summary, methods_meta)))) |>
        ggplot2::ggplot(ggplot2::aes(x = .data$Type, y = .data$MethodLabel, colour = .data$Type)) +
        ggplot2::geom_point(
            ggplot2::aes(alpha = -.data$MeanRank, size = .data$SDRank),
            shape = "square",
            na.rm = TRUE
        ) +
        ggplot2::scale_colour_manual(values = type_palette_spatial) +
        ggplot2::scale_size_continuous(
            trans = "reverse",
            limits = c(max(ranking$SDRank, na.rm = TRUE), 0),
            range = c(0.2, 3)
        ) +
        ggplot2::labs(x = NULL, y = NULL, size = "Rank SD", alpha = "Mean rank", colour = NULL) +
        theme_features_pub() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
            panel.grid = ggplot2::element_blank(),
            legend.position = "right"
        )
}

plot_task_panels <- function(metrics_scaled, methods_meta) {
    task_data <- metrics_scaled |>
        dplyr::filter(.data$Metric %in% c("dASW", "dLISI", "ILL", "bASW", "iLISI", "GC",
                                          "ari", "nmi", "CHAOS", "PAS",
                                          "Accuracy", "Ratio",
                                          "slice_repr_ARI", "slice_repr_NMI")) |>
        dplyr::mutate(
            MethodLabel = methods_meta$Name[match(.data$Method, methods_meta$Method)],
            MetricLabel = dplyr::case_when(
                .data$Metric == "ari" ~ "ARI",
                .data$Metric == "nmi" ~ "NMI",
                .data$Metric == "slice_repr_ARI" ~ "ARI",
                .data$Metric == "slice_repr_NMI" ~ "NMI",
                TRUE ~ .data$Metric
            ),
            TaskFacet = factor(.data$Type, levels = c("Integration", "Clustering", "Alignment", "SliceRepresentation"))
        )

    summaries <- task_data |>
        dplyr::group_by(.data$MethodLabel, .data$IntegrationLabel, .data$MetricLabel, .data$TaskFacet) |>
        dplyr::summarise(
            MeanValue = mean(.data$ScaledValue, na.rm = TRUE),
            SDValue = sd(.data$ScaledValue, na.rm = TRUE),
            .groups = "drop"
        )

    ggplot2::ggplot(task_data, ggplot2::aes(x = .data$ScaledValue, y = .data$MethodLabel, colour = .data$IntegrationLabel)) +
        ggplot2::geom_jitter(alpha = 0.25, height = 0.18, width = 0) +
        ggplot2::geom_linerange(
            data = summaries,
            ggplot2::aes(
                x = .data$MeanValue,
                y = .data$MethodLabel,
                xmin = .data$MeanValue - .data$SDValue,
                xmax = .data$MeanValue + .data$SDValue,
                colour = .data$IntegrationLabel
            ),
            linewidth = 0.7,
            position = ggplot2::position_dodge(width = 0.3),
            inherit.aes = FALSE
        ) +
        ggplot2::geom_point(
            data = summaries,
            ggplot2::aes(x = .data$MeanValue, y = .data$MethodLabel, colour = .data$IntegrationLabel),
            size = 1.8,
            position = ggplot2::position_dodge(width = 0.3),
            inherit.aes = FALSE
        ) +
        ggplot2::facet_grid(.data$TaskFacet ~ .data$MetricLabel, scales = "free_x", space = "free_x") +
        ggplot2::labs(x = "Scaled value", y = NULL, colour = "Integration") +
        theme_features_pub() +
        ggplot2::theme(
            panel.grid.major.y = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)
        )
}

plot_num_features_summary <- function(metrics_summary, methods_meta, datasets_meta) {
    num_data <- metrics_summary |>
        dplyr::group_by(.data$Dataset, .data$MethodBase, .data$SelFeatures) |>
        dplyr::summarise(
            Integration = mean(.data$Integration, na.rm = TRUE),
            Clustering = mean(.data$Clustering, na.rm = TRUE),
            Alignment = mean(.data$Alignment, na.rm = TRUE),
            SliceRepresentation = mean(.data$SliceRepresentation, na.rm = TRUE),
            Overall = mean(.data$Overall, na.rm = TRUE),
            .groups = "drop"
        ) |>
        tidyr::pivot_longer(
            cols = c("Integration", "Clustering", "Alignment", "SliceRepresentation", "Overall"),
            names_to = "Type",
            values_to = "Value"
        ) |>
        dplyr::mutate(
            MethodLabel = methods_meta$Name[match(.data$MethodBase, methods_meta$Method)],
            DatasetLabel = datasets_meta$Name[match(.data$Dataset, datasets_meta$Dataset)],
            Type = factor(.data$Type, levels = names(type_names_spatial), labels = type_names_spatial)
        )

    top_counts <- num_data |>
        dplyr::group_by(.data$Dataset, .data$MethodLabel, .data$Type) |>
        dplyr::slice_max(.data$Value, n = 1, with_ties = FALSE) |>
        dplyr::group_by(.data$Type, .data$SelFeatures) |>
        dplyr::count(name = "Count")

    counts_plot <- ggplot2::ggplot(top_counts, ggplot2::aes(x = factor(.data$SelFeatures), y = .data$Count, fill = .data$Type)) +
        ggplot2::geom_col() +
        ggplot2::facet_wrap(~ .data$Type, nrow = 1) +
        ggplot2::scale_fill_manual(values = type_palette_spatial) +
        ggplot2::labs(x = "Number of selected features", y = "Top-value\ncounts", fill = NULL) +
        theme_features_pub() +
        ggplot2::theme(
            legend.position = "bottom",
            axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)
        )

    heatmap_plot <- ggplot2::ggplot(
        num_data,
        ggplot2::aes(x = factor(.data$SelFeatures), y = .data$MethodLabel, fill = .data$Value)
    ) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_gradient2(low = "#762a83", mid = "white", high = "#1b7837", midpoint = 0.5) +
        ggplot2::facet_grid(.data$DatasetLabel ~ .data$Type) +
        ggplot2::labs(x = "Number of selected features", y = NULL, fill = "Mean scaled\nvalue") +
        theme_features_pub() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
            panel.grid = ggplot2::element_blank(),
            strip.text.x = ggplot2::element_blank()
        )

    patchwork::wrap_plots(counts_plot, heatmap_plot, ncol = 1, heights = c(1, 6), guides = "collect")
}

plot_stability_summary <- function(metrics_summary, datasets_meta, methods_meta) {
    n_datasets <- length(unique(metrics_summary$Dataset))
    n_platforms <- length(unique(metrics_summary$Platform))

    stability_data <- metrics_summary |>
        dplyr::group_by(.data$Dataset, .data$Platform, .data$Method) |>
        dplyr::summarise(
            OverallMean = mean(.data$Overall, na.rm = TRUE),
            OverallSD = sd(.data$Overall, na.rm = TRUE),
            .groups = "drop"
        ) |>
        dplyr::mutate(
            DatasetLabel = datasets_meta$Name[match(.data$Dataset, datasets_meta$Dataset)],
            MethodLabel = methods_meta$Name[match(.data$Method, methods_meta$Method)]
        )

    if (n_datasets > 1 || n_platforms > 1) {
        ggplot2::ggplot(
            stability_data,
            ggplot2::aes(x = .data$DatasetLabel, y = .data$MethodLabel, colour = .data$OverallMean, size = .data$OverallSD)
        ) +
            ggplot2::geom_point(shape = "square") +
            ggplot2::scale_colour_gradient2(low = "#762a83", mid = "white", high = "#1b7837", midpoint = 0.5) +
            ggplot2::scale_size_continuous(trans = "reverse", range = c(0.2, 3)) +
            ggplot2::labs(
                x = NULL,
                y = NULL,
                colour = "Mean overall\nscaled value",
                size = "Overall SD",
                title = "Dataset/platform stability"
            ) +
            theme_features_pub() +
            ggplot2::theme(
                axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)
            )
    } else {
        integration_stability <- metrics_summary |>
            dplyr::group_by(.data$IntegrationLabel, .data$Method) |>
            dplyr::summarise(
                OverallMean = mean(.data$Overall, na.rm = TRUE),
                OverallSD = sd(.data$Overall, na.rm = TRUE),
                .groups = "drop"
            ) |>
            dplyr::mutate(
                MethodLabel = methods_meta$Name[match(.data$Method, methods_meta$Method)]
            )
        ggplot2::ggplot(
            integration_stability,
            ggplot2::aes(x = .data$IntegrationLabel, y = .data$MethodLabel, colour = .data$OverallMean, size = .data$OverallSD)
        ) +
            ggplot2::geom_point(shape = "square") +
            ggplot2::scale_colour_gradient2(low = "#762a83", mid = "white", high = "#1b7837", midpoint = 0.5) +
            ggplot2::scale_size_continuous(trans = "reverse", range = c(0.2, 3)) +
            ggplot2::labs(
                x = NULL,
                y = NULL,
                colour = "Mean overall\nscaled value",
                size = "Overall SD",
                title = "Stability summary",
                subtitle = "Current run contains one dataset/platform, so stability is shown across integration methods"
            ) +
            theme_features_pub() +
            ggplot2::theme(
                axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)
            )
    }
}
