scale_spatial_metrics <- function(metrics, baseline_ranges) {
    metrics |>
        dplyr::left_join(
            baseline_ranges,
            by = c("Dataset", "Metric", "Type")
        ) |>
        dplyr::mutate(
            ScaledValue = (.data$Value - .data$Lower) / .data$Range
        ) |>
        dplyr::mutate(
            ScaledValue = pmax(pmin(.data$ScaledValue, 1), 0)
        )
}

summarise_spatial_metrics <- function(metrics, baseline_ranges,
                                      type_weights = c(
                                          "Integration" = 0.25,
                                          "Clustering" = 0.25,
                                          "Alignment" = 0.25,
                                          "SliceRepresentation" = 0.25
                                      ),
                                      require_types_for_overall = names(type_weights)) {
    scaled <- scale_spatial_metrics(metrics, baseline_ranges)

    type_means <- scaled |>
        dplyr::group_by(
            .data$Dataset,
            .data$Platform,
            .data$NSlices,
            .data$Method,
            .data$MethodBase,
            .data$SelFeatures,
            IntegrationMethod = .data$Integration,
            .data$IntegrationLabel,
            .data$Type
        ) |>
        dplyr::summarise(
            TypeMean = mean(.data$ScaledValue, na.rm = TRUE),
            .groups = "drop"
        )

    type_wide <- type_means |>
        tidyr::pivot_wider(
            names_from = "Type",
            values_from = "TypeMean"
        )

    available_types <- intersect(names(type_weights), colnames(type_wide))
    required_types <- intersect(require_types_for_overall, colnames(type_wide))
    weighted_sum <- Reduce(
        `+`,
        purrr::map(
            available_types,
            function(.type) {
                tidyr::replace_na(type_wide[[.type]], 0) * type_weights[[.type]]
            }
        )
    )
    required_complete <- Reduce(
        `&`,
        purrr::map(
            required_types,
            function(.type) {
                !is.na(type_wide[[.type]])
            }
        ),
        init = rep(TRUE, nrow(type_wide))
    )
    fixed_weight_total <- sum(type_weights[required_types])

    type_wide |>
        dplyr::mutate(
            Overall = dplyr::if_else(
                required_complete & fixed_weight_total > 0,
                weighted_sum / fixed_weight_total,
                NA_real_
            )
        )
}
