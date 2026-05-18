#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    if (!requireNamespace("remotes", quietly = TRUE)) {
        stop("Package 'remotes' is required before running this installer.", call. = FALSE)
    }
})

install_if_missing <- function(pkg, expr) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        message("Installing missing R package: ", pkg)
        force(expr)
    } else {
        message("R package already available: ", pkg)
    }
}

install_if_missing(
    "singleCellHaystack",
    remotes::install_version(
        "singleCellHaystack",
        version = "0.3.4",
        repos = "https://cloud.r-project.org",
        dependencies = FALSE
    )
)

install_if_missing(
    "DUBStepR",
    remotes::install_github(
        "prabhakarlab/DUBStepR@76aa39485742d5f5bcfb86346a8a1784ee08f6b9",
        dependencies = FALSE
    )
)

install_if_missing(
    "akmedoids",
    remotes::install_version(
        "akmedoids",
        version = "1.3.0",
        repos = "https://cloud.r-project.org",
        dependencies = FALSE
    )
)

install_if_missing(
    "scPNMF",
    remotes::install_github(
        "JSB-UCLA/scPNMF@47d5b10cb09450255aea9b53ace555a95ab69502",
        dependencies = FALSE
    )
)
