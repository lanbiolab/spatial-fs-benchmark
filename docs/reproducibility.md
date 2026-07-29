# Reproducibility Notes

## Workflow

The benchmark uses versioned Python/R scripts and YAML configuration files; it
does not use Snakemake or Nextflow. The auditable sequence is:

1. convert each public dataset to a local AnnData input;
2. validate the non-negative integer-count contract;
3. run feature selection and the downstream representation models;
4. reconstruct metric rows from immutable task records;
5. audit embedding, selected-feature, task, and metric keys;
6. freeze one Dataset-by-Metric scaling table from the primary matrix;
7. run sensitivity and validation analyses without redefining those ranges;
8. generate manuscript source tables and figures with R.

## Frozen Dataset Set

Seven count-validated datasets enter the downstream benchmark: DLPFC,
MouseBrainSerialSections, STOmics0212, STOmics0218, STOmics0224, E8p5Embryo,
and E9p5Embryo. STOmicsVisium5Samples is excluded because its available
`raw_count` layer is transformed rather than integer-valued. E8.5 and E9.5 lack
benchmark labels, so Wilcoxon and label-agreement metrics are missing for those
datasets. Alignment is evaluated only for declared DLPFC and Mouse Brain pairs.

## Integrators And Scores

scVI and CellCharter form the full primary matrix. GraphST is a canonical
third-integrator sensitivity analysis for the 29 competitive methods at seeds
0, 1, and 2. GraphST values use the primary frozen ranges and do not redefine
them.

Thirteen metrics are oriented so larger values indicate better performance.
Integration gives equal weight to batch mixing and biological conservation;
Clustering gives equal weight to label-independent and label-agreement
components when both are available. `CoreOverall` is the equal-weight mean of
Integration and Clustering. Alignment is reported separately. Missing values
are never replaced by zero.

## Representative Settings And Validation

Representative figures use all available genes for `all_features`, 500 genes
for `random` and `scsegindex`, and 2,000 genes for other selectors. Feature-number
analyses request 100, 200, 500, 1,000, 2,000, 5,000, and 10,000 genes.

Validation includes semi-synthetic spatial truth, six leave-one-slice-out folds,
a cross-statistic nnSVG held-out reference, alternative scaling and weighting,
dataset bootstrap intervals, exploratory HVG-SVG union/intersection controls,
and a standardized feature-selector resource profile.

See `docs/rebuild_v1_protocol.md`, `docs/rebuild_v1_scoring.md`, and the manifest
under `results/v2/frozen_scores/` for exact definitions and checksums.
