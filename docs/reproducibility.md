# Reproducibility Notes

## Workflow

The benchmark uses versioned Python/R scripts and YAML configuration files:

1. convert public datasets to local h5ad inputs;
2. validate the integer-count contract;
3. run feature selection, scVI/CellCharter, and task evaluation;
4. reconstruct metric rows from immutable task records;
5. audit embedding and metric keys;
6. build one frozen Dataset-by-Metric scaling table;
7. generate manuscript source-data tables and figures with R.

No Snakemake or Nextflow workflow is used.

## Frozen Dataset Set

Seven count-validated datasets enter the rebuilt downstream benchmark:

- DLPFC
- MouseBrainSerialSections
- STOmics0212
- STOmics0218
- STOmics0224
- E8p5Embryo
- E9p5Embryo

STOmicsVisium5Samples is excluded because its available `raw_count` layer is
transformed rather than integer-valued. E8.5 and E9.5 lack benchmark labels;
Wilcoxon, label-agreement metrics, and alignment are therefore missing for
these datasets. Alignment is evaluated only for explicit DLPFC and Mouse Brain
section pairs.

## Scores

Thirteen metrics are oriented so larger values indicate better performance and
scaled with frozen observed ranges within each Dataset-by-Metric stratum.
Integration gives equal weight to batch mixing and biological conservation;
Clustering gives equal weight to label-independent and label-agreement
components when both are available.

`CoreOverall` is the equal-weight mean of Integration and Clustering and is
defined for all seven datasets. Alignment is reported separately. Missing
values are never replaced by zero.

## Representative Settings

Representative ranking figures use:

- `all_features`: all available genes;
- `random`: 500 genes;
- `scsegindex`: 500 genes;
- other selectors: 2,000 genes where returned.

Feature-number analyses use setting-level results at 100, 200, 500, 1,000,
2,000, 5,000, and 10,000 requested genes. Canonical settings use seeds 0, 1,
and 2; noncanonical feature-number settings use seed 0.

See `docs/rebuild_v1_protocol.md` and the manifest under
`results/v2/frozen_scores/` for exact aggregation and checksums.
