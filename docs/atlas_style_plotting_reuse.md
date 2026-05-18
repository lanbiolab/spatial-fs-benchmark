# Atlas-Style Plotting Reuse Plan

This document records how the plotting system from
`atlas-feature-selection-benchmark` is reused for the spatial benchmark.

## 1. Reuse audit

### Directly reusable

- `external/atlas-feature-selection-benchmark/analysis/R/plotting.R`
  - `theme_features()`
  - `theme_features_pub()`
  - `save_figure_files()`
- `external/atlas-feature-selection-benchmark/analysis/R/summarisation.R`
  - `scale_metrics()` logic
- `external/atlas-feature-selection-benchmark/analysis/04-benchmark.Rmd`
  - overall score layout
  - ranking heatmap / ranking bar-chart logic
- `external/atlas-feature-selection-benchmark/analysis/03-num-features.Rmd`
  - number-of-features ablation structure
  - top-value count bar chart + heatmap stack
- `external/atlas-feature-selection-benchmark/reports/functions.R`
  - method-vs-metric dotplot summary logic

### Reusable with small modifications

- `analysis/04-benchmark.Rmd`
  - replace scRNA task groups with spatial task groups
  - replace `Method` interpretation from feature-selection setting to spatial FS setting
  - replace `Integration == "scVI-1"` filtering with spatial integration methods
- `analysis/03-num-features.Rmd`
  - keep the same feature-number narrative
  - map `SelFeatures` to `n_features`
  - change metric families from scRNA categories to spatial task categories
- `reports/functions.R`
  - keep dotplot / variance structure
  - adapt faceting from scRNA integration families to spatial integration methods

### Must be rewritten

- query/reference mapping figures that are specific to scRNA transfer tasks
- lineage / HLCA subset analyses that are tightly tied to scRNA biological subsets
- any UMAP-based case-study notebook that depends on atlas-specific AnnData outputs

## 2. Input compatibility layer

The spatial benchmark results are adapted into an atlas-style long table:

- `results/atlas_style/data/spatial-benchmark-long.tsv`
- `results/atlas_style/data/benchmark.tsv`
- `results/atlas_style/data/num-features.tsv`

These files preserve the first-paper plotting expectations while injecting
spatial content.

### Core compatibility fields

- `Dataset`
- `Platform`
- `NSlices`
- `Method`
- `MethodBase`
- `SelFeatures`
- `Integration`
- `IntegrationLabel`
- `Task`
- `TaskLabel`
- `Type`
- `Metric`
- `MetricName`
- `Value`
- `ValueRaw`
- `HigherBetter`
- `Seed`
- `Runtime`
- `Notes`

## 3. Field mapping

### Benchmark core fields

- spatial `dataset` -> atlas-style `Dataset`
- spatial `platform` -> atlas-style `Platform`
- spatial `n_slices` -> atlas-style `NSlices`
- spatial `fs_method` -> atlas-style `MethodBase`
- spatial `fs_method + n_features` -> atlas-style `Method`
- spatial `n_features` -> atlas-style `SelFeatures`
- spatial `integration_method` -> atlas-style `Integration`
- spatial `task` -> atlas-style `Task`
- spatial grouped task family -> atlas-style `Type`
- spatial `metric_name` -> atlas-style `Metric`
- spatial `metric_value` -> atlas-style `ValueRaw`
- direction-adjusted metric value -> atlas-style `Value`
- spatial `random_seed` -> atlas-style `Seed`
- spatial `runtime` -> atlas-style `Runtime`
- spatial `notes` -> atlas-style `Notes`

### Task group mapping

- `integration_eval` -> `Integration`
- `clustering_eval` -> `Clustering`
- `alignment_eval` -> `Alignment`
- `slice_representation_eval` -> `SliceRepresentation`

### Metric direction mapping

Higher-is-better after compatibility adjustment:

- `dASW`, `dLISI`, `ILL`, `bASW`, `iLISI`, `GC`
- `ARI`, `NMI`
- `Accuracy`

Lower-is-better before adjustment:

- `CHAOS`
- `PAS`
- `Ratio`
- `slice_repr_distance_mean`

## 4. Plotting modules to keep / adapt / rewrite

### Keep

- theme and export helpers from the first paper
- score scaling and type-level aggregation pattern
- overall ranking layout
- num-features ablation layout

### Adapt

- overall summary weighting:
  - first paper: integration + query/classification/unseen
  - spatial version: integration + clustering + alignment + slice representation
- metric metadata:
  - first paper uses scRNA metrics metadata files
  - spatial version uses spatial task metadata TSVs
- method metadata:
  - first paper labels each FS configuration as a method
  - spatial version keeps both `MethodBase` and `Method`

### Rewrite

- spatial case study
- slice-level spatial domain map visualisation
- alignment visualisation

## 5. Script organisation

### Adapter layer

- `scripts/adapt_atlas_style_results.py`
  - converts benchmark CSV to atlas-style long tables
  - writes metric, method, dataset, integration metadata
  - writes scaling ranges

### R plotting layer

- `analysis/atlas_style/R/spatial_summarisation.R`
  - spatial equivalent of atlas summary aggregation
- `analysis/atlas_style/R/spatial_plotting.R`
  - sources atlas plotting helpers directly
  - defines spatial figure constructors
- `scripts/generate_atlas_style_figures.R`
  - generates Figures 1-4

### Python case-study layer

- `scripts/generate_spatial_case_study.py`
  - generates Figure 5

## 6. Recommended pipeline structure

```text
analysis/
  atlas_style/
    R/
      spatial_plotting.R
      spatial_summarisation.R
docs/
  atlas_style_plotting_reuse.md
results/
  atlas_style/
    data/
    output/
    figures/
scripts/
  adapt_atlas_style_results.py
  generate_atlas_style_figures.R
  generate_spatial_case_study.py
```

## 7. Final figure narrative

- Figure 1: overview heatmap + overall ranking
- Figure 2: per-task metric panels
- Figure 3: feature-number ablation
- Figure 4: dataset/platform stability
- Figure 5: spatial case study

This preserves the first paper's high-level storytelling order while replacing
the content with spatial multi-slice integration and downstream tasks.
