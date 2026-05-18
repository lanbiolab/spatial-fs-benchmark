# Reproducibility Notes

## Workflow

The benchmark is implemented as a Python/R script workflow:

1. convert/download datasets into h5ad files under `data/raw/`
2. run feature selection, integration, and task evaluation with
   `scripts/run_benchmark.py`
3. merge metric records into `results/current_rank/data/benchmark.tsv`
4. generate scaling ranges and summary tables
5. render manuscript figures with R scripts

No Snakemake or Nextflow workflow is used in this implementation.

## Main Dataset Set

The final benchmark summary tables include eight datasets:

- DLPFC
- MouseBrainSerialSections
- STOmics0212
- STOmics0218
- STOmics0224
- STOmicsVisium5Samples
- E8p5Embryo
- E9p5Embryo

Alignment metrics require labels and are unavailable for E8p5Embryo and
E9p5Embryo in the current benchmark.

## Main Evaluation Tasks

The benchmark evaluates three task groups:

- Integration
- Clustering
- Alignment

Task-level scores are computed from oriented and scaled metric values. Overall
scores are computed as the equal-weight mean of the three task groups when all
three are available.

## Representative Settings

Family-level ranking figures use representative feature-number settings:

- `all_features`: `N=all`
- `random`: `N=500`
- other feature-selection methods: preferentially `N=2000`, with fallback where
  needed

Feature-number sensitivity figures use setting-level results rather than a
single representative setting.
