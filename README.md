# Spatial Feature-Selection Benchmark

This repository contains the core code and summary outputs for benchmarking how
feature selection affects multi-slice spatial transcriptomics integration and
downstream analysis.

The project adapts the benchmark design philosophy of the scRNA-seq feature
selection benchmark to spatial transcriptomics. It evaluates feature selectors
through a common pipeline:

```text
spatial dataset -> feature selection -> integration -> evaluation tasks -> ranking/figures
```

## Directory Structure

```text
analysis/          R helper functions for result summarisation and plotting
configs/           Benchmark, dataset, and method configuration files
data/              Small benchmark resources only; raw h5ad files are not tracked
docs/              Notes on methods, metrics, plotting reuse, and reproducibility
external/          Minimal vendored plotting helper from the reference benchmark
results/           Small summary tables and final plotted outputs
scripts/           Data conversion, benchmark execution, summarisation, and plotting scripts
src/               Python package implementing the benchmark framework
tests/             Smoke tests for core registries and benchmark components
```

Large raw datasets, intermediate embeddings, model checkpoints, conda
environments, temporary CellCharter directories, and full run logs are excluded
from git.

## Environment

Create the main environment with:

```bash
conda env create -f environment.yml
conda activate spatial-fs-benchmark
pip install -e .[dev]
```

Some integration methods, especially CellCharter, GPSA, STAligner, and scVI, may
require additional GPU-specific dependencies depending on the server.

## Core Benchmark Run

The benchmark is implemented with Python and R scripts rather than Snakemake or
Nextflow. Dataset-level runs are configured in `configs/benchmark/`.

Example:

```bash
python scripts/run_benchmark.py --config configs/benchmark/dlpfc_spatial_main_native.yaml
```

The final main benchmark used eight datasets and four integration strategies:

- `scVI`
- `CellCharter`
- `GPSA`
- `STAligner`

The main ranking figures use representative feature-number settings. Feature
number sensitivity and method-family analyses are generated separately from the
saved benchmark tables.

## Results and Figures

Small summary tables needed for manuscript figures are included under
`results/`. Top-level figure PNG exports are not tracked in git; figures can be
regenerated from the saved summary tables and plotting scripts.

Regenerate representative figures with scripts such as:

```bash
Rscript scripts/generate_metric_overview_figure.R
Rscript scripts/generate_baselines_figure.R
Rscript scripts/generate_num_features_benchmark_figure.R
Rscript scripts/generate_fig4a_spatial_benchmark.R
Rscript scripts/generate_fig4b_spatial_benchmark.R
Rscript scripts/generate_fig4e_spatial_benchmark.R
Rscript scripts/generate_fig5abc_spatial_lineages.R
Rscript scripts/generate_fig5d_spatial_lineages.R
Rscript scripts/generate_fig6_spatial_integration_comparison.R
```

## Data Availability

Raw and processed h5ad datasets are not included in this repository because of
size and licensing constraints. Dataset metadata and benchmark summary tables are
included under `results/current_rank/`.

See `data/README.md` and `docs/reproducibility.md` for the expected local data
layout and known dataset sources.

## Reference

This repository reuses plotting helper code from:

Zappia et al. _Feature selection methods affect the performance of scRNA-seq
data integration and querying._ Nature Methods, 2025.

The vendored helper is kept under `external/atlas-feature-selection-benchmark/`
with its original license.
