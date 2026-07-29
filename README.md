# Spatial Feature-Selection Benchmark

This repository contains the code, frozen configurations, audited summary
tables, and figure-generation scripts for benchmarking feature selection in
multi-slice spatial transcriptomics.

The rebuilt benchmark separates two questions:

1. Does a selector recover reproducible spatial expression signal?
2. Does the selected representation support integration, clustering, and
   cross-section correspondence?

These objectives are evaluated separately because spatial signal recovery does
not imply uniformly better downstream performance.

## Frozen v2 Design

The manuscript analysis uses:

- seven public count-validated datasets: DLPFC, MouseBrainSerialSections,
  STOmics0212, STOmics0218, STOmics0224, E8p5Embryo, and E9p5Embryo;
- 34 representative methods or references, including 29 competitive
  unsupervised selectors and five coordinate-aware spatial selectors;
- seven requested feature budgets: 100, 200, 500, 1,000, 2,000, 5,000, and
  10,000 genes;
- two downstream representation strategies: scVI and CellCharter;
- 13 integration, clustering, and alignment metrics;
- three seeds for canonical 2,000-gene settings;
- semi-synthetic spatial truth, six held-out-slice folds, biological-context
  sensitivity, score-construction sensitivity, and dataset bootstrap analyses.

The primary score is `CoreOverall`, the equal-weight mean of Integration and
Clustering. Alignment is reported separately because it is defined only for
explicit DLPFC and Mouse Brain section pairs. Missing metrics are stored as
missing and are never replaced by zero.

The frozen scoring protocol is documented in
[`docs/rebuild_v1_protocol.md`](docs/rebuild_v1_protocol.md). Machine-readable
score tables and the input checksum are under
[`results/v2/frozen_scores/`](results/v2/frozen_scores/).

## Method Classification

The coordinate-aware selectors are Moran's I, SPARK-X, nnSVG, SpatialDE, and
SOMDE. Hotspot is classified as expression-driven in this benchmark because its
30-nearest-neighbor graph is constructed in PCA expression space rather than
from tissue coordinates. Wilcoxon uses benchmark labels and is retained only as
a label-informed oracle; it is missing for E8.5 and E9.5, which have no
benchmark labels.

The complete method, input, version, and classification table is available at
[`results/v2/source_data/methods_inventory.tsv`](results/v2/source_data/methods_inventory.tsv).

## Repository Structure

```text
analysis/       R helpers for result summarization and plotting
configs/        dataset and frozen rebuild configurations
data/           small metadata/resources; raw h5ad files are not tracked
docs/           protocol, metric, and reproducibility documentation
external/       provenance and minimal referenced helper code
results/v2/     audited frozen scores and manuscript source-data tables
scripts/        conversion, execution, audit, scoring, and R plotting scripts
src/            spatial_fs_benchmark Python package
tests/          unit and smoke tests
```

Raw datasets, model embeddings, environments, temporary directories, and full
run logs are excluded because of size or source-data licensing constraints.

## Environment

```bash
conda env create -f environment.yml
conda activate spatial-fs-benchmark
pip install -e .[dev]
pytest -q
```

The frozen environment used Python 3.10, Scanpy 1.10.4, scvi-tools
1.1.6.post2, and CellCharter 0.3.5. R and method-specific versions are recorded
in the method inventory and SVG provenance table.

## Reproducing the Rebuild

Dataset-level execution is controlled by YAML rather than Snakemake or
Nextflow. Example:

```bash
spatial-fs-run --config configs/rebuild_v1/downstream/canonical/dlpfc.yaml
```

After all dataset records are available:

```bash
python scripts/rebuild_results_from_records.py \
  --results-root results/spatial_svg_rebuild_v1 \
  --include-glob '*' \
  --merged-output results/spatial_svg_rebuild_v1/merged_results.csv

python scripts/audit_svg_downstream_results.py
python scripts/build_frozen_scores.py
python scripts/build_submission_robustness_tables.py
```

Official SVG wrapper provenance is recorded in
[`external/svg-method-references/PROVENANCE.tsv`](external/svg-method-references/PROVENANCE.tsv).

## Results and Figures

The repository includes saved source-data tables rather than large
intermediate embeddings. Main and supplemental plots are regenerated from
these tables with R scripts. The score-construction and bootstrap robustness
analysis is generated with:

```bash
Rscript scripts/generate_supplemental_robustness_figure.R \
  results/submission_robustness_v1 \
  results/submission_robustness_v1/figure
```

## Data Availability

Dataset accessions and benchmark-level spot/gene counts are listed in
[`results/v2/source_data/datasets_inventory.tsv`](results/v2/source_data/datasets_inventory.tsv).
Raw data should be downloaded from the original public repositories and placed
according to [`data/README.md`](data/README.md).
