# Spatial Feature-Selection Benchmark

This repository contains the audited workflow for benchmarking feature-selection objectives in spatial transcriptomics. It compares expression-driven selectors, spatially informed selectors, controls, and a label-informed oracle across integration, clustering, alignment, spatial-truth recovery, and held-out-slice validation.

## Benchmark design

- Seven count-validated public datasets: DLPFC, Mouse Brain Serial Sections, STOmics0212, STOmics0218, STOmics0224, E8.5 embryo, and E9.5 embryo.
- Thirty-four representative methods or reference settings, including five spatially informed selectors.
- Feature budgets of 100, 200, 500, 1,000, 2,000, 5,000, and 10,000 genes where supported.
- scVI and CellCharter as primary downstream representations, with a canonical GraphST sensitivity analysis.
- Thirteen metrics grouped into integration, clustering, and alignment endpoints.
- Semi-synthetic spatial truth, six held-out-slice folds, score sensitivity, dataset bootstrap, and exploratory HVG-SVG controls.

The frozen protocol and scoring definitions are documented in [`docs/rebuild_v1_protocol.md`](docs/rebuild_v1_protocol.md) and [`docs/rebuild_v1_scoring.md`](docs/rebuild_v1_scoring.md). Missing metrics are retained as missing and are never replaced by zero.

## Repository layout

```text
configs/                 Dataset and frozen experiment configurations
docs/                    Protocol, metric, scoring, and figure contracts
external/                Provenance records and vendored method interfaces
results/v2/              Frozen score summaries and manuscript source data
scripts/                 Experiment, audit, summary, and R figure scripts
src/spatial_fs_benchmark Reusable Python package
tests/                   Unit and regression tests
```

Raw data and full embedding outputs are not tracked in Git. Dataset accessions and preprocessing contracts are recorded in the dataset configs and manuscript source-data tables.

## Installation

```bash
conda env create -f environment.yml
conda activate spatial-fs-benchmark
pip install -e .
pytest -q
```

Several selectors use isolated R or method-specific environments. Their provenance and versions are listed in `external/svg-method-references/PROVENANCE.tsv` and the environment files at repository root.

## Core workflow

Generate and run the rebuilt feature-selection and downstream configurations:

```bash
python scripts/generate_svg_rebuild_configs.py
python scripts/generate_svg_downstream_configs.py
python scripts/generate_nonspatial_downstream_configs.py
spatial-fs-run --config configs/rebuild_v1/downstream/canonical/dlpfc.yaml
```

Rebuild metric tables from immutable task records and construct the frozen scores:

```bash
python scripts/rebuild_results_from_records.py \
  --results-root results/spatial_svg_rebuild_v1 \
  --include-glob '*' \
  --merged-output results/spatial_svg_rebuild_v1/merged_results.csv

python scripts/build_frozen_scores.py
```

Validation and sensitivity entry points include:

```bash
python scripts/run_heldout_slice_validation.py --help
python scripts/build_submission_robustness_tables.py
python scripts/summarize_m3_m7_validations.py
python scripts/summarize_graphst_canonical.py
```

Publication figures are generated from saved TSV tables using the R scripts under `scripts/`. Each plotting script exports editable PDF/SVG, 600-dpi TIFF, and source-data tables.

Dataset accessions and benchmark-level spot/gene counts are listed in
[`results/v2/source_data/datasets_inventory.tsv`](results/v2/source_data/datasets_inventory.tsv).
Raw data should be downloaded from the original public repositories and placed
according to [`data/README.md`](data/README.md).

## Citation

Citation metadata are provided in [`CITATION.cff`](CITATION.cff). A versioned archival DOI will be added after the corresponding GitHub release is deposited with Zenodo.

## License

MIT. Third-party methods retain their original licenses; consult their repositories and the provenance table before redistribution.
