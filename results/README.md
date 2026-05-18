# Results Directory

This directory contains small summary tables and final figure outputs needed to
reproduce the manuscript-level visualisations.

Included:

- `current_rank/data/`: combined benchmark metric tables and metadata
- `current_rank/output/baseline-ranges.tsv`: scaling ranges used by ranking scripts
- `fig4a_spatial_benchmark/figures/`: method ranking summaries
- `fig4b_spatial_benchmark/figures/`: feature-set overlap summaries
- `fig4e_spatial_benchmark/figures/`: batch-aware selector summaries
- `fig5*_spatial_lineages/figures/`: lineage subset summaries
- `fig6_spatial_integration/figures/`: integration-strategy comparison summaries

Excluded:

- per-run embeddings
- selected-feature JSON files for every dataset/method/feature number
- model checkpoints
- full intermediate result directories

Those excluded files are large and should be regenerated from `configs/` and
`scripts/` when a full rerun is required.
