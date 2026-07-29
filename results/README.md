# Tracked Result Summaries

The public repository tracks only compact manuscript-level tables. Raw matrices,
selected-feature caches, embeddings, model checkpoints, and per-run task payloads
are intentionally excluded.

`v2/frozen_scores/` contains the complete frozen primary score tables for scVI
and CellCharter. `v2/graphst_summary/` contains the audited canonical GraphST
sensitivity tables. `v2/source_data/` contains the plotted values, rank and
bootstrap summaries, method/dataset inventories, held-out cross-statistic
validation, hybrid-control comparisons, and the standardized feature-selector
resource profile supplied with the manuscript.

The manifest in `v2/frozen_scores/frozen_score_manifest.json` records the source
checksum and scoring contract. Missing values remain missing and are never
encoded as numerical zeros.
