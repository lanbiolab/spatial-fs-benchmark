# Spatial feature-selection benchmark rebuild v1

## Primary question

The benchmark compares expression-driven, spatially informed, label-informed, predefined, and control feature sets. The primary question is whether spatially informed selectors provide task-specific advantages over commonly used scRNA-seq feature-selection procedures in spatial transcriptomics workflows.

## Method groups

- Spatially informed: Moran's I, SPARK-X, nnSVG, SpatialDE, SOMDE.
- Expression-driven: Scanpy/Seurat HVG variants and the existing unsupervised selector panel.
- Label-informed oracle: Wilcoxon. This method is not treated as an unsupervised competitor.
- Controls: all features, random features, transcription factors, and stable-gene controls.
- scGCO: implementation pilot only. It enters the main benchmark only if throughput and validation checks pass.

All spatial methods operate on each slice independently. Slice-level gene rankings are converted to percentile ranks and averaged across slices. Coordinates from unregistered slices are never placed in one spatial graph.

## Input contracts

- Count-based methods receive non-negative integer values from `layers["counts"]`.
- If `counts` is absent, `raw_count` is preferred over `X`.
- scVI, CellCharter, SPARK-X, SOMDE, scGCO, Pearson residuals, and Seurat v3 fail explicitly when integer counts are unavailable.
- Expression-based spatial methods receive the dataset's normalized/log-transformed expression representation.
- The count source and count validation result are recorded in every dataset signature.

STOmicsVisium5Samples is temporarily excluded from count-dependent methods because the available `raw_count` layer is transformed rather than integer-valued. It is restored only after genuine raw counts are obtained.

## Feature-number design

- General fixed settings: 100, 200, 500, 1,000, 2,000, 5,000, and 10,000.
- Moran's I and SPARK-X: all fixed settings.
- nnSVG and SpatialDE: up to 2,000 after an explicitly recorded candidate-gene prefilter.
- SOMDE: up to 5,000.
- Canonical method comparison: 2,000 features when available; method-specific fallback is reported, never silent.
- Random control: 500 features.
- All-features control: all available genes.

Methods that cannot produce a requested number are recorded as unavailable. Unranked genes are not used to pad a feature set.

## Repetition

- Integration and clustering: seeds 0, 1, and 2 for canonical settings.
- Feature-number screening: seed 0, with repeated canonical settings used to quantify stochastic variation.
- Random feature sets and selectors that subsample observations use seed-specific feature sets.
- Deterministic rankings are computed once and reused across integration seeds.

## Evaluation tasks

### Integration

- Batch mixing: bASW and iLISI, available with or without biological labels.
- Biological conservation: dASW, dLISI, ILL, and GC, evaluated only where benchmark labels exist.
- Mixing and conservation are summarized separately before receiving equal weight in a combined integration score.

### Clustering

- Label-independent: Silhouette, CHAOS, and PAS.
- Label agreement: ARI and NMI where benchmark labels exist.
- The two clustering components are summarized separately before equal-weight combination.

### Spatial alignment

- Alignment is evaluated only for explicitly declared serial-section pairs.
- DLPFC pairs: 151507--151508 and 151673--151674.
- MouseBrainSerialSections pair: section1--section2.
- Independent samples or conditions are not treated as adjacent slices.
- Embedding-neighbor correspondence is described as such and is not presented as native coordinate registration.

## Aggregate scores

- `CoreOverall`: equal-weight mean of Integration and Clustering, available across the main dataset collection.
- `Alignment`: reported separately on alignment-eligible datasets.
- `AlignmentEligibleOverall`: equal-weight mean of Integration, Clustering, and Alignment, reported only as a secondary sensitivity summary on explicitly paired serial sections.
- Missing metrics are never replaced by zero.

## Scaling

Metric direction is oriented before scaling. A single frozen Dataset x Metric range table is calculated from seed-averaged setting-level results after the run matrix is complete. The same table is reused by all main and supplemental figures. Figure-specific baseline or min-max ranges are prohibited. Feature-number effect standardization, when used, is explicitly labeled as a secondary within-method analysis.

## Independent validation

- Semi-synthetic spatial datasets contain known domain, gradient, focal, and periodic SVG patterns across effect sizes and prevalence levels.
- SVG recovery is evaluated by AUROC, AUPRC, precision/recall at fixed N, and pattern-specific recall.
- Held-out-slice analysis selects features on training slices and evaluates spatial reproducibility and downstream performance on held-out slices.
- Marker analyses use external or original-study marker resources. Dataset-derived Wilcoxon markers are retained only as internal-consistency analyses.

## Execution phases

1. Generate and validate spatial feature rankings.
2. Recompute count-affected feature rankings and scVI results.
3. Run canonical settings with three seeds.
4. Run seed-0 feature-number screening.
5. Run workflow sensitivity on a prespecified representative selector panel.
6. Run semi-synthetic and held-out-slice validation.
7. Freeze ranges, build audited summary tables, and regenerate every figure.
8. Rewrite Methods and Results from machine-readable coverage and provenance tables.

Each phase writes to `results/spatial_svg_rebuild_v1` and must pass coverage, count-source, feature-hash, missingness, and duplicate-row audits before entering the next phase.
