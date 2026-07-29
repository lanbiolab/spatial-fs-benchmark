# Manuscript figure contract: rebuild v2

## Scientific thesis

Feature selection is a consequential modelling decision in spatial transcriptomics. Its utility depends jointly on feature budget, downstream task and integration strategy; spatially informed selectors provide context-dependent rather than universal gains.

## Shared visual system

- Backend: R only (`ggplot2`, `patchwork`, `ComplexHeatmap` where needed).
- Final width: 183 mm for all main figures; panel heights are set by evidence density.
- Typeface: Liberation Sans; 7.0 pt body, 8.0 pt axis titles, 9.0 pt panel labels.
- Background: white; no decorative grids, shadows or gradients.
- Method families: expression-based steel blue, spatial warm coral, controls neutral gray, supervised oracle muted violet.
- Tasks: Integration blue, Clustering vermilion, Alignment green, CoreOverall charcoal.
- Integration strategies: scVI navy and CellCharter ochre.
- Primary export: SVG; companion PDF, 600 dpi TIFF and 300 dpi PNG preview.
- Every figure receives a source-data TSV bundle and a rendered-preview inspection.

## Figure 1: benchmark framework

**Claim:** Feature selection changes the representation entering integration and must be evaluated across budgets, tasks and integration strategies.

**Structure:** A dominant left-to-right benchmark schematic, a compact dataset/coverage matrix and a score-construction inset. The schematic starts with multi-slice count matrices and coordinates, branches into expression-based, spatially informed and control selectors, passes selected features into scVI or CellCharter, and ends in Integration, Clustering and Alignment evidence. The score inset separates CoreOverall from the secondary alignment analysis.

**Reviewer risks:** Do not imply that all datasets support alignment; do not present supervised Wilcoxon as a competitive method; do not imply that missing metrics are zeros.

## Figure 2: global benchmark result

**Claim:** Feature selection materially changes performance, but no method wins across every task and integration strategy.

**Structure:** The hero panel shows competitive mean ranks for scVI and CellCharter. Supporting panels show task-specific performance, control/reference behaviour, and the cross-integrator rank association. Methods are ordered once using the primary scVI CoreOverall result.

**Reviewer risks:** Distinguish score from rank, show rank 1 as best, label controls separately and use the frozen representative settings only.

## Figure 3: feature-budget effect

**Claim:** Intermediate feature budgets generally outperform very small and very large sets, with dataset-, task- and method-specific optima.

**Structure:** The hero panel shows aggregate CoreOverall trajectories across feature numbers. Supporting panels split Integration and Clustering, summarize dataset-level effects and show method-specific optimal budgets.

**Reviewer risks:** Use setting-level results rather than representative settings; avoid treating repeated seeds as independent datasets; retain the logarithmic feature-number spacing.

## Figure 4: spatial selector context

**Claim:** Spatially informed selectors can be highly competitive in specific datasets, but their gains depend on the integration strategy and their selected sets are not interchangeable with expression-based sets.

**Structure:** Dataset-level SVG ranks, scVI-to-CellCharter rank shifts, a representative-setting Jaccard matrix and family-pair overlap distributions.

**Reviewer risks:** Avoid claiming universal SVG superiority; compute overlap within datasets before summarizing; keep all-feature controls out of pairwise Jaccard summaries.

## Figure 5: independent validation

**Claim:** Reserved for semi-synthetic truth recovery, held-out-slice validation and independently sourced marker evidence.

**Status:** Not drawn until those result tables are complete. No proxy benchmark result will be relabelled as biological truth.

## Figure 6: evidence-backed guidance

**Claim:** Method choice should follow task, feature budget and integration strategy rather than a single global winner.

**Structure:** A compact decision map backed only by audited benchmark summaries, accompanied by sensitivity and caveat annotations. It is a synthesis figure, not a new ranking.

## Extended Data

- Metric behaviour and correlation overview.
- Frozen baseline-range scaling and score construction.
- Complete method-by-dataset rank matrices.
- Complete representative-setting feature overlap.
- Coverage, missingness and sensitivity audits.
