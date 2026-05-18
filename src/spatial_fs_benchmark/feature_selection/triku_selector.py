from __future__ import annotations

import numpy as np
import scanpy as sc
import triku as tk

from spatial_fs_benchmark.data.spatial_object import SpatialDataset
from spatial_fs_benchmark.feature_selection.base import FeatureSelectionResult, FeatureSelector


class TrikuSelector(FeatureSelector):
    name = "triku"

    def __init__(self, min_genes: int = 50, min_cells: int = 10, max_cells: int | None = None) -> None:
        self.min_genes = min_genes
        self.min_cells = min_cells
        self.max_cells = max_cells

    def select(self, dataset: SpatialDataset, n_features: int, random_seed: int = 0) -> FeatureSelectionResult:
        gene_names = dataset.adata.var_names.to_numpy()
        score_map = {gene: 0.0 for gene in gene_names}

        adata = self._maybe_subsample_obs(self._adata_from_counts(dataset), self.max_cells, random_seed)
        sc.pp.filter_cells(adata, min_genes=self.min_genes)
        sc.pp.filter_genes(adata, min_cells=self.min_cells)
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        sc.pp.pca(adata)
        sc.pp.neighbors(adata, metric="cosine", n_neighbors=int(0.5 * max(2, len(adata)) ** 0.5))
        tk.tl.triku(adata, verbose="error")

        if "triku_distance" in adata.var:
            raw_scores = np.asarray(adata.var["triku_distance"], dtype=float)
        else:
            raw_scores = np.asarray(adata.var["highly_variable"].astype(float), dtype=float)
        raw_scores = np.nan_to_num(raw_scores, nan=0.0, posinf=0.0, neginf=0.0)
        for gene_name, score in zip(adata.var_names, raw_scores, strict=True):
            score_map[str(gene_name)] = float(score)

        scores = np.asarray([score_map[gene] for gene in gene_names], dtype=float)
        return self._build_result(
            method_name=self.name,
            gene_names=gene_names,
            scores=scores,
            n_features=n_features,
            metadata={"max_cells": self.max_cells, "min_genes": self.min_genes, "min_cells": self.min_cells},
        )
