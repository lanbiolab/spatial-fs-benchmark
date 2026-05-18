from __future__ import annotations

import numpy as np
from scipy import sparse
from sklearn.neighbors import kneighbors_graph

from spatial_fs_benchmark.data.spatial_object import SpatialDataset
from spatial_fs_benchmark.feature_selection.base import FeatureSelectionResult, FeatureSelector


class SVGSelector(FeatureSelector):
    name = "svg"

    def __init__(self, n_neighbors: int = 8, chunk_size: int = 512, max_cells: int | None = None) -> None:
        self.n_neighbors = n_neighbors
        self.chunk_size = chunk_size
        self.max_cells = max_cells

    def select(self, dataset: SpatialDataset, n_features: int, random_seed: int = 0) -> FeatureSelectionResult:
        if self.max_cells is not None and dataset.n_obs > self.max_cells:
            rng = np.random.default_rng(random_seed)
            keep = np.sort(rng.choice(dataset.n_obs, size=self.max_cells, replace=False))
            coords = dataset.coords[keep]
            matrix = dataset.adata.X[keep]
        else:
            coords = dataset.coords
            matrix = dataset.adata.X
        graph = kneighbors_graph(
            coords,
            n_neighbors=min(self.n_neighbors, max(1, coords.shape[0] - 1)),
            mode="connectivity",
            include_self=False,
        ).tocsr()
        total_weight = float(graph.sum())
        n_obs = coords.shape[0]
        scores = np.zeros(dataset.n_vars, dtype=float)
        for start in range(0, dataset.n_vars, self.chunk_size):
            end = min(start + self.chunk_size, dataset.n_vars)
            chunk = matrix[:, start:end]
            if sparse.issparse(chunk):
                chunk = chunk.toarray()
            else:
                chunk = np.asarray(chunk, dtype=np.float32)
            centered = chunk - chunk.mean(axis=0, keepdims=True)
            numerator = np.asarray(centered * (graph @ centered)).sum(axis=0)
            denominator = np.square(centered).sum(axis=0) + 1e-12
            scores[start:end] = (n_obs / total_weight) * (numerator / denominator)
        scores = np.nan_to_num(scores, nan=0.0, posinf=0.0, neginf=0.0)
        return self._build_result(
            method_name=self.name,
            gene_names=np.asarray(dataset.gene_names),
            scores=scores,
            n_features=n_features,
            metadata={"n_neighbors": self.n_neighbors, "max_cells": self.max_cells},
        )
