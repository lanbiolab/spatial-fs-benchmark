from __future__ import annotations

import numpy as np
from scipy import sparse
from scipy.stats import rankdata
from sklearn.neighbors import kneighbors_graph

from spatial_fs_benchmark.data.spatial_object import SpatialDataset
from spatial_fs_benchmark.feature_selection.base import FeatureSelectionResult, FeatureSelector


class SVGSelector(FeatureSelector):
    """Rank genes by slice-wise Moran's I and aggregate ranks across slices."""

    name = "morans_i"
    implementation_version = "v2_slice_wise_rank_aggregation"

    def __init__(self, n_neighbors: int = 8, chunk_size: int = 512, max_cells: int | None = None) -> None:
        self.n_neighbors = n_neighbors
        self.chunk_size = chunk_size
        self.max_cells = max_cells
        self.stochastic_selection = max_cells is not None
        self._score_cache: dict[tuple[str, int], np.ndarray] = {}

    @staticmethod
    def _symmetric_spatial_graph(coords: np.ndarray, n_neighbors: int) -> sparse.csr_matrix:
        graph = kneighbors_graph(
            coords,
            n_neighbors=min(n_neighbors, max(1, coords.shape[0] - 1)),
            mode="connectivity",
            include_self=False,
        ).tocsr()
        graph = graph.maximum(graph.transpose()).tocsr()
        graph.setdiag(0)
        graph.eliminate_zeros()
        return graph

    def _moran_scores(self, matrix, coords: np.ndarray) -> np.ndarray:
        graph = self._symmetric_spatial_graph(coords, self.n_neighbors)
        total_weight = float(graph.sum())
        n_obs, n_vars = matrix.shape
        scores = np.full(n_vars, -np.inf, dtype=float)
        if total_weight <= 0 or n_obs < 3:
            return scores

        for start in range(0, n_vars, self.chunk_size):
            end = min(start + self.chunk_size, n_vars)
            chunk = matrix[:, start:end]
            if sparse.issparse(chunk):
                chunk = chunk.toarray()
            else:
                chunk = np.asarray(chunk)
            chunk = chunk.astype(np.float32, copy=False)
            centered = chunk - chunk.mean(axis=0, keepdims=True)
            denominator = np.square(centered).sum(axis=0)
            valid = denominator > 1e-12
            if np.any(valid):
                numerator = np.asarray(centered * (graph @ centered)).sum(axis=0)
                values = (n_obs / total_weight) * (numerator / np.maximum(denominator, 1e-12))
                scores[start:end][valid] = values[valid]
        return scores

    @staticmethod
    def _percentile_ranks(scores: np.ndarray) -> np.ndarray:
        finite = np.isfinite(scores)
        percentiles = np.zeros(scores.shape[0], dtype=float)
        if np.any(finite):
            percentiles[finite] = rankdata(scores[finite], method="average") / finite.sum()
        return percentiles

    def _select_slice_indices(
        self,
        indices: np.ndarray,
        random_seed: int,
        slice_number: int,
        n_slices: int,
    ) -> np.ndarray:
        if self.max_cells is None:
            return indices
        per_slice = max(3, int(np.ceil(self.max_cells / max(1, n_slices))))
        if len(indices) <= per_slice:
            return indices
        rng = np.random.default_rng(random_seed + slice_number * 1009)
        return np.sort(rng.choice(indices, size=per_slice, replace=False))

    def _aggregate_scores(self, dataset: SpatialDataset, random_seed: int) -> np.ndarray:
        slice_ids = dataset.slice_ids
        unique_slices = np.unique(slice_ids)
        slice_ranks = []
        for slice_number, slice_id in enumerate(unique_slices):
            indices = np.flatnonzero(slice_ids == slice_id)
            indices = self._select_slice_indices(indices, random_seed, slice_number, len(unique_slices))
            if len(indices) < 3:
                continue
            matrix = dataset.adata.X[indices]
            coords = dataset.coords[indices]
            slice_ranks.append(self._percentile_ranks(self._moran_scores(matrix, coords)))
        if not slice_ranks:
            raise ValueError(f"Dataset '{dataset.name}' has no slice with at least three spatial observations.")
        return np.mean(np.vstack(slice_ranks), axis=0)

    def select(self, dataset: SpatialDataset, n_features: int, random_seed: int = 0) -> FeatureSelectionResult:
        cache_key = (dataset.source_path or dataset.name, random_seed)
        if cache_key not in self._score_cache:
            self._score_cache[cache_key] = self._aggregate_scores(dataset, random_seed)
        scores = self._score_cache[cache_key]
        return self._build_result(
            method_name=self.name,
            gene_names=np.asarray(dataset.gene_names),
            scores=scores,
            n_features=n_features,
            metadata={
                "n_neighbors": self.n_neighbors,
                "max_cells": self.max_cells,
                "aggregation": "mean within-slice percentile rank",
                "n_slices": int(len(np.unique(dataset.slice_ids))),
            },
        )
