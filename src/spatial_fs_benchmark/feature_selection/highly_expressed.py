from __future__ import annotations

import numpy as np
from scipy import sparse

from spatial_fs_benchmark.data.spatial_object import SpatialDataset
from spatial_fs_benchmark.feature_selection.base import FeatureSelectionResult, FeatureSelector


class HighlyExpressedSelector(FeatureSelector):
    name = "highly_expressed"

    def select(self, dataset: SpatialDataset, n_features: int, random_seed: int = 0) -> FeatureSelectionResult:
        matrix = dataset.adata.X
        if sparse.issparse(matrix):
            scores = np.asarray(matrix.mean(axis=0)).ravel()
        else:
            scores = np.asarray(matrix, dtype=np.float32).mean(axis=0)
        return self._build_result(
            method_name=self.name,
            gene_names=dataset.adata.var_names.to_numpy(),
            scores=np.asarray(scores, dtype=float),
            n_features=n_features,
        )
