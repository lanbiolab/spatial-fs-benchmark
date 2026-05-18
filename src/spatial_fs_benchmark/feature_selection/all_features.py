from __future__ import annotations

import numpy as np

from spatial_fs_benchmark.data.spatial_object import SpatialDataset
from spatial_fs_benchmark.feature_selection.base import FeatureSelectionResult, FeatureSelector


class AllFeaturesSelector(FeatureSelector):
    name = "all_features"

    def select(self, dataset: SpatialDataset, n_features: int, random_seed: int = 0) -> FeatureSelectionResult:
        scores = np.arange(dataset.n_vars, 0, -1, dtype=float)
        return self._build_result(
            method_name=self.name,
            gene_names=dataset.adata.var_names.to_numpy(),
            scores=scores,
            n_features=dataset.n_vars,
            metadata={"requested_n_features_ignored": int(n_features)},
        )
