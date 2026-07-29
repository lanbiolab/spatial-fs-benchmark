from __future__ import annotations

import numpy as np

from spatial_fs_benchmark.data.spatial_object import SpatialDataset
from spatial_fs_benchmark.feature_selection.base import FeatureSelectionResult, FeatureSelector


class RandomSelector(FeatureSelector):
    name = "random"
    stochastic_selection = True

    def select(self, dataset: SpatialDataset, n_features: int, random_seed: int = 0) -> FeatureSelectionResult:
        rng = np.random.default_rng(random_seed)
        scores = rng.random(dataset.n_vars)
        return self._build_result(
            method_name=self.name,
            gene_names=dataset.adata.var_names.to_numpy(),
            scores=scores,
            n_features=n_features,
            metadata={"seed": random_seed},
        )
