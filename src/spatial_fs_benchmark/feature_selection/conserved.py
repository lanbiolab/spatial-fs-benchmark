from __future__ import annotations

import numpy as np
import pandas as pd

from spatial_fs_benchmark.data.spatial_object import SpatialDataset
from spatial_fs_benchmark.feature_selection.base import FeatureSelectionResult, FeatureSelector


class CrossSliceConservedSelector(FeatureSelector):
    name = "conserved"

    def select(self, dataset: SpatialDataset, n_features: int, random_seed: int = 0) -> FeatureSelectionResult:
        matrix = dataset.to_dense_matrix()
        slice_ids = dataset.slice_ids
        frame = pd.DataFrame(matrix, columns=dataset.gene_names)
        frame["slice_id"] = slice_ids
        slice_means = frame.groupby("slice_id", observed=True).mean(numeric_only=True)
        mean_expression = slice_means.mean(axis=0).to_numpy()
        std_expression = slice_means.std(axis=0).to_numpy()
        cv = std_expression / (mean_expression + 1e-8)
        scores = mean_expression / (1.0 + cv)
        return self._build_result(
            method_name=self.name,
            gene_names=np.asarray(dataset.gene_names),
            scores=scores,
            n_features=n_features,
        )
