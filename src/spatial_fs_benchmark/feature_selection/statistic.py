from __future__ import annotations

import numpy as np
import scanpy as sc
from scipy import sparse

from spatial_fs_benchmark.data.spatial_object import SpatialDataset
from spatial_fs_benchmark.feature_selection.base import FeatureSelectionResult, FeatureSelector


class StatisticSelector(FeatureSelector):
    name = "statistic"

    def __init__(self, statistic: str = "mean") -> None:
        if statistic not in {"mean", "variance"}:
            raise ValueError(f"Unsupported statistic: {statistic}")
        self.statistic = statistic

    def select(self, dataset: SpatialDataset, n_features: int, random_seed: int = 0) -> FeatureSelectionResult:
        adata = self._adata_from_counts(dataset)
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

        if self.statistic == "mean":
            scores = np.asarray(adata.X.mean(axis=0)).ravel()
        else:
            matrix = adata.X
            if sparse.issparse(matrix):
                squared = matrix.copy()
                squared.data **= 2
                scores = np.asarray(squared.mean(axis=0) - np.square(matrix.mean(axis=0))).ravel()
            else:
                scores = np.asarray(np.var(np.asarray(matrix), axis=0), dtype=float)

        return self._build_result(
            method_name=self.name,
            gene_names=adata.var_names.to_numpy(),
            scores=np.asarray(scores, dtype=float),
            n_features=n_features,
            metadata={"statistic": self.statistic},
        )
