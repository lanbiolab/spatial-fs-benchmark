from __future__ import annotations

import numpy as np
import scanpy as sc

from spatial_fs_benchmark.data.spatial_object import SpatialDataset
from spatial_fs_benchmark.feature_selection.base import FeatureSelectionResult, FeatureSelector


class HVGSelector(FeatureSelector):
    name = "hvg"

    def __init__(self, flavor: str = "seurat_v3") -> None:
        self.flavor = flavor

    def select(self, dataset: SpatialDataset, n_features: int, random_seed: int = 0) -> FeatureSelectionResult:
        adata = dataset.adata.copy()
        sc.pp.highly_variable_genes(
            adata,
            flavor=self.flavor,
            n_top_genes=min(n_features, adata.n_vars),
            batch_key=dataset.slice_key,
            subset=False,
        )
        scores = adata.var.get("dispersions_norm")
        if scores is None:
            scores = adata.var["highly_variable_rank"].max() - adata.var["highly_variable_rank"]
        score_array = np.nan_to_num(np.asarray(scores, dtype=float), nan=0.0)
        return self._build_result(
            method_name=self.name,
            gene_names=adata.var_names.to_numpy(),
            scores=score_array,
            n_features=n_features,
            metadata={"flavor": self.flavor},
        )
