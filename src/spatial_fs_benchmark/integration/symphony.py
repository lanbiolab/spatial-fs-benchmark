from __future__ import annotations

from hashlib import md5

import numpy as np
import scanpy as sc
from symphonypy.preprocessing import harmony_integrate

from spatial_fs_benchmark.data.spatial_object import SpatialDataset
from spatial_fs_benchmark.feature_selection.base import FeatureSelectionResult
from spatial_fs_benchmark.integration.base import IntegrationResult, SpatialIntegrator


class SymphonyIntegrator(SpatialIntegrator):
    name = "symphony"

    def __init__(self, n_components: int = 30) -> None:
        self.n_components = n_components

    def fit_transform(
        self,
        dataset: SpatialDataset,
        selected_features: FeatureSelectionResult,
        random_seed: int = 0,
    ) -> IntegrationResult:
        subset = dataset.subset_features(selected_features.feature_names)
        adata = subset.adata.copy()
        adata.obs["Batch"] = adata.obs[subset.slice_key].astype(str)

        matrix = adata.X
        if hasattr(matrix, "toarray"):
            min_value = float(matrix.min())
        else:
            min_value = float(np.min(np.asarray(matrix)))
        if min_value >= 0:
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
        else:
            adata.uns.pop("log1p", None)

        sc.pp.scale(adata, max_value=10)
        sc.tl.pca(adata, n_comps=min(self.n_components, adata.n_vars, adata.n_obs - 1), zero_center=False)
        harmony_integrate(adata, key="Batch", ref_basis_adjusted="X_emb", verbose=False, random_seed=random_seed)
        embedding = np.asarray(adata.obsm["X_emb"], dtype=float)

        return IntegrationResult(
            method_name=self.name,
            embedding=embedding,
            slice_embedding=self._build_slice_embedding(dataset, embedding),
            metadata={
                "n_input_features": int(adata.n_vars),
                "input_feature_hash": md5("\n".join(selected_features.feature_names).encode("utf-8")).hexdigest(),
                "n_components": int(embedding.shape[1]),
                "input_mode": "normalized_log1p" if min_value >= 0 else "pretransformed",
            },
        )
