from __future__ import annotations

from hashlib import md5

import numpy as np
import scvi

from spatial_fs_benchmark.data.spatial_object import SpatialDataset
from spatial_fs_benchmark.feature_selection.base import FeatureSelectionResult
from spatial_fs_benchmark.integration.base import IntegrationResult, SpatialIntegrator


class SCVIIntegrator(SpatialIntegrator):
    name = "scvi"

    def __init__(
        self,
        n_latent: int = 30,
        max_epochs: int = 100,
        batch_size: int = 2048,
    ) -> None:
        self.n_latent = n_latent
        self.max_epochs = max_epochs
        self.batch_size = batch_size

    def fit_transform(
        self,
        dataset: SpatialDataset,
        selected_features: FeatureSelectionResult,
        random_seed: int = 0,
    ) -> IntegrationResult:
        subset = dataset.subset_features(selected_features.feature_names)
        adata = subset.adata.copy()
        if "counts" in adata.layers:
            adata.X = adata.layers["counts"].copy()
            adata.uns.pop("log1p", None)

        matrix = adata.X
        if hasattr(matrix, "toarray"):
            min_value = float(matrix.min())
        else:
            min_value = float(np.min(np.asarray(matrix)))
        if min_value < 0:
            raise ValueError(
                f"scVI requires non-negative count-like input, but dataset '{dataset.name}' has negative values."
            )

        adata.obs["Batch"] = adata.obs[subset.slice_key].astype(str)
        scvi.settings.seed = random_seed
        scvi.model.SCVI.setup_anndata(adata, batch_key="Batch")
        model = scvi.model.SCVI(
            adata,
            n_latent=min(self.n_latent, max(2, adata.n_vars)),
            use_layer_norm="both",
            use_batch_norm="none",
            encode_covariates=True,
            dropout_rate=0.2,
            n_layers=2,
        )
        model.train(max_epochs=self.max_epochs, batch_size=self.batch_size, check_val_every_n_epoch=None)
        embedding = model.get_latent_representation()
        return IntegrationResult(
            method_name=self.name,
            embedding=embedding,
            slice_embedding=self._build_slice_embedding(dataset, embedding),
            metadata={
                "n_input_features": int(adata.n_vars),
                "input_feature_hash": md5("\n".join(selected_features.feature_names).encode("utf-8")).hexdigest(),
                "n_latent": int(embedding.shape[1]),
                "max_epochs": self.max_epochs,
                "batch_size": self.batch_size,
            },
        )
