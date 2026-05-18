from __future__ import annotations

from hashlib import md5

import numpy as np
import torch
from sklearn.decomposition import PCA

from spatial_fs_benchmark.data.spatial_object import SpatialDataset
from spatial_fs_benchmark.feature_selection.base import FeatureSelectionResult
from spatial_fs_benchmark.integration.base import IntegrationResult, SpatialIntegrator


class GPSAIntegrator(SpatialIntegrator):
    name = "gpsa"
    implementation_version = "v2_use_log_pca_input"

    def __init__(
        self,
        n_input_dims: int = 20,
        m_x_per_view: int = 50,
        m_g: int = 50,
        max_epochs: int = 200,
        learning_rate: float = 1e-2,
        fixed_view_idx: int = 0,
        device: str = "cpu",
    ) -> None:
        self.n_input_dims = n_input_dims
        self.m_x_per_view = m_x_per_view
        self.m_g = m_g
        self.max_epochs = max_epochs
        self.learning_rate = learning_rate
        self.fixed_view_idx = fixed_view_idx
        self.device = device

    def fit_transform(
        self,
        dataset: SpatialDataset,
        selected_features: FeatureSelectionResult,
        random_seed: int = 0,
    ) -> IntegrationResult:
        import gpsa.models.gpsa as gpsa_base
        import gpsa.models.vgpsa as gpsa_variational
        from gpsa import VariationalGPSA

        subset = dataset.subset_features(selected_features.feature_names)
        adata = subset.adata.copy()
        # Use the dataset's preprocessed expression matrix as GPSA input.
        # In this benchmark, X already reflects the configured normalize_total/log1p
        # pipeline, which preserves biological structure better than raw-count PCA.
        matrix = adata.X
        if hasattr(matrix, "toarray"):
            matrix = matrix.toarray()
        matrix = np.asarray(matrix, dtype=np.float32)
        coords = np.asarray(adata.obsm[subset.coord_key], dtype=np.float32)
        slice_ids = adata.obs[subset.slice_key].astype(str).to_numpy()
        slice_names = np.unique(slice_ids).tolist()
        n_samples_list = [int(np.sum(slice_ids == slice_name)) for slice_name in slice_names]

        if matrix.min() < 0:
            matrix = matrix - float(matrix.min())
        if matrix.shape[1] > self.n_input_dims:
            matrix = PCA(
                n_components=min(self.n_input_dims, matrix.shape[0] - 1, matrix.shape[1]),
                random_state=random_seed,
            ).fit_transform(matrix).astype(np.float32)

        # GPSA expects one modality with per-view sample counts.
        device = torch.device(self.device if self.device != "auto" else ("cuda" if torch.cuda.is_available() else "cpu"))
        # GPSA stores the default device in module-level globals. Force it to match this run.
        gpsa_base.device = str(device)
        gpsa_variational.device = str(device)
        data_dict = {
            "expression": {
                "spatial_coords": torch.tensor(coords, dtype=torch.float32, device=device),
                "outputs": torch.tensor(matrix, dtype=torch.float32, device=device),
                "n_samples_list": n_samples_list,
            }
        }
        model = VariationalGPSA(
            data_dict,
            m_X_per_view=min(self.m_x_per_view, min(n_samples_list)),
            m_G=min(self.m_g, coords.shape[0]),
            n_spatial_dims=coords.shape[1],
            fixed_view_idx=self.fixed_view_idx,
            n_latent_gps={"expression": None},
        ).to(device)
        optimizer = torch.optim.Adam(model.parameters(), lr=self.learning_rate)
        X_spatial = {"expression": data_dict["expression"]["spatial_coords"]}

        torch.manual_seed(random_seed)
        np.random.seed(random_seed)
        model.train()
        for _ in range(self.max_epochs):
            optimizer.zero_grad()
            G_means, _, _, F_obs = model.forward(X_spatial, model.view_idx, model.Ns, S=1)
            loss = model.loss_fn(data_dict, F_obs)
            loss.backward()
            optimizer.step()

        model.eval()
        with torch.no_grad():
            G_means, _, _, F_obs = model.forward(X_spatial, model.view_idx, model.Ns, S=1)
        aligned_coords = G_means["expression"].detach().cpu().numpy().astype(np.float32)
        latent_outputs = F_obs["expression"].squeeze(0).detach().cpu().numpy().astype(np.float32)
        embedding = latent_outputs
        return IntegrationResult(
            method_name=self.name,
            embedding=embedding,
            slice_embedding=self._build_slice_embedding(dataset, embedding),
            metadata={
                "aligned_coords": aligned_coords,
                "n_input_features": int(adata.n_vars),
                "input_feature_hash": md5("\n".join(selected_features.feature_names).encode("utf-8")).hexdigest(),
                "n_input_dims": matrix.shape[1],
                "input_representation": "preprocessed_X_pca",
                "m_x_per_view": min(self.m_x_per_view, min(n_samples_list)),
                "m_g": min(self.m_g, coords.shape[0]),
                "max_epochs": self.max_epochs,
                "learning_rate": self.learning_rate,
                "fixed_view_idx": self.fixed_view_idx,
                "device": str(device),
                "representation_type": "latent",
            },
        )
