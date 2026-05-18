from __future__ import annotations

import sys
from hashlib import md5
from pathlib import Path

import anndata as ad
import numpy as np
import scipy.linalg
import scipy.sparse as sp
import torch

from spatial_fs_benchmark.data.spatial_object import SpatialDataset
from spatial_fs_benchmark.feature_selection.base import FeatureSelectionResult
from spatial_fs_benchmark.integration.base import IntegrationResult, SpatialIntegrator


def _staligner_module_path() -> Path:
    return Path(__file__).resolve().parents[3] / "external" / "STAligner"


class STAlignerIntegrator(SpatialIntegrator):
    name = "staligner"
    implementation_version = "v2_hnsw_safe"

    def __init__(
        self,
        hidden_dims: list[int] | None = None,
        n_epochs: int = 600,
        lr: float = 1e-3,
        knn_neigh: int = 50,
        radius_cutoff: float = 50.0,
        device: str = "cpu",
    ) -> None:
        self.hidden_dims = hidden_dims or [128, 30]
        self.n_epochs = n_epochs
        self.lr = lr
        self.knn_neigh = knn_neigh
        self.radius_cutoff = radius_cutoff
        self.device = device

    def fit_transform(
        self,
        dataset: SpatialDataset,
        selected_features: FeatureSelectionResult,
        random_seed: int = 0,
    ) -> IntegrationResult:
        module_root = _staligner_module_path()
        if str(module_root) not in sys.path:
            sys.path.insert(0, str(module_root))
        try:
            import STAligner
        except Exception as exc:
            raise ImportError(
                "STAligner dependencies are still incomplete in the current environment. "
                "The remaining blocker is the PyG sparse stack (for example torch_sparse)."
            ) from exc

        subset = dataset.subset_features(selected_features.feature_names)
        adata = subset.adata.copy()
        if "counts" in adata.layers:
            adata.X = adata.layers["counts"].copy()
            adata.uns.pop("log1p", None)
        if hasattr(adata.X, "toarray"):
            adata.X = sp.csr_matrix(adata.X)
        else:
            adata.X = sp.csr_matrix(np.asarray(adata.X))

        batches = np.unique(subset.slice_ids).tolist()
        batch_list = []
        adj_list = []
        original_obs_names = adata.obs_names.astype(str).to_numpy()
        original_slice_ids = adata.obs[subset.slice_key].astype(str).to_numpy()
        target_obs_names = np.array(
            [f"{obs}_{slice_id}" for obs, slice_id in zip(original_obs_names, original_slice_ids, strict=True)],
            dtype=object,
        )
        for batch in batches:
            batch_adata = adata[adata.obs[subset.slice_key].astype(str) == batch].copy()
            batch_adata.obs_names = [f"{obs}_{batch}" for obs in batch_adata.obs_names]
            batch_adata.obsm["spatial"] = np.asarray(batch_adata.obsm[subset.coord_key], dtype=np.float32)
            max_neigh = max(1, min(50, batch_adata.n_obs - 1))
            STAligner.Cal_Spatial_Net(
                batch_adata,
                rad_cutoff=self.radius_cutoff,
                max_neigh=max_neigh,
                model="Radius",
                verbose=False,
            )
            batch_list.append(batch_adata)
            adj_list.append(batch_adata.uns["adj"])

        adata_concat = ad.concat(batch_list, label="slice_name", keys=batches)
        adata_concat.obs["batch_name"] = adata_concat.obs["slice_name"].astype("category")
        adj_concat = np.asarray(adj_list[0].todense())
        for idx in range(1, len(adj_list)):
            adj_concat = scipy.linalg.block_diag(adj_concat, np.asarray(adj_list[idx].todense()))
        adata_concat.uns["edgeList"] = np.nonzero(adj_concat)

        device = torch.device(self.device if self.device != "auto" else ("cuda" if torch.cuda.is_available() else "cpu"))
        trained = STAligner.train_STAligner(
            adata_concat,
            hidden_dims=self.hidden_dims,
            n_epochs=self.n_epochs,
            lr=self.lr,
            key_added="STAligner",
            random_seed=random_seed,
            knn_neigh=self.knn_neigh,
            device=device,
        )
        trained_index = {obs_name: idx for idx, obs_name in enumerate(trained.obs_names.astype(str))}
        order = np.array([trained_index[name] for name in target_obs_names], dtype=int)
        embedding = np.asarray(trained.obsm["STAligner"][order], dtype=np.float32)
        return IntegrationResult(
            method_name=self.name,
            embedding=embedding,
            slice_embedding=self._build_slice_embedding(dataset, embedding),
            metadata={
                "n_input_features": int(adata.n_vars),
                "input_feature_hash": md5("\n".join(selected_features.feature_names).encode("utf-8")).hexdigest(),
                "hidden_dims": list(self.hidden_dims),
                "n_epochs": self.n_epochs,
                "lr": self.lr,
                "knn_neigh": self.knn_neigh,
                "radius_cutoff": self.radius_cutoff,
                "device": str(device),
                "representation_type": "latent",
                "implementation_version": self.implementation_version,
            },
        )
