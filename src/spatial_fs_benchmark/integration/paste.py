from __future__ import annotations

from hashlib import md5

import numpy as np
import ot
import scanpy as sc
from ot.backend import NumpyBackend
from paste.helper import extract_data_matrix, to_dense_array

from spatial_fs_benchmark.data.spatial_object import SpatialDataset
from spatial_fs_benchmark.feature_selection.base import FeatureSelectionResult
from spatial_fs_benchmark.integration.base import IntegrationResult, SpatialIntegrator


class PASTEIntegrator(SpatialIntegrator):
    name = "paste"

    def __init__(
        self,
        alpha: float = 0.1,
        dissimilarity: str = "kl",
        use_rep: str | None = None,
    ) -> None:
        self.alpha = alpha
        self.dissimilarity = dissimilarity
        self.use_rep = use_rep

    def fit_transform(
        self,
        dataset: SpatialDataset,
        selected_features: FeatureSelectionResult,
        random_seed: int = 0,
    ) -> IntegrationResult:
        import paste as pst

        subset = dataset.subset_features(selected_features.feature_names)
        adata = subset.adata.copy()
        if "counts" in adata.layers:
            adata.X = adata.layers["counts"].copy()
            adata.uns.pop("log1p", None)
        elif hasattr(adata.X, "min") and float(adata.X.min()) < 0:
            raise ValueError(
                f"PASTE requires non-negative expression-like input, but dataset '{dataset.name}' contains negative values."
            )

        slice_names = np.unique(subset.slice_ids).tolist()
        slices = [adata[adata.obs[subset.slice_key].astype(str) == slice_name].copy() for slice_name in slice_names]
        for slice_adata in slices:
            slice_adata.var_names_make_unique()
            if self.use_rep == "pca":
                sc.pp.pca(slice_adata, n_comps=min(30, max(2, slice_adata.n_vars - 1)))

        pis = []
        for idx in range(len(slices) - 1):
            slice_a = slices[idx]
            slice_b = slices[idx + 1]
            coords_a = np.asarray(slice_a.obsm[subset.coord_key], dtype=np.float64)
            coords_b = np.asarray(slice_b.obsm[subset.coord_key], dtype=np.float64)
            d_a = ot.dist(coords_a, coords_a, metric="euclidean")
            d_b = ot.dist(coords_b, coords_b, metric="euclidean")
            if np.any(d_a > 0):
                d_a = d_a / np.min(d_a[d_a > 0])
            if np.any(d_b > 0):
                d_b = d_b / np.min(d_b[d_b > 0])

            x_a = np.asarray(to_dense_array(extract_data_matrix(slice_a, self.use_rep)), dtype=np.float64)
            x_b = np.asarray(to_dense_array(extract_data_matrix(slice_b, self.use_rep)), dtype=np.float64)
            if self.dissimilarity.lower() in {"euclidean", "euc"}:
                m = ot.dist(x_a, x_b, metric="euclidean")
            else:
                x_a = x_a + 1e-2
                x_b = x_b + 1e-2
                x_a = x_a / x_a.sum(axis=1, keepdims=True)
                x_b = x_b / x_b.sum(axis=1, keepdims=True)
                m = np.sum(x_a[:, None, :] * np.log(x_a[:, None, :] / x_b[None, :, :]), axis=2)
            p = np.full(slice_a.n_obs, 1.0 / slice_a.n_obs, dtype=np.float64)
            q = np.full(slice_b.n_obs, 1.0 / slice_b.n_obs, dtype=np.float64)
            pi0 = pst.match_spots_using_spatial_heuristic(coords_a, coords_b, use_ot=True)
            pi = ot.gromov.fused_gromov_wasserstein(
                m,
                d_a,
                d_b,
                p=p,
                q=q,
                loss_fun="square_loss",
                alpha=self.alpha,
                G0=pi0 / pi0.sum(),
                max_iter=200,
                tol_rel=1e-7,
                tol_abs=1e-7,
            )
            pis.append(pi)

        aligned_slices = pst.stack_slices_pairwise(slices, pis)
        aligned = sc.concat(aligned_slices, label=subset.slice_key, keys=slice_names, join="outer")
        aligned = aligned[adata.obs_names].copy()
        embedding = np.asarray(aligned.obsm[subset.coord_key], dtype=np.float32)
        return IntegrationResult(
            method_name=self.name,
            embedding=embedding,
            slice_embedding=self._build_slice_embedding(dataset, embedding),
            metadata={
                "aligned_representation": embedding,
                "n_input_features": int(adata.n_vars),
                "input_feature_hash": md5("\n".join(selected_features.feature_names).encode("utf-8")).hexdigest(),
                "alpha": self.alpha,
                "dissimilarity": self.dissimilarity,
                "use_rep": self.use_rep or "",
            },
        )
