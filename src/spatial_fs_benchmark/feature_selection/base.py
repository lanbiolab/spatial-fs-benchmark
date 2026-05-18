from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import anndata as ad
import numpy as np
from hashlib import md5

from spatial_fs_benchmark.data.spatial_object import SpatialDataset


@dataclass(slots=True)
class FeatureSelectionResult:
    method_name: str
    feature_names: list[str]
    feature_indices: list[int]
    scores: list[float]
    metadata: dict[str, Any] = field(default_factory=dict)


class FeatureSelector:
    name = "base"

    def select(self, dataset: SpatialDataset, n_features: int, random_seed: int = 0) -> FeatureSelectionResult:
        raise NotImplementedError

    @staticmethod
    def _adata_from_counts(dataset: SpatialDataset) -> ad.AnnData:
        adata = dataset.adata.copy()
        if "counts" in adata.layers:
            adata.X = adata.layers["counts"].copy()
            adata.uns.pop("log1p", None)
        return adata

    @staticmethod
    def _maybe_subsample_obs(adata: ad.AnnData, max_cells: int | None, random_seed: int) -> ad.AnnData:
        if max_cells is None or adata.n_obs <= max_cells:
            return adata
        rng = np.random.default_rng(random_seed)
        keep = np.sort(rng.choice(adata.n_obs, size=max_cells, replace=False))
        return adata[keep].copy()

    @staticmethod
    def _build_result(
        method_name: str,
        gene_names: np.ndarray,
        scores: np.ndarray,
        n_features: int,
        metadata: dict[str, Any] | None = None,
    ) -> FeatureSelectionResult:
        order = np.argsort(scores)[::-1][:n_features]
        selected_names = gene_names[order].tolist()
        payload = metadata.copy() if metadata is not None else {}
        payload.update(
            {
                "requested_n_features": int(n_features),
                "effective_n_features": int(len(order)),
                "feature_name_hash": md5("\n".join(selected_names).encode("utf-8")).hexdigest(),
            }
        )
        return FeatureSelectionResult(
            method_name=method_name,
            feature_names=selected_names,
            feature_indices=order.tolist(),
            scores=scores[order].tolist(),
            metadata=payload,
        )

    @staticmethod
    def _build_ranked_result(
        method_name: str,
        gene_names: np.ndarray,
        ranked_features: list[str],
        ranked_scores: list[float] | np.ndarray,
        n_features: int,
        metadata: dict[str, Any] | None = None,
    ) -> FeatureSelectionResult:
        index_by_gene = {str(gene): idx for idx, gene in enumerate(gene_names.tolist())}
        seen: set[str] = set()
        selected_names: list[str] = []
        selected_indices: list[int] = []
        selected_scores: list[float] = []
        n_take = min(int(n_features), len(ranked_features))
        for feature_name, score in zip(ranked_features, ranked_scores, strict=False):
            feature_name = str(feature_name)
            if feature_name in seen or feature_name not in index_by_gene:
                continue
            seen.add(feature_name)
            selected_names.append(feature_name)
            selected_indices.append(index_by_gene[feature_name])
            selected_scores.append(float(score))
            if len(selected_names) >= n_take:
                break
        payload = metadata.copy() if metadata is not None else {}
        payload.update(
            {
                "requested_n_features": int(n_features),
                "effective_n_features": int(len(selected_names)),
                "feature_name_hash": md5("\n".join(selected_names).encode("utf-8")).hexdigest(),
            }
        )
        return FeatureSelectionResult(
            method_name=method_name,
            feature_names=selected_names,
            feature_indices=selected_indices,
            scores=selected_scores,
            metadata=payload,
        )
