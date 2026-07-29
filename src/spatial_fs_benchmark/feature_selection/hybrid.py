from __future__ import annotations

import numpy as np

from spatial_fs_benchmark.data.spatial_object import SpatialDataset
from spatial_fs_benchmark.feature_selection.base import FeatureSelectionResult, FeatureSelector
from spatial_fs_benchmark.feature_selection.scanpy_hvg import ScanpyHVGSelector
from spatial_fs_benchmark.feature_selection.svg import SVGSelector


class HybridHVGSVGSelector(FeatureSelector):
    """Combine batch-aware Seurat-v3 HVGs with slice-wise Moran-ranked genes."""

    implementation_version = "v1_scanpy_seurat_v3_batch_plus_moran8"

    def __init__(self, mode: str = "balanced_union", n_neighbors: int = 8) -> None:
        if mode not in {"balanced_union", "intersection"}:
            raise ValueError(f"Unsupported hybrid mode: {mode}")
        self.mode = mode
        self.n_neighbors = n_neighbors
        self.name = f"hybrid_hvg_svg_{mode}"
        self.stochastic_selection = False
        self.hvg_selector = ScanpyHVGSelector(flavor="seurat_v3", batch=True)
        self.svg_selector = SVGSelector(n_neighbors=n_neighbors)

    @staticmethod
    def _balanced_union(hvg: list[str], svg: list[str], target: int) -> list[str]:
        half = target // 2
        selected: list[str] = []
        seen: set[str] = set()

        def append_unique(values: list[str]) -> None:
            for value in values:
                if value not in seen:
                    seen.add(value)
                    selected.append(value)

        append_unique(hvg[:half])
        append_unique(svg[: target - half])
        depth = half
        while len(selected) < target and (depth < len(hvg) or depth < len(svg)):
            if depth < len(hvg):
                append_unique([hvg[depth]])
            if len(selected) >= target:
                break
            if depth < len(svg):
                append_unique([svg[depth]])
            depth += 1
        return selected[:target]

    @staticmethod
    def _intersection(hvg: list[str], svg: list[str], target: int) -> list[str]:
        hvg_top = hvg[:target]
        svg_top = svg[:target]
        hvg_rank = {gene: index for index, gene in enumerate(hvg_top)}
        svg_rank = {gene: index for index, gene in enumerate(svg_top)}
        shared = set(hvg_rank) & set(svg_rank)
        return sorted(shared, key=lambda gene: (hvg_rank[gene] + svg_rank[gene], gene))

    def select(
        self,
        dataset: SpatialDataset,
        n_features: int,
        random_seed: int = 0,
    ) -> FeatureSelectionResult:
        ranking_depth = min(dataset.n_vars, max(int(n_features) * 2, int(n_features)))
        hvg = self.hvg_selector.select(dataset, ranking_depth, random_seed=random_seed)
        svg = self.svg_selector.select(dataset, ranking_depth, random_seed=random_seed)
        if self.mode == "balanced_union":
            ranked = self._balanced_union(hvg.feature_names, svg.feature_names, int(n_features))
        else:
            ranked = self._intersection(hvg.feature_names, svg.feature_names, int(n_features))
        scores = np.arange(len(ranked), 0, -1, dtype=float)
        return self._build_ranked_result(
            method_name=self.name,
            gene_names=np.asarray(dataset.gene_names),
            ranked_features=ranked,
            ranked_scores=scores,
            n_features=n_features,
            metadata={
                "mode": self.mode,
                "hvg_component": "scanpy_seurat_v3_batch",
                "svg_component": "morans_i",
                "moran_n_neighbors": self.n_neighbors,
                "hvg_ranking_depth": len(hvg.feature_names),
                "svg_ranking_depth": len(svg.feature_names),
                "shared_top_n": len(set(hvg.feature_names[:n_features]) & set(svg.feature_names[:n_features])),
            },
        )
