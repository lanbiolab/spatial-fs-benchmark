from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import rankdata

from spatial_fs_benchmark.data.preprocess import is_nonnegative_integer_matrix
from spatial_fs_benchmark.data.spatial_object import SpatialDataset
from spatial_fs_benchmark.feature_selection.base import FeatureSelectionResult, FeatureSelector


def _legacy_scipy_compatibility() -> None:
    import scipy
    import scipy.misc

    def derivative(function, x0, dx=1.0, n=1, args=(), order=3):
        if n == 1:
            return (function(x0 + dx, *args) - function(x0 - dx, *args)) / (2 * dx)
        if n == 2:
            return (function(x0 + dx, *args) - 2 * function(x0, *args) + function(x0 - dx, *args)) / (dx**2)
        raise ValueError("Only first and second derivatives are supported")

    scipy.misc.derivative = derivative
    for name in ("arange", "argsort", "array", "zeros_like"):
        if not hasattr(scipy, name):
            setattr(scipy, name, getattr(np, name))


def _percentile_scores_from_pvalues(result: pd.DataFrame) -> dict[str, float]:
    if result.empty:
        raise ValueError("The SVG method returned no genes that could be ranked.")
    pvalues = result["pval"].to_numpy(dtype=float)
    finite = np.isfinite(pvalues)
    scores = np.zeros(len(result), dtype=float)
    if finite.any():
        statistic = -np.log10(np.maximum(pvalues[finite], np.finfo(float).tiny))
        scores[finite] = rankdata(statistic, method="average") / len(statistic)
    return dict(zip(result["g"].astype(str), scores, strict=True))


def _sanitize_pvalues(values) -> np.ndarray:
    pvalues = np.asarray(values, dtype=float)
    return np.clip(np.where(np.isfinite(pvalues), pvalues, 1.0), 0.0, 1.0)


class _SliceWisePythonSVGSelector(FeatureSelector):
    max_cells_per_slice: int | None
    max_genes: int | None

    def __init__(
        self,
        max_cells_per_slice: int | None = None,
        max_genes: int | None = None,
    ) -> None:
        self.max_cells_per_slice = max_cells_per_slice
        self.max_genes = max_genes
        self.stochastic_selection = max_cells_per_slice is not None
        self._score_cache: dict[tuple[str, int], np.ndarray] = {}

    def _matrix(self, dataset: SpatialDataset):
        return dataset.adata.X

    def _score_slice(self, matrix, coords: np.ndarray, genes: np.ndarray) -> dict[str, float]:
        raise NotImplementedError

    def _method_metadata(self) -> dict[str, object]:
        return {}

    def _minimum_cells_per_slice(self) -> int:
        return 10

    @staticmethod
    def _row_variances(matrix) -> np.ndarray:
        if sparse.issparse(matrix):
            means = np.asarray(matrix.mean(axis=1)).ravel()
            means_sq = np.asarray(matrix.multiply(matrix).mean(axis=1)).ravel()
            return np.maximum(means_sq - means**2, 0)
        return np.var(np.asarray(matrix), axis=1)

    def _aggregate_scores(self, dataset: SpatialDataset, random_seed: int) -> np.ndarray:
        matrix = self._matrix(dataset)
        slice_ids = dataset.slice_ids
        genes = dataset.adata.var_names.to_numpy().astype(str)
        score_rows = []
        for slice_number, slice_id in enumerate(np.unique(slice_ids)):
            indices = np.flatnonzero(slice_ids == slice_id)
            if self.max_cells_per_slice is not None and len(indices) > self.max_cells_per_slice:
                rng = np.random.default_rng(random_seed + slice_number * 1009)
                indices = np.sort(rng.choice(indices, self.max_cells_per_slice, replace=False))
            if len(indices) < self._minimum_cells_per_slice():
                continue
            slice_matrix = matrix[indices].T
            variances = self._row_variances(slice_matrix)
            candidates = np.flatnonzero(np.isfinite(variances) & (variances > 0))
            if self.max_genes is not None and len(candidates) > self.max_genes:
                candidates = candidates[np.argsort(variances[candidates])[::-1][: self.max_genes]]
            score_map = self._score_slice(slice_matrix[candidates], dataset.coords[indices], genes[candidates])
            if not score_map:
                continue
            row = np.zeros(len(genes), dtype=float)
            for candidate in candidates:
                row[candidate] = score_map.get(genes[candidate], 0.0)
            score_rows.append(row)
        if not score_rows:
            raise ValueError(f"No slice produced a valid {self.name} ranking for dataset '{dataset.name}'.")
        return np.mean(np.vstack(score_rows), axis=0)

    def select(self, dataset: SpatialDataset, n_features: int, random_seed: int = 0) -> FeatureSelectionResult:
        cache_key = (dataset.source_path or dataset.name, random_seed)
        if cache_key not in self._score_cache:
            self._score_cache[cache_key] = self._aggregate_scores(dataset, random_seed)
        scores = self._score_cache[cache_key]
        ranked_indices = np.flatnonzero(scores > 0)
        ranked_indices = ranked_indices[np.argsort(scores[ranked_indices])[::-1]]
        return self._build_ranked_result(
            method_name=self.name,
            gene_names=dataset.adata.var_names.to_numpy(),
            ranked_features=dataset.adata.var_names.to_numpy()[ranked_indices].tolist(),
            ranked_scores=scores[ranked_indices],
            n_features=n_features,
            metadata={
                "max_cells_per_slice": self.max_cells_per_slice,
                "max_genes": self.max_genes,
                "aggregation": "mean within-slice percentile rank",
                **self._method_metadata(),
            },
        )


class SpatialDESelector(_SliceWisePythonSVGSelector):
    name = "spatialde"
    implementation_version = "v2_spatialde_1.1.3_official_naivede_slice_wise"

    def _matrix(self, dataset: SpatialDataset):
        counts = dataset.adata.layers["counts"]
        if not is_nonnegative_integer_matrix(counts):
            source = dataset.adata.uns.get("spatial_fs_benchmark", {}).get("counts_source", "unknown")
            raise ValueError(f"SpatialDE requires integer counts, but dataset '{dataset.name}' uses '{source}'.")
        return counts

    def _method_metadata(self) -> dict[str, object]:
        return {
            "input_assay": "counts",
            "official_preprocessing": "NaiveDE stabilize plus log-total-count regression",
        }

    def _score_slice(self, matrix, coords: np.ndarray, genes: np.ndarray) -> dict[str, float]:
        _legacy_scipy_compatibility()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            import NaiveDE
            import SpatialDE

        expression = matrix.T.toarray() if sparse.issparse(matrix) else np.asarray(matrix.T)
        nonzero_spots = expression.sum(axis=1) > 0
        expression_frame = pd.DataFrame(expression[nonzero_spots], columns=genes)
        sample_info = pd.DataFrame(
            {
                "x": coords[nonzero_spots, 0],
                "y": coords[nonzero_spots, 1],
                "total_counts": expression_frame.sum(axis=1).to_numpy(),
            },
            index=expression_frame.index,
        )
        normalized = NaiveDE.stabilize(expression_frame.T).T
        residual = NaiveDE.regress_out(
            sample_info,
            normalized.T,
            "np.log(total_counts)",
        ).T
        result = SpatialDE.run(sample_info[["x", "y"]].to_numpy(dtype=float), residual)
        return _percentile_scores_from_pvalues(result)


class SOMDESelector(_SliceWisePythonSVGSelector):
    name = "somde"
    implementation_version = "v4_somde_0.1.8_finite_pvalues_slice_wise"

    def __init__(self, spots_per_node: int = 20, **kwargs) -> None:
        super().__init__(**kwargs)
        self.spots_per_node = spots_per_node

    def _matrix(self, dataset: SpatialDataset):
        counts = dataset.adata.layers["counts"]
        if not is_nonnegative_integer_matrix(counts):
            source = dataset.adata.uns.get("spatial_fs_benchmark", {}).get("counts_source", "unknown")
            raise ValueError(f"SOMDE requires integer counts, but dataset '{dataset.name}' uses '{source}'.")
        return counts

    def _method_metadata(self) -> dict[str, object]:
        return {
            "input_assay": "counts",
            "official_preprocessing": "SomNode spatial aggregation plus SOMDE norm",
            "spots_per_node": self.spots_per_node,
            "minimum_cells_per_slice": self._minimum_cells_per_slice(),
        }

    def _minimum_cells_per_slice(self) -> int:
        # SomNode uses floor(sqrt(n_spots / spots_per_node)) nodes per axis.
        # Require a 2 x 2 grid rather than allowing an invalid 0 x 0 grid.
        return max(10, 4 * self.spots_per_node)

    def _score_slice(self, matrix, coords: np.ndarray, genes: np.ndarray) -> dict[str, float]:
        _legacy_scipy_compatibility()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            import somde

        expression = matrix.toarray() if sparse.issparse(matrix) else np.asarray(matrix)
        frame = pd.DataFrame(expression, index=genes)
        model = somde.SomNode(coords.astype(np.float32), self.spots_per_node)
        model.mtx(frame)
        model.norm()
        # Degenerate fits can produce NaN p-values. SOMDE's qvalue helper
        # asserts on them, so conservatively treat them as non-significant.
        import somde.som as som_module

        original_qvalue = som_module.qvalue

        def finite_qvalue(pvalues, pi0=None):
            sanitized = _sanitize_pvalues(pvalues)
            if sanitized.size == 0:
                return sanitized
            return original_qvalue(sanitized, pi0=pi0)

        som_module.qvalue = finite_qvalue
        try:
            result, _ = model.run()
        finally:
            som_module.qvalue = original_qvalue
        if result.empty:
            return {}
        return _percentile_scores_from_pvalues(result)
