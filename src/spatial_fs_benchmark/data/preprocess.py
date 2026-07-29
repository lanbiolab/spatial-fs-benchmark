from __future__ import annotations

import numpy as np
import scanpy as sc
from scipy import sparse

from spatial_fs_benchmark.data.spatial_object import SpatialDataset


PREPROCESS_IMPLEMENTATION_VERSION = "v2_count_layer_contract"


def is_nonnegative_integer_matrix(matrix, tolerance: float = 1e-6) -> bool:
    values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix).ravel()
    if values.size == 0:
        return True
    return bool(
        np.all(np.isfinite(values))
        and np.min(values) >= -tolerance
        and np.all(np.abs(values - np.rint(values)) <= tolerance)
    )


def _resolve_count_source(adata, requested: object | None) -> tuple[str, object]:
    if requested is not None:
        source = str(requested)
        if source == "X":
            return source, adata.X
        if source not in adata.layers:
            raise ValueError(f"Configured count layer '{source}' is not present in the dataset.")
        return source, adata.layers[source]
    for source in ("counts", "raw_count"):
        if source in adata.layers:
            return source, adata.layers[source]
    return "X", adata.X


def preprocess_dataset(dataset: SpatialDataset, options: dict[str, object]) -> SpatialDataset:
    processed = dataset.copy()
    adata = processed.adata
    count_source, count_matrix = _resolve_count_source(adata, options.get("count_layer"))
    adata.layers["counts"] = count_matrix.copy()
    count_like = is_nonnegative_integer_matrix(adata.layers["counts"])
    adata.uns.setdefault("spatial_fs_benchmark", {}).update(
        {
            "preprocess_implementation_version": PREPROCESS_IMPLEMENTATION_VERSION,
            "counts_source": count_source,
            "counts_are_nonnegative_integers": count_like,
        }
    )
    if bool(options.get("require_integer_counts", False)) and not count_like:
        raise ValueError(
            f"Dataset '{dataset.name}' does not provide non-negative integer counts "
            f"through source '{count_source}'."
        )
    min_cells_per_gene = int(options.get("min_cells_per_gene", 0) or 0)
    if min_cells_per_gene > 0:
        sc.pp.filter_genes(adata, min_cells=min_cells_per_gene)
    if bool(options.get("normalize_total", False)):
        sc.pp.normalize_total(adata, target_sum=float(options.get("target_sum", 1e4)))
    if bool(options.get("log1p", False)):
        sc.pp.log1p(adata)
    if bool(options.get("scale", False)):
        sc.pp.scale(adata, zero_center=True, max_value=10.0)
    return processed
