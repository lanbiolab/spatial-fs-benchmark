from __future__ import annotations

import scanpy as sc

from spatial_fs_benchmark.data.spatial_object import SpatialDataset


def preprocess_dataset(dataset: SpatialDataset, options: dict[str, object]) -> SpatialDataset:
    processed = dataset.copy()
    adata = processed.adata
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()
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
