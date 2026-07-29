from __future__ import annotations

from pathlib import Path

import anndata as ad

from spatial_fs_benchmark.config import DatasetConfig
from spatial_fs_benchmark.data.spatial_object import SpatialDataset


def load_dataset(config: DatasetConfig) -> SpatialDataset:
    path = Path(config.path)
    adata = ad.read_h5ad(path)
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    if config.slice_key in adata.obs:
        adata.obs[config.slice_key] = adata.obs[config.slice_key].astype("category")
    if config.label_key is not None and config.label_key in adata.obs:
        adata.obs[config.label_key] = adata.obs[config.label_key].astype("category")
    if config.slice_class_key is not None and config.slice_class_key in adata.obs:
        adata.obs[config.slice_class_key] = adata.obs[config.slice_class_key].astype("category")
    return SpatialDataset(
        name=config.name,
        adata=adata,
        slice_key=config.slice_key,
        coord_key=config.coord_key,
        label_key=config.label_key,
        slice_class_key=config.slice_class_key,
        platform=config.platform,
        species=config.species,
        source_path=str(path),
        alignment_pairs=[(str(left), str(right)) for left, right in config.alignment_pairs],
    )
