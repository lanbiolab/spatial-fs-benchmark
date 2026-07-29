"""Independent validation utilities for the spatial feature benchmark."""

from spatial_fs_benchmark.validation.semi_synthetic import (
    SemiSyntheticSpec,
    build_semi_synthetic_dataset,
    score_feature_ranking,
)
from spatial_fs_benchmark.validation.heldout_slice import (
    HeldoutFold,
    heldout_metrics_frame,
    heldout_moran_reference,
    score_heldout_clustering,
    score_spatial_reproducibility,
    split_heldout_slice,
)

__all__ = [
    "SemiSyntheticSpec",
    "build_semi_synthetic_dataset",
    "score_feature_ranking",
    "HeldoutFold",
    "heldout_metrics_frame",
    "heldout_moran_reference",
    "score_heldout_clustering",
    "score_spatial_reproducibility",
    "split_heldout_slice",
]
