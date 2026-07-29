from __future__ import annotations

import numpy as np

from spatial_fs_benchmark.feature_selection import build_feature_selector
from spatial_fs_benchmark.validation import (
    SemiSyntheticSpec,
    build_semi_synthetic_dataset,
    heldout_moran_reference,
    score_heldout_clustering,
    score_spatial_reproducibility,
    split_heldout_slice,
)


def _dataset():
    return build_semi_synthetic_dataset(
        SemiSyntheticSpec(prevalence=0.12, seed=19, n_genes=300, grid_size=9, n_slices=3)
    )


def test_heldout_split_is_disjoint_and_preserves_count_layers() -> None:
    dataset = _dataset()
    training, heldout = split_heldout_slice(dataset, "slice3")
    assert set(training.slice_ids) == {"slice1", "slice2"}
    assert set(heldout.slice_ids) == {"slice3"}
    assert set(training.adata.obs_names).isdisjoint(set(heldout.adata.obs_names))
    assert "counts" in training.adata.layers
    assert "counts" in heldout.adata.layers


def test_training_moran_ranking_reproduces_better_than_random_on_heldout_slice() -> None:
    dataset = _dataset()
    training, heldout = split_heldout_slice(dataset, "slice3")
    moran = build_feature_selector("morans_i", n_neighbors=6).select(
        training, n_features=80, random_seed=19
    )
    random = build_feature_selector("random").select(training, n_features=80, random_seed=19)
    reference = heldout_moran_reference(heldout, random_seed=19, n_neighbors=6)
    moran_scores = score_spatial_reproducibility(heldout, moran, reference, n_features=80)
    random_scores = score_spatial_reproducibility(heldout, random, reference, n_features=80)
    assert moran_scores["heldout_moran_percentile"] > random_scores["heldout_moran_percentile"]


def test_heldout_clustering_returns_finite_metrics() -> None:
    dataset = _dataset()
    training, heldout = split_heldout_slice(dataset, "slice3")
    selection = build_feature_selector("morans_i", n_neighbors=6).select(
        training, n_features=80, random_seed=19
    )
    scores = score_heldout_clustering(heldout, selection, n_features=80, random_seed=19)
    assert {"heldout_ari", "heldout_nmi", "heldout_silhouette", "heldout_CHAOS", "heldout_PAS"} <= set(scores)
    assert all(np.isfinite(list(scores.values())))
