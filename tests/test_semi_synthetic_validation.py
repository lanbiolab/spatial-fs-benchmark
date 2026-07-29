from __future__ import annotations

import numpy as np

from spatial_fs_benchmark.feature_selection import build_feature_selector
from spatial_fs_benchmark.validation import (
    SemiSyntheticSpec,
    build_semi_synthetic_dataset,
    score_feature_ranking,
)


def test_semi_synthetic_generation_is_deterministic_and_count_valid() -> None:
    spec = SemiSyntheticSpec(prevalence=0.1, seed=7, n_genes=240, grid_size=8, n_slices=2)
    first = build_semi_synthetic_dataset(spec)
    second = build_semi_synthetic_dataset(spec)
    assert first.adata.shape == (128, 240)
    assert np.array_equal(first.adata.layers["counts"].toarray(), second.adata.layers["counts"].toarray())
    assert int(first.adata.var["is_svg"].sum()) == spec.n_truth
    assert set(first.adata.var.loc[first.adata.var["is_svg"], "pattern"]) == {
        "domain",
        "gradient",
        "focal",
        "periodic",
    }


def test_morans_i_recovers_more_truth_than_random_at_fixed_cutoff() -> None:
    dataset = build_semi_synthetic_dataset(
        SemiSyntheticSpec(prevalence=0.1, seed=11, n_genes=360, grid_size=10, n_slices=2)
    )
    moran = build_feature_selector("morans_i", n_neighbors=6).select(
        dataset, n_features=dataset.n_vars, random_seed=11
    )
    random = build_feature_selector("random").select(
        dataset, n_features=dataset.n_vars, random_seed=11
    )
    moran_scores = score_feature_ranking(dataset, moran, cutoffs=(50,))
    random_scores = score_feature_ranking(dataset, random, cutoffs=(50,))
    moran_auprc = moran_scores.loc[moran_scores["metric"].eq("AUPRC"), "value"].iloc[0]
    random_auprc = random_scores.loc[random_scores["metric"].eq("AUPRC"), "value"].iloc[0]
    assert moran_auprc > random_auprc
