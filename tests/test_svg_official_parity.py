from __future__ import annotations

import importlib.metadata
import subprocess
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from scipy.stats import rankdata

from spatial_fs_benchmark.data.preprocess import preprocess_dataset
from spatial_fs_benchmark.data.spatial_object import SpatialDataset
from spatial_fs_benchmark.feature_selection import build_feature_selector
from spatial_fs_benchmark.feature_selection.spatial_python_methods import (
    SpatialDESelector,
    _legacy_scipy_compatibility,
)


def _synthetic_dataset(n_spots: int = 100, n_genes: int = 24) -> SpatialDataset:
    rng = np.random.default_rng(12)
    side = int(np.sqrt(n_spots))
    coords = np.column_stack([np.arange(n_spots) % side, np.arange(n_spots) // side]).astype(float)
    means = np.full((n_spots, n_genes), 2.0)
    means[:, 0] = 0.3 + 8 * coords[:, 0] / max(coords[:, 0].max(), 1)
    means[:, 1] = 0.3 + 8 * coords[:, 1] / max(coords[:, 1].max(), 1)
    counts = rng.poisson(means).astype(np.int32)
    adata = ad.AnnData(X=sparse.csr_matrix(counts))
    adata.var_names = [f"g{i}" for i in range(n_genes)]
    adata.obs["slice"] = "s1"
    adata.obsm["spatial"] = coords
    return preprocess_dataset(SpatialDataset("synthetic", adata, "slice", "spatial"), {})


def test_spatialde_wrapper_matches_official_tutorial_calls() -> None:
    for package in ("SpatialDE", "NaiveDE"):
        try:
            importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            pytest.skip(f"{package} is not installed")

    dataset = _synthetic_dataset()
    genes = dataset.adata.var_names.to_numpy().astype(str)
    counts = dataset.adata.layers["counts"].T
    wrapper_scores = SpatialDESelector()._score_slice(counts, dataset.coords, genes)

    _legacy_scipy_compatibility()
    import NaiveDE
    import SpatialDE

    expression = pd.DataFrame(counts.T.toarray(), columns=genes)
    sample_info = pd.DataFrame(
        {
            "x": dataset.coords[:, 0],
            "y": dataset.coords[:, 1],
            "total_counts": expression.sum(axis=1).to_numpy(),
        },
        index=expression.index,
    )
    normalized = NaiveDE.stabilize(expression.T).T
    residual = NaiveDE.regress_out(sample_info, normalized.T, "np.log(total_counts)").T
    direct = SpatialDE.run(sample_info[["x", "y"]].to_numpy(), residual)
    statistic = -np.log10(np.maximum(direct["pval"].to_numpy(), np.finfo(float).tiny))
    expected = dict(zip(direct["g"].astype(str), rankdata(statistic) / len(statistic), strict=True))

    assert wrapper_scores.keys() == expected.keys()
    np.testing.assert_allclose(list(wrapper_scores.values()), list(expected.values()))
    assert list(sorted(wrapper_scores, key=wrapper_scores.get, reverse=True))[:2] == ["g0", "g1"]


@pytest.mark.parametrize("method", ["sparkx", "nnsvg"])
def test_r_svg_wrappers_recover_spatial_gradients_from_counts(method: str) -> None:
    rscript = Path(sys.executable).with_name("Rscript")
    if not rscript.exists():
        pytest.skip("Rscript is not installed")
    package = "SPARK" if method == "sparkx" else "nnSVG"
    package_check = subprocess.run(
        [str(rscript), "-e", f"quit(status=ifelse(requireNamespace('{package}', quietly=TRUE), 0, 1))"],
        check=False,
    )
    if package_check.returncode != 0:
        pytest.skip(f"R package {package} is not installed")

    result = build_feature_selector(method, n_threads=2).select(_synthetic_dataset(), n_features=5)

    assert {"g0", "g1"}.issubset(result.feature_names)
    assert result.metadata["input_assay"] == "counts"


def test_somde_official_workflow_recovers_spatial_gradients() -> None:
    try:
        importlib.metadata.version("SOMDE")
    except importlib.metadata.PackageNotFoundError:
        pytest.skip("SOMDE is not installed")

    rng = np.random.default_rng(22)
    n_spots, n_genes, side = 400, 80, 20
    coords = np.column_stack([np.arange(n_spots) % side, np.arange(n_spots) // side]).astype(float)
    means = np.full((n_spots, n_genes), 2.0)
    means[:, 0] = 0.2 + 12 * coords[:, 0] / 19
    means[:, 1] = 0.2 + 12 * coords[:, 1] / 19
    counts = rng.negative_binomial(1, 1 / (1 + means)).astype(np.int32)
    adata = ad.AnnData(X=sparse.csr_matrix(counts))
    adata.var_names = [f"g{i}" for i in range(n_genes)]
    adata.obs["slice"] = "s1"
    adata.obsm["spatial"] = coords
    dataset = preprocess_dataset(SpatialDataset("synthetic_nb", adata, "slice", "spatial"), {})

    result = build_feature_selector("somde", spots_per_node=5).select(dataset, n_features=5)

    assert {"g0", "g1"}.issubset(result.feature_names)
    assert result.metadata["official_preprocessing"] == "SomNode spatial aggregation plus SOMDE norm"


def test_official_svg_wrappers_reject_transformed_pseudo_counts() -> None:
    dataset = _synthetic_dataset()
    dataset.adata.layers["counts"] = dataset.adata.layers["counts"].astype(float)
    dataset.adata.layers["counts"].data = np.log1p(dataset.adata.layers["counts"].data)
    dataset.adata.uns["spatial_fs_benchmark"]["counts_source"] = "transformed"

    for method in ("sparkx", "nnsvg", "spatialde", "somde"):
        with pytest.raises(ValueError, match="requires integer counts"):
            build_feature_selector(method).select(dataset, n_features=5)
