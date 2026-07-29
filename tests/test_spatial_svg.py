import anndata as ad
import numpy as np

from spatial_fs_benchmark.data.spatial_object import SpatialDataset
from spatial_fs_benchmark.feature_selection.spatial_python_methods import _sanitize_pvalues
from spatial_fs_benchmark.feature_selection.svg import SVGSelector


def test_svg_pvalue_sanitization_is_conservative() -> None:
    observed = _sanitize_pvalues(np.asarray([-1e-12, 0.2, 1 + 1e-12, np.nan, np.inf]))
    assert np.array_equal(observed, np.asarray([0.0, 0.2, 1.0, 1.0, 1.0]))


def test_morans_i_ranks_reproducible_spatial_gradient_first() -> None:
    rng = np.random.default_rng(7)
    coords_one = np.column_stack([np.arange(30), np.zeros(30)])
    coords = np.vstack([coords_one, coords_one])
    gradient = np.tile(np.linspace(0, 10, 30), 2)
    noise = rng.normal(size=60)
    constant = np.ones(60)
    matrix = np.column_stack([gradient, noise, constant]).astype(np.float32)
    adata = ad.AnnData(X=matrix)
    adata.var_names = ["gradient", "noise", "constant"]
    adata.obs["slice"] = ["s1"] * 30 + ["s2"] * 30
    adata.obsm["spatial"] = coords
    dataset = SpatialDataset("synthetic", adata, slice_key="slice", coord_key="spatial")

    result = SVGSelector(n_neighbors=4, chunk_size=2).select(dataset, n_features=2)

    assert result.feature_names[0] == "gradient"
    assert result.metadata["n_slices"] == 2
