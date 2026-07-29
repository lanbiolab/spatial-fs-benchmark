import anndata as ad
import numpy as np
from scipy import sparse

from spatial_fs_benchmark.data.preprocess import preprocess_dataset
from spatial_fs_benchmark.data.spatial_object import SpatialDataset


def test_preprocess_prefers_raw_count_layer() -> None:
    counts = sparse.csr_matrix([[1, 0], [2, 3]], dtype=np.float32)
    adata = ad.AnnData(X=counts.copy())
    adata.layers["raw_count"] = counts.copy()
    adata.X = counts.copy()
    adata.X.data = np.log1p(adata.X.data)
    adata.obs["slice"] = ["a", "b"]
    adata.obsm["spatial"] = np.array([[0, 0], [1, 0]], dtype=float)
    dataset = SpatialDataset("test", adata, slice_key="slice", coord_key="spatial")

    processed = preprocess_dataset(dataset, {})

    np.testing.assert_array_equal(processed.adata.layers["counts"].toarray(), counts.toarray())
    assert processed.adata.uns["spatial_fs_benchmark"]["counts_source"] == "raw_count"
    assert processed.adata.uns["spatial_fs_benchmark"]["counts_are_nonnegative_integers"] is True


def test_preprocess_rejects_transformed_counts_when_required() -> None:
    adata = ad.AnnData(X=np.array([[0.1, 0.2], [0.3, 0.4]], dtype=np.float32))
    adata.obs["slice"] = ["a", "b"]
    adata.obsm["spatial"] = np.array([[0, 0], [1, 0]], dtype=float)
    dataset = SpatialDataset("test", adata, slice_key="slice", coord_key="spatial")

    try:
        preprocess_dataset(dataset, {"require_integer_counts": True})
    except ValueError as exc:
        assert "non-negative integer counts" in str(exc)
    else:
        raise AssertionError("Transformed expression was accepted as integer counts")
