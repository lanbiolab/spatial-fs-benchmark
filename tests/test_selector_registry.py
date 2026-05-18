from __future__ import annotations

from spatial_fs_benchmark.feature_selection import build_feature_selector


def test_atlas_reuse_selectors_build() -> None:
    selector_names = [
        "TFs",
        "scanpy_seurat",
        "scanpy_seurat_batch",
        "scanpy_seurat_v3",
        "scanpy_seurat_v3_batch",
        "scanpy_cell_ranger",
        "scanpy_cell_ranger_batch",
        "scanpy_pearson",
        "scanpy_pearson_batch",
        "seurat_vst",
        "seurat_mvp",
        "seurat_disp",
        "seurat_sct",
        "scsegindex",
        "dubstepr",
        "nbumi",
        "osca",
        "scry",
        "singleCellHaystack",
        "Brennecke",
        "scPNMF",
        "triku",
        "hotspot",
        "anticor",
        "statistic_mean",
        "statistic_variance",
        "wilcoxon",
        "all",
    ]
    for name in selector_names:
        selector = build_feature_selector(name)
        assert selector is not None
