"""Feature selection registry."""

from spatial_fs_benchmark.feature_selection.anticor import AnticorSelector
from spatial_fs_benchmark.feature_selection.all_features import AllFeaturesSelector
from spatial_fs_benchmark.feature_selection.base import FeatureSelectionResult, FeatureSelector
from spatial_fs_benchmark.feature_selection.conserved import CrossSliceConservedSelector
from spatial_fs_benchmark.feature_selection.highly_expressed import HighlyExpressedSelector
from spatial_fs_benchmark.feature_selection.hotspot_selector import HotspotSelector
from spatial_fs_benchmark.feature_selection.hybrid import HybridHVGSVGSelector
from spatial_fs_benchmark.feature_selection.hvg import HVGSelector
from spatial_fs_benchmark.feature_selection.random_baseline import RandomSelector
from spatial_fs_benchmark.feature_selection.r_method import RMethodSelector
from spatial_fs_benchmark.feature_selection.scanpy_hvg import ScanpyHVGSelector
from spatial_fs_benchmark.feature_selection.statistic import StatisticSelector
from spatial_fs_benchmark.feature_selection.tfs import TFSelector
from spatial_fs_benchmark.feature_selection.triku_selector import TrikuSelector
from spatial_fs_benchmark.feature_selection.svg import SVGSelector
from spatial_fs_benchmark.feature_selection.spatial_r_method import SpatialRMethodSelector
from spatial_fs_benchmark.feature_selection.spatial_python_methods import SOMDESelector, SpatialDESelector
from spatial_fs_benchmark.feature_selection.scgco_selector import ScGCOSelector
from spatial_fs_benchmark.feature_selection.wilcoxon import WilcoxonSelector

SELECTOR_REGISTRY = {
    "anticor": AnticorSelector,
    "all": AllFeaturesSelector,
    "all_features": AllFeaturesSelector,
    "conserved": CrossSliceConservedSelector,
    "cross_slice_conserved": CrossSliceConservedSelector,
    "highly_expressed": HighlyExpressedSelector,
    "hotspot": HotspotSelector,
    "hybrid_hvg_svg": HybridHVGSVGSelector,
    "hvg": HVGSelector,
    "random": RandomSelector,
    "r_method": RMethodSelector,
    "scanpy_hvg": ScanpyHVGSelector,
    "statistic": StatisticSelector,
    "TFs": TFSelector,
    "tfs": TFSelector,
    "triku": TrikuSelector,
    "svg": SVGSelector,
    "morans_i": SVGSelector,
    "spatial_r_method": SpatialRMethodSelector,
    "spatialde": SpatialDESelector,
    "somde": SOMDESelector,
    "scgco": ScGCOSelector,
    "wilcoxon": WilcoxonSelector,
}


def _selector_alias(name: str, params: dict[str, object]) -> tuple[str, dict[str, object]]:
    resolved = dict(params)
    scanpy_aliases = {
        "scanpy_seurat": {"flavor": "seurat", "batch": False},
        "scanpy_seurat_batch": {"flavor": "seurat", "batch": True},
        "scanpy_cell_ranger": {"flavor": "cell_ranger", "batch": False},
        "scanpy_cell_ranger_batch": {"flavor": "cell_ranger", "batch": True},
        "scanpy_pearson": {"flavor": "pearson", "batch": False},
        "scanpy_pearson_batch": {"flavor": "pearson", "batch": True},
        "scanpy_seurat_v3": {"flavor": "seurat_v3", "batch": False},
        "scanpy_seurat_v3_batch": {"flavor": "seurat_v3", "batch": True},
    }
    statistic_aliases = {
        "statistic_mean": {"statistic": "mean"},
        "statistic_variance": {"statistic": "variance"},
    }
    r_method_aliases = {
        "seurat_vst": {"method": "seurat_vst"},
        "seurat_mvp": {"method": "seurat_mvp"},
        "seurat_disp": {"method": "seurat_disp"},
        "seurat_sct": {"method": "seurat_sct"},
        "scsegindex": {"method": "scsegindex"},
        "dubstepr": {"method": "dubstepr"},
        "nbumi": {"method": "nbumi"},
        "osca": {"method": "osca"},
        "scry": {"method": "scry"},
        "singleCellHaystack": {"method": "singleCellHaystack"},
        "Brennecke": {"method": "Brennecke"},
        "brennecke": {"method": "Brennecke"},
        "scPNMF": {"method": "scPNMF"},
        "scpnmf": {"method": "scPNMF"},
    }
    spatial_r_method_aliases = {
        "sparkx": {"method": "sparkx"},
        "nnsvg": {"method": "nnsvg"},
    }
    hybrid_aliases = {
        "hybrid_hvg_svg_union": {"mode": "balanced_union"},
        "hybrid_hvg_svg_intersection": {"mode": "intersection"},
    }
    if name in scanpy_aliases:
        resolved.update(scanpy_aliases[name])
        return "scanpy_hvg", resolved
    if name in statistic_aliases:
        resolved.update(statistic_aliases[name])
        return "statistic", resolved
    if name in r_method_aliases:
        resolved.update(r_method_aliases[name])
        return "r_method", resolved
    if name in spatial_r_method_aliases:
        resolved.update(spatial_r_method_aliases[name])
        return "spatial_r_method", resolved
    if name in hybrid_aliases:
        resolved.update(hybrid_aliases[name])
        return "hybrid_hvg_svg", resolved
    return name, resolved


def build_feature_selector(name: str, **params: object) -> FeatureSelector:
    name, params = _selector_alias(name, params)
    selector_cls = SELECTOR_REGISTRY[name]
    return selector_cls(**params)


__all__ = [
    "FeatureSelectionResult",
    "FeatureSelector",
    "build_feature_selector",
]
