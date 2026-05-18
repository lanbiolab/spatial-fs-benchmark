"""Integration method registry."""

from __future__ import annotations

from importlib import import_module

from spatial_fs_benchmark.integration.base import IntegrationResult, SpatialIntegrator

INTEGRATOR_REGISTRY = {
    "cellcharter": "spatial_fs_benchmark.integration.cellcharter:CellCharterIntegrator",
    "combat": "spatial_fs_benchmark.integration.combat:CombatIntegrator",
    "gpsa": "spatial_fs_benchmark.integration.gpsa:GPSAIntegrator",
    "graphst": "spatial_fs_benchmark.integration.graphst:GraphSTIntegrator",
    "paste": "spatial_fs_benchmark.integration.paste:PASTEIntegrator",
    "pca": "spatial_fs_benchmark.integration.pca:PCAIntegrator",
    "scanorama": "spatial_fs_benchmark.integration.scanorama:ScanoramaIntegrator",
    "scvi": "spatial_fs_benchmark.integration.scvi:SCVIIntegrator",
    "staligner": "spatial_fs_benchmark.integration.staligner:STAlignerIntegrator",
    "symphony": "spatial_fs_benchmark.integration.symphony:SymphonyIntegrator",
}


def build_integrator(name: str, **params: object) -> SpatialIntegrator:
    module_path, class_name = INTEGRATOR_REGISTRY[name].split(":")
    module = import_module(module_path)
    integrator_cls = getattr(module, class_name)
    return integrator_cls(**params)


__all__ = [
    "IntegrationResult",
    "SpatialIntegrator",
    "build_integrator",
]
