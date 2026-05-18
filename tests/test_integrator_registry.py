from __future__ import annotations

from spatial_fs_benchmark.integration import build_integrator


def test_spatial_integrators_build() -> None:
    names = [
        "cellcharter",
        "combat",
        "gpsa",
        "graphst",
        "paste",
        "pca",
        "scanorama",
        "scvi",
        "staligner",
        "symphony",
    ]
    for name in names:
        integrator = build_integrator(name)
        assert integrator is not None
