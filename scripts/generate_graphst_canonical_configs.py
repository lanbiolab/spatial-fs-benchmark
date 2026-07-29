from __future__ import annotations

from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = ROOT / "configs/rebuild_v1/downstream/graphst_canonical"
DATASETS = (
    "dlpfc",
    "mouse_brain_serial_sections",
    "stomics_0212",
    "stomics_0218",
    "stomics_0224",
    "e8p5_embryo",
    "e9p5_embryo",
)
COMPETITIVE_METHODS = {
    "Brennecke",
    "anticor",
    "dubstepr",
    "hotspot",
    "morans_i",
    "nbumi",
    "nnsvg",
    "osca",
    "scPNMF",
    "scanpy_cell_ranger",
    "scanpy_cell_ranger_batch",
    "scanpy_pearson",
    "scanpy_pearson_batch",
    "scanpy_seurat",
    "scanpy_seurat_batch",
    "scanpy_seurat_v3",
    "scanpy_seurat_v3_batch",
    "scry",
    "seurat_disp",
    "seurat_mvp",
    "seurat_sct",
    "seurat_vst",
    "singleCellHaystack",
    "somde",
    "sparkx",
    "spatialde",
    "statistic_mean",
    "statistic_variance",
    "triku",
}


def read_yaml(path: Path) -> dict:
    return yaml.safe_load(path.read_text(encoding="utf-8"))


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    for dataset in DATASETS:
        nonspatial = read_yaml(
            ROOT / f"configs/rebuild_v1/non_spatial_downstream/canonical/{dataset}.yaml"
        )
        spatial = read_yaml(ROOT / f"configs/rebuild_v1/downstream/canonical/{dataset}.yaml")
        candidates = {
            method["name"]: method
            for method in nonspatial["feature_selection_methods"]
            + spatial["feature_selection_methods"]
        }
        missing = COMPETITIVE_METHODS.difference(candidates)
        if missing:
            raise ValueError(f"{dataset} is missing competitive methods: {sorted(missing)}")
        methods = [candidates[name] for name in sorted(COMPETITIVE_METHODS)]
        payload = {
            "name": f"{dataset}_graphst_canonical_v1",
            # Reuse the audited feature-selection cache in the primary results tree.
            "output_dir": nonspatial["output_dir"],
            "datasets": nonspatial["datasets"],
            "feature_selection_methods": methods,
            "integration_methods": [
                {
                    "name": "graphst",
                    "params": {
                        "epochs": 600,
                        "dim_output": 32,
                        "pca_components": 20,
                        "n_neighbors": 3,
                        "device": "auto",
                    },
                }
            ],
            "tasks": nonspatial["tasks"],
            "n_features": [2000],
            "seeds": [0, 1, 2],
            "save_embeddings": True,
        }
        destination = OUTPUT_DIR / f"{dataset}.yaml"
        destination.write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")


if __name__ == "__main__":
    main()
