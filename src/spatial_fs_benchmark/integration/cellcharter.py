from __future__ import annotations

import json
import subprocess
import tempfile
from hashlib import md5
from pathlib import Path

import numpy as np

from spatial_fs_benchmark.data.spatial_object import SpatialDataset
from spatial_fs_benchmark.feature_selection.base import FeatureSelectionResult
from spatial_fs_benchmark.integration.base import IntegrationResult, SpatialIntegrator


def _default_cellcharter_python() -> Path:
    return Path(__file__).resolve().parents[3] / ".conda-env-cellcharter-clone" / "bin" / "python"


class CellCharterIntegrator(SpatialIntegrator):
    name = "cellcharter"

    def __init__(
        self,
        n_latent: int = 30,
        nhood_layers: int = 4,
        spatial_neighbors: int | None = None,
        max_epochs: int = 100,
        cellcharter_python: str | None = None,
    ) -> None:
        self.n_latent = n_latent
        self.nhood_layers = nhood_layers
        self.spatial_neighbors = spatial_neighbors
        self.max_epochs = max_epochs
        self.cellcharter_python = Path(cellcharter_python) if cellcharter_python is not None else _default_cellcharter_python()

    def fit_transform(
        self,
        dataset: SpatialDataset,
        selected_features: FeatureSelectionResult,
        random_seed: int = 0,
    ) -> IntegrationResult:
        subset = dataset.subset_features(selected_features.feature_names)
        if not self.cellcharter_python.exists():
            raise FileNotFoundError(
                f"CellCharter native environment python was not found at '{self.cellcharter_python}'."
            )

        with tempfile.TemporaryDirectory(prefix="cellcharter_", dir=Path.cwd()) as tmpdir:
            tmpdir_path = Path(tmpdir)
            input_path = tmpdir_path / "input.h5ad"
            output_path = tmpdir_path / "embedding.npy"
            metadata_path = tmpdir_path / "metadata.json"

            adata = subset.adata.copy()
            if "counts" not in adata.layers:
                adata.layers["counts"] = adata.X.copy()
            adata.write_h5ad(input_path)

            runner = Path(__file__).resolve().parents[3] / "scripts" / "run_cellcharter_native.py"
            command = [
                str(self.cellcharter_python),
                str(runner),
                "--input",
                str(input_path),
                "--output",
                str(output_path),
                "--metadata-output",
                str(metadata_path),
                "--slice-key",
                subset.slice_key,
                "--coord-key",
                subset.coord_key,
                "--n-latent",
                str(self.n_latent),
                "--nhood-layers",
                str(self.nhood_layers),
                "--max-epochs",
                str(self.max_epochs),
                "--seed",
                str(random_seed),
            ]
            completed = subprocess.run(command, capture_output=True, text=True, check=False)
            if completed.returncode != 0:
                raise RuntimeError(
                    "Native CellCharter execution failed.\n"
                    f"stdout:\n{completed.stdout}\n"
                    f"stderr:\n{completed.stderr}"
                )

            embedding = np.load(output_path).astype(np.float32)
            metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
            metadata.update(
                {
                    "n_input_features": int(adata.n_vars),
                    "input_feature_hash": md5("\n".join(selected_features.feature_names).encode("utf-8")).hexdigest(),
                    "cellcharter_python": str(self.cellcharter_python),
                    "spatial_neighbors": self.spatial_neighbors,
                }
            )
            return IntegrationResult(
                method_name=self.name,
                embedding=embedding,
                slice_embedding=self._build_slice_embedding(dataset, embedding),
                metadata=metadata,
            )
