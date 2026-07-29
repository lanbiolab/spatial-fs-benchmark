from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import io as spio
from scipy import sparse

from spatial_fs_benchmark.data.preprocess import is_nonnegative_integer_matrix
from spatial_fs_benchmark.data.spatial_object import SpatialDataset
from spatial_fs_benchmark.feature_selection.base import FeatureSelectionResult, FeatureSelector


class ScGCOSelector(FeatureSelector):
    name = "scgco"
    implementation_version = "v1_scgco_1.1.3_isolated"

    def __init__(
        self,
        max_cells_per_slice: int | None = 3000,
        max_genes: int = 2000,
        n_threads: int = 8,
    ) -> None:
        self.max_cells_per_slice = max_cells_per_slice
        self.max_genes = max_genes
        self.n_threads = n_threads
        self.stochastic_selection = max_cells_per_slice is not None
        self._score_cache: dict[tuple[str, int], np.ndarray] = {}

    @staticmethod
    def _project_root() -> Path:
        return Path(__file__).resolve().parents[3]

    def _compute_scores(self, dataset: SpatialDataset, random_seed: int) -> np.ndarray:
        counts = dataset.adata.layers["counts"]
        if not is_nonnegative_integer_matrix(counts):
            source = dataset.adata.uns.get("spatial_fs_benchmark", {}).get("counts_source", "unknown")
            raise ValueError(f"scGCO requires integer counts, but dataset '{dataset.name}' uses '{source}'.")

        root = self._project_root()
        python = root / ".venv-scgco" / "bin" / "python"
        if not python.exists():
            raise FileNotFoundError(f"The scGCO environment is missing: {python}")
        with tempfile.TemporaryDirectory(prefix="scgco_", dir=root / ".cache") as tmpdir:
            input_dir = Path(tmpdir)
            output_path = input_dir / "ranked_features.tsv"
            matrix = counts.T
            spio.mmwrite(input_dir / "counts.mtx", matrix.tocoo() if sparse.issparse(matrix) else sparse.coo_matrix(matrix))
            pd.Series(dataset.adata.var_names.astype(str)).to_csv(
                input_dir / "genes.tsv", sep="\t", index=False, header=False
            )
            coords = dataset.coords
            pd.DataFrame(
                {"slice": dataset.slice_ids, "x": coords[:, 0], "y": coords[:, 1]}
            ).to_csv(input_dir / "observations.tsv", sep="\t", index=False)
            command = [
                str(python),
                str(root / "scripts" / "run_scgco_selector.py"),
                "--input-dir",
                str(input_dir),
                "--output",
                str(output_path),
                "--max-cells-per-slice",
                str(int(self.max_cells_per_slice or 0)),
                "--max-genes",
                str(self.max_genes),
                "--n-threads",
                str(self.n_threads),
                "--seed",
                str(random_seed),
            ]
            completed = subprocess.run(command, cwd=root, check=False, capture_output=True, text=True)
            if completed.returncode != 0:
                raise RuntimeError(
                    f"scGCO failed for dataset '{dataset.name}':\n"
                    f"STDOUT:\n{completed.stdout}\nSTDERR:\n{completed.stderr}"
                )
            ranked = pd.read_csv(output_path, sep="\t")
        score_map = dict(zip(ranked["Feature"].astype(str), ranked["Score"].astype(float), strict=True))
        return np.asarray([score_map.get(str(gene), 0.0) for gene in dataset.adata.var_names], dtype=float)

    def select(self, dataset: SpatialDataset, n_features: int, random_seed: int = 0) -> FeatureSelectionResult:
        cache_key = (dataset.source_path or dataset.name, random_seed)
        if cache_key not in self._score_cache:
            self._score_cache[cache_key] = self._compute_scores(dataset, random_seed)
        scores = self._score_cache[cache_key]
        ranked_indices = np.flatnonzero(scores > 0)
        ranked_indices = ranked_indices[np.argsort(scores[ranked_indices])[::-1]]
        return self._build_ranked_result(
            method_name=self.name,
            gene_names=dataset.adata.var_names.to_numpy(),
            ranked_features=dataset.adata.var_names.to_numpy()[ranked_indices].tolist(),
            ranked_scores=scores[ranked_indices],
            n_features=n_features,
            metadata={
                "max_cells_per_slice": self.max_cells_per_slice,
                "max_genes": self.max_genes,
                "n_threads": self.n_threads,
                "aggregation": "mean within-slice percentile rank",
            },
        )
