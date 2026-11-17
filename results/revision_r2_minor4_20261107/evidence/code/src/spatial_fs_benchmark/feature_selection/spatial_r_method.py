from __future__ import annotations

import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import io as spio
from scipy import sparse

from spatial_fs_benchmark.data.preprocess import is_nonnegative_integer_matrix
from spatial_fs_benchmark.data.spatial_object import SpatialDataset
from spatial_fs_benchmark.feature_selection.base import FeatureSelectionResult, FeatureSelector


class SpatialRMethodSelector(FeatureSelector):
    implementation_version = "v2_official_preprocessing_slice_wise"
    supported_methods = {"sparkx", "nnsvg"}

    def __init__(
        self,
        method: str,
        max_cells_per_slice: int | None = None,
        max_genes: int | None = None,
        n_threads: int = 8,
    ) -> None:
        method = method.lower()
        if method not in self.supported_methods:
            raise ValueError(f"Unsupported spatial R selector: {method}")
        self.method = method
        self.name = method
        self.max_cells_per_slice = max_cells_per_slice
        self.max_genes = max_genes
        self.n_threads = n_threads
        self.stochastic_selection = max_cells_per_slice is not None
        self._score_cache: dict[tuple[str, int], np.ndarray] = {}

    @staticmethod
    def _project_root() -> Path:
        return Path(__file__).resolve().parents[3]

    @staticmethod
    def _write_matrix(path: Path, matrix) -> None:
        transposed = matrix.T
        if not sparse.issparse(transposed):
            transposed = sparse.coo_matrix(np.asarray(transposed))
        else:
            transposed = transposed.tocoo()
        spio.mmwrite(path, transposed)

    def _compute_scores(self, dataset: SpatialDataset, random_seed: int) -> np.ndarray:
        adata = dataset.adata
        if not is_nonnegative_integer_matrix(adata.layers["counts"]):
            source = adata.uns.get("spatial_fs_benchmark", {}).get("counts_source", "unknown")
            raise ValueError(
                f"{self.method} requires integer counts, but dataset '{dataset.name}' uses '{source}'."
            )

        project_root = self._project_root()
        cache_dir = project_root / ".cache"
        cache_dir.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix=f"{self.method}_", dir=cache_dir) as tmpdir:
            input_dir = Path(tmpdir)
            output_path = input_dir / "ranked_features.tsv"
            matrix = adata.layers["counts"]
            self._write_matrix(input_dir / "matrix.mtx", matrix)
            pd.Series(adata.var_names.astype(str)).to_csv(
                input_dir / "genes.tsv", sep="\t", index=False, header=False
            )
            coords = dataset.coords
            pd.DataFrame(
                {
                    "slice": dataset.slice_ids,
                    "x": coords[:, 0],
                    "y": coords[:, 1],
                }
            ).to_csv(input_dir / "observations.tsv", sep="\t", index=False)

            rscript = Path(sys.executable).with_name("Rscript")
            runner = project_root / "scripts" / "run_spatial_r_feature_selector.R"
            command = [
                str(rscript),
                str(runner),
                str(input_dir),
                str(output_path),
                self.method,
                str(int(self.max_cells_per_slice or 0)),
                str(int(self.max_genes or 0)),
                str(int(self.n_threads)),
                str(int(random_seed)),
            ]
            env = os.environ.copy()
            env["OMP_NUM_THREADS"] = str(self.n_threads)
            env["OPENBLAS_NUM_THREADS"] = str(self.n_threads)
            completed = subprocess.run(
                command,
                cwd=project_root,
                env=env,
                check=False,
                capture_output=True,
                text=True,
            )
            if completed.returncode != 0:
                raise RuntimeError(
                    f"{self.method} failed for dataset '{dataset.name}':\n"
                    f"STDOUT:\n{completed.stdout}\nSTDERR:\n{completed.stderr}"
                )
            ranked = pd.read_csv(output_path, sep="\t")

        score_by_gene = dict(zip(ranked["Feature"].astype(str), ranked["Score"].astype(float), strict=True))
        return np.asarray([score_by_gene.get(str(gene), 0.0) for gene in adata.var_names], dtype=float)

    def select(self, dataset: SpatialDataset, n_features: int, random_seed: int = 0) -> FeatureSelectionResult:
        cache_key = (dataset.source_path or dataset.name, random_seed)
        if cache_key not in self._score_cache:
            self._score_cache[cache_key] = self._compute_scores(dataset, random_seed)
        scores = self._score_cache[cache_key]
        ranked_indices = np.flatnonzero(scores > 0)
        ranked_indices = ranked_indices[np.argsort(scores[ranked_indices])[::-1]]
        ranked_features = dataset.adata.var_names.to_numpy()[ranked_indices].tolist()
        ranked_scores = scores[ranked_indices]
        return self._build_ranked_result(
            method_name=self.method,
            gene_names=dataset.adata.var_names.to_numpy(),
            ranked_features=ranked_features,
            ranked_scores=ranked_scores,
            n_features=n_features,
            metadata={
                "max_cells_per_slice": self.max_cells_per_slice,
                "max_genes": self.max_genes,
                "n_threads": self.n_threads,
                "aggregation": "mean within-slice percentile rank",
                "input_assay": "counts",
                "official_preprocessing": (
                    "SPARK-X mixture kernels with mitochondrial genes removed"
                    if self.method == "sparkx"
                    else "nnSVG filter_genes plus library-size normalization and logcounts"
                ),
            },
        )
