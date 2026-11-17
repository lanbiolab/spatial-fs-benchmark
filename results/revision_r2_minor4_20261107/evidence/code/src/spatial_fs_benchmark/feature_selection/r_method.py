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

from spatial_fs_benchmark.data.spatial_object import SpatialDataset
from spatial_fs_benchmark.feature_selection.base import FeatureSelectionResult, FeatureSelector


class RMethodSelector(FeatureSelector):
    name = "r_method"
    implementation_version = "v3_count_matrix_io_seeded_subsampling"

    def __init__(self, method: str, max_cells: int | None = None) -> None:
        self.method = method
        self.max_cells = max_cells
        self.stochastic_selection = max_cells is not None

    @staticmethod
    def _project_root() -> Path:
        return Path(__file__).resolve().parents[3]

    def _prepare_input(self, dataset: SpatialDataset, random_seed: int) -> tuple[np.ndarray, object]:
        adata = self._maybe_subsample_obs(dataset.adata.copy(), self.max_cells, random_seed)
        if "counts" not in adata.layers:
            adata.layers["counts"] = adata.X.copy()
        adata.obs["Batch"] = adata.obs[dataset.slice_key].astype(str).to_numpy()
        if dataset.species:
            adata.uns["Species"] = str(dataset.species).title()
        return adata.var_names.to_numpy(), adata

    def select(self, dataset: SpatialDataset, n_features: int, random_seed: int = 0) -> FeatureSelectionResult:
        gene_names, adata = self._prepare_input(dataset, random_seed)
        project_root = self._project_root()
        runner = project_root / "scripts" / "run_r_feature_selector.R"
        rscript = Path(sys.executable).with_name("Rscript")
        if not rscript.exists():
            raise FileNotFoundError(f"Rscript not found next to Python executable: {rscript}")
        cache_dir = project_root / ".cache"
        cache_dir.mkdir(parents=True, exist_ok=True)

        with tempfile.TemporaryDirectory(prefix=f"{self.method}_", dir=cache_dir) as tmpdir:
            tmpdir_path = Path(tmpdir)
            input_dir = tmpdir_path / "input"
            output_path = tmpdir_path / "features.tsv"
            input_dir.mkdir(parents=True, exist_ok=True)
            counts = adata.layers["counts"] if "counts" in adata.layers else adata.X
            counts = counts.T
            if not sparse.issparse(counts):
                counts = sparse.coo_matrix(np.asarray(counts))
            else:
                counts = counts.tocoo()
            spio.mmwrite(input_dir / "counts.mtx", counts)
            pd.Series(adata.var_names.astype(str)).to_csv(
                input_dir / "genes.tsv",
                sep="\t",
                index=False,
                header=False,
            )
            pd.Series(adata.obs_names.astype(str)).to_csv(
                input_dir / "cells.tsv",
                sep="\t",
                index=False,
                header=False,
            )
            adata.obs.to_csv(input_dir / "obs.tsv", sep="\t", index=True)

            env = os.environ.copy()
            env["RETICULATE_PYTHON"] = sys.executable
            env["R_FUTURE_PLAN"] = "sequential"
            env["R_FUTURE_GLOBALS_MAXSIZE"] = str(8 * 1024**3)
            command = [
                str(rscript),
                str(runner),
                "--input-dir",
                str(input_dir),
                "--output",
                str(output_path),
                "--method",
                self.method,
                "--n-features",
                str(int(n_features)),
            ]
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
                    f"R selector '{self.method}' failed with exit code {completed.returncode}\n"
                    f"STDOUT:\n{completed.stdout}\nSTDERR:\n{completed.stderr}"
                )
            ranking = pd.read_csv(output_path, sep="\t")

        if "Feature" not in ranking.columns:
            raise ValueError(f"R selector '{self.method}' did not write a Feature column.")
        if "Score" in ranking.columns:
            scores = ranking["Score"].astype(float).tolist()
        else:
            scores = np.linspace(len(ranking), 1, num=len(ranking), dtype=float).tolist()
        return self._build_ranked_result(
            method_name=self.name,
            gene_names=gene_names,
            ranked_features=ranking["Feature"].astype(str).tolist(),
            ranked_scores=scores,
            n_features=n_features,
            metadata={
                "method": self.method,
                "max_cells": self.max_cells,
                "implementation_version": self.implementation_version,
            },
        )
