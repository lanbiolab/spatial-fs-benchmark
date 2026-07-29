from __future__ import annotations

import numpy as np
import scanpy as sc
from scipy import sparse

from spatial_fs_benchmark.data.spatial_object import SpatialDataset
from spatial_fs_benchmark.data.preprocess import is_nonnegative_integer_matrix
from spatial_fs_benchmark.feature_selection.base import FeatureSelectionResult, FeatureSelector


class ScanpyHVGSelector(FeatureSelector):
    name = "scanpy_hvg"
    implementation_version = "v3_robust_batch_loess"

    def __init__(self, flavor: str = "seurat", batch: bool = False) -> None:
        self.flavor = flavor
        self.batch = batch

    def select(self, dataset: SpatialDataset, n_features: int, random_seed: int = 0) -> FeatureSelectionResult:
        full_adata = self._adata_from_counts(dataset)
        full_gene_names = full_adata.var_names.to_numpy()
        adata = full_adata
        batch_key = dataset.slice_key if self.batch else None
        used_bin_jitter = False
        applied_low_count_filter = False
        batch_fallback_mode = "none"
        excluded_batches: list[str] = []
        batch_fallback_reason: str | None = None
        effective_seurat_v3_span = 0.3

        if self.flavor in {"pearson", "seurat_v3"} and not is_nonnegative_integer_matrix(adata.X):
            source = adata.uns.get("spatial_fs_benchmark", {}).get("counts_source", "unknown")
            raise ValueError(
                f"Scanpy flavor '{self.flavor}' requires integer counts, but dataset '{dataset.name}' "
                f"uses non-count input from '{source}'."
            )

        min_value = float(adata.X.min()) if not sparse.issparse(adata.X) else float(adata.X.min())
        if np.isfinite(min_value) and min_value < 0:
            # Some imaging datasets ship transformed continuous values instead of raw counts.
            # Shift them into a non-negative range so Scanpy HVG routines do not fail on log bins.
            adata.X = adata.X - min_value

        if self.flavor == "cell_ranger":
            matrix = adata.X
            if sparse.issparse(matrix):
                detected = np.asarray((matrix > 0).sum(axis=0)).ravel()
                total = np.asarray(matrix.sum(axis=0)).ravel()
            else:
                dense = np.asarray(matrix)
                detected = np.asarray((dense > 0).sum(axis=0)).ravel()
                total = np.asarray(dense.sum(axis=0)).ravel()
            keep = (total > 0) & (detected > 10)
            if keep.sum() >= min(1000, adata.n_vars):
                adata = adata[:, keep].copy()
                applied_low_count_filter = True

        if self.flavor in {"seurat", "cell_ranger"}:
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)

        if self.flavor == "pearson":
            sc.experimental.pp.highly_variable_genes(
                adata,
                flavor="pearson_residuals",
                n_top_genes=min(n_features, adata.n_vars),
                batch_key=batch_key,
            )
            raw_scores = adata.var.get("variances_norm", adata.var.get("highly_variable_rank"))
        else:
            try:
                sc.pp.highly_variable_genes(
                    adata,
                    flavor=self.flavor,
                    n_top_genes=min(n_features, adata.n_vars),
                    batch_key=batch_key,
                    subset=False,
                )
            except ValueError as exc:
                if self.flavor == "seurat_v3" and batch_key is not None:
                    batch_fallback_reason = str(exc)
                    fitted = False
                    # Preserve all slices whenever possible by increasing the LOESS
                    # smoothing span before considering any slice exclusion.
                    for span in (0.5, 0.75, 1.0):
                        candidate = full_adata.copy()
                        try:
                            sc.pp.highly_variable_genes(
                                candidate,
                                flavor=self.flavor,
                                n_top_genes=min(n_features, candidate.n_vars),
                                batch_key=batch_key,
                                subset=False,
                                span=span,
                            )
                        except ValueError:
                            continue
                        adata = candidate
                        effective_seurat_v3_span = span
                        batch_fallback_mode = "expanded_loess_span"
                        fitted = True
                        break
                    # Very small slices can make the per-batch LOESS fit singular. Drop
                    # the smallest slices one at a time, retaining batch-aware ranking
                    # whenever at least two slices remain.
                    if not fitted:
                        counts = adata.obs[batch_key].astype(str).value_counts().sort_values()
                        ordered_batches = counts.index.tolist()
                        for n_drop in range(1, max(1, len(ordered_batches) - 1)):
                            excluded = ordered_batches[:n_drop]
                            keep = ~adata.obs[batch_key].astype(str).isin(excluded)
                            candidate = adata[keep].copy()
                            try:
                                sc.pp.highly_variable_genes(
                                    candidate,
                                    flavor=self.flavor,
                                    n_top_genes=min(n_features, candidate.n_vars),
                                    batch_key=batch_key,
                                    subset=False,
                                )
                            except ValueError:
                                continue
                            adata = candidate
                            excluded_batches = [str(value) for value in excluded]
                            batch_fallback_mode = "exclude_small_slices"
                            fitted = True
                            break
                    if not fitted:
                        adata = full_adata.copy()
                        sc.pp.highly_variable_genes(
                            adata,
                            flavor=self.flavor,
                            n_top_genes=min(n_features, adata.n_vars),
                            batch_key=None,
                            subset=False,
                        )
                        batch_fallback_mode = "global_seurat_v3"
                elif self.flavor == "cell_ranger" and "Bin edges must be unique" in str(exc):
                    used_bin_jitter = True
                    # Cell Ranger binning can fail on datasets with many tied post-log means.
                    # Apply an infinitesimal gene-wise scaling to break ties without changing ranking materially.
                    scale = 1.0 + np.arange(adata.n_vars, dtype=float) * 1e-12
                    if sparse.issparse(adata.X):
                        adata.X = adata.X.multiply(scale).tocsr()
                    else:
                        adata.X = np.asarray(adata.X, dtype=float) * scale
                    sc.pp.highly_variable_genes(
                        adata,
                        flavor=self.flavor,
                        n_top_genes=min(n_features, adata.n_vars),
                        batch_key=batch_key,
                        subset=False,
                    )
                else:
                    raise
            raw_scores = (
                adata.var.get("dispersions_norm")
                if self.flavor in {"seurat", "cell_ranger"}
                else adata.var.get("highly_variable_rank")
            )

        if raw_scores is None:
            raw_scores = adata.var["highly_variable"].astype(float)
        score_array = np.nan_to_num(np.asarray(raw_scores, dtype=float), nan=0.0, posinf=0.0, neginf=0.0)
        if "highly_variable_rank" in adata.var:
            rank = np.asarray(adata.var["highly_variable_rank"], dtype=float)
            score_array = np.where(np.isfinite(rank), np.nanmax(rank) - rank, score_array)

        if adata.n_vars != len(full_gene_names):
            full_scores = np.zeros(len(full_gene_names), dtype=float)
            selected_index = full_adata.var_names.get_indexer(adata.var_names)
            valid = selected_index >= 0
            full_scores[selected_index[valid]] = score_array[valid]
            score_array = full_scores

        return self._build_result(
            method_name=self.name,
            gene_names=full_gene_names,
            scores=score_array,
            n_features=n_features,
            metadata={
                "flavor": self.flavor,
                "batch": self.batch,
                "used_bin_jitter": used_bin_jitter,
                "applied_low_count_filter": applied_low_count_filter,
                "batch_fallback_mode": batch_fallback_mode,
                "excluded_batches": excluded_batches,
                "batch_fallback_reason": batch_fallback_reason,
                "effective_seurat_v3_span": effective_seurat_v3_span,
            },
        )
