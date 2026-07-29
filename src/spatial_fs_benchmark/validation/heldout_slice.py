from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import rankdata, spearmanr
from sklearn.decomposition import TruncatedSVD

from spatial_fs_benchmark.data.spatial_object import SpatialDataset
from spatial_fs_benchmark.feature_selection.base import FeatureSelectionResult
from spatial_fs_benchmark.feature_selection.svg import SVGSelector
from spatial_fs_benchmark.metrics.clustering_metrics import cluster_embedding, clustering_scores
from spatial_fs_benchmark.metrics.spatial_clustering import chaos_score, pas_score


@dataclass(frozen=True, slots=True)
class HeldoutFold:
    dataset: str
    heldout_slice: str

    @property
    def name(self) -> str:
        safe_slice = self.heldout_slice.replace("/", "_").replace(" ", "_")
        return f"{self.dataset}__heldout_{safe_slice}"


def subset_slices(dataset: SpatialDataset, slices: list[str], name: str) -> SpatialDataset:
    keep = np.isin(dataset.slice_ids, np.asarray(slices, dtype=str))
    if not np.any(keep):
        raise ValueError(f"No observations matched slices {slices!r} in dataset '{dataset.name}'.")
    subset = dataset.adata[keep].copy()
    return SpatialDataset(
        name=name,
        adata=subset,
        slice_key=dataset.slice_key,
        coord_key=dataset.coord_key,
        label_key=dataset.label_key,
        slice_class_key=dataset.slice_class_key,
        platform=dataset.platform,
        species=dataset.species,
        source_path=f"{dataset.source_path}::{name}",
    )


def split_heldout_slice(
    dataset: SpatialDataset,
    heldout_slice: str,
) -> tuple[SpatialDataset, SpatialDataset]:
    available = sorted(np.unique(dataset.slice_ids).tolist())
    if heldout_slice not in available:
        raise ValueError(f"Unknown held-out slice '{heldout_slice}'; available={available}")
    training_slices = [value for value in available if value != heldout_slice]
    if not training_slices:
        raise ValueError("Held-out validation requires at least two slices.")
    training = subset_slices(dataset, training_slices, f"{dataset.name}_training_without_{heldout_slice}")
    heldout = subset_slices(dataset, [heldout_slice], f"{dataset.name}_heldout_{heldout_slice}")
    return training, heldout


def heldout_moran_reference(
    heldout: SpatialDataset,
    random_seed: int = 0,
    n_neighbors: int = 8,
) -> FeatureSelectionResult:
    selector = SVGSelector(n_neighbors=n_neighbors)
    return selector.select(heldout, n_features=heldout.n_vars, random_seed=random_seed)


def _rank_vector(gene_names: list[str], ranking: FeatureSelectionResult) -> np.ndarray:
    score_by_gene = {
        str(gene): float(len(ranking.feature_names) - index)
        for index, gene in enumerate(ranking.feature_names)
    }
    return np.asarray([score_by_gene.get(str(gene), 0.0) for gene in gene_names], dtype=float)


def score_spatial_reproducibility(
    heldout: SpatialDataset,
    training_selection: FeatureSelectionResult,
    heldout_reference: FeatureSelectionResult,
    n_features: int,
) -> dict[str, float]:
    if training_selection.method_name == "all_features":
        return {
            "heldout_topn_jaccard": float("nan"),
            "heldout_moran_percentile": float("nan"),
            "heldout_rank_spearman": float("nan"),
        }
    cutoff = min(int(n_features), heldout.n_vars)
    training_top = [
        str(gene) for gene in training_selection.feature_names if str(gene) in heldout.adata.var_names
    ][:cutoff]
    reference_top = [str(gene) for gene in heldout_reference.feature_names][:cutoff]
    training_set = set(training_top)
    reference_set = set(reference_top)
    union = training_set | reference_set

    heldout_rank = _rank_vector(heldout.gene_names, heldout_reference)
    heldout_percentile = rankdata(heldout_rank, method="average") / len(heldout_rank)
    index_by_gene = {gene: index for index, gene in enumerate(heldout.gene_names)}
    selected_indices = [index_by_gene[gene] for gene in training_top]
    mean_percentile = (
        float(np.mean(heldout_percentile[selected_indices])) if selected_indices else float("nan")
    )

    training_rank = _rank_vector(heldout.gene_names, training_selection)
    correlation = spearmanr(training_rank, heldout_rank).statistic
    return {
        "heldout_topn_jaccard": len(training_set & reference_set) / max(len(union), 1),
        "heldout_moran_percentile": mean_percentile,
        "heldout_rank_spearman": float(correlation),
    }


def _normalized_log_matrix(dataset: SpatialDataset, features: list[str]):
    adata = dataset.adata[:, features]
    matrix = adata.layers["counts"] if "counts" in adata.layers else adata.X
    matrix = sparse.csr_matrix(matrix, dtype=np.float32)
    library_sizes = np.asarray(matrix.sum(axis=1)).ravel()
    scale = 1e4 / np.maximum(library_sizes, 1.0)
    matrix = matrix.multiply(scale[:, None]).tocsr()
    matrix.data = np.log1p(matrix.data)
    return matrix


def score_heldout_clustering(
    heldout: SpatialDataset,
    training_selection: FeatureSelectionResult,
    n_features: int,
    random_seed: int = 0,
    n_components: int = 30,
) -> dict[str, float]:
    available = set(heldout.gene_names)
    selected = [str(gene) for gene in training_selection.feature_names if str(gene) in available]
    if training_selection.method_name != "all_features":
        selected = selected[: min(int(n_features), len(selected))]
    if len(selected) < 2:
        raise ValueError("At least two selected features are required for held-out clustering.")
    matrix = _normalized_log_matrix(heldout, selected)
    n_components = min(int(n_components), matrix.shape[0] - 1, matrix.shape[1] - 1)
    if n_components < 2:
        raise ValueError("Held-out matrix is too small for dimensionality reduction.")
    embedding = TruncatedSVD(n_components=n_components, random_state=random_seed).fit_transform(matrix)

    labels = heldout.labels
    if labels is None:
        raise ValueError(f"Dataset '{heldout.name}' does not provide benchmark labels.")
    valid_labels = np.asarray([str(value) for value in labels])
    n_clusters = len(np.unique(valid_labels))
    predicted = cluster_embedding(embedding, n_clusters=n_clusters, random_seed=random_seed)
    scores = clustering_scores(
        embedding,
        predicted_labels=predicted,
        true_labels=valid_labels,
        random_seed=random_seed,
    )
    scores["CHAOS"] = chaos_score(predicted, heldout.coords)
    scores["PAS"] = pas_score(predicted, heldout.coords, n_neighbors=8)
    return {f"heldout_{metric}": float(value) for metric, value in scores.items()}


def heldout_metrics_frame(
    fold: HeldoutFold,
    method: str,
    spatial_scores: dict[str, float],
    clustering_scores_: dict[str, float],
    n_features: int,
    seed: int,
    runtime_seconds: float,
) -> pd.DataFrame:
    rows = []
    for role, scores in (
        ("spatial_reproducibility", spatial_scores),
        ("heldout_clustering", clustering_scores_),
    ):
        rows.extend(
            {
                "dataset": fold.dataset,
                "heldout_slice": fold.heldout_slice,
                "method": method,
                "n_features": int(n_features),
                "seed": int(seed),
                "evaluation_role": role,
                "metric": metric,
                "value": value,
                "runtime_seconds": float(runtime_seconds),
            }
            for metric, value in scores.items()
        )
    return pd.DataFrame(rows)
