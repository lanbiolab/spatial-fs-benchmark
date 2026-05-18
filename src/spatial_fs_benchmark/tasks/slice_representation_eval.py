from __future__ import annotations

import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

from spatial_fs_benchmark.data.spatial_object import SpatialDataset
from spatial_fs_benchmark.feature_selection.base import FeatureSelectionResult
from spatial_fs_benchmark.integration.base import IntegrationResult
from spatial_fs_benchmark.metrics.clustering_metrics import cluster_embedding
from spatial_fs_benchmark.tasks.base import BenchmarkTask, TaskOutput


def _slice_abundance_matrix(slice_ids: np.ndarray, domain_labels: np.ndarray) -> tuple[np.ndarray, list[str], list[str]]:
    slices = np.unique(slice_ids)
    domains = np.unique(domain_labels)
    abundance = np.zeros((len(slices), len(domains)), dtype=float)
    for i, slice_id in enumerate(slices):
        mask = slice_ids == slice_id
        present, counts = np.unique(domain_labels[mask], return_counts=True)
        domain_to_index = {domain: idx for idx, domain in enumerate(domains)}
        total = counts.sum()
        for label, count in zip(present, counts, strict=True):
            abundance[i, domain_to_index[label]] = count / total
    return abundance, slices.tolist(), domains.tolist()


class SliceRepresentationEvaluationTask(BenchmarkTask):
    name = "slice_representation_eval"

    def __init__(
        self,
        n_clusters: int = 6,
        domain_n_clusters: int | None = None,
        slice_n_clusters: int | None = None,
        cluster_method: str = "kmeans",
        distance: str = "euclidean",
        use_true_domains_for_abundance: bool = False,
    ) -> None:
        self.n_clusters = n_clusters
        self.domain_n_clusters = domain_n_clusters if domain_n_clusters is not None else n_clusters
        self.slice_n_clusters = slice_n_clusters if slice_n_clusters is not None else n_clusters
        self.cluster_method = cluster_method
        self.distance = distance
        self.use_true_domains_for_abundance = use_true_domains_for_abundance

    def evaluate(
        self,
        dataset: SpatialDataset,
        selected_features: FeatureSelectionResult,
        integration_result: IntegrationResult,
        random_seed: int = 0,
    ) -> TaskOutput:
        if self.use_true_domains_for_abundance and dataset.labels is not None:
            domains = dataset.labels
            domain_source = "ground_truth"
        else:
            domains = cluster_embedding(
                integration_result.embedding,
                self.domain_n_clusters,
                random_seed=random_seed,
            ).astype(str)
            domain_source = "predicted_from_embedding"
        abundance, slice_order, domain_order = _slice_abundance_matrix(dataset.slice_ids, domains)
        if abundance.shape[0] <= 1:
            predicted = np.ones(abundance.shape[0], dtype=int)
        elif self.cluster_method == "hclust":
            linkage_matrix = linkage(abundance, method="complete", metric=self.distance)
            predicted = fcluster(
                linkage_matrix,
                t=min(self.slice_n_clusters, abundance.shape[0]),
                criterion="maxclust",
            )
        else:
            predicted = KMeans(
                n_clusters=min(self.slice_n_clusters, abundance.shape[0]),
                n_init=20,
                random_state=random_seed,
            ).fit_predict(abundance)
        metrics = {
            "slice_repr_ARI": float("nan"),
            "slice_repr_NMI": float("nan"),
        }
        if dataset.slice_classes is not None:
            slice_truth = []
            for slice_id in slice_order:
                mask = dataset.slice_ids == slice_id
                slice_truth.append(dataset.slice_classes[mask][0])
            metrics["slice_repr_ARI"] = float(adjusted_rand_score(slice_truth, predicted))
            metrics["slice_repr_NMI"] = float(normalized_mutual_info_score(slice_truth, predicted))
        else:
            distances = pdist(abundance, metric=self.distance) if abundance.shape[0] > 1 else np.array([0.0])
            metrics["slice_repr_distance_mean"] = float(np.mean(distances))
        artifacts = {
            "slice_order": slice_order,
            "domain_order": domain_order,
            "domain_source": domain_source,
            "domain_n_clusters": self.domain_n_clusters,
            "slice_n_clusters": self.slice_n_clusters,
            "abundance_matrix": abundance.tolist(),
            "spot_domain_labels": domains.tolist(),
            "predicted_slice_clusters": predicted.tolist(),
        }
        return TaskOutput(metrics=metrics, artifacts=artifacts)
