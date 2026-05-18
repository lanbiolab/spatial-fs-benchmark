from __future__ import annotations

import numpy as np
from sklearn.cluster import KMeans, MiniBatchKMeans
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, silhouette_score


def cluster_embedding(embedding: np.ndarray, n_clusters: int, random_seed: int = 0) -> np.ndarray:
    if embedding.shape[0] > 10000:
        model = MiniBatchKMeans(
            n_clusters=n_clusters,
            n_init=5,
            batch_size=2048,
            random_state=random_seed,
        )
    else:
        model = KMeans(n_clusters=n_clusters, n_init=20, random_state=random_seed)
    return model.fit_predict(embedding)


def clustering_scores(
    embedding: np.ndarray,
    predicted_labels: np.ndarray,
    true_labels: np.ndarray | None = None,
    random_seed: int = 0,
    silhouette_sample_size: int | None = 10000,
) -> dict[str, float]:
    metrics: dict[str, float] = {}
    if len(np.unique(predicted_labels)) > 1:
        sample_size = None
        if silhouette_sample_size is not None and embedding.shape[0] > silhouette_sample_size:
            sample_size = silhouette_sample_size
        metrics["silhouette"] = float(
            silhouette_score(
                embedding,
                predicted_labels,
                sample_size=sample_size,
                random_state=random_seed,
            )
        )
    else:
        metrics["silhouette"] = 0.0
    if true_labels is not None:
        metrics["ari"] = float(adjusted_rand_score(true_labels, predicted_labels))
        metrics["nmi"] = float(normalized_mutual_info_score(true_labels, predicted_labels))
    return metrics
