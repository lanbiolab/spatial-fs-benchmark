from __future__ import annotations

import numpy as np
from sklearn.metrics import silhouette_score
from sklearn.neighbors import NearestNeighbors


def neighbor_indices(matrix: np.ndarray, n_neighbors: int) -> np.ndarray:
    n_obs = matrix.shape[0]
    effective_neighbors = min(max(2, n_neighbors + 1), n_obs)
    model = NearestNeighbors(n_neighbors=effective_neighbors)
    model.fit(matrix)
    indices = model.kneighbors(return_distance=False)
    return indices[:, 1:]


def slice_mixing_entropy(embedding: np.ndarray, slice_ids: np.ndarray, n_neighbors: int = 15) -> float:
    neighbors = neighbor_indices(embedding, n_neighbors)
    entropies = []
    for row in neighbors:
        neighbor_labels = slice_ids[row]
        values, counts = np.unique(neighbor_labels, return_counts=True)
        probabilities = counts / counts.sum()
        entropy = -(probabilities * np.log(probabilities + 1e-12)).sum()
        if len(values) > 1:
            entropy /= np.log(len(values))
        entropies.append(float(entropy))
    return float(np.mean(entropies))


def inverse_slice_silhouette(embedding: np.ndarray, slice_ids: np.ndarray) -> float:
    if len(np.unique(slice_ids)) < 2:
        return 0.0
    raw = silhouette_score(embedding, slice_ids)
    return float(1.0 - ((raw + 1.0) / 2.0))
