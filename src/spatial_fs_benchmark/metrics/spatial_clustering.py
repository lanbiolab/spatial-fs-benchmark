from __future__ import annotations

import numpy as np
from sklearn.neighbors import NearestNeighbors


def chaos_score(cluster_labels: np.ndarray, coords: np.ndarray) -> float:
    labels = np.asarray(cluster_labels)
    valid = np.array(
        [label is not None and not (isinstance(label, float) and np.isnan(label)) for label in labels],
        dtype=bool,
    )
    labels = labels[valid]
    coords = coords[valid]
    if coords.shape[0] <= 1:
        return 0.0
    scaled = (coords - coords.mean(axis=0, keepdims=True)) / (coords.std(axis=0, keepdims=True) + 1e-12)
    distance_total = 0.0
    for label in np.unique(labels):
        mask = labels == label
        cluster_coords = scaled[mask]
        if cluster_coords.shape[0] <= 2:
            continue
        model = NearestNeighbors(n_neighbors=2)
        model.fit(cluster_coords)
        distances, _ = model.kneighbors(cluster_coords)
        distance_total += float(distances[:, 1].sum())
    return float(distance_total / len(labels))


def pas_score(cluster_labels: np.ndarray, coords: np.ndarray, n_neighbors: int = 10) -> float:
    if coords.shape[0] <= 1:
        return 0.0
    effective = min(max(2, n_neighbors + 1), coords.shape[0])
    model = NearestNeighbors(n_neighbors=effective)
    model.fit(coords)
    neighbors = model.kneighbors(return_distance=False)[:, 1:]
    labels = np.asarray(cluster_labels)
    disagreements = []
    for idx, row in enumerate(neighbors):
        disagree = np.sum(labels[row] != labels[idx]) > (len(row) / 2.0)
        disagreements.append(float(disagree))
    return float(np.mean(disagreements))
