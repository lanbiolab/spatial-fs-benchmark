from __future__ import annotations

import numpy as np
from sklearn.neighbors import kneighbors_graph


def spatial_coherence_score(labels: np.ndarray, coords: np.ndarray, n_neighbors: int = 8) -> float:
    graph = kneighbors_graph(
        coords,
        n_neighbors=min(max(1, n_neighbors), max(1, coords.shape[0] - 1)),
        mode="connectivity",
        include_self=False,
    )
    adjacency = graph.tocoo()
    if adjacency.nnz == 0:
        return 0.0
    matches = labels[adjacency.row] == labels[adjacency.col]
    return float(np.mean(matches))


def morans_i(values: np.ndarray, coords: np.ndarray, n_neighbors: int = 8) -> float:
    graph = kneighbors_graph(
        coords,
        n_neighbors=min(max(1, n_neighbors), max(1, coords.shape[0] - 1)),
        mode="connectivity",
        include_self=False,
    ).tocsr()
    centered = values - values.mean()
    denominator = float(np.square(centered).sum()) + 1e-12
    numerator = float(np.sum(centered * (graph @ centered)))
    return float((len(values) / float(graph.sum())) * (numerator / denominator))
