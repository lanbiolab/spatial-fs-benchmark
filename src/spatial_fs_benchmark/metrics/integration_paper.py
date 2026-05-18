from __future__ import annotations

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from sklearn.metrics import silhouette_score
from sklearn.neighbors import NearestNeighbors


def _neighbors(matrix: np.ndarray, n_neighbors: int) -> np.ndarray:
    effective = min(max(2, n_neighbors + 1), matrix.shape[0])
    model = NearestNeighbors(n_neighbors=effective)
    model.fit(matrix)
    return model.kneighbors(return_distance=False)[:, 1:]


def scaled_silhouette_labels(embedding: np.ndarray, labels: np.ndarray) -> float:
    if len(np.unique(labels)) < 2:
        return float("nan")
    sample_size = min(10000, embedding.shape[0])
    kwargs = {}
    if embedding.shape[0] > sample_size:
        kwargs["sample_size"] = sample_size
        kwargs["random_state"] = 0
    return float((silhouette_score(embedding, labels, **kwargs) + 1.0) / 2.0)


def scaled_silhouette_batch(embedding: np.ndarray, batch_ids: np.ndarray, labels: np.ndarray) -> float:
    unique_labels = np.unique(labels)
    scores = []
    for label in unique_labels:
        mask = labels == label
        if mask.sum() < 3:
            continue
        label_batches = batch_ids[mask]
        if len(np.unique(label_batches)) < 2:
            continue
        label_embedding = embedding[mask]
        sample_size = min(10000, label_embedding.shape[0])
        kwargs = {}
        if label_embedding.shape[0] > sample_size:
            kwargs["sample_size"] = sample_size
            kwargs["random_state"] = 0
        score = silhouette_score(label_embedding, label_batches, **kwargs)
        scores.append(1.0 - abs(float(score)))
    if not scores:
        return float("nan")
    return float(np.mean(scores))


def _inverse_simpson(neighbor_labels: np.ndarray) -> float:
    _, counts = np.unique(neighbor_labels, return_counts=True)
    probabilities = counts / counts.sum()
    return float(1.0 / np.square(probabilities).sum())


def scaled_lisi(
    embedding: np.ndarray,
    values: np.ndarray,
    n_neighbors: int,
    mode: str,
) -> float:
    neighbors = _neighbors(embedding, n_neighbors)
    scores = np.array([_inverse_simpson(values[row]) for row in neighbors], dtype=float)
    n_classes = max(1, len(np.unique(values)))
    if n_classes == 1:
        return 1.0
    if mode == "batch":
        scaled = (scores - 1.0) / (n_classes - 1.0)
    elif mode == "label":
        scaled = 1.0 - ((scores - 1.0) / (n_classes - 1.0))
    else:
        raise ValueError(f"Unsupported LISI mode: {mode}")
    return float(np.clip(np.mean(scaled), 0.0, 1.0))


def isolated_label_f1(
    embedding: np.ndarray,
    labels: np.ndarray,
    batch_ids: np.ndarray,
    n_neighbors: int,
) -> float:
    label_batches = {label: np.unique(batch_ids[labels == label]) for label in np.unique(labels)}
    isolated = [label for label, batches in label_batches.items() if len(batches) == 1]
    target_labels = isolated if isolated else np.unique(labels).tolist()
    neighbors = _neighbors(embedding, n_neighbors)
    f1_scores = []
    for label in target_labels:
        positives = labels == label
        predicted = np.array([np.any(labels[row] == label) for row in neighbors], dtype=bool)
        tp = int(np.sum(predicted & positives))
        fp = int(np.sum(predicted & ~positives))
        fn = int(np.sum(~predicted & positives))
        precision = tp / (tp + fp + 1e-12)
        recall = tp / (tp + fn + 1e-12)
        f1_scores.append((2.0 * precision * recall) / (precision + recall + 1e-12))
    return float(np.mean(f1_scores))


def graph_connectivity_score(
    embedding: np.ndarray,
    labels: np.ndarray,
    n_neighbors: int,
) -> float:
    graph = NearestNeighbors(n_neighbors=min(max(2, n_neighbors + 1), embedding.shape[0]))
    graph.fit(embedding)
    neighbor_ids = graph.kneighbors(return_distance=False)[:, 1:]
    n_obs = embedding.shape[0]
    row_index = np.repeat(np.arange(n_obs), neighbor_ids.shape[1])
    col_index = neighbor_ids.reshape(-1)
    adjacency = csr_matrix((np.ones_like(row_index), (row_index, col_index)), shape=(n_obs, n_obs))
    adjacency = adjacency.maximum(adjacency.transpose())
    scores = []
    for label in np.unique(labels):
        mask = labels == label
        if mask.sum() <= 1:
            scores.append(1.0)
            continue
        subgraph = adjacency[mask][:, mask]
        _, components = connected_components(subgraph, directed=False)
        largest = np.bincount(components).max()
        scores.append(largest / mask.sum())
    return float(np.mean(scores))
