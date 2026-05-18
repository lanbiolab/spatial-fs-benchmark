from __future__ import annotations

import numpy as np

from spatial_fs_benchmark.metrics.batch_mixing import neighbor_indices


def neighbor_overlap_score(reference: np.ndarray, embedding: np.ndarray, n_neighbors: int = 15) -> float:
    reference_neighbors = neighbor_indices(reference, n_neighbors)
    embedding_neighbors = neighbor_indices(embedding, n_neighbors)
    overlaps = []
    for ref_row, emb_row in zip(reference_neighbors, embedding_neighbors, strict=True):
        ref_set = set(ref_row.tolist())
        emb_set = set(emb_row.tolist())
        overlaps.append(len(ref_set & emb_set) / max(1, len(ref_set | emb_set)))
    return float(np.mean(overlaps))


def label_neighbor_purity(embedding: np.ndarray, labels: np.ndarray, n_neighbors: int = 15) -> float:
    neighbors = neighbor_indices(embedding, n_neighbors)
    purity = []
    for idx, row in enumerate(neighbors):
        purity.append(float(np.mean(labels[row] == labels[idx])))
    return float(np.mean(purity))
