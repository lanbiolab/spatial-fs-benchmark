from __future__ import annotations

import numpy as np
from scipy.spatial.distance import cdist


def alignment_accuracy_ratio(
    representation: np.ndarray,
    labels: np.ndarray,
    slice_ids: np.ndarray,
    alignment_pairs: list[tuple[str, str]] | None = None,
    max_points_per_slice: int | None = None,
    random_seed: int = 0,
) -> tuple[float, float]:
    pairs = list(alignment_pairs or [])
    if not pairs:
        return float("nan"), float("nan")
    available_slices = set(np.unique(slice_ids).astype(str))
    missing = sorted({slice_id for pair in pairs for slice_id in pair} - available_slices)
    if missing:
        raise ValueError(f"Configured alignment slices are absent from the dataset: {missing}")
    rng = np.random.default_rng(random_seed)
    accuracies = []
    ratios = []
    slice_ids = slice_ids.astype(str)
    for left, right in pairs:
        left_mask = slice_ids == left
        right_mask = slice_ids == right
        left_repr = representation[left_mask]
        right_repr = representation[right_mask]
        left_labels = labels[left_mask]
        right_labels = labels[right_mask]
        if max_points_per_slice is not None and len(left_repr) > max_points_per_slice:
            keep = rng.choice(len(left_repr), size=max_points_per_slice, replace=False)
            left_repr = left_repr[keep]
            left_labels = left_labels[keep]
        if max_points_per_slice is not None and len(right_repr) > max_points_per_slice:
            keep = rng.choice(len(right_repr), size=max_points_per_slice, replace=False)
            right_repr = right_repr[keep]
            right_labels = right_labels[keep]
        distance_matrix = cdist(left_repr, right_repr)
        left_to_right = np.argmin(distance_matrix, axis=1)
        right_to_left = np.argmin(distance_matrix, axis=0)
        same_rate_lr = np.mean(left_labels == right_labels[left_to_right])
        same_rate_rl = np.mean(left_labels[right_to_left] == right_labels)
        accuracies.append(float(np.mean([same_rate_lr, same_rate_rl])))
        denominator = min(len(left_labels), len(right_labels))
        ratio_lr = abs(np.log2(denominator / max(1, len(np.unique(left_to_right)))))
        ratio_rl = abs(np.log2(denominator / max(1, len(np.unique(right_to_left)))))
        ratios.append(float(np.mean([ratio_lr, ratio_rl])))
    return float(np.mean(accuracies)), float(np.mean(ratios))
