import numpy as np

from spatial_fs_benchmark.metrics.alignment_metrics import alignment_accuracy_ratio
from spatial_fs_benchmark.metrics.integration_paper import scaled_silhouette_batch


def test_alignment_uses_only_explicit_pairs() -> None:
    representation = np.array([[0.0], [10.0], [100.0], [110.0], [0.1], [10.1]])
    labels = np.array(["a", "b", "x", "y", "a", "b"])
    slices = np.array(["left", "left", "unrelated", "unrelated", "right", "right"])

    accuracy, ratio = alignment_accuracy_ratio(
        representation,
        labels,
        slices,
        alignment_pairs=[("left", "right")],
    )

    assert accuracy == 1.0
    assert ratio == 0.0


def test_unlabelled_batch_asw_is_finite() -> None:
    embedding = np.array([[0.0], [0.2], [0.1], [0.3]])
    batches = np.array(["a", "a", "b", "b"])

    score = scaled_silhouette_batch(embedding, batches, labels=None)

    assert np.isfinite(score)
