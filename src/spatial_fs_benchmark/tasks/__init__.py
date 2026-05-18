"""Task registry."""

from spatial_fs_benchmark.tasks.alignment_eval import AlignmentEvaluationTask
from spatial_fs_benchmark.tasks.base import BenchmarkTask
from spatial_fs_benchmark.tasks.clustering_eval import ClusteringEvaluationTask
from spatial_fs_benchmark.tasks.integration_eval import IntegrationEvaluationTask
from spatial_fs_benchmark.tasks.slice_representation_eval import SliceRepresentationEvaluationTask

TASK_REGISTRY = {
    "alignment_eval": AlignmentEvaluationTask,
    "clustering_eval": ClusteringEvaluationTask,
    "integration_eval": IntegrationEvaluationTask,
    "slice_representation_eval": SliceRepresentationEvaluationTask,
}


def build_task(name: str, **params: object) -> BenchmarkTask:
    task_cls = TASK_REGISTRY[name]
    return task_cls(**params)


__all__ = ["BenchmarkTask", "build_task"]
