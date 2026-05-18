from __future__ import annotations

import json
import resource
import time
from hashlib import md5
from pathlib import Path

import numpy as np

from spatial_fs_benchmark.benchmark.result_schema import MetricRecord, records_to_frame
from spatial_fs_benchmark.data.spatial_object import SpatialDataset
from spatial_fs_benchmark.feature_selection.base import FeatureSelectionResult
from spatial_fs_benchmark.integration.base import IntegrationResult
from spatial_fs_benchmark.tasks.base import TaskOutput


def max_memory_mb() -> float:
    return float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0)


def run_timed_block(func):
    start = time.perf_counter()
    output = func()
    runtime = time.perf_counter() - start
    return output, runtime


def cache_signature(payload: dict) -> str:
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return md5(canonical.encode("utf-8")).hexdigest()


def _to_jsonable(value):
    if isinstance(value, dict):
        return {str(key): _to_jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_to_jsonable(item) for item in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    return value


def write_feature_selection(path: str | Path, selection: FeatureSelectionResult) -> None:
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "method_name": selection.method_name,
        "feature_names": selection.feature_names,
        "feature_indices": selection.feature_indices,
        "scores": selection.scores,
        "metadata": selection.metadata,
    }
    output_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def read_feature_selection(path: str | Path) -> FeatureSelectionResult:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    return FeatureSelectionResult(
        method_name=payload["method_name"],
        feature_names=list(payload["feature_names"]),
        feature_indices=list(payload["feature_indices"]),
        scores=list(payload["scores"]),
        metadata=dict(payload.get("metadata", {})),
    )


def write_integration_result(path: str | Path, dataset: SpatialDataset, result: IntegrationResult) -> None:
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    metadata = _to_jsonable(result.metadata)
    np.savez_compressed(
        output_path,
        embedding=result.embedding,
        corrected_matrix=result.corrected_matrix if result.corrected_matrix is not None else np.array([], dtype=np.float32),
        metadata_json=json.dumps(metadata, sort_keys=True),
        obs_names=np.asarray(dataset.adata.obs_names.astype(str)),
        slice_ids=dataset.slice_ids,
        labels=dataset.labels if dataset.labels is not None else np.array([], dtype=str),
    )


def read_integration_result(path: str | Path, method_name: str) -> IntegrationResult:
    payload = np.load(Path(path), allow_pickle=False)
    corrected = payload["corrected_matrix"]
    metadata_json = payload["metadata_json"].item() if payload["metadata_json"].shape == () else str(payload["metadata_json"])
    return IntegrationResult(
        method_name=method_name,
        embedding=payload["embedding"],
        corrected_matrix=None if corrected.size == 0 else corrected,
        metadata=json.loads(metadata_json) if metadata_json else {},
    )


def write_metric_records(path: str | Path, records: list[MetricRecord]) -> None:
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(records_to_frame(records).to_json(orient="records", indent=2), encoding="utf-8")


def read_metric_records(path: str | Path) -> list[MetricRecord]:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    return [MetricRecord(**row) for row in payload]


def write_json(path: str | Path, payload: dict) -> None:
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def read_json(path: str | Path) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def metric_records_from_task_output(
    dataset: SpatialDataset,
    fs_method: str,
    n_features: int,
    effective_n_features: int | None,
    integration_method: str,
    task_name: str,
    task_output: TaskOutput,
    random_seed: int,
    runtime: float,
    notes: str = "",
) -> list[MetricRecord]:
    memory_usage = max_memory_mb()
    return [
        MetricRecord(
            dataset=dataset.name,
            platform=dataset.platform,
            n_slices=len(np.unique(dataset.slice_ids)),
            fs_method=fs_method,
            n_features=n_features,
            integration_method=integration_method,
            task=task_name,
            metric_name=metric_name,
            metric_value=float(metric_value),
            random_seed=random_seed,
            runtime=runtime,
            memory_usage=memory_usage,
            effective_n_features=effective_n_features,
            notes=notes,
        )
        for metric_name, metric_value in task_output.metrics.items()
    ]
