from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml


@dataclass(slots=True)
class MethodConfig:
    name: str
    params: dict[str, Any] = field(default_factory=dict)


@dataclass(slots=True)
class DatasetConfig:
    name: str
    path: str
    sample_data_dir: str | None = None
    platform: str | None = None
    species: str | None = None
    slice_key: str = "slice_id"
    label_key: str | None = None
    slice_class_key: str | None = None
    coord_key: str = "spatial"
    preprocess: dict[str, Any] = field(default_factory=dict)


@dataclass(slots=True)
class BenchmarkConfig:
    name: str
    output_dir: str
    datasets: list[DatasetConfig]
    feature_selection_methods: list[MethodConfig]
    integration_methods: list[MethodConfig]
    tasks: list[MethodConfig]
    n_features: list[int]
    seeds: list[int]
    save_embeddings: bool = False


def _read_yaml(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def load_dataset_config(path: str | Path) -> DatasetConfig:
    config_path = Path(path)
    data = _read_yaml(config_path)
    for key in ("path", "sample_data_dir"):
        value = data.get(key)
        if value is not None and not Path(value).is_absolute():
            relative_to_config = (config_path.parent / value).resolve()
            relative_to_cwd = Path(value).resolve()
            data[key] = str(relative_to_config if relative_to_config.exists() else relative_to_cwd)
    return DatasetConfig(**data)


def load_benchmark_config(path: str | Path) -> BenchmarkConfig:
    config_path = Path(path)
    data = _read_yaml(config_path)
    base_dir = config_path.parent
    dataset_configs = []
    for dataset_entry in data["datasets"]:
        dataset_path = Path(dataset_entry)
        if not dataset_path.is_absolute():
            candidate = (base_dir / dataset_path).resolve()
            dataset_path = candidate if candidate.exists() else Path(dataset_entry).resolve()
        dataset_configs.append(load_dataset_config(dataset_path))
    return BenchmarkConfig(
        name=data["name"],
        output_dir=data["output_dir"],
        datasets=dataset_configs,
        feature_selection_methods=[MethodConfig(**item) for item in data["feature_selection_methods"]],
        integration_methods=[MethodConfig(**item) for item in data["integration_methods"]],
        tasks=[MethodConfig(**item) for item in data["tasks"]],
        n_features=list(data["n_features"]),
        seeds=list(data["seeds"]),
        save_embeddings=bool(data.get("save_embeddings", False)),
    )
