from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any

import pandas as pd


@dataclass(slots=True)
class MetricRecord:
    dataset: str
    platform: str | None
    n_slices: int
    fs_method: str
    n_features: int
    integration_method: str
    task: str
    metric_name: str
    metric_value: float
    random_seed: int
    runtime: float
    memory_usage: float
    effective_n_features: int | None = None
    notes: str = ""

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def records_to_frame(records: list[MetricRecord]) -> pd.DataFrame:
    return pd.DataFrame([record.to_dict() for record in records])
