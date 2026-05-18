from __future__ import annotations

import pandas as pd


def add_within_task_ranks(frame: pd.DataFrame) -> pd.DataFrame:
    ranked = frame.copy()
    ranked["rank"] = ranked.groupby(["dataset", "task", "metric_name"])["metric_value"].rank(
        ascending=False,
        method="average",
    )
    return ranked
