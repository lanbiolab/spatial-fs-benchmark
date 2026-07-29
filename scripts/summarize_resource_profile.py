from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results/resource_profile_v1"


def parse_time(path: Path) -> dict[str, float]:
    text = path.read_text(encoding="utf-8")
    rss = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", text)
    elapsed = re.search(r"Elapsed \(wall clock\) time.*\):\s*([^\n]+)", text)
    if not rss or not elapsed:
        raise ValueError(f"Incomplete GNU time output: {path}")
    fields = [float(value) for value in elapsed.group(1).strip().split(":")]
    wall_seconds = sum(value * 60**power for power, value in enumerate(reversed(fields)))
    return {"wall_seconds_process": wall_seconds, "peak_rss_mib": int(rss.group(1)) / 1024}


def main() -> None:
    rows = []
    failures = []
    methods = sorted(path.stem for path in (RESULTS / "time").glob("*.txt"))
    for method in methods:
        json_path = RESULTS / "json" / f"{method}.json"
        if not json_path.exists():
            failures.append({"method": method, "log": str(RESULTS / f"{method}.log")})
            continue
        row = json.loads(json_path.read_text(encoding="utf-8"))
        row.update(parse_time(RESULTS / "time" / f"{method}.txt"))
        rows.append(row)
    frame = pd.DataFrame(rows).sort_values("wall_seconds_process")
    frame.to_csv(RESULTS / "feature_selection_resource_profile.tsv", sep="\t", index=False)
    pd.DataFrame(failures, columns=["method", "log"]).to_csv(
        RESULTS / "failures.tsv", sep="\t", index=False
    )
    print(f"completed={len(frame)} failures={len(failures)}")


if __name__ == "__main__":
    main()
