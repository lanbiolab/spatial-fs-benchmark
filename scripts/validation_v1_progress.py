#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]
CONFIGS = [
    ROOT / "configs/rebuild_v1/validation/semi_synthetic.yaml",
    ROOT / "configs/rebuild_v1/validation/heldout_slice.yaml",
]


def _status_counts(root: Path) -> tuple[int, int]:
    complete = 0
    failed = 0
    for path in root.glob("*/status/*.json"):
        status = json.loads(path.read_text(encoding="utf-8")).get("status")
        complete += status == "complete"
        failed += status == "failed"
    return complete, failed


def main() -> None:
    for config_path in CONFIGS:
        config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
        scenarios = (
            len(config["prevalence"]) * len(config["seeds"])
            if config["kind"] == "semi_synthetic"
            else len(config["folds"])
        )
        expected = scenarios * len(config["methods"])
        complete, failed = _status_counts(ROOT / config["output_root"])
        finished = complete + failed
        percentage = 100.0 * finished / expected
        width = 30
        filled = round(width * finished / expected)
        bar = "#" * filled + "." * (width - filled)
        print(
            f"{config['kind']:<18} [{bar}] {finished:>3}/{expected:<3} "
            f"({percentage:5.1f}%) complete={complete} failed={failed}"
        )


if __name__ == "__main__":
    main()
