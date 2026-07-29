from __future__ import annotations

import json
from pathlib import Path

from spatial_fs_benchmark.config import load_benchmark_config
from spatial_fs_benchmark.feature_selection import build_feature_selector


ROOT = Path(__file__).resolve().parents[1]
CONFIG_DIR = ROOT / "configs" / "rebuild_v1"


def progress_bar(done: int, total: int, width: int = 24) -> str:
    filled = round(width * done / total) if total else width
    return "[" + "#" * filled + "." * (width - filled) + "]"


def main() -> None:
    total_done = 0
    total_expected = 0
    for config_path in sorted(CONFIG_DIR.glob("*_svg_features.yaml")):
        config = load_benchmark_config(config_path)
        dataset = config.datasets[0]
        dataset_dir = ROOT / config.output_dir / dataset.name.lower()
        done = 0
        expected = 0
        for method in config.feature_selection_methods:
            selector = build_feature_selector(method.name, **method.params)
            implementation = getattr(selector, "implementation_version", None)
            seeds = config.seeds if getattr(selector, "stochastic_selection", False) else [0]
            for n_features in method.n_features or config.n_features:
                for seed in seeds:
                    expected += 1
                    metadata_path = (
                        dataset_dir
                        / "feature_selection"
                        / method.name
                        / f"n{n_features}"
                        / f"seed{seed}"
                        / "selected_features.meta.json"
                    )
                    if not metadata_path.exists():
                        continue
                    try:
                        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
                    except (json.JSONDecodeError, OSError):
                        continue
                    if metadata.get("selector_implementation_version") == implementation:
                        done += 1
        total_done += done
        total_expected += expected
        percent = 100 * done / expected if expected else 100
        print(f"{dataset.name:28s} {progress_bar(done, expected)} {done:3d}/{expected:<3d} {percent:6.1f}%")
    percent = 100 * total_done / total_expected if total_expected else 100
    print(f"{'TOTAL':28s} {progress_bar(total_done, total_expected)} {total_done:3d}/{total_expected:<3d} {percent:6.1f}%")


if __name__ == "__main__":
    main()
