from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

from spatial_fs_benchmark.config import load_benchmark_config


ROOT = Path(__file__).resolve().parents[1]
CONFIG_ROOT = ROOT / "configs" / "rebuild_v1" / "downstream"
INTEGRATOR_VERSIONS = {
    "scvi": "v2_integer_counts_batch_covariate",
    "cellcharter": "v2_native_scvi_spatial_aggregation_counts",
}


def bar(done: int, total: int, width: int = 24) -> str:
    filled = round(width * done / total) if total else width
    return "[" + "#" * filled + "." * (width - filled) + "]"


def main() -> None:
    expected: dict[str, set[tuple[Path, str]]] = defaultdict(set)
    for config_path in sorted(CONFIG_ROOT.glob("*/*.yaml")):
        config = load_benchmark_config(config_path)
        dataset = config.datasets[0]
        dataset_dir = ROOT / config.output_dir / dataset.name.lower()
        for method in config.feature_selection_methods:
            for n_features in method.n_features or config.n_features:
                for seed in config.seeds:
                    for integrator in config.integration_methods:
                        path = (
                            dataset_dir
                            / integrator.name
                            / method.name
                            / f"n{n_features}"
                            / f"seed{seed}"
                            / "embedding.meta.json"
                        )
                        expected[dataset.name].add((path, integrator.name))

    total_done = 0
    total_expected = 0
    for dataset_name in sorted(expected):
        done = 0
        for path, integrator_name in expected[dataset_name]:
            if not path.exists():
                continue
            try:
                payload = json.loads(path.read_text(encoding="utf-8"))
            except (json.JSONDecodeError, OSError):
                continue
            if payload.get("integrator_implementation_version") == INTEGRATOR_VERSIONS[integrator_name]:
                done += 1
        count = len(expected[dataset_name])
        total_done += done
        total_expected += count
        pct = 100 * done / count if count else 100
        print(f"{dataset_name:28s} {bar(done, count)} {done:3d}/{count:<3d} {pct:6.1f}%")
    pct = 100 * total_done / total_expected if total_expected else 100
    print(f"{'TOTAL':28s} {bar(total_done, total_expected)} {total_done:3d}/{total_expected:<3d} {pct:6.1f}%")


if __name__ == "__main__":
    main()
