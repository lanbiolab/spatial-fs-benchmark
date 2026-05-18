from __future__ import annotations

import argparse
import traceback

from spatial_fs_benchmark.benchmark.runner import BenchmarkRunner
from spatial_fs_benchmark.config import load_benchmark_config


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run the spatial feature-selection benchmark.")
    parser.add_argument("--config", required=True, help="Path to the benchmark YAML config.")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    config = load_benchmark_config(args.config)
    runner = BenchmarkRunner(config)
    try:
        runner.run()
    except Exception as exc:  # noqa: BLE001
        runner.logger.exception("Benchmark failed: %s", exc)
        traceback.print_exc()
        raise


if __name__ == "__main__":
    main()
