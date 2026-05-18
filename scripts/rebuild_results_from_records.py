from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

ROOT_DIR = Path(__file__).resolve().parents[1]
SRC_DIR = ROOT_DIR / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from spatial_fs_benchmark.benchmark.experiment import read_metric_records
from spatial_fs_benchmark.benchmark.result_schema import records_to_frame


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Rebuild benchmark results.csv files from *_records.json payloads."
    )
    parser.add_argument(
        "--results-root",
        required=True,
        help="Root directory containing per-benchmark output directories.",
    )
    parser.add_argument(
        "--include-glob",
        default="*_spatial_main",
        help="Glob used to select per-benchmark result directories under results-root.",
    )
    parser.add_argument(
        "--merged-output",
        default="",
        help="Optional merged CSV path built from every rebuilt results.csv.",
    )
    return parser


def rebuild_one(results_dir: Path) -> tuple[Path | None, int]:
    records = []
    for path in sorted(results_dir.rglob("*_records.json")):
        records.extend(read_metric_records(path))
    if not records:
        return None, 0
    frame = records_to_frame(records)
    mask = (frame["fs_method"] == "all_features") & (frame["n_features"].astype(str) == "1")
    frame.loc[mask, "n_features"] = "all"
    output = results_dir / "results.csv"
    output.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output, index=False)
    return output, len(frame)


def main() -> None:
    args = build_parser().parse_args()
    results_root = Path(args.results_root)
    result_dirs = sorted(p for p in results_root.glob(args.include_glob) if p.is_dir())

    rebuilt_outputs: list[Path] = []
    for results_dir in result_dirs:
        output, row_count = rebuild_one(results_dir)
        if output is None:
            print(f"skipped\t{results_dir}\t0")
            continue
        rebuilt_outputs.append(output)
        print(f"rebuilt\t{output}\t{row_count}")

    if args.merged_output:
        frames = [pd.read_csv(path) for path in rebuilt_outputs if path.exists()]
        merged = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
        merged_output = Path(args.merged_output)
        merged_output.parent.mkdir(parents=True, exist_ok=True)
        merged.to_csv(merged_output, index=False)
        print(f"merged\t{merged_output}\t{len(merged)}")


if __name__ == "__main__":
    main()
