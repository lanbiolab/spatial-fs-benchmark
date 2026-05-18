from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Merge multiple benchmark result CSV files.")
    parser.add_argument("--inputs", nargs="+", required=True, help="Input results.csv files.")
    parser.add_argument("--output", required=True, help="Output merged CSV path.")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    frames = [pd.read_csv(path) for path in args.inputs]
    merged = pd.concat(frames, ignore_index=True)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(output, index=False)
    print(f"Wrote merged results to {output} with {len(merged)} rows")


if __name__ == "__main__":
    main()
