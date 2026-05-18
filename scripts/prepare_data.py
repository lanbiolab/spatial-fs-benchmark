from __future__ import annotations

import argparse
import shutil
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Prepare benchmark datasets in the local data directory.")
    parser.add_argument("--dataset", required=True, choices=["baristaseq"], help="Dataset identifier.")
    parser.add_argument("--source", required=True, help="Source directory to copy from.")
    parser.add_argument("--destination-root", default="data/raw", help="Destination root directory.")
    return parser


def copy_dataset(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        shutil.rmtree(destination)
    shutil.copytree(source, destination)


def main() -> None:
    args = build_parser().parse_args()
    source = Path(args.source).resolve()
    if not source.exists():
        raise FileNotFoundError(f"Source dataset directory does not exist: {source}")
    destination = Path(args.destination_root).resolve() / args.dataset
    copy_dataset(source, destination)
    print(f"Copied {args.dataset} data into {destination}")


if __name__ == "__main__":
    main()
