from __future__ import annotations

import argparse
from pathlib import Path


SUBSETS = {
    "stomics0212immune": Path("stomics_0212_immune_subset/stomics0212immune"),
    "stomics0212epithelial": Path("stomics_0212_epithelial_subset/stomics0212epithelial"),
}
EXPECTED_EMBEDDINGS = 2 * 34 * 3
EXPECTED_TASK_RECORDS = EXPECTED_EMBEDDINGS * 2


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--results-root",
        default="results/biological_subsets_rebuild_v1",
        type=Path,
    )
    args = parser.parse_args()

    embeddings = list(args.results_root.glob("**/embedding.npz"))
    records = list(args.results_root.glob("**/*_records.json"))
    pct = 100 * len(embeddings) / EXPECTED_EMBEDDINGS
    width = 30
    complete = min(width, round(width * pct / 100))
    bar = "#" * complete + "-" * (width - complete)
    print(f"Embeddings [{bar}] {len(embeddings)}/{EXPECTED_EMBEDDINGS} ({pct:.1f}%)")
    print(f"Task records: {len(records)}/{EXPECTED_TASK_RECORDS}")
    for subset, relative_path in SUBSETS.items():
        subset_dir = args.results_root / relative_path
        subset_embeddings = len(list(subset_dir.glob("**/embedding.npz")))
        subset_records = len(list(subset_dir.glob("**/*_records.json")))
        print(f"{subset}: {subset_embeddings}/102 embeddings, {subset_records}/204 records")


if __name__ == "__main__":
    main()
