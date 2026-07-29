from __future__ import annotations

import argparse
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    parser = argparse.ArgumentParser(description="Create a disjoint-seed GraphST accelerator config.")
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--seeds", required=True, help="Comma-separated integer seeds")
    args = parser.parse_args()

    source = ROOT / f"configs/rebuild_v1/downstream/graphst_canonical/{args.dataset}.yaml"
    payload = yaml.safe_load(source.read_text(encoding="utf-8"))
    seeds = [int(value) for value in args.seeds.split(",")]
    payload["name"] = f"{args.dataset}_graphst_seed_{'_'.join(map(str, seeds))}_accelerator"
    payload["seeds"] = seeds
    output = (
        ROOT
        / "configs/rebuild_v1/downstream/graphst_seed_accelerators"
        / f"{args.dataset}_seed_{'_'.join(map(str, seeds))}.yaml"
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")
    print(output)


if __name__ == "__main__":
    main()
