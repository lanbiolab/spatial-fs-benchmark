from __future__ import annotations

import argparse
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]


def load(path: Path) -> dict:
    return yaml.safe_load(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate a primary-integrator rerun for a selector whose feature hash changed."
    )
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--method", required=True)
    args = parser.parse_args()

    graphst_path = ROOT / f"configs/rebuild_v1/downstream/graphst_canonical/{args.dataset}.yaml"
    primary_path = ROOT / f"configs/rebuild_v1/non_spatial_downstream/canonical/{args.dataset}.yaml"
    graphst = load(graphst_path)
    primary = load(primary_path)
    selected = [item for item in graphst["feature_selection_methods"] if item["name"] == args.method]
    if len(selected) != 1:
        raise ValueError(f"Expected one {args.method} setting in {graphst_path}, found {len(selected)}")

    payload = {
        "name": f"{args.dataset}_{args.method}_feature_hash_primary_integrator_repair",
        "output_dir": primary["output_dir"],
        "datasets": primary["datasets"],
        "feature_selection_methods": selected,
        "integration_methods": primary["integration_methods"],
        "tasks": primary["tasks"],
        "n_features": primary["n_features"],
        "seeds": [0, 1, 2],
        "save_embeddings": True,
    }
    output = (
        ROOT
        / "configs/rebuild_v1/downstream/selector_hash_repairs"
        / f"{args.dataset}_{args.method}.yaml"
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")
    print(output)


if __name__ == "__main__":
    main()
