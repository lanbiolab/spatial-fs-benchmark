from __future__ import annotations

from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]
CONFIG_DIR = ROOT / "configs" / "benchmark"
OUTPUT_ROOT = "results/spatial_main_native_seed0_fix2"

SOURCE_CONFIGS = [
    "dlpfc_spatial_main.yaml",
    "mouse_brain_serial_sections_spatial_main.yaml",
    "sagittal_region_spatial_main.yaml",
    "stomics_visium_5samples_spatial_main.yaml",
    "e8p5_embryo_spatial_main.yaml",
    "e9p5_embryo_spatial_main.yaml",
]


def main() -> None:
    for source_name in SOURCE_CONFIGS:
        source_path = CONFIG_DIR / source_name
        payload = yaml.safe_load(source_path.read_text(encoding="utf-8"))
        payload["name"] = str(payload["name"]).replace("_spatial_main", "_spatial_main_native")
        output_dir = str(payload["output_dir"])
        dataset_slug = output_dir.split("/")[-1]
        payload["output_dir"] = f"{OUTPUT_ROOT}/{dataset_slug}"
        payload["seeds"] = [0]
        destination_name = source_name.replace("_spatial_main", "_spatial_main_native")
        destination_path = CONFIG_DIR / destination_name
        destination_path.write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")
        print(destination_path.relative_to(ROOT))


if __name__ == "__main__":
    main()
