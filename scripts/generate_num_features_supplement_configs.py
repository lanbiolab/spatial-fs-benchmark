#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path
import copy
import yaml


ROOT = Path(__file__).resolve().parents[1]
CONFIG_DIR = ROOT / "configs" / "benchmark"

MAIN_CONFIGS = [
    "dlpfc_spatial_main_native.yaml",
    "mouse_brain_serial_sections_spatial_main_native.yaml",
    "stomics_visium_5samples_spatial_main_native.yaml",
    "e8p5_embryo_spatial_main_native.yaml",
    "e9p5_embryo_spatial_main_native.yaml",
]

WILCOXON_CONFIGS = [
    "stomics_0218_wilcoxon_spatial_main_native.yaml",
    "stomics_0224_wilcoxon_spatial_main_native.yaml",
    "stomics_0212_wilcoxon_spatial_main_native.yaml",
]

EXTRA_N_FEATURES = [100, 200, 5000, 10000]
ALL_FEATURES_DUMMY_N = [1]


def load_yaml(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def dump_yaml(path: Path, payload: dict) -> None:
    with path.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(payload, handle, sort_keys=False, allow_unicode=False)


def build_extra_config(src_name: str) -> None:
    src = CONFIG_DIR / src_name
    cfg = load_yaml(src)
    out = copy.deepcopy(cfg)
    out["name"] = f'{cfg["name"]}_num_features_extra'
    out["n_features"] = EXTRA_N_FEATURES
    dst = CONFIG_DIR / f'{src.stem}_num_features_extra.yaml'
    dump_yaml(dst, out)


def build_all_config(src_name: str) -> None:
    src = CONFIG_DIR / src_name
    cfg = load_yaml(src)
    out = copy.deepcopy(cfg)
    out["name"] = f'{cfg["name"]}_all_features'
    out["feature_selection_methods"] = [{"name": "all_features"}]
    out["n_features"] = ALL_FEATURES_DUMMY_N
    dst = CONFIG_DIR / f'{src.stem}_all_features.yaml'
    dump_yaml(dst, out)


def main() -> None:
    for name in MAIN_CONFIGS + WILCOXON_CONFIGS:
        build_extra_config(name)
        build_all_config(name)


if __name__ == "__main__":
    main()
