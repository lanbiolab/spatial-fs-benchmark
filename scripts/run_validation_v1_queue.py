#!/usr/bin/env python3

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

import yaml


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser()
    result.add_argument("config")
    return result


def _run(command: list[str], threads: int) -> None:
    env = os.environ.copy()
    for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        env[variable] = str(threads)
    env["R_FUTURE_PLAN"] = "sequential"
    print("RUN " + " ".join(command), flush=True)
    subprocess.run(command, check=True, env=env)


def main() -> None:
    args = parser().parse_args()
    config_path = Path(args.config)
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    python = str(Path(sys.executable))
    methods = ",".join(config["methods"])
    threads = int(config.get("threads", 4))

    if config["kind"] == "semi_synthetic":
        for prevalence in config["prevalence"]:
            for seed in config["seeds"]:
                _run(
                    [
                        python,
                        "scripts/run_semi_synthetic_validation.py",
                        "--prevalence",
                        str(prevalence),
                        "--seed",
                        str(seed),
                        "--n-genes",
                        str(config["n_genes"]),
                        "--grid-size",
                        str(config["grid_size"]),
                        "--n-slices",
                        str(config["n_slices"]),
                        "--methods",
                        methods,
                        "--output-root",
                        str(config["output_root"]),
                    ],
                    threads,
                )
    elif config["kind"] == "heldout_slice":
        for fold in config["folds"]:
            _run(
                [
                    python,
                    "scripts/run_heldout_slice_validation.py",
                    "--dataset-config",
                    str(fold["dataset_config"]),
                    "--heldout-slice",
                    str(fold["heldout_slice"]),
                    "--n-features",
                    str(config["n_features"]),
                    "--seed",
                    str(config["seed"]),
                    "--methods",
                    methods,
                    "--output-root",
                    str(config["output_root"]),
                ],
                threads,
            )
    else:
        raise ValueError(f"Unsupported validation kind: {config['kind']}")


if __name__ == "__main__":
    main()
