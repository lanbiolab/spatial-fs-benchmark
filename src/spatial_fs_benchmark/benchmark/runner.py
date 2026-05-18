from __future__ import annotations

import json
from pathlib import Path

from spatial_fs_benchmark.benchmark.experiment import (
    cache_signature,
    metric_records_from_task_output,
    read_json,
    read_feature_selection,
    read_integration_result,
    read_metric_records,
    run_timed_block,
    write_json,
    write_feature_selection,
    write_integration_result,
    write_metric_records,
)
from spatial_fs_benchmark.benchmark.result_schema import MetricRecord, records_to_frame
from spatial_fs_benchmark.config import BenchmarkConfig
from spatial_fs_benchmark.data.io import load_dataset
from spatial_fs_benchmark.data.preprocess import preprocess_dataset
from spatial_fs_benchmark.feature_selection import build_feature_selector
from spatial_fs_benchmark.integration import build_integrator
from spatial_fs_benchmark.tasks import build_task
from spatial_fs_benchmark.utils.logging import configure_logging


class BenchmarkRunner:
    def __init__(self, config: BenchmarkConfig) -> None:
        self.config = config
        self.output_dir = Path(config.output_dir)
        self.logger = configure_logging(self.output_dir / "logs" / "benchmark.log")

    @staticmethod
    def _dataset_signature(dataset_config, dataset) -> dict:
        preprocess = json.loads(json.dumps(dataset_config.preprocess, sort_keys=True))
        return {
            "dataset_name": dataset.name,
            "dataset_path": dataset_config.path,
            "platform": dataset.platform,
            "n_obs": int(dataset.n_obs),
            "n_vars": int(dataset.n_vars),
            "slice_key": dataset.slice_key,
            "label_key": dataset.label_key,
            "coord_key": dataset.coord_key,
            "preprocess": preprocess,
        }

    def run(self) -> None:
        self.output_dir.mkdir(parents=True, exist_ok=True)
        all_records: list[MetricRecord] = []
        for dataset_config in self.config.datasets:
            dataset = preprocess_dataset(load_dataset(dataset_config), dataset_config.preprocess)
            dataset_sig = self._dataset_signature(dataset_config, dataset)
            self.logger.info("Loaded dataset %s with %s spots and %s genes", dataset.name, dataset.n_obs, dataset.n_vars)
            dataset_dir = self.output_dir / dataset.name.lower()
            dataset_dir.mkdir(parents=True, exist_ok=True)
            for seed in self.config.seeds:
                for selector_config in self.config.feature_selection_methods:
                    try:
                        selector = build_feature_selector(selector_config.name, **selector_config.params)
                    except Exception as exc:
                        self.logger.exception(
                            "Failed to build selector for dataset=%s selector=%s seed=%s: %s",
                            dataset.name,
                            selector_config.name,
                            seed,
                            exc,
                        )
                        continue
                    for n_features in self.config.n_features:
                        selection_dir = dataset_dir / "feature_selection" / selector_config.name / f"n{n_features}" / f"seed{seed}"
                        selection_path = selection_dir / "selected_features.json"
                        selection_meta_path = selection_dir / "selected_features.meta.json"
                        selection_payload = {
                            "stage": "feature_selection",
                            "dataset": dataset_sig,
                            "selector_name": selector_config.name,
                            "selector_params": selector_config.params,
                            "n_features": int(n_features),
                            "seed": int(seed),
                        }
                        selection_sig = cache_signature(selection_payload)
                        try:
                            selection_cache_valid = False
                            if selection_path.exists() and selection_meta_path.exists():
                                selection_meta = read_json(selection_meta_path)
                                selection_cache_valid = selection_meta.get("cache_signature") == selection_sig
                            if selection_cache_valid:
                                selection = read_feature_selection(selection_path)
                                selection_runtime = 0.0
                            else:
                                selection, selection_runtime = run_timed_block(
                                    lambda selector=selector, n_features=n_features, seed=seed: selector.select(
                                        dataset,
                                        n_features=n_features,
                                        random_seed=seed,
                                    )
                                )
                                effective_n_features = int(selection.metadata.get("effective_n_features", len(selection.feature_names)))
                                if effective_n_features < n_features:
                                    self.logger.warning(
                                        "Requested %s features for %s on %s, but only %s genes are available after preprocessing",
                                        n_features,
                                        selector_config.name,
                                        dataset.name,
                                        effective_n_features,
                                    )
                                selection.metadata["cache_signature"] = selection_sig
                                write_feature_selection(selection_path, selection)
                                write_json(
                                    selection_meta_path,
                                    {
                                        **selection_payload,
                                        "cache_signature": selection_sig,
                                        "effective_n_features": effective_n_features,
                                    },
                                )
                        except Exception as exc:
                            self.logger.exception(
                                "Feature selection failed for dataset=%s selector=%s n_features=%s seed=%s: %s",
                                dataset.name,
                                selector_config.name,
                                n_features,
                                seed,
                                exc,
                            )
                            continue
                        for integrator_config in self.config.integration_methods:
                            try:
                                integrator = build_integrator(integrator_config.name, **integrator_config.params)
                            except Exception as exc:
                                self.logger.exception(
                                    "Failed to build integrator for dataset=%s selector=%s n_features=%s integrator=%s seed=%s: %s",
                                    dataset.name,
                                    selector_config.name,
                                    n_features,
                                    integrator_config.name,
                                    seed,
                                    exc,
                                )
                                continue
                            combination_dir = (
                                dataset_dir
                                / integrator_config.name
                                / selector_config.name
                                / f"n{n_features}"
                                / f"seed{seed}"
                            )
                            embedding_path = combination_dir / "embedding.npz"
                            integration_meta_path = combination_dir / "embedding.meta.json"
                            integration_payload = {
                                "stage": "integration",
                                "dataset": dataset_sig,
                                "selector_name": selector_config.name,
                                "selector_params": selector_config.params,
                                "n_features": int(n_features),
                                "seed": int(seed),
                                "selection_feature_hash": selection.metadata.get("feature_name_hash", ""),
                                "integrator_name": integrator_config.name,
                                "integrator_params": integrator_config.params,
                                "integrator_implementation_version": getattr(
                                    integrator,
                                    "implementation_version",
                                    None,
                                ),
                            }
                            integration_sig = cache_signature(integration_payload)
                            try:
                                integration_cache_valid = False
                                if self.config.save_embeddings and embedding_path.exists() and integration_meta_path.exists():
                                    integration_meta = read_json(integration_meta_path)
                                    integration_cache_valid = integration_meta.get("cache_signature") == integration_sig
                                if integration_cache_valid:
                                    integration_result = read_integration_result(embedding_path, integrator_config.name)
                                    integration_runtime = 0.0
                                else:
                                    integration_result, integration_runtime = run_timed_block(
                                        lambda integrator=integrator, selection=selection, seed=seed: integrator.fit_transform(
                                            dataset,
                                            selection,
                                            random_seed=seed,
                                        )
                                    )
                                    integration_result.metadata["cache_signature"] = integration_sig
                                    if self.config.save_embeddings:
                                        write_integration_result(embedding_path, dataset, integration_result)
                                        write_json(
                                            integration_meta_path,
                                            {
                                                **integration_payload,
                                                "cache_signature": integration_sig,
                                            },
                                        )
                            except Exception as exc:
                                self.logger.exception(
                                    "Integration failed for dataset=%s selector=%s n_features=%s integrator=%s seed=%s: %s",
                                    dataset.name,
                                    selector_config.name,
                                    n_features,
                                    integrator_config.name,
                                    seed,
                                    exc,
                                )
                                continue
                            for task_config in self.config.tasks:
                                task_records_path = combination_dir / f"{task_config.name}_records.json"
                                task_meta_path = combination_dir / f"{task_config.name}_records.meta.json"
                                task_payload = {
                                    "stage": "task",
                                    "dataset": dataset_sig,
                                    "selector_name": selector_config.name,
                                    "selector_params": selector_config.params,
                                    "n_features": int(n_features),
                                    "seed": int(seed),
                                    "integrator_name": integrator_config.name,
                                    "integrator_params": integrator_config.params,
                                    "task_name": task_config.name,
                                    "task_params": task_config.params,
                                    "selection_feature_hash": selection.metadata.get("feature_name_hash", ""),
                                    "integration_cache_signature": integration_sig,
                                }
                                task_sig = cache_signature(task_payload)
                                task_cache_valid = False
                                if task_records_path.exists() and task_meta_path.exists():
                                    task_meta = read_json(task_meta_path)
                                    task_cache_valid = task_meta.get("cache_signature") == task_sig
                                if task_cache_valid:
                                    all_records.extend(read_metric_records(task_records_path))
                                    continue
                                try:
                                    task = build_task(task_config.name, **task_config.params)
                                except Exception as exc:
                                    self.logger.exception(
                                        "Failed to build task for dataset=%s selector=%s n_features=%s integrator=%s task=%s seed=%s: %s",
                                        dataset.name,
                                        selector_config.name,
                                        n_features,
                                        integrator_config.name,
                                        task_config.name,
                                        seed,
                                        exc,
                                    )
                                    continue
                                try:
                                    task_output, task_runtime = run_timed_block(
                                        lambda task=task, selection=selection, integration_result=integration_result, seed=seed: task.evaluate(
                                            dataset,
                                            selection,
                                            integration_result,
                                            random_seed=seed,
                                        )
                                    )
                                    task_records = metric_records_from_task_output(
                                        dataset=dataset,
                                        fs_method=selector_config.name,
                                        n_features=n_features,
                                        effective_n_features=int(
                                            selection.metadata.get("effective_n_features", len(selection.feature_names))
                                        ),
                                        integration_method=integrator_config.name,
                                        task_name=task_config.name,
                                        task_output=task_output,
                                        random_seed=seed,
                                        runtime=selection_runtime + integration_runtime + task_runtime,
                                        notes=f"cache_signature={task_sig}",
                                    )
                                    write_metric_records(task_records_path, task_records)
                                    write_json(
                                        task_meta_path,
                                        {
                                            **task_payload,
                                            "cache_signature": task_sig,
                                        },
                                    )
                                    all_records.extend(task_records)
                                except Exception as exc:
                                    self.logger.exception(
                                        "Task failed for dataset=%s selector=%s n_features=%s integrator=%s task=%s seed=%s: %s",
                                        dataset.name,
                                        selector_config.name,
                                        n_features,
                                        integrator_config.name,
                                        task_config.name,
                                        seed,
                                        exc,
                                    )
                                    continue
        results = records_to_frame(all_records)
        results_path = self.output_dir / "results.csv"
        results.to_csv(results_path, index=False)
        self.logger.info("Wrote %s metric rows to %s", len(results), results_path)
