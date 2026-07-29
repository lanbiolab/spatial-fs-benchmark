from __future__ import annotations

from dataclasses import asdict, dataclass

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.metrics import average_precision_score, roc_auc_score

from spatial_fs_benchmark.data.spatial_object import SpatialDataset
from spatial_fs_benchmark.feature_selection.base import FeatureSelectionResult


PATTERNS = ("domain", "gradient", "focal", "periodic")
EFFECTS = ("weak", "moderate", "strong")
EFFECT_STRENGTH = {"weak": 0.35, "moderate": 0.70, "strong": 1.10}


@dataclass(frozen=True, slots=True)
class SemiSyntheticSpec:
    prevalence: float
    seed: int
    n_genes: int = 3000
    grid_size: int = 15
    n_slices: int = 3
    target_sum: float = 1e4
    dispersion: float = 12.0

    @property
    def n_spots_per_slice(self) -> int:
        return self.grid_size**2

    @property
    def n_obs(self) -> int:
        return self.n_spots_per_slice * self.n_slices

    @property
    def n_truth(self) -> int:
        return max(len(PATTERNS) * len(EFFECTS), round(self.n_genes * self.prevalence))


def _standardize(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    return (values - values.mean()) / (values.std() + 1e-12)


def _pattern_profiles(coords: np.ndarray, slice_number: int) -> dict[str, np.ndarray]:
    x = coords[:, 0]
    y = coords[:, 1]
    focal_x = 0.32 + 0.10 * (slice_number % 3)
    focal_y = 0.68 - 0.08 * (slice_number % 3)
    phase = slice_number * np.pi / 7
    return {
        "domain": _standardize((x + 0.15 * y) > 0.58),
        "gradient": _standardize(0.75 * x + 0.25 * y),
        "focal": _standardize(np.exp(-((x - focal_x) ** 2 + (y - focal_y) ** 2) / 0.035)),
        "periodic": _standardize(np.sin(4 * np.pi * x + phase) + 0.35 * np.cos(2 * np.pi * y)),
    }


def _allocate_truth_strata(n_truth: int) -> list[tuple[str, str]]:
    strata = [(pattern, effect) for pattern in PATTERNS for effect in EFFECTS]
    return [strata[index % len(strata)] for index in range(n_truth)]


def build_semi_synthetic_dataset(spec: SemiSyntheticSpec) -> SpatialDataset:
    if not 0 < spec.prevalence < 1:
        raise ValueError("prevalence must be between zero and one")
    if spec.n_truth >= spec.n_genes:
        raise ValueError("prevalence leaves no null genes")

    rng = np.random.default_rng(spec.seed)
    axis = np.linspace(0, 1, spec.grid_size)
    xx, yy = np.meshgrid(axis, axis, indexing="xy")
    base_coords = np.column_stack([xx.ravel(), yy.ravel()])
    coords = np.vstack([base_coords for _ in range(spec.n_slices)]).astype(np.float32)
    slice_ids = np.repeat([f"slice{index + 1}" for index in range(spec.n_slices)], spec.n_spots_per_slice)
    labels = np.where(
        coords[:, 0] < 0.5,
        np.where(coords[:, 1] < 0.5, "lower-left", "upper-left"),
        np.where(coords[:, 1] < 0.5, "lower-right", "upper-right"),
    )

    gene_names = np.asarray([f"G{index:05d}" for index in range(spec.n_genes)])
    truth_indices = rng.choice(spec.n_genes, size=spec.n_truth, replace=False)
    truth_strata = _allocate_truth_strata(spec.n_truth)
    rng.shuffle(truth_strata)
    truth_pattern = np.full(spec.n_genes, "null", dtype=object)
    truth_effect = np.full(spec.n_genes, "null", dtype=object)
    truth_sign = np.ones(spec.n_genes, dtype=float)
    for gene_index, (pattern, effect) in zip(truth_indices, truth_strata, strict=True):
        truth_pattern[gene_index] = pattern
        truth_effect[gene_index] = effect
        truth_sign[gene_index] = rng.choice((-1.0, 1.0))

    base_means = np.exp(rng.normal(np.log(1.2), 0.75, size=spec.n_genes))
    library_factors = np.exp(rng.normal(0, 0.28, size=spec.n_obs))
    means = library_factors[:, None] * base_means[None, :]

    for slice_number in range(spec.n_slices):
        start = slice_number * spec.n_spots_per_slice
        end = start + spec.n_spots_per_slice
        profiles = _pattern_profiles(base_coords, slice_number)
        for gene_index in truth_indices:
            pattern = str(truth_pattern[gene_index])
            effect = str(truth_effect[gene_index])
            log_multiplier = (
                EFFECT_STRENGTH[effect]
                * truth_sign[gene_index]
                * profiles[pattern]
            )
            means[start:end, gene_index] *= np.exp(np.clip(log_multiplier, -2.2, 2.2))

    theta = float(spec.dispersion)
    probabilities = theta / (theta + means)
    counts = rng.negative_binomial(theta, probabilities).astype(np.int32)
    counts_sparse = sparse.csr_matrix(counts)
    library_sizes = np.asarray(counts_sparse.sum(axis=1)).ravel()
    scale = spec.target_sum / np.maximum(library_sizes, 1)
    normalized = counts_sparse.multiply(scale[:, None]).tocsr().astype(np.float32)
    normalized.data = np.log1p(normalized.data)

    obs = pd.DataFrame(
        {
            "slice_id": pd.Categorical(slice_ids),
            "domain_label": pd.Categorical(labels),
        },
        index=[f"spot{index:05d}" for index in range(spec.n_obs)],
    )
    var = pd.DataFrame(
        {
            "is_svg": np.isin(np.arange(spec.n_genes), truth_indices),
            "pattern": truth_pattern,
            "effect": truth_effect,
        },
        index=gene_names,
    )
    adata = ad.AnnData(X=normalized, obs=obs, var=var)
    adata.layers["counts"] = counts_sparse
    adata.obsm["spatial"] = coords
    adata.uns["semi_synthetic_spec"] = asdict(spec)
    adata.uns["spatial_fs_benchmark"] = {
        "counts_source": "generated_negative_binomial",
        "counts_are_nonnegative_integers": True,
    }
    scenario = f"semi_synthetic_p{spec.prevalence:.2f}_seed{spec.seed}"
    return SpatialDataset(
        name=scenario,
        adata=adata,
        slice_key="slice_id",
        coord_key="spatial",
        label_key="domain_label",
        platform="semi-synthetic grid",
        species="synthetic",
        source_path=scenario,
    )


def ranking_scores(
    gene_names: np.ndarray,
    selection: FeatureSelectionResult,
) -> np.ndarray:
    index = {str(gene): position for position, gene in enumerate(gene_names)}
    scores = np.zeros(len(gene_names), dtype=float)
    n_ranked = len(selection.feature_names)
    for rank, gene in enumerate(selection.feature_names):
        position = index.get(str(gene))
        if position is not None:
            scores[position] = n_ranked - rank
    return scores


def score_feature_ranking(
    dataset: SpatialDataset,
    selection: FeatureSelectionResult,
    cutoffs: tuple[int, ...] = (100, 200, 500),
) -> pd.DataFrame:
    truth = dataset.adata.var["is_svg"].to_numpy(dtype=bool)
    if selection.method_name == "all_features":
        prevalence = float(np.mean(truth))
        return pd.DataFrame(
            [
                {
                    "summary": "overall",
                    "stratum": "all_svg",
                    "cutoff": pd.NA,
                    "metric": "AUROC",
                    "value": 0.5,
                },
                {
                    "summary": "overall",
                    "stratum": "all_svg",
                    "cutoff": pd.NA,
                    "metric": "AUPRC",
                    "value": prevalence,
                },
                {
                    "summary": "cutoff",
                    "stratum": "all_svg",
                    "cutoff": dataset.n_vars,
                    "metric": "Precision",
                    "value": prevalence,
                },
                {
                    "summary": "cutoff",
                    "stratum": "all_svg",
                    "cutoff": dataset.n_vars,
                    "metric": "Recall",
                    "value": 1.0,
                },
            ]
        )
    scores = ranking_scores(dataset.adata.var_names.to_numpy(), selection)
    order = np.argsort(scores, kind="stable")[::-1]
    rows: list[dict[str, object]] = [
        {
            "summary": "overall",
            "stratum": "all_svg",
            "cutoff": pd.NA,
            "metric": "AUROC",
            "value": float(roc_auc_score(truth, scores)),
        },
        {
            "summary": "overall",
            "stratum": "all_svg",
            "cutoff": pd.NA,
            "metric": "AUPRC",
            "value": float(average_precision_score(truth, scores)),
        },
    ]
    for cutoff in cutoffs:
        effective = min(int(cutoff), len(order))
        selected_mask = np.zeros(len(order), dtype=bool)
        selected_mask[order[:effective]] = True
        true_positive = int(np.sum(selected_mask & truth))
        rows.extend(
            [
                {
                    "summary": "cutoff",
                    "stratum": "all_svg",
                    "cutoff": effective,
                    "metric": "Precision",
                    "value": true_positive / max(effective, 1),
                },
                {
                    "summary": "cutoff",
                    "stratum": "all_svg",
                    "cutoff": effective,
                    "metric": "Recall",
                    "value": true_positive / max(int(truth.sum()), 1),
                },
            ]
        )

    pattern = dataset.adata.var["pattern"].astype(str).to_numpy()
    effect = dataset.adata.var["effect"].astype(str).to_numpy()
    validation_cutoff = min(500, len(order))
    top_mask = np.zeros(len(order), dtype=bool)
    top_mask[order[:validation_cutoff]] = True
    for summary, values, levels in (
        ("pattern", pattern, PATTERNS),
        ("effect", effect, EFFECTS),
    ):
        for level in levels:
            stratum_truth = truth & (values == level)
            rows.append(
                {
                    "summary": summary,
                    "stratum": level,
                    "cutoff": validation_cutoff,
                    "metric": "Recall",
                    "value": float(np.sum(top_mask & stratum_truth) / max(int(stratum_truth.sum()), 1)),
                }
            )
    return pd.DataFrame(rows)
