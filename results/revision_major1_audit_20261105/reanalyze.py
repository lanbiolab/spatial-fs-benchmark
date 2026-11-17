"""Reproduce the R1-major-1 descriptive reanalysis from the frozen tables.

Uses the existing dataset/metric scaling bounds and seed-averaged task scores;
no embedding, feature set, or frozen range is refit. Run from any directory.
Requires numpy and pandas. Rank bootstrap samples datasets (not seeds).
"""
from pathlib import Path
import numpy as np
import pandas as pd

def rankdata(a, method="average"):
    return pd.Series(np.asarray(a)).rank(method=method).to_numpy()

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
OUT = ROOT / "results"
OUT.mkdir(exist_ok=True)
SEED = 20261105
DRAWS = 5000
LABELLED = ["DLPFC", "MouseBrainSerialSections", "STOmics0212", "STOmics0218", "STOmics0224"]
LOWER_CONFIDENCE = ["MouseBrainSerialSections", "STOmics0224"]
ALL = LABELLED + ["E8p5Embryo", "E9p5Embryo"]
SCENARIOS = {
    "7_datasets_original": ALL,
    "5_labelled": LABELLED,
    "3_label_provenance_sensitivity": [d for d in LABELLED if d not in LOWER_CONFIDENCE],
}
GROUPS = ["expression_driven", "spatially_informed"]

def write(df, name):
    df.to_csv(OUT / name, sep="\t", index=False, float_format="%.8f", na_rep="NA")

primary = pd.read_csv(DATA / "representative_task_scores.tsv", sep="\t", dtype={"n_features": str})
graphst = pd.read_csv(DATA / "graphst_setting_task_scores.tsv", sep="\t", dtype={"n_features": str})
comp = pd.concat([primary[primary.MethodGroup.isin(GROUPS)], graphst[graphst.MethodGroup.isin(GROUPS)]], ignore_index=True)
keys = ["dataset", "integration_method", "fs_method"]
assert not comp.duplicated(keys).any(), "duplicate competitive method/dataset/integrator keys"
assert comp.groupby("integration_method").size().to_dict() == {"cellcharter": 203, "graphst": 203, "scvi": 203}
assert comp.CoreOverallMean.notna().all()
assert (comp.groupby(["dataset", "integration_method"]).size() == 29).all()
for column in ["BatchMixingMean", "BiologicalConservationMean", "LabelIndependentClusteringMean", "LabelAgreementMean"]:
    assert comp.loc[comp.dataset.isin(LABELLED), column].notna().all(), ("annotated metric missing", column)

dataset_ranks = comp[keys + ["n_features", "MethodGroup", "BatchMixingMean", "BiologicalConservationMean", "IntegrationMean", "LabelIndependentClusteringMean", "LabelAgreementMean", "ClusteringMean", "CoreOverallMean"]].copy()
dataset_ranks["DatasetRank"] = dataset_ranks.groupby(["dataset", "integration_method"])["CoreOverallMean"].rank(ascending=False, method="average")
write(dataset_ranks.sort_values(["integration_method", "dataset", "DatasetRank"]), "dataset_scores_and_ranks.tsv")

rows = []
for scenario, datasets in SCENARIOS.items():
    for integrator, g in dataset_ranks[dataset_ranks.dataset.isin(datasets)].groupby("integration_method"):
        rank_matrix = g.pivot(index="dataset", columns="fs_method", values="DatasetRank").loc[datasets]
        assert rank_matrix.notna().all().all() and len(rank_matrix.columns) == 29
        scores = g.pivot(index="dataset", columns="fs_method", values="CoreOverallMean").loc[datasets]
        means = rank_matrix.mean(axis=0)
        global_rank = pd.Series(rankdata(means, method="average"), index=rank_matrix.columns)
        for method in rank_matrix.columns:
            rows.append(dict(scenario=scenario, integration_method=integrator, fs_method=method,
                             n_datasets=len(datasets), mean_dataset_rank=means[method],
                             global_rank=global_rank[method], mean_core_score=scores[method].mean()))
rank_table = pd.DataFrame(rows)
write(rank_table.sort_values(["scenario", "integration_method", "global_rank", "fs_method"]), "global_ranks_by_scenario.tsv")

prior = pd.read_csv(DATA / "prior_sensitivity_ranks.tsv", sep="\t")
prior = prior[(prior.Scenario == "labelled_datasets_only") & (prior.integration_method.isin(["scvi", "cellcharter"]))]
merged = rank_table[rank_table.scenario == "5_labelled"].merge(prior, on=["fs_method", "integration_method"], validate="one_to_one")
assert len(merged) == 58
assert np.allclose(merged.global_rank, merged.GlobalRank)
assert np.allclose(merged.mean_dataset_rank, merged.MeanDatasetRank)
assert np.allclose(merged.mean_core_score, merged.MeanScore, atol=1e-8)
prior_seven = pd.read_csv(DATA / "prior_seven_global_ranks.tsv", sep="\t")
original = rank_table[rank_table.scenario == "7_datasets_original"].merge(
    prior_seven, on=["fs_method", "integration_method"], validate="one_to_one")
assert len(original) == 58
assert np.allclose(original.global_rank, original.GlobalCoreOverallRank)
assert np.allclose(original.mean_dataset_rank, original.MeanCoreOverallRank)
assert np.allclose(original.mean_core_score, original.MeanCoreOverallScore, atol=1e-8)
prior_graphst = pd.read_csv(DATA / "prior_three_integrator_global_ranks.tsv", sep="\t")
prior_graphst = prior_graphst[prior_graphst.integration_method == "graphst"]
graphcheck = rank_table[(rank_table.scenario == "7_datasets_original") & (rank_table.integration_method == "graphst")].merge(
    prior_graphst, on=["fs_method", "integration_method"], validate="one_to_one")
assert len(graphcheck) == 29
assert np.allclose(graphcheck.global_rank, graphcheck.GlobalCoreOverallRank)
assert np.allclose(graphcheck.mean_dataset_rank, graphcheck.MeanCoreOverallRank)
assert np.allclose(graphcheck.mean_core_score, graphcheck.MeanCoreOverallScore, atol=1e-8)

comparison = rank_table.pivot(index=["integration_method", "fs_method"], columns="scenario", values="global_rank").reset_index()
comparison["five_minus_seven"] = comparison["5_labelled"] - comparison["7_datasets_original"]
comparison["three_minus_five"] = comparison["3_label_provenance_sensitivity"] - comparison["5_labelled"]
write(comparison.sort_values(["integration_method", "5_labelled"]), "rank_changes.tsv")

rng = np.random.default_rng(SEED)
bootrows = []
for scenario in ["7_datasets_original", "5_labelled"]:
    datasets = SCENARIOS[scenario]
    for integrator, g in dataset_ranks[dataset_ranks.dataset.isin(datasets)].groupby("integration_method"):
        pivot = g.pivot(index="dataset", columns="fs_method", values="DatasetRank").loc[datasets]
        a = pivot.to_numpy()
        sample_indices = rng.integers(0, len(datasets), size=(DRAWS, len(datasets)))
        means = a[sample_indices].mean(axis=1)
        ranks = np.stack([rankdata(row, method="average") for row in means])
        obs = rankdata(a.mean(axis=0), method="average")
        for j, method in enumerate(pivot.columns):
            bootrows.append(dict(scenario=scenario, integration_method=integrator, fs_method=method,
                n_datasets=len(datasets), replicates=DRAWS, rng_seed=SEED, observed_global_rank=obs[j],
                bootstrap_rank_median=np.median(ranks[:,j]), rank_lower_2_5=np.quantile(ranks[:,j],.025),
                rank_upper_97_5=np.quantile(ranks[:,j],.975), fraction_top_5=(ranks[:,j]<=5).mean()))
write(pd.DataFrame(bootrows).sort_values(["scenario", "integration_method", "observed_global_rank"]), "dataset_bootstrap_rank_intervals.tsv")

control = primary[primary.fs_method.isin(["random", "all_features"])].copy()
assert len(control) == 28 and set(control.loc[control.fs_method == "random", "n_features"]) == {"500"}
diagrows = []
for (ds, integrator), g in primary.groupby(["dataset", "integration_method"]):
    competitive = g[g.MethodGroup.isin(GROUPS)]
    controls = g[g.fs_method.isin(["random", "all_features"])]
    assert len(competitive) == 29 and len(controls) == 2
    for metric in ["BatchMixingMean", "BiologicalConservationMean", "LabelIndependentClusteringMean", "LabelAgreementMean", "IntegrationMean", "ClusteringMean", "CoreOverallMean"]:
        v = pd.to_numeric(competitive[metric], errors="coerce").dropna()
        for _, control_row in controls.iterrows():
            x = control_row[metric]
            if pd.isna(x) or v.empty:
                position = np.nan; percentile = np.nan; delta = np.nan
            else:
                position = 1 + (v > x).sum()
                percentile = ((v < x).sum() + 0.5*(v == x).sum()) / len(v)
                delta = x - v.median()
            diagrows.append(dict(dataset=ds, integration_method=integrator, control=control_row.fs_method,
                control_budget=control_row.n_features, component=metric, control_value=x,
                competitive_n=len(v), competitive_median=v.median() if len(v) else np.nan,
                difference_vs_competitive_median=delta, percentile_vs_competitive=percentile,
                position_among_29_competitors_plus_this_control=position))
diagnostics = pd.DataFrame(diagrows)
write(diagnostics.sort_values(["dataset", "integration_method", "control", "component"]), "control_component_diagnostics.tsv")

# Independently rebuild batch mixing from the two frozen seed-level metric rows.
seed_metrics = pd.read_csv(DATA / "batch_mixing_metric_seed.tsv", sep="\t", dtype={"n_features": str})
seed_keys = ["dataset", "integration_method", "fs_method", "n_features", "random_seed", "MethodGroup"]
assert len(seed_metrics) == 7 * 2 * 31 * 3 * 2
assert not seed_metrics.duplicated(seed_keys + ["metric_name"]).any()
seed_wide = seed_metrics.pivot(index=seed_keys, columns="metric_name", values="ScaledValue").reset_index()
assert seed_wide[["bASW", "iLISI"]].notna().all().all()
seed_wide["BatchMixing"] = seed_wide[["bASW", "iLISI"]].mean(axis=1)
recomputed = seed_wide.groupby(["dataset", "integration_method", "fs_method"], as_index=False).BatchMixing.mean()
verified = recomputed.merge(primary[keys + ["BatchMixingMean"]], on=keys, validate="one_to_one")
assert len(verified) == 7 * 2 * 31
assert np.allclose(verified.BatchMixing, verified.BatchMixingMean, atol=1e-8)
seed_results = []
for (dataset, integrator, seed), g in seed_wide.groupby(["dataset", "integration_method", "random_seed"]):
    competitors = g[g.MethodGroup.isin(GROUPS)].BatchMixing
    assert len(competitors) == 29
    for _, control_row in g[g.fs_method.isin(["random", "all_features"])].iterrows():
        x = control_row.BatchMixing
        seed_results.append(dict(dataset=dataset, integration_method=integrator, seed=seed,
            control=control_row.fs_method, control_budget=control_row.n_features,
            batch_mixing=x, competitive_median=competitors.median(),
            position_among_29_plus_control=1 + (competitors > x).sum()))
write(pd.DataFrame(seed_results).sort_values(["dataset", "integration_method", "control", "seed"]),
      "control_batch_mixing_by_seed.tsv")

summaries=[]
for integrator, g in comparison.groupby("integration_method"):
    rho = g["7_datasets_original"].rank().corr(g["5_labelled"].rank())
    summaries.append(dict(integration_method=integrator, rank_spearman_5_vs_7=rho,
                          n_changed=(g.five_minus_seven != 0).sum(),
                          median_absolute_rank_change=g.five_minus_seven.abs().median(),
                          max_absolute_rank_change=g.five_minus_seven.abs().max()))
write(pd.DataFrame(summaries), "rank_comparison_summary.tsv")
print("Validated prior five-dataset rankings (58/58) and seven-dataset rankings (87/87, including GraphST).")
print("Recomputed 434 seed-averaged batch-mixing scores from 2,604 seed/metric rows; all match.")
print("Datasets: 7, 5, and exploratory 3; integrators: scVI, CellCharter, GraphST.")
print("Dataset bootstrap: 5,000 replicates per integrator/scenario; seed",SEED)
print("All outputs:", OUT)
