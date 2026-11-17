"""R1 minor 1/3: provenance-aware resampling and feature-budget uncertainty."""
import hashlib
import importlib.util
import itertools
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import rankdata
import yaml

ROOT = Path('/home/xd/y/Benchmark')
OUT = ROOT / 'results/revision_minor1_3_20261106'
FIVE = ROOT / 'revision_r1_major1_five_20261105/frozen_five'
OLD_ANALYSIS = FIVE.parent / 'analysis'
ANALYSIS = OUT / 'analysis'
ANALYSIS.mkdir(parents=True, exist_ok=True)
GROUPS = {'expression_driven', 'spatially_informed'}
BUDGETS = [100, 200, 500, 1000, 2000, 5000, 10000]
DATASETS = ['DLPFC', 'MouseBrainSerialSections', 'STOmics0212', 'STOmics0218', 'STOmics0224']
SEED = 20261106
NBOOT = 5000
input_paths = [FIVE / n for n in ['setting_task_scores.tsv', 'combined_integrator_scores.tsv', 'setting_metric_seed_summary.tsv', 'frozen_metric_ranges.tsv']]
input_paths += [ROOT / 'scripts/generate_num_features_benchmark_figure.R', ROOT / 'external/atlas-feature-selection-benchmark/analysis/R/plotting.R',
                OLD_ANALYSIS / 'five_dataset_bootstrap_ranks.tsv', OLD_ANALYSIS / 'five_dataset_budget_summary.tsv']
input_hashes = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in input_paths}
SOURCES = [
    dict(dataset='DLPFC', config='dlpfc', study='Maynard2021_spatialLIBD', species='Human', tissue='DLPFC',
         source='Maynard et al., Nature Neuroscience (2021); spatialLIBD/HumanPilot',
         identifier='10.1038/s41593-020-00787-0', url='https://doi.org/10.1038/s41593-020-00787-0', role='primary',
         within_study_dependence='4 sections from donors Br5292 and Br8100; one study-level benchmark score'),
    dict(dataset='MouseBrainSerialSections', config='mouse_brain_serial_sections', study='10x2019_sagittal_anterior_serial', species='Mouse', tissue='Brain',
         source='10x Genomics sagittal-anterior serial-section collection', identifier='10x serial sections 1/2, Space Ranger 1.0.0',
         url='https://www.10xgenomics.com/datasets/mouse-brain-serial-section-1-sagittal-anterior-1-standard-1-0-0', role='primary',
         within_study_dependence='2 serial sections combined as one collection unit, not 2 bootstrap units'),
    dict(dataset='STOmics0212', config='stomics_0212_wilcoxon', study='GSE197023', species='Human', tissue='Skin',
         source='Mitamura et al.; original atopic-dermatitis skin study', identifier='GSE197023 / STDS0000212',
         url='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE197023', role='primary',
         within_study_dependence='10 selected biopsies, including paired lesion/non-lesion donors; one study unit'),
    dict(dataset='STOmics0218', config='stomics_0218_wilcoxon', study='GSE217058', species='Mouse', tissue='Brain',
         source='Ballester Roig et al.; rhynchophylline/sleep spatial transcriptomics study', identifier='GSE217058 / STDS0000218',
         url='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217058', role='primary',
         within_study_dependence='7 selected sex/time/treatment samples; one original experimental study'),
    dict(dataset='STOmics0224', config='stomics_0224_wilcoxon', study='GSE184510', species='Human', tissue='Orbitofrontal/temporal cortex',
         source='Mohammadi et al.; consecutive-section spatial transcriptomics study', identifier='GSE184510 / STDS0000224',
         url='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184510', role='primary',
         within_study_dependence='8 sections, donor IDs 2339/2248 and 2 cortical regions; one study unit'),
    dict(dataset='E8p5Embryo', config='e8p5_embryo', study='GSE197353', species='Mouse', tissue='E8.5 embryo',
         source='Whole-mouse-embryo Slide-seq study', identifier='GSE197353',
         url='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE197353', role='separate_label_free',
         within_study_dependence='Shares original study with E9.5; excluded from primary bootstrap'),
    dict(dataset='E9p5Embryo', config='e9p5_embryo', study='GSE197353', species='Mouse', tissue='E9.5 embryo',
         source='Whole-mouse-embryo Slide-seq study', identifier='GSE197353',
         url='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE197353', role='separate_label_free',
         within_study_dependence='Shares original study with E8.5; excluded from primary bootstrap'),
]
sources = pd.DataFrame(SOURCES)
sources.to_csv(ANALYSIS / 'dataset_original_study_inventory.tsv', sep='\t', index=False)
primary = sources[sources.role.eq('primary')]
assert len(primary) == primary.study.nunique() == 5
sample_rows = []
for r in primary.itertuples():
    cfg = yaml.safe_load((ROOT / 'configs/datasets' / f'{r.config}.yaml').read_text())
    data = ad.read_h5ad(ROOT / cfg['path'], backed='r')
    for sid, part in data.obs.groupby(cfg['slice_key'], observed=True):
        sample_rows.append({'dataset': r.dataset, 'study': r.study, 'sample_id': str(sid), 'n_spots': len(part),
                            'subject_ids': ';'.join(sorted(part.subject.astype(str).unique())) if 'subject' in part else ''})
    data.file.close()
pd.DataFrame(sample_rows).to_csv(ANALYSIS / 'sample_to_study_audit.tsv', sep='\t', index=False)

spec = importlib.util.spec_from_file_location('original_robustness', ROOT / 'scripts/build_submission_robustness_tables.py')
robust = importlib.util.module_from_spec(spec)
spec.loader.exec_module(robust)
_, selected = robust.representative_metrics(ROOT / 'results/spatial_svg_rebuild_v1/frozen_scores')
selected = selected[selected.dataset.isin(DATASETS) & selected.MethodGroup.isin(GROUPS)]
base = robust.task_scores(selected, 'ScaledValue')
assert len(base) == 2 * 5 * 29 and set(base.integration_method) == {'scvi', 'cellcharter'}
mc_reproduced = robust.bootstrap_ranks(base, NBOOT, np.random.default_rng(20261105))
old_mc = pd.read_csv(OLD_ANALYSIS / 'five_dataset_bootstrap_ranks.tsv', sep='\t')
keys = ['integration_method', 'fs_method']
a = mc_reproduced.set_index(keys).sort_index()
b = old_mc.set_index(keys).sort_index()
assert a.index.equals(b.index)
np.testing.assert_allclose(a.to_numpy(dtype=float), b.to_numpy(dtype=float), atol=1e-12, rtol=1e-12)
mc_reproduced.to_csv(ANALYSIS / 'original_primary_bootstrap_reproduction.tsv', sep='\t', index=False)

combined = pd.read_csv(FIVE / 'combined_integrator_scores.tsv', sep='\t')
graph = combined[combined.MethodGroup.isin(GROUPS) & combined.integration_method.eq('graphst')]
assert len(graph) == 5 * 29
graph = graph[['dataset', 'fs_method', 'integration_method', 'CoreOverallMean']].rename(columns={'CoreOverallMean': 'CoreOverall'})
canonical = pd.concat([base[['dataset', 'fs_method', 'integration_method', 'CoreOverall']], graph], ignore_index=True)
exact_samples = np.asarray(list(itertools.product(range(5), repeat=5)))
exact_rows, loo_rows, bootstrap_identity = [], [], []
for integ, part in canonical.groupby('integration_method', observed=True):
    part = part.copy()
    part['DatasetRank'] = part.groupby('dataset').CoreOverall.rank(ascending=False, method='average')
    pivot = part.pivot(index='dataset', columns='fs_method', values='DatasetRank').loc[DATASETS]
    values = pivot.to_numpy()
    assert values.shape == (5, 29) and np.isfinite(values).all()
    observed_mean = values.mean(axis=0)
    observed_rank = rankdata(observed_mean, method='average')
    # Five original studies map one-to-one onto five aggregated dataset scores.
    study_values = pivot.copy()
    study_values.index = [primary.set_index('dataset').loc[x, 'study'] for x in DATASETS]
    np.testing.assert_array_equal(study_values.to_numpy(), values)
    dataset_draws = rankdata(values[exact_samples].mean(axis=1), axis=1, method='average')
    study_draws = rankdata(study_values.to_numpy()[exact_samples].mean(axis=1), axis=1, method='average')
    np.testing.assert_array_equal(dataset_draws, study_draws)
    bootstrap_identity.append({'integration_method': integ, 'original_studies': 5, 'dataset_units': 5,
                               'all_ordered_bootstrap_samples': len(exact_samples), 'dataset_study_draws_identical': True})
    leave_out = np.stack([rankdata(values[np.arange(5) != k].mean(axis=0), method='average') for k in range(5)])
    for j, method in enumerate(pivot.columns):
        exact_rows.append({'integration_method': integ, 'fs_method': method,
                           'analysis_role': 'secondary_canonical' if integ == 'graphst' else 'primary',
                           'ObservedMeanDatasetRank': observed_mean[j], 'ObservedGlobalRank': observed_rank[j],
                           'BootstrapMedian': np.median(study_draws[:, j]),
                           'Lower95': np.quantile(study_draws[:, j], .025), 'Upper95': np.quantile(study_draws[:, j], .975),
                           'Top5Frequency': np.mean(study_draws[:, j] <= 5),
                           'LOSO_MinRank': leave_out[:, j].min(), 'LOSO_MaxRank': leave_out[:, j].max(),
                           'NStudies': 5, 'NExactResamples': len(exact_samples)})
        for k, dataset in enumerate(DATASETS):
            loo_rows.append({'integration_method': integ, 'fs_method': method, 'omitted_dataset': dataset,
                             'omitted_study': primary.set_index('dataset').loc[dataset, 'study'], 'GlobalRank': leave_out[k, j]})
pd.DataFrame(exact_rows).to_csv(ANALYSIS / 'study_bootstrap_and_leave_one_study_out.tsv', sep='\t', index=False)
pd.DataFrame(loo_rows).to_csv(ANALYSIS / 'leave_one_study_out_ranks.tsv', sep='\t', index=False)
pd.DataFrame(bootstrap_identity).to_csv(ANALYSIS / 'study_dataset_bootstrap_identity_audit.tsv', sep='\t', index=False)

scores = pd.read_csv(FIVE / 'setting_task_scores.tsv', sep='\t')
scores = scores[scores.MethodGroup.isin(GROUPS) & scores.n_features.ne('all')].copy()
scores['N'] = pd.to_numeric(scores.n_features).astype(int)
assert len(scores) == 5 * 29 * 7 * 2
assert not scores.duplicated(['dataset', 'fs_method', 'N', 'integration_method']).any()
assert np.isfinite(scores.CoreOverallMean).all()
methods = sorted(scores.fs_method.unique())
composition = []
for (ds, integ, n), part in scores.groupby(['dataset', 'integration_method', 'N']):
    actual_methods = sorted(part.fs_method.tolist())
    assert actual_methods == methods
    composition.append({'dataset': ds, 'integration_method': integ, 'N': n, 'NMethods': len(actual_methods),
                        'method_names': ';'.join(actual_methods),
                        'composition_sha256': hashlib.sha256('\n'.join(actual_methods).encode()).hexdigest()})
pd.DataFrame(composition).to_csv(ANALYSIS / 'budget_method_composition_audit.tsv', sep='\t', index=False)
metadata = pd.read_csv(FIVE / 'setting_metric_seed_summary.tsv', sep='\t')
metadata = metadata[metadata.fs_method.isin(methods) & metadata.n_features.ne('all')].copy()
metadata['N'] = pd.to_numeric(metadata.n_features).astype(int)
effective = metadata.groupby(['dataset', 'fs_method', 'N'], as_index=False).agg(
    EffectiveFeaturesMin=('EffectiveFeaturesMin', 'min'), EffectiveFeaturesMax=('EffectiveFeaturesMax', 'max'),
    MaxEvaluationSeeds=('NSeeds', 'max'))
assert len(effective) == 5 * 29 * 7
configs = dict(zip(primary.dataset, primary.config.str.replace('_wilcoxon', '', regex=False)))
selection_audit = []
for row in effective.itertuples():
    folder = ROOT / 'results/spatial_svg_rebuild_v1' / configs[row.dataset] / row.dataset.lower() / 'feature_selection' / row.fs_method / f'n{row.N}'
    sizes = []
    for seed in range(int(row.MaxEvaluationSeeds)):
        p = folder / f'seed{seed}' / 'selected_features.json'
        if not p.exists():
            p = folder / 'seed0/selected_features.json'
        assert p.exists(), p
        data = json.loads(p.read_text())
        names = data['feature_names']
        assert len(names) == len(set(names)) <= row.N
        sizes.append(len(names))
        selection_audit.append({'dataset': row.dataset, 'method': row.fs_method, 'N': row.N, 'evaluation_seed': seed,
                                'selection_file': str(p.relative_to(ROOT)), 'actual_features': len(names),
                                'sha256': hashlib.sha256(p.read_bytes()).hexdigest()})
    assert min(sizes) == row.EffectiveFeaturesMin and max(sizes) == row.EffectiveFeaturesMax, (row, sizes)
effective['ExactBudget'] = effective.EffectiveFeaturesMin.eq(effective.N) & effective.EffectiveFeaturesMax.eq(effective.N)
effective.to_csv(ANALYSIS / 'requested_vs_effective_feature_counts.tsv', sep='\t', index=False)
pd.DataFrame(selection_audit).to_csv(ANALYSIS / 'selected_feature_count_hash_audit.tsv', sep='\t', index=False)
eligibility = effective.groupby('fs_method', as_index=False).agg(ExactAtAllBudgetsAndDatasets=('ExactBudget', 'all'),
                                                               NExactSettings=('ExactBudget', 'sum'), NSettings=('ExactBudget', 'size'))
common = eligibility.loc[eligibility.ExactAtAllBudgetsAndDatasets, 'fs_method'].tolist()
assert len(common) > 10
eligibility.to_csv(ANALYSIS / 'fixed_exact_budget_panel_methods.tsv', sep='\t', index=False)

rng = np.random.default_rng(SEED)
resamples = rng.integers(0, 5, size=(NBOOT, 5))
all_curves, contrast_rows = [], []
for panel, allowed in [('all_29_requested_budget', methods), ('fixed_exact_budget_panel', common)]:
    sub = scores[scores.fs_method.isin(allowed)]
    bystudy = sub.groupby(['dataset', 'integration_method', 'N'], as_index=False).agg(CoreOverall=('CoreOverallMean', 'mean'), NMethods=('fs_method', 'nunique'))
    assert bystudy.NMethods.eq(len(allowed)).all()
    bystudy['panel'] = panel
    bystudy.to_csv(ANALYSIS / f'{panel}_study_means.tsv', sep='\t', index=False)
    for integ, part in bystudy.groupby('integration_method', observed=True):
        matrix = part.pivot(index='dataset', columns='N', values='CoreOverall').loc[DATASETS, BUDGETS].to_numpy()
        boot = matrix[resamples].mean(axis=1)
        for j, n in enumerate(BUDGETS):
            all_curves.append({'panel': panel, 'integration_method': integ, 'N': n, 'NMethods': len(allowed), 'NStudies': 5,
                               'Mean': matrix[:, j].mean(), 'SDStudies': matrix[:, j].std(ddof=1),
                               'Lower95': np.quantile(boot[:, j], .025), 'Upper95': np.quantile(boot[:, j], .975),
                               'NBootstrap': NBOOT, 'Seed': SEED})
        for n1, n2 in [(500, 1000), (500, 2000), (1000, 2000), (2000, 10000), (500, 10000)]:
            delta = boot[:, BUDGETS.index(n2)] - boot[:, BUDGETS.index(n1)]
            contrast_rows.append({'panel': panel, 'integration_method': integ, 'N_from': n1, 'N_to': n2,
                                  'MeanDifference': matrix[:, BUDGETS.index(n2)].mean() - matrix[:, BUDGETS.index(n1)].mean(),
                                  'Lower95': np.quantile(delta, .025), 'Upper95': np.quantile(delta, .975),
                                  'NStudies': 5, 'NMethods': len(allowed), 'paired_study_resampling': True})
pd.DataFrame(all_curves).to_csv(ANALYSIS / 'budget_curves_with_study_uncertainty.tsv', sep='\t', index=False)
pd.DataFrame(contrast_rows).to_csv(ANALYSIS / 'paired_budget_contrasts.tsv', sep='\t', index=False)
old_budget = pd.read_csv(OLD_ANALYSIS / 'five_dataset_budget_summary.tsv', sep='\t')
new_budget = pd.DataFrame(all_curves).query("panel == 'all_29_requested_budget'")
check = new_budget.merge(old_budget[['integration_method', 'N', 'core_mean']], on=['integration_method', 'N'])
assert len(check) == 14
np.testing.assert_allclose(check.Mean, check.core_mean, atol=1e-12, rtol=1e-12)
# Reproduce the original Figure 3 panel-a standardisation, then resample whole
# study profiles. These bars are on that plot's standardised scale, not raw score.
metric_columns = {'IntegrationMean': 'Integration', 'ClusteringMean': 'Clustering',
                  'AlignmentMean': 'Alignment', 'CoreOverallMean': 'Overall'}
long = scores.melt(id_vars=['dataset', 'fs_method', 'integration_method', 'N'],
                   value_vars=list(metric_columns), var_name='Type', value_name='Value').dropna(subset=['Value'])
long.Type = long.Type.map(metric_columns)
groups = long.groupby(['dataset', 'fs_method', 'integration_method', 'Type']).Value
long['StandardValue'] = (long.Value - groups.transform('mean')) / groups.transform('std')
assert long.StandardValue.notna().all()
plot_study = long.groupby(['dataset', 'Type', 'N'], as_index=False).agg(StudyMean=('StandardValue', 'mean'))
plot_study.to_csv(ANALYSIS / 'figure3_standardized_study_profiles.tsv', sep='\t', index=False)
plot_ci = []
for t, part in plot_study.groupby('Type'):
    matrix = part.pivot(index='dataset', columns='N', values='StudyMean').loc[:, BUDGETS].to_numpy()
    assert np.isfinite(matrix).all()
    bootstrap_indices = np.random.default_rng(SEED).integers(0, len(matrix), size=(NBOOT, len(matrix)))
    draws = matrix[bootstrap_indices].mean(axis=1)
    for j, n in enumerate(BUDGETS):
        plot_ci.append({'Type': t, 'N': n, 'MeanStandardValue': matrix[:, j].mean(),
                        'Lower95': np.quantile(draws[:, j], .025), 'Upper95': np.quantile(draws[:, j], .975),
                        'NStudies': len(matrix), 'NBootstrap': NBOOT})
pd.DataFrame(plot_ci).to_csv(ANALYSIS / 'figure3_panelA_standardized_study_ci.tsv', sep='\t', index=False)
scores.groupby(['integration_method', 'N'], as_index=False).agg(
    MinDownstreamSeeds=('CoreOverallNSeeds', 'min'), MaxDownstreamSeeds=('CoreOverallNSeeds', 'max')).to_csv(
        ANALYSIS / 'budget_downstream_seed_coverage.tsv', sep='\t', index=False)
old7 = pd.read_csv(ROOT / 'results/spatial_svg_rebuild_v1/frozen_scores/setting_task_scores.tsv', sep='\t')
old7 = old7[old7.MethodGroup.isin(GROUPS) & old7.n_features.ne('all')]
old7.groupby(['integration_method', 'n_features'], as_index=False).agg(OriginalSevenDatasetMean=('CoreOverallMean', 'mean'),
    NMethodDatasetSettings=('dataset', 'size'), NDatasets=('dataset', 'nunique')).to_csv(ANALYSIS / 'historical_seven_dataset_curve_reproduction.tsv', sep='\t', index=False)

protocol = {'reviewer_comments': ['R1-m1', 'R1-m3'], 'primary_datasets': DATASETS,
            'resampling_unit': 'Original study or 10x collection; all slices aggregated within each dataset/study',
            'primary_original_studies': primary.study.tolist(), 'primary_study_to_dataset_mapping': 'one-to-one',
            'bootstrap_rank': 'All 5^5=3125 ordered study resamples; original 5000-draw primary Monte Carlo bootstrap also reproduced',
            'curve_bootstrap': NBOOT, 'curve_seed': SEED, 'paired_budget_resampling': True,
            'feature_budgets': BUDGETS, 'primary_curve_methods': len(methods), 'fixed_exact_budget_methods': common,
            'canonical_seed_caveat': 'N=2000 has three downstream seeds; most other budgets have one. Error intervals resample studies, not seeds or methods.',
            'original_rank_reproduction_seed': 20261105,
            'no_new_feature_selection_or_model_fitting': True,
            'scope': 'Original zero-filled aggregation and original frozen scaling; M3 missing-rank sensitivity is not substituted for primary figures'}
(OUT / 'analysis_protocol.json').write_text(json.dumps(protocol, indent=2))
for path, checksum in input_hashes.items():
    assert hashlib.sha256(Path(path).read_bytes()).hexdigest() == checksum
pd.DataFrame([{'path': path, 'sha256': checksum} for path, checksum in input_hashes.items()]).to_csv(
    OUT / 'original_input_hashes.tsv', sep='\t', index=False)
audit = {'all_checks_passed': True, 'original_primary_bootstrap_rows_reproduced': len(mc_reproduced),
         'study_bootstrap_methods_integrators': len(exact_rows), 'exact_resamples_per_integrator': len(exact_samples),
         'study_dataset_bootstrap_draws_identical': True, 'budget_rows': len(scores),
         'budget_method_composition_constant': True, 'NMethodsEveryBudget': len(methods),
         'requested_effective_settings_checked': len(effective), 'selected_feature_records_checked': len(selection_audit),
         'non_exact_budget_settings': int((~effective.ExactBudget).sum()), 'fixed_exact_budget_methods': len(common),
         'original_five_dataset_budget_means_reproduced': True, 'original_inputs_unchanged': True,
         'figure3_standardized_ci_rows': len(plot_ci)}
(ANALYSIS / 'final_audit.json').write_text(json.dumps(audit, indent=2))
print(json.dumps(audit, indent=2), flush=True)
print('Exact-budget methods excluded:', eligibility[~eligibility.ExactAtAllBudgetsAndDatasets].fs_method.tolist(), flush=True)
print(pd.DataFrame(all_curves).query("N in [500,1000,2000,10000]").to_string(index=False), flush=True)
print('Contrasts:', pd.DataFrame(contrast_rows).to_string(index=False), flush=True)
