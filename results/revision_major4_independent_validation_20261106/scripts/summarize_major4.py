import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

ROOT = Path('/home/xd/y/Benchmark')
OUT = ROOT / 'results/revision_major4_independent_validation_20261106'
ANALYSIS = OUT / 'analysis'
ANALYSIS.mkdir(exist_ok=True)
spec = json.loads((OUT / 'analysis_protocol.json').read_text())
rows, fold_audits = [], []
for fold in spec['folds']:
    job = OUT / 'folds' / fold
    assert (job / 'DONE').exists(), fold
    fold_audits.append(json.loads((job / 'fold_audit.json').read_text()))
    for method in spec['methods']:
        for seed in spec['evaluation_seeds']:
            r = json.loads((job / f'{method}_seed{seed}.json').read_text())
            with np.load(job / f'{method}_seed{seed}.npz') as z:
                mask = z['valid_annotation']
                for key, value in {
                    'legacy_ari': adjusted_rand_score(z['true_labels'], z['original_predictions']),
                    'legacy_nmi': normalized_mutual_info_score(z['true_labels'], z['original_predictions']),
                    'valid_annotation_ari': adjusted_rand_score(z['true_labels'][mask], z['valid_annotation_predictions'][mask]),
                    'valid_annotation_nmi': normalized_mutual_info_score(z['true_labels'][mask], z['valid_annotation_predictions'][mask]),
                }.items():
                    assert np.isfinite(value) and abs(value - r[key]) < 1e-12
                assert len(z['spot_ids']) == len(set(z['spot_ids'])) == len(mask)
            rows.append(r)
d = pd.DataFrame(rows)
assert len(d) == 252 and not d.duplicated(['fold', 'method', 'evaluation_seed']).any()
old = d[d.evaluation_seed.eq(0)]
assert len(old) == 84 and old.seed0_absolute_ari_difference.max() < 1e-8, old.sort_values('seed0_absolute_ari_difference').tail()
d.to_csv(ANALYSIS / 'heldout_clustering_all_seed_results.tsv', sep='\t', index=False)
pd.DataFrame(fold_audits).to_csv(ANALYSIS / 'heldout_fold_annotation_audit.tsv', sep='\t', index=False)
summaries = {}
for label, column in [('original_protocol', 'legacy_ari'), ('valid_annotation', 'valid_annotation_ari')]:
    for seeds, data in [('seed0', old), ('three_seeds', d)]:
        fold_means = data.groupby(['dataset', 'fold', 'method', 'method_role'], as_index=False).agg(
            ARI=(column, 'mean'), NSeeds=('evaluation_seed', 'nunique'))
        fold_means['WithinFoldCompetitiveRank'] = np.nan
        mask = fold_means.method_role.eq('competitive')
        fold_means.loc[mask, 'WithinFoldCompetitiveRank'] = fold_means[mask].groupby('fold').ARI.rank(ascending=False, method='average')
        summary = fold_means.groupby(['method', 'method_role'], as_index=False).agg(
            MeanARI=('ARI', 'mean'), SDFoldARI=('ARI', 'std'), MinFoldARI=('ARI', 'min'), MaxFoldARI=('ARI', 'max'),
            NFolds=('fold', 'nunique'), MeanWithinFoldCompetitiveRank=('WithinFoldCompetitiveRank', 'mean'))
        assert summary.NFolds.eq(6).all()
        summary['RankByMeanARI'] = np.nan
        mask = summary.method_role.eq('competitive')
        summary.loc[mask, 'RankByMeanARI'] = summary.loc[mask, 'MeanARI'].rank(ascending=False, method='average')
        summary = summary.sort_values('MeanARI', ascending=False)
        summary.to_csv(ANALYSIS / f'{label}_{seeds}_method_summary.tsv', sep='\t', index=False)
        fold_means.to_csv(ANALYSIS / f'{label}_{seeds}_fold_results.tsv', sep='\t', index=False)
        fold_means.groupby(['dataset', 'method', 'method_role'], as_index=False).agg(MeanARI=('ARI', 'mean'), NFolds=('fold', 'nunique')).to_csv(
            ANALYSIS / f'{label}_{seeds}_dataset_summary.tsv', sep='\t', index=False)
        summaries[(label, seeds)] = summary

marker_frames = [pd.read_csv(OUT / 'folds' / fold / 'external_marker_recovery.tsv', sep='\t')
                 for fold in spec['folds'] if fold.startswith('DLPFC__')]
markers = pd.concat(marker_frames, ignore_index=True)
assert len(markers) == 2 * 14 * 9 and markers.background_genes.nunique() == 1
assert (markers.overlap <= markers.available_markers).all() and (markers.overlap <= markers.selected_genes).all()
assert markers[markers.marker_set.eq('All human-supported layers')].available_markers.nunique() == 1
markers.to_csv(ANALYSIS / 'dlpfc_external_marker_recovery_by_fold.tsv', sep='\t', index=False)
marker_summary = markers.groupby(['method', 'method_role', 'marker_set'], as_index=False).agg(
    NFolds=('fold', 'nunique'), AvailableMarkers=('available_markers', 'first'), MeanOverlap=('overlap', 'mean'),
    MeanRecall=('marker_recall', 'mean'), MinRecall=('marker_recall', 'min'), MaxRecall=('marker_recall', 'max'),
    MeanUniformEnrichment=('uniform_fold_enrichment', 'mean'), MeanMatchedEnrichment=('matched_fold_enrichment', 'mean'),
    MeanMatchedExpectedOverlap=('matched_expected_overlap', 'mean'))
assert marker_summary.NFolds.eq(2).all()
marker_summary.to_csv(ANALYSIS / 'dlpfc_external_marker_summary.tsv', sep='\t', index=False)
mapping = pd.read_csv(OUT / 'folds/DLPFC__heldout_151507/external_marker_mapping.tsv', sep='\t')
other = pd.read_csv(OUT / 'folds/DLPFC__heldout_151673/external_marker_mapping.tsv', sep='\t')
assert mapping.equals(other)
human = mapping[mapping.human_supported & mapping.dataset_gene_id.notna()]
for fold, fold_rows in markers.groupby('fold'):
    background_bins = pd.read_csv(OUT / 'folds' / fold / 'training_background_bins.tsv', sep='\t')
    universe = set(background_bins.gene)
    bin_sizes = background_bins.groupby('joint_bin').size()
    for r in fold_rows.itertuples():
        saved = json.loads((ROOT / 'results/validation_v1/heldout_slice' / fold / 'rankings' / f'{r.method}.json').read_text())
        selected = [str(g) for g in saved['feature_names'] if str(g) in universe]
        if r.method != 'all_features':
            selected = selected[:2000]
        selected = set(selected)
        if r.marker_set == 'All human-supported layers':
            expected_markers = set(human.dataset_gene_id)
        else:
            layer = r.marker_set.removeprefix('L')
            expected_markers = set(human.loc[human[layer].eq(1), 'dataset_gene_id'])
        hits = selected & expected_markers
        assert len(hits) == r.overlap and len(expected_markers) == r.available_markers
        assert set(str(r.recovered_marker_ids).split(';')) - {'nan', ''} == hits
        assert abs(r.marker_recall - len(hits) / len(expected_markers)) < 1e-12
        selected_counts = background_bins[background_bins.gene.isin(selected)].groupby('joint_bin').size()
        marker_counts = background_bins[background_bins.gene.isin(expected_markers)].groupby('joint_bin').size()
        expectation = (selected_counts * marker_counts / bin_sizes).sum()
        assert abs(expectation - r.matched_expected_overlap) < 1e-9
        assert r.matched_null_lower <= r.matched_null_upper
mapping.to_csv(ANALYSIS / 'external_marker_mapping_audit.tsv', sep='\t', index=False)
for r in pd.read_csv(OUT / 'original_input_hashes.tsv', sep='\t').itertuples():
    assert hashlib.sha256(Path(r.path).read_bytes()).hexdigest() == r.sha256, r.path
source = json.loads((OUT / 'sources/source_manifest.json').read_text())
for r in source['files']:
    assert hashlib.sha256((OUT / 'sources' / Path(r['file']).name).read_bytes()).hexdigest() == r['sha256']
audit = {'all_checks_passed': True, 'folds': 6, 'methods_including_references': 14,
         'competitive_selectors': 11, 'evaluation_seeds': 3, 'selection_seed': 0,
         'clustering_jobs': len(d), 'prediction_archives_independently_scored': len(d),
         'original_seed0_scores_reproduced': len(old), 'max_original_ari_difference': float(old.seed0_absolute_ari_difference.max()),
         'original_inputs_unchanged': True, 'marker_method_fold_set_rows': len(markers),
         'marker_overlap_and_matched_expectation_independently_verified': True,
         'source_marker_rows': len(mapping), 'human_supported_marker_rows': int(mapping.human_supported.sum()),
         'mapped_human_supported_markers': int((mapping.human_supported & mapping.dataset_gene_id.notna()).sum()),
         'source_symbol_id_disagreements_not_used_for_mapping': int((mapping.dataset_gene_id.notna() & ~mapping.source_id_agrees).sum())}
(ANALYSIS / 'final_audit.json').write_text(json.dumps(audit, indent=2))
names = {'morans_i': "Moran’s I", 'nnsvg': 'nnSVG', 'somde': 'SOMDE', 'spatialde': 'SpatialDE', 'sparkx': 'SPARK-X',
         'scanpy_seurat_v3_batch': 'Scanpy Seurat v3, batch', 'scanpy_cell_ranger_batch': 'Scanpy Cell Ranger, batch',
         'singleCellHaystack': 'singleCellHaystack', 'seurat_vst': 'Seurat VST', 'scsegindex': 'scSEGIndex', 'triku': 'triku',
         'random': 'Random (2,000)', 'all_features': 'All features', 'wilcoxon': 'Wilcoxon (label-informed)'}
legacy = summaries[('original_protocol', 'seed0')].set_index('method')
new = summaries[('valid_annotation', 'three_seeds')].set_index('method')
marker_all = marker_summary[marker_summary.marker_set.eq('All human-supported layers')].set_index('method')
parts = ['# Reviewer 1 Major 4：独立参照补充实验', '',
         '已完成六个留出切片的直接 ARI 汇总、三个聚类随机种子的稳定性检查，以及 DLPFC 外部皮层分层标记基因恢复分析。本轮不修改主稿、主图或回复信。', '',
         '## 1. 留出切片 ARI：直接评价标注一致性', '',
         '沿用原来的六个留出切片及训练切片上选出的基因集合：DLPFC 151507/151673、小鼠脑 section1/section2、STOmics0212 的 AD_4_LS/HE_3。11 个竞争选择器与 Random、All features、label-informed Wilcoxon 参照均完整保留。除 All features 使用全部基因外，预算均为 2,000。', '',
         '原始 seed 0 的 ARI 全部重新计算核对。新增 seed 1、2 只检验固定基因集合下的 SVD/聚类稳定性，没有重选基因。另一套有效标注分析仅用非缺失参考标注确定聚类数并计算 ARI；聚类输入仍包括留出切片的全部点，点的真实标签不作为聚类输入。', '',
         '下表两列计算口径不同，不应用它们的差值估计方法改进：左列为原流程、seed 0；右列为有效标注分析、三个聚类种子平均。各列均对六折等权平均。完整的同口径比较、逐折及逐种子结果另存 TSV。ARI 越大表示与参考标注越一致。', '',
         '| 方法 | 原流程 ARI，seed 0 | 有效标注 ARI，3 种子 | 用途 |', '| --- | ---: | ---: | --- |']
for method in new.sort_values('MeanARI', ascending=False).index:
    parts.append(f'| {names[method]} | {legacy.loc[method, "MeanARI"]:.4f} | {new.loc[method, "MeanARI"]:.4f} | {"竞争方法" if new.loc[method, "method_role"] == "competitive" else "参照"} |')
parts += ['', '缺失标注处理明细：', '', '| 留出切片 | 点数 | 无有效标注点数 | 原聚类数 | 有效标注聚类数 |', '| --- | ---: | ---: | ---: | ---: |']
for r in fold_audits:
    parts.append(f'| {r["dataset"]} / {r["heldout_slice"]} | {r["n_heldout_spots"]} | {r["excluded_unannotated_spots"]} | {r["legacy_k"]} | {r["valid_annotation_k"]} |')
parts += ['', '## 2. 外部标记恢复：独立于空间自相关统计量', '',
          '参考集采用 Zeng et al., Cell (2012) 的皮层分层标记，由 Maynard et al. 原作者 HumanPilot 仓库提供。原始研究在独立的人类皮层组织中进行原位杂交；其名单早于本研究，并非在本次 DLPFC 留出切片上重新做差异分析得到。该名单是外部皮层参考，不能称为在本次 DLPFC 上独立发现的标记或全基因组生物学真值。', '',
          f'原表共有 {len(mapping)} 行，其中 H/MH 人类支持条目 {audit["human_supported_marker_rows"]} 个。按数据中的 gene_name 做唯一、不区分大小写的精确匹配后，{audit["mapped_human_supported_markers"]} 个进入可测标记分母。未匹配条目和不一致的原始 Ensembl 注释均在 mapping audit 中保留；没有猜测别名，也没有使用不一致的源 Ensembl ID 强行映射。', '',
          '分别用两个 DLPFC 留出折的原始训练基因集合计算：整体和逐层标记召回、相对均匀随机背景的富集倍数，以及按训练数据表达量和检出率匹配背景的富集倍数。匹配背景联合使用两个变量各 10 个分位组，并固定 10,000 次条件超几何抽样。没有按结果筛选层、基因或方法。', '',
          '| 方法 | 两折平均命中标记数 | 标记召回 | 均匀背景富集倍数 | 表达/检出匹配富集倍数 |', '| --- | ---: | ---: | ---: | ---: |']
for method in marker_all.sort_values('MeanRecall', ascending=False).index:
    r = marker_all.loc[method]
    parts.append(f'| {names[method]} | {r.MeanOverlap:.1f}/{r.AvailableMarkers} | {r.MeanRecall:.1%} | {r.MeanUniformEnrichment:.2f} | {r.MeanMatchedEnrichment:.2f} |')
parts += ['', '## 3. 这条意见如何落实', '',
          '将留出标注 ARI 作为与空间自相关重复性并列的评价结果，而不是只报告二者相关性；外部标记恢复提供另一项不使用 Moran’s I/nnSVG 统计量定义参考集的证据。方法优劣按实际结果分别报告，不预设空间方法在所有端点上领先。两项结果分别呈现，不合成为新的总分。摘要和主稿的最终措辞待写回复及统一修改主稿时处理。', '',
          '原 Figure 2 及其他主图未修改。新增图表使用原绘图脚本中的主题和方法名称，作为独立补充候选图保存。', '',
          '## 4. 核验与来源', '',
          f'- 252 项聚类评估完成，预测结果均独立重算 ARI/NMI 通过；84 个原 seed-0 ARI 的最大绝对误差为 {audit["max_original_ari_difference"]:.3g}。',
          '- 原始选择、指标与状态文件的哈希保持一致；六个折的留出切片均未进入对应训练切片名单。',
          '- 14 方法 × 2 DLPFC 折 × 9 个标记集合的 252 行结果已齐全，包括总集合和 L1–L6、L6b、WM。',
          '- 主要表格见 analysis/；逐点预测、原标注、匹配背景和日志完整保留在服务器 folds/ 中。', '',
          '[Zeng et al., Cell (2012)](https://doi.org/10.1016/j.cell.2012.02.052)；[Maynard et al., Nature Neuroscience (2021)](https://doi.org/10.1038/s41593-020-00787-0)。', '',
          '冻结的原作者名单：' + source['files'][0]['url'], '', '版本：`' + source['commit_or_ref'] + '`。下载校验值见 sources/source_manifest.json。', '']
(OUT / 'Major4_补充实验结果.md').write_text('\n'.join(parts), encoding='utf-8-sig')
(OUT / 'COMPLETE').write_text('all independent checks passed\n')
print(json.dumps(audit, indent=2), flush=True)
print(summaries[('original_protocol', 'seed0')][['method', 'MeanARI']].to_string(index=False), flush=True)
print(marker_all[['MeanOverlap', 'MeanRecall', 'MeanMatchedEnrichment']].to_string(), flush=True)
