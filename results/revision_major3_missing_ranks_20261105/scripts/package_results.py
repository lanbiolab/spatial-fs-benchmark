import json
import sys
from pathlib import Path

import pandas as pd

BASE = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parents[1] / 'outputs/GR_R1_major3_sensitivity'
OUT = Path(sys.argv[2]) if len(sys.argv) > 2 else BASE / 'final'
NAMES = {'morans_i': 'Moran’s I', 'somde': 'SOMDE', 'spatialde': 'SpatialDE', 'sparkx': 'SPARK-X', 'nnsvg': 'nnSVG'}
METHODS = list(NAMES)
INTEGRATORS = {'scvi': 'scVI', 'cellcharter': 'CellCharter', 'graphst': 'GraphST'}
u = pd.read_csv(OUT / 'method_gene_universe_by_dataset.tsv', sep='\t')
s = pd.read_csv(OUT / 'method_gene_universe_by_slice.tsv', sep='\t')
changes = pd.read_csv(OUT / 'five_spatial_methods_rank_changes.tsv', sep='\t')
paired = pd.read_csv(OUT / 'five_spatial_methods_paired_refit_rank_changes.tsv', sep='\t')
audit = json.loads((OUT / 'final_audit_summary.json').read_text())
assert audit['all_checks_passed']

universe_rows = []
for (ds, method), g in u.groupby(['dataset', 'method'], sort=False):
    per_slice = s[(s.dataset == ds) & (s.method == method)]
    universe_rows.append({
        '数据集': ds, '分析用途': '胚胎单独报告' if 'Embryo' in ds else '五数据集主分析',
        '方法': NAMES[method], '输入基因数': int(g.total_genes.iloc[0]),
        '切片数': int(g.slices_used.iloc[0]), '筛选种子数': len(g),
        '每片有效排名基因数_最少': int(per_slice.ranked_genes.min()),
        '每片有效排名基因数_最多': int(per_slice.ranked_genes.max()),
        '至少一片有效排名基因数_最少': int(g.union_ranked_genes.min()),
        '至少一片有效排名基因数_最多': int(g.union_ranked_genes.max()),
        '所有切片均有排名基因数_最少': int(g.intersection_ranked_genes.min()),
        '所有切片均有排名基因数_最多': int(g.intersection_ranked_genes.max()),
        '基因切片缺失比例_最少': float((1 - g.ranked_gene_slice_cells / g.total_gene_slice_cells).min()),
        '基因切片缺失比例_最多': float((1 - g.ranked_gene_slice_cells / g.total_gene_slice_cells).max()),
        '新旧2000基因重合数_最少': int(g.new_vs_original_top2000_overlap.min()),
        '新旧2000基因重合数_最多': int(g.new_vs_original_top2000_overlap.max()),
    })
uni = pd.DataFrame(universe_rows)
uni.to_csv(OUT / '基因范围汇总_五方法七数据集.tsv', sep='\t', index=False, encoding='utf-8-sig')
rank_rows = []
for integ in INTEGRATORS:
    for method in METHODS:
        r = changes[(changes.integration_method == integ) & (changes.fs_method == method)].iloc[0]
        p = paired[(paired.integration_method == integ) & (paired.fs_method == method)].iloc[0]
        rank_rows.append({'整合方法': INTEGRATORS[integ], '方法': NAMES[method],
                         '原结果_五数据集平均名次': r.MeanDatasetRank_zero,
                         '旧规则配对基线_五数据集平均名次': p.MeanDatasetRank_zero_refit,
                         '缺失排除_五数据集平均名次': r.MeanDatasetRank_missing,
                         '原结果_全局名次': r.GlobalRank_zero,
                         '旧规则配对基线_全局名次': p.GlobalRank_zero_refit,
                         '缺失排除_全局名次': r.GlobalRank_missing,
                         '原结果_平均综合分': r.MeanCoreOverall_zero,
                         '旧规则配对基线_平均综合分': p.MeanCoreOverall_zero_refit,
                         '缺失排除_平均综合分': r.MeanCoreOverall_missing})
pd.DataFrame(rank_rows).to_csv(OUT / '排名变化汇总_三种整合方法.tsv', sep='\t', index=False, encoding='utf-8-sig')

parts = [
    '# Reviewer 1 Major 3：补充实验最终结果',
    '',
    '夜间计算于 2026-11-06 04:25（北京时间）完成，本轮独立复核全部通过。本文档仅整理实验材料；本轮未修改主稿、主图或 Major 3 回复正文。',
    '',
    '## 审稿人要求与落实',
    '',
    '审稿人要求补充一项敏感性分析：方法在某张切片未给出有效排名的基因不再按百分位 0 计入平均，而只在有有效排名的切片上求平均；展示五种空间方法的排名变化，并报告每个方法在每个数据集中的基因范围。',
    '',
    '已对 Moran’s I、SOMDE、SpatialDE、SPARK-X、nnSVG 完成两种规则的比较。两种规则使用同一份逐切片百分位矩阵，唯一的聚合区别是：旧规则用 0 补缺并除以全部有效切片数；新规则除以该基因实际有有效排名的切片数。所有切片均未排名的基因不进入新规则候选集合。每组均选取 2,000 个基因。',
    '',
    '沿用原数据预处理、筛选参数、整合模型参数和三个下游随机种子。五个有标签数据集为 MouseBrainSerialSections、DLPFC、STOmics0212、STOmics0218、STOmics0224；分别使用 scVI、CellCharter、GraphST 评估。E8.5、E9.5 的基因范围和基因集合变化单独报告，不进入主排名。',
    '',
    '## 五数据集平均名次：主要报告口径',
    '',
    '每个数据集先在 29 个表达驱动/空间方法中按综合分排序，再对五个数据集取平均；数值越小越好。下表列出原冻结结果、旧规则配对基线和新规则。原规则与新规则均使用冻结的指标归一化范围。',
    '',
]
for integ, title in INTEGRATORS.items():
    parts += [f'### {title}', '', '| 方法 | 原冻结结果 | 旧规则配对基线 | 新规则：只对可用切片平均 |', '| --- | ---: | ---: | ---: |']
    for method in METHODS:
        r = changes[(changes.integration_method == integ) & (changes.fs_method == method)].iloc[0]
        p = paired[(paired.integration_method == integ) & (paired.fs_method == method)].iloc[0]
        parts.append(f'| {NAMES[method]} | {r.MeanDatasetRank_zero:.1f} | {p.MeanDatasetRank_zero_refit:.1f} | {r.MeanDatasetRank_missing:.1f} |')
    parts.append('')
parts += [
    '“平均名次”与“全局名次”不同：全局名次是把上述平均名次再次排序。两种数值均保存在 `排名变化汇总_三种整合方法.tsv`；不能把全局名次写成平均名次。审稿人引述的旧稿 20.1、9.4 来自原先分析口径，以上结果统一按目前约定的五数据集口径计算。',
    '',
    '旧规则配对基线的处理：49 个筛选任务都准确复现原来的前 2,000 个基因集合，其中 25 个顺序也完全相同；其余 24 个任务集合相同但顺序不同。对主分析中顺序不同的情况，用重建的旧规则基因顺序重新运行原下游流程；顺序完全相同的情况复用原结果。两个基线均完整保留，避免把原结果与重跑之间的变化全部归因于缺失处理。',
    '',
    '## 可直接用于回复的事实',
    '',
    '- 新规则没有一致地提高或降低五个方法的表现，变化取决于方法和整合器。',
    '- scVI 下，SPARK-X 的五数据集平均名次从 24.4 变为 24.6，SOMDE 从 12.2 变为 11.8；配对旧规则基线分别为 24.6、12.2。因此，该敏感性分析没有消除二者的排名差距。',
    '- nnSVG 的平均名次在 scVI 下从 14.2 变为 12.6，CellCharter 下从 14.4 变为 11.2，GraphST 下从 13.4 变为 15.4。应分别报告，不概括为所有整合器下均改善。',
    '- SpatialDE 在 scVI 下原冻结结果为 15.2，新规则为 13.4；旧规则配对基线为 13.2。因此，不能把相对于原冻结结果的改善全部归因于缺失排除。',
    '',
    '## 基因范围',
    '',
    '完整数据见 `基因范围汇总_五方法七数据集.tsv`；逐切片、逐筛选种子明细见 `method_gene_universe_by_slice.tsv`。有效排名范围按实际输出统计，包含内部过滤和其他导致无有效分数的情况，不能把所有缺失都解释为同一种生物学过滤。',
    '',
    '下表以 nnSVG 为例：“并集”为至少一张切片有排名的基因数；“交集”为所有切片均有排名的基因数；“缺失比例”为未排名的基因×切片位置占全部基因×切片位置的比例。',
    '',
    '| 数据集 | 输入基因数 | 并集 | 交集 | 基因×切片缺失比例 | 新旧2000基因重合数 |',
    '| --- | ---: | ---: | ---: | ---: | ---: |',
]
for r in u[u.method == 'nnsvg'].itertuples():
    parts.append(f'| {r.dataset} | {r.total_genes:,} | {r.union_ranked_genes:,} | {r.intersection_ranked_genes:,} | {1-r.ranked_gene_slice_cells/r.total_gene_slice_cells:.1%} | {r.new_vs_original_top2000_overlap:,} |')
parts += [
    '',
    '## 完整性核验',
    '',
    f'- 49/49 筛选任务通过：逐切片有效基因数量、两种均值公式、选中基因分数、前 2,000 基因及原基因集合均重新核对。',
    f'- 225/225 新规则下游组合完成（5 数据集 × 5 方法 × 3 种子 × 3 整合器）；另有 66 项旧规则基线，共 291 项结果或原结果复用检查，其中实际重算 237 项、完全相同顺序的原结果复用 54 项。',
    f'- 75 个原冻结综合分重新计算全部一致，最大绝对误差 {audit["max_frozen_score_difference"]:.3g}。',
    '- 未发现未完成下游任务或 FAILED 标记。检查明细保存在 `final_selection_audit.tsv`、`final_downstream_audit.tsv` 和 `final_audit_summary.json`。',
    '',
    '服务器归档目录：`/home/xd/y/Benchmark/results/revision_major3_missing_ranks_20261105`。完整原始逐切片矩阵、两种基因列表、评估记录和日志均保留在服务器；本地 `final/` 保存最终汇总、指标明细及复核记录。',
    '',
]
(BASE / 'Major3_补充实验最终结果.md').write_text('\n'.join(parts), encoding='utf-8-sig')
print('Packaged verified results:', BASE / 'Major3_补充实验最终结果.md')
