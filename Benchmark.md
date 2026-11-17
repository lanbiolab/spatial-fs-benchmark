



# 空间组学中“特征选择影响多切片整合与下游任务”的 Benchmark 重新实施说明

## 1. 项目概述

本项目是一个**新的空间转录组 benchmark**，灵感来自两篇文章：

1. **Feature selection methods affect the performance of scRNA-seq data integration and querying** 

   项目代码：https://github.com/theislab/atlas-feature-selection-benchmark

2. **Benchmarking multi-slice integration and downstream applications in spatial transcriptomics data analysis**

   Python and R were used for data preprocessing, format conversion, and evaluation metric calculation. The relevant code is available at “https://github.com/bm2-lab/iSTBench” [41]. R scripts for result analysis and visualization can be accessed through “https://github.com/bm2-lab/iSTBench/tree/main/Analysis.” Our benchmarking workflow is provided as a reproducible pipeline at “https://github.com/bm2-lab/iSTBench/tree/main/Benchmark.” We also uploaded the relevant code to Zenodo which is available at https://zenodo.org/records/17089390 [44].

项目的核心思想是：

> 把第一篇 scRNA-seq 文章提出的核心问题迁移到空间组学中：  
> **不同的特征选择（feature selection）策略，会如何影响空间转录组多切片整合以及下游任务表现？**

这个项目**不是主要做一个新的模型**，而是一个：

- **benchmark / evaluation framework**
- 聚焦于**前处理（preprocessing）**
- 重点研究**feature selection** 这一步对后续分析的影响

重点关注的任务包括：

- 多切片整合（multi-slice integration）
- 空间聚类 / 空间域识别（spatial clustering / domain identification）
- 切片对齐（slice alignment）
- 切片表征（slice representation）
- 可选：reference-query mapping / label transfer / unseen domain detection

---

## 2. 核心研究问题

项目要回答的问题是：

### 总问题
**在空间转录组数据分析中，feature selection 是否会显著影响 multi-slice integration 以及下游任务的性能？**

### 细化问题
1. 不同的 feature selection 方法是否会显著改变 integration 效果？
2. 不同 integration 方法对 feature selection 的敏感性是否不同？
3. 哪类 feature selection 更适合哪类空间任务？
4. 最优 feature selection 是否依赖于：
   - 平台（Visium / Stereo-seq / MERFISH 等）
   - 切片数量
   - 分辨率（spot-level / cell-level）
   - 组织复杂度
   - 下游任务类型
5. scRNA-seq 中常见的 HVG 是否在空间组学中仍然最优，还是 spatially variable genes / cross-slice conserved genes 更优？

---

## 3. 项目定位

这个项目应被定位为：

### 3.1 一个空间组学 benchmark
核心是系统评估，而不是发明单一新方法。

### 3.2 一个 preprocessing-aware benchmark
重点不是只比较 integration model，而是研究：
**feature selection → integration → downstream tasks**
这条链路上的影响。

### 3.3 一个 data-centric 计算组学项目
该项目强调：模型性能不仅取决于模型本身，也取决于输入特征如何构建。

---

## 4. 与两篇论文的关系

### 4.1 第一篇论文提供什么
第一篇论文提供的是**主问题和 benchmark 哲学**：

- 把 feature selection 作为主要自变量
- 控制 integration / query 方法
- 比较多个下游任务
- 输出不同场景下的推荐策略

因此，第一篇适合作为：

- 实验设计模板
- benchmark 结构模板
- 配置和运行逻辑模板
- 结果汇总与统计模板

### 4.2 第二篇论文提供什么
第二篇论文提供的是**空间组学任务层**：

- multi-slice integration
- spatial clustering
- slice alignment
- slice representation
- 不同空间平台和数据集下的评估思路

因此，第二篇适合作为：

- 空间任务定义来源
- 空间数据集组织方式参考
- 空间任务指标参考
- multi-slice spatial benchmark 范式来源

### 4.3 本项目的正确组合方式
最合理的思路是：

> **以第一篇的 benchmark 骨架为基础，加入第二篇的空间任务层，重构为新的 spatial benchmark。**

换句话说：

- **复用第一篇的“benchmark 框架思想”**
- **采用第二篇的“空间任务设计”**
- **不要原样照搬第一篇的 scRNA-specific 实现**

---

## 5. 是否直接复用第一篇代码

结论：

> **可以优先以第一篇代码框架为基础，但只能复用“实验骨架”和模块化设计，不能直接照搬其任务实现。**

### 值得复用的部分
如果第一篇代码写得比较规范，优先复用以下内容：

- 数据预处理总流程
- feature selection 模块接口
- benchmark 配置文件结构
- 批量跑实验的 runner
- 结果记录和汇总逻辑
- 表格生成和画图脚本
- ranking / stability analysis 框架

### 必须重写的部分
下面这些部分几乎一定要按空间组学重写：

- 数据对象定义（加入 spatial coordinates / slice ID）
- multi-slice integration 接口
- spatial clustering 任务接口
- alignment 任务接口
- slice representation 任务接口
- spatial-aware metrics
- reference-query mapping 在空间场景下的定义
- 部分数据预处理步骤（如空间图构建）

### 实际原则
**复用“外壳”，重写“内核”。**

---

## 6. 建议的 benchmark 总体结构

建议将项目拆成以下层次：

### 6.1 输入层：数据集
输入应为多个空间转录组数据集，每个数据集包含多个切片，至少包括：

- expression matrix
- gene list
- slice ID
- spatial coordinates
- optional annotation（domain / cell type / region）
- optional registration / alignment truth
- platform metadata

### 6.2 自变量 1：feature selection
候选 feature selection 方法可以包括：

- HVG（Highly Variable Genes）
- SVG（Spatially Variable Genes）
- 高表达基因（Highly Expressed Genes）
- marker genes
- cross-slice conserved genes
- random genes（作为 baseline）
- union / intersection 策略
- supervised feature selection（如果后续扩展）

### 6.3 自变量 2：integration methods
选取若干空间组学 multi-slice integration 方法进行组合测试。

建议初期先选 2–4 个代表方法，而不是一开始铺太大。

### 6.4 输出层：下游任务
每组 feature selection × integration method 组合都评估若干任务：

- integration quality
- spatial clustering
- alignment
- slice representation
- optional query mapping / label transfer

---

## 7. 最小可实施版本（MVP）

建议先做一个最小 benchmark，不要一开始就全铺开。

### 7.1 最小实验矩阵
建议起步配置：

- **feature selection 方法：3–4 个**
- **integration 方法：2–3 个**
- **数据集：2–3 个**
- **任务：3 个**
- **每个设置重复运行若干次**

### 7.2 推荐最小 feature selection 集合
建议第一版先用：

1. HVG
2. SVG
3. 高表达基因
4. random genes（负对照）

如果资源允许，可以再加：

5. cross-slice conserved genes

### 7.3 推荐最小任务集合
建议第一版先评估：

1. multi-slice integration quality
2. spatial clustering / domain identification
3. slice alignment

如果资源允许，再加：

4. slice representation
5. spatial query / label transfer

---

## 8. 建议代码架构

建议使用如下 repo 结构：

```text
spatial-fs-benchmark/
│
├── configs/
│   ├── datasets/
│   ├── feature_selection/
│   ├── integration/
│   ├── tasks/
│   └── benchmark/
│
├── data/
│   ├── raw/
│   ├── processed/
│   └── metadata/
│
├── src/
│   ├── data/
│   │   ├── io.py
│   │   ├── preprocess.py
│   │   └── spatial_object.py
│   │
│   ├── feature_selection/
│   │   ├── base.py
│   │   ├── hvg.py
│   │   ├── svg.py
│   │   ├── highly_expressed.py
│   │   ├── conserved.py
│   │   └── random_baseline.py
│   │
│   ├── integration/
│   │   ├── base.py
│   │   ├── method_a.py
│   │   ├── method_b.py
│   │   └── wrappers.py
│   │
│   ├── tasks/
│   │   ├── integration_eval.py
│   │   ├── clustering_eval.py
│   │   ├── alignment_eval.py
│   │   ├── slice_representation_eval.py
│   │   └── query_mapping_eval.py
│   │
│   ├── metrics/
│   │   ├── biological_conservation.py
│   │   ├── batch_mixing.py
│   │   ├── spatial_coherence.py
│   │   ├── clustering_metrics.py
│   │   ├── alignment_metrics.py
│   │   └── ranking.py
│   │
│   ├── benchmark/
│   │   ├── runner.py
│   │   ├── experiment.py
│   │   ├── logging.py
│   │   └── result_schema.py
│   │
│   └── plotting/
│       ├── summary_plots.py
│       ├── heatmaps.py
│       ├── rankings.py
│       └── ablations.py
│
├── notebooks/
│   ├── exploratory/
│   └── figures/
│
├── scripts/
│   ├── run_benchmark.py
│   ├── prepare_data.py
│   ├── summarize_results.py
│   └── generate_figures.py
│
├── results/
│   ├── intermediate/
│   ├── tables/
│   ├── plots/
│   └── logs/
│
├── README.md
└── requirements.txt
```

## 9. 建议统一的数据对象

建议定义一个统一的 spatial data object，至少包含：

- `X`: expression matrix
- `genes`: gene names
- `obs`: spot / cell metadata
- `coords`: spatial coordinates
- `slice_id`: 每个 spot/cell 所属切片
- `labels`: 可选真值标签（domain / cell type）
- `platform`: 平台信息
- `adjacency_graph`: 可选空间邻接图
- `image`: 可选组织图像

### 示例接口

```
class SpatialDataset:
    def __init__(
        self,
        X,
        genes,
        obs,
        coords,
        slice_id,
        labels=None,
        platform=None,
        adjacency_graph=None,
        image=None,
    ):
        self.X = X
        self.genes = genes
        self.obs = obs
        self.coords = coords
        self.slice_id = slice_id
        self.labels = labels
        self.platform = platform
        self.adjacency_graph = adjacency_graph
        self.image = image
```

### 统一要求

所有 feature selection、integration、tasks 都应围绕这个统一对象开发，避免方法之间各自定义输入格式。

------

## 10. Feature selection 模块设计建议

建议定义统一父类：

```
class FeatureSelector:
    name = "base"

    def fit(self, dataset):
        raise NotImplementedError

    def select(self, dataset, n_features):
        raise NotImplementedError
```

### 各方法建议实现

- `HVGSelector`
- `SVGSelector`
- `HighlyExpressedSelector`
- `RandomSelector`
- `CrossSliceConservedSelector`

### 输出格式建议

统一输出为：

```
selected_features = {
    "feature_names": [...],
    "feature_indices": [...],
    "scores": [...],  # optional
    "metadata": {...}
}
```

### 关键注意事项

1. 所有 selector 要保证输出特征数量可控，比如 `n_features=1000/2000/3000`
2. 所有 selector 需要支持多切片输入
3. 需要明确：
   - 是先按切片选，再合并
   - 还是全局选
4. 对于跨切片策略，要保留 union / intersection 的可配置性

------

## 11. Integration 模块设计建议

定义统一接口：

```
class SpatialIntegrator:
    name = "base"

    def fit_transform(self, dataset, selected_features):
        raise NotImplementedError
```

### 输出建议

每个 integration 方法至少输出：

- integrated embedding
- corrected expression（如果有）
- slice-level embedding（如果有）
- method metadata

例如：

```
integration_output = {
    "cell_embedding": ...,
    "corrected_X": ...,
    "slice_embedding": ...,
    "metadata": {...}
}
```

### 设计原则

1. 不要在 integration 模块里隐式重新做 feature selection
2. 输入特征必须由外部显式传入
3. 所有 integration 方法都应统一输出 embedding，便于后续下游任务复用

------

## 12. 下游任务设计

### 12.1 任务一：integration quality evaluation

目标：评估整合后是否实现跨切片融合，同时保留生物学和空间结构。

可考虑指标：

- batch mixing / slice mixing
- biological conservation
- neighborhood preservation
- domain consistency

### 12.2 任务二：spatial clustering / domain identification

目标：评估 embedding 是否有利于空间域识别。

可考虑指标：

- ARI
- NMI
- domain purity
- spatial coherence
- boundary consistency

### 12.3 任务三：slice alignment

目标：评估不同切片是否被正确对齐。

可考虑指标：

- alignment accuracy
- registration error
- matched neighborhood consistency
- landmark consistency（如果有真值）

### 12.4 任务四：slice representation

目标：评估方法是否学到合理的切片级表示。

可考虑指标：

- slice retrieval accuracy
- 同组织/同条件切片是否更接近
- 低维可分性
- slice-level clustering quality

### 12.5 可选任务五：reference-query mapping

如果要真正把第一篇的问题迁移过来，可以加入：

- reference slice / reference dataset
- query slice 映射到 reference latent space
- label transfer
- unknown domain detection / unseen population detection

------

## 13. 指标体系建议

指标需要分层组织，否则后期难以解释结果。

### 13.1 integration 层

- batch / slice mixing
- biological conservation
- neighborhood preservation
- latent structure stability

### 13.2 spatial 层

- spatial autocorrelation preservation
- spatial coherence
- adjacency consistency
- local smoothness / boundary sharpness

### 13.3 task 层

- clustering: ARI / NMI / silhouette
- alignment: registration / matching-related metrics
- slice representation: retrieval / ranking metrics
- query mapping: transfer accuracy / OOD detection

### 13.4 robustness 层

- seed stability
- dataset stability
- platform transferability
- feature number sensitivity

------

## 14. 推荐的实验设计

建议使用如下实验主循环：

```
for dataset in datasets:
    load dataset
    preprocess dataset

    for fs_method in feature_selection_methods:
        selected_features = fs_method.select(dataset, n_features=k)

        for integration_method in integration_methods:
            embedding = integration_method.fit_transform(dataset, selected_features)

            for task in tasks:
                score = evaluate(task, dataset, embedding)

                save result
```

### 建议控制变量

必须尽量控制以下因素：

- 相同数据预处理流程
- 相同特征数
- 相同随机种子集合
- 相同下游聚类参数或调参策略
- 相同评价脚本

### 建议重复实验

每个设置至少跑多个随机种子，避免偶然性。

------

## 15. 结果保存格式建议

建议每次实验保存标准化结果表，至少包含：

- dataset
- platform
- n_slices
- fs_method
- n_features
- integration_method
- task
- metric_name
- metric_value
- random_seed
- runtime
- memory_usage
- notes

例如 csv 表头：

```
dataset,platform,n_slices,fs_method,n_features,integration_method,task,metric_name,metric_value,random_seed,runtime,memory_usage,notes
```

### 为什么必须标准化

因为 benchmark 后期最耗时的不是跑实验，而是：

- 汇总
- 画图
- 排名
- 写结论

如果结果格式不统一，后期基本会崩。

------

## 16. 推荐画图与结果展示形式

建议统一生成以下图：

### 16.1 总体 performance heatmap

行：feature selection
 列：integration method × task
 值：标准化 score

### 16.2 ranking 图

展示不同 feature selection 的平均排名和稳定性。

### 16.3 feature number ablation

横轴：特征数
 纵轴：任务性能
 比较不同 feature selection 方法在不同特征数下的表现。

### 16.4 dataset-wise 分图

展示不同数据集上的方法差异。

### 16.5 robustness 图

展示不同随机种子、不同平台下的波动。

### 16.6 case study 可视化

选择 1–2 个代表数据集画：

- integration embedding
- spatial domains
- alignment before/after
- 错误案例

------

## 17. 需要优先解决的关键实现问题

### 17.1 feature selection 的空间化定义

必须明确：

- SVG 如何定义和实现
- cross-slice conserved genes 如何定义
- 多切片时选特征是按切片还是全局

### 17.2 不同方法对输入格式的兼容

很多空间 integration 方法输入格式不一致，需要 wrapper 层统一。

### 17.3 空间指标的选取

不能只用 scRNA 的 integration 指标，否则会丢失空间问题本质。

### 17.4 query mapping 是否进入第一版

建议第一版先不做太复杂的 query mapping，先把主 benchmark 跑通。
 如果资源充足，再加入 spatial query 子任务。

------

## 18. 实施顺序建议

### 阶段 1：搭框架

目标：把 benchmark 主干跑通。

任务：

1. 确定 repo 结构
2. 实现统一数据对象
3. 实现 2–3 个 feature selection
4. 实现 1–2 个 integration wrapper
5. 实现 2 个任务和基本评估
6. 实现 runner 和结果保存

### 阶段 2：做 MVP

目标：得到首批可解释结果。

任务：

1. 选 2–3 个数据集
2. 跑 3–4 个 feature selection
3. 跑 2–3 个 integration 方法
4. 生成首批 heatmap / ranking
5. 检查结论是否合理

### 阶段 3：扩展任务与方法

目标：补充 benchmark 完整性。

任务：

1. 增加 alignment / slice representation
2. 增加 conserved genes 等 feature selection
3. 增加平台和数据集覆盖
4. 做 feature number ablation
5. 做 robustness analysis

### 阶段 4：论文化整理

目标：形成文章叙事。

任务：

1. 提炼主要发现
2. 写推荐结论
3. 生成标准图表
4. 补案例分析
5. 写方法和补充材料