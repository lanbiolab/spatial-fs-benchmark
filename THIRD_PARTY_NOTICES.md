# Third-party code, attribution and distribution terms

Attribution audit prepared 2026-11-06; incorporated into the 2026-11-17 repository revision. Audit base commit `98feb759a666d02828a0a781354aa1238c0b7b23`.
The original audit and patch are retained under provenance/R1_M5/. Applying these
attribution and licensing changes does not change numerical implementations.

## Distribution scope

The combined source distribution is prepared on an AGPL-3.0 basis, following the
explicit GraphST upstream LICENSE.md. The top-level LICENSE contains that complete
text. This distribution follows that explicit license text, not a claim that
the conflicting upstream package metadata has been clarified by its authors.

Previously MIT-licensed original benchmark contributions retain their existing MIT
grant and notice in LICENSES/MIT-original-benchmark.txt. Earlier MIT grants are not
revoked. The original MIT permission does not purport to cover third-party code.
Third-party copyrights, source-specific licenses and notices remain applicable.
The retained iSTBench GPL-3.0 terms also apply to its covered portions; those files
are not relicensed as MIT or stripped of their existing terms. GPLv3/AGPLv3 section
13 describes the special terms for a combined work. Complete source and the local
modifications are included in this source distribution.

## Zappia et al.: feature-selection benchmark infrastructure

Reference: Zappia L et al. Feature selection methods affect the performance of
scRNA-seq data integration and querying. Nature Methods 22, 834-844 (2025).
https://doi.org/10.1038/s41592-025-02624-3

Repository: https://github.com/theislab/atlas-feature-selection-benchmark
Audited version: `b89fc0f66747062e6e1b4b35bd392b27ad035295`. Upstream license: MIT;
Copyright (c) 2022 Theis Lab. The complete upstream license is retained at
external/atlas-feature-selection-benchmark/LICENSE.

The R plotting helper is directly sourced by the benchmark figure scripts.
theme_features(), theme_features_pub() and the export helper are reused/adapted.
Local changes specify Arial-compatible typography, font/export handling and the
ggplot2 linewidth interface. The plotting helper is synchronized with the version used by the
experiments/figures and carries an attribution comment. Historical hashes and the
modification diff are retained in provenance/R1_M5/.

Observed-range metric scaling and grouped aggregation follow the upstream
summarisation.R framework, with a separate Python implementation in
scripts/build_frozen_scores.py and spatial R aggregation. Ranking/heatmap,
feature-budget, and method-versus-metric dotplot designs are adapted in the
figure-specific benchmark scripts. This is distinct from claiming that every
upstream plotting script is called unchanged. The spatial task definitions and
frozen ranges/weighting are documented separately.

## iSTBench and GraphST implementation provenance

iSTBench: https://github.com/bm2-lab/iSTBench , audited commit `e51f60e6127dbbdd953a6fb43351b8864d9e646a`;
its complete GPL-3.0 text and README are retained at external/iSTBench/.
Associated study: Dong K et al. Benchmarking multi-slice integration and downstream
applications in spatial transcriptomics data analysis. Genome Biology 26, 318
(2025), https://doi.org/10.1186/s13059-025-03796-z .
This is distinct from Hu et al.'s BenchmarkST study (Genome Biology 25, 212, 2024;
https://doi.org/10.1186/s13059-024-03361-0), which is related prior work but is not
the source of the bm2-lab/iSTBench implementation.

The GraphST wrapper imports the vendored GraphST package through
external/iSTBench/Benchmark/RunModel/GraphST/ . Attribution therefore includes
both the intermediary implementation source iSTBench and the original method:
Long Y et al. Nature Communications 14, 1155 (2023),
https://doi.org/10.1038/s41467-023-36796-3 .

Original GraphST repository: https://github.com/JinmiaoChenLab/GraphST , commit
`d62b0b7b6cd38ee285f3ac8cd67b7341a10bcc74`.
All five Python source files match the corresponding iSTBench committed originals
before the documented benchmark adaptations. LICENSE.md at that exact GraphST
commit contains AGPL-3.0, whereas setup.py says MIT. Both facts are recorded in
docs/revision_major5_source_audit.json. This candidate retains the actual AGPL
text beside the GraphST source and follows that text rather than treating the
setup.py field as authority to redistribute these files under MIT.

GraphST.py and model.py have documented local sparse adaptations. preprocess.py,
utils.py and __init__.py are unchanged from the audited upstream files. The sparse
adaptation is described in the adjacent MODIFICATIONS.md. No iSTBench or GraphST
copyright is assigned to the benchmark authors by these notices.

## Other method software

External installed software and additional method resources retain their own
licenses. See external/svg-method-references/PROVENANCE.tsv and the environment
files. This notice does not replace the license of an independently installed
dependency, and it does not assert ownership of data or upstream publications.

License texts and source records:
- GraphST original LICENSE.md: https://raw.githubusercontent.com/JinmiaoChenLab/GraphST/d62b0b7b6cd38ee285f3ac8cd67b7341a10bcc74/LICENSE.md
- GraphST setup.py: https://raw.githubusercontent.com/JinmiaoChenLab/GraphST/d62b0b7b6cd38ee285f3ac8cd67b7341a10bcc74/setup.py
- GPL/AGPL combined works: https://www.gnu.org/licenses/gpl-faq.en.html#GPLv3WithAGPLv3
