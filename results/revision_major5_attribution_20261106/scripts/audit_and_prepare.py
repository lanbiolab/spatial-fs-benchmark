import ast
import datetime
import difflib
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

ROOT = Path('/home/xd/y/Benchmark')
RELEASE = ROOT / '_github_release/spatial-fs-benchmark'
OUT = ROOT / 'results/revision_major5_attribution_20261106'
UPSTREAM = OUT / 'upstream'
CANDIDATE = OUT / 'release_candidate'
GRAPHST = Path('external/iSTBench/Benchmark/RunModel/GraphST/GraphST')
ATLAS = Path('external/atlas-feature-selection-benchmark')
DATE = '2026-11-06'

def digest(b):
    return hashlib.sha256(b).hexdigest()

def git(repo, *args):
    return subprocess.check_output(['git', '-C', str(repo), *args])

OUT.mkdir(exist_ok=True)
assert git(RELEASE, 'status', '--porcelain') == b'', 'Release checkout is not clean'
head = git(RELEASE, 'rev-parse', 'HEAD').decode().strip()
atlas_head = git(ROOT / ATLAS, 'rev-parse', 'HEAD').decode().strip()
ist_head = git(ROOT / 'external/iSTBench', 'rev-parse', 'HEAD').decode().strip()
if not CANDIDATE.exists():
    subprocess.run(['git', '-C', str(RELEASE), 'worktree', 'add', '--detach', str(CANDIDATE), head], check=True)
assert git(CANDIDATE, 'rev-parse', 'HEAD').decode().strip() == head

rows = []
for name in ['GraphST.py', 'model.py', 'preprocess.py', 'utils.py', '__init__.py']:
    official = (UPSTREAM / 'GraphST' / name).read_bytes()
    baseline = git(ROOT / 'external/iSTBench', 'show', 'HEAD:Benchmark/RunModel/GraphST/GraphST/' + name)
    running = (ROOT / GRAPHST / name).read_bytes()
    released = (RELEASE / GRAPHST / name).read_bytes()
    assert official == baseline, name
    assert running == released, f'Release/runtime differ: {name}'
    before = official.decode().replace('\r\n', '\n').replace('\r', '\n')
    after = running.decode().replace('\r\n', '\n').replace('\r', '\n')
    (OUT / f'GraphST_{name}.diff').write_text(''.join(difflib.unified_diff(before.splitlines(True), after.splitlines(True), fromfile='official/' + name, tofile='benchmark/' + name)))
    rows.append({'file': name, 'official_sha256': digest(official), 'istbench_commit_sha256': digest(baseline),
                 'runtime_sha256': digest(running), 'release_sha256': digest(released),
                 'official_equals_istbench_before_adaptation': official == baseline,
                 'runtime_equals_release': running == released, 'local_adaptation': official != running})
import pandas as pd
pd.DataFrame(rows).to_csv(OUT / 'GraphST_source_hash_audit.tsv', sep='\t', index=False)
theme = Path('analysis/R/plotting.R')
theme_original = git(ROOT / ATLAS, 'show', 'HEAD:' + str(theme))
theme_current = (ROOT / ATLAS / theme).read_bytes()
theme_released = (RELEASE / ATLAS / theme).read_bytes()
(OUT / 'atlas_plotting_modifications.diff').write_bytes(git(ROOT / ATLAS, 'diff', '--', str(theme)))

license_text = (UPSTREAM / 'LICENSE.md').read_bytes()
setup_text = (UPSTREAM / 'setup.py').read_text()
assert b'GNU AFFERO GENERAL PUBLIC LICENSE' in license_text and '__license__ = "MIT"' in setup_text
mit = (RELEASE / 'LICENSE').read_bytes()
gpl = (RELEASE / 'external/iSTBench/LICENSE').read_bytes()
assert b'MIT License' in mit and b'GNU GENERAL PUBLIC LICENSE' in gpl
(CANDIDATE / 'LICENSES').mkdir(exist_ok=True)
(CANDIDATE / 'LICENSES/MIT-original-benchmark.txt').write_bytes(mit)
(CANDIDATE / 'LICENSE').write_bytes(license_text)
(CANDIDATE / GRAPHST / 'LICENSE.md').write_bytes(license_text)

for filename in ['GraphST.py', 'model.py']:
    original = (RELEASE / GRAPHST / filename).read_bytes()
    header = f'''# GraphST by Yahui Long and colleagues; obtained through iSTBench.
# Upstream GraphST commit: d62b0b7b6cd38ee285f3ac8cd67b7341a10bcc74.
# Benchmark adaptation: sparse graph/neighborhood operations; see MODIFICATIONS.md.
# Modification notice added {DATE}; computational body retained verbatim.
# This candidate follows the upstream AGPL-3.0 license text in LICENSE.md.
'''.encode()
    (CANDIDATE / GRAPHST / filename).write_bytes(header + original)
    assert ast.dump(ast.parse((header + original).decode())) == ast.dump(ast.parse(original.decode()))
theme_header = f'''# Adapted from Zappia et al. (2025), theislab/atlas-feature-selection-benchmark.
# Upstream version: {atlas_head}; upstream MIT notice is retained.
# Local font/export and ggplot2 compatibility changes documented in THIRD_PARTY_NOTICES.md.
# Attribution notice added {DATE}; existing plotting body retained verbatim.
'''.encode()
(CANDIDATE / ATLAS / theme).write_bytes(theme_header + theme_current)

for filename in ['pyproject.toml', 'CITATION.cff']:
    text = (RELEASE / filename).read_text()
    if filename == 'pyproject.toml':
        assert text.count('license = { text = "MIT" }') == 1
        text = text.replace('license = { text = "MIT" }', 'license = { text = "AGPL-3.0-only" }')
    else:
        assert text.count('license: MIT') == 1
        text = text.replace('license: MIT', 'license: AGPL-3.0-only')
    (CANDIDATE / filename).write_text(text)

notices = f'''# Third-party code, attribution and distribution terms

Revision candidate prepared {DATE}, based on release commit `{head}`.
This is a local, unpublished release candidate. No experiment output or numerical
implementation was changed to prepare this candidate.

## Distribution scope

The combined source distribution is prepared on an AGPL-3.0 basis, following the
explicit GraphST upstream LICENSE.md. The top-level LICENSE contains that complete
text. This is a conservative release decision for this candidate, not a claim that
the conflicting upstream package metadata has been clarified by its authors.

Previously MIT-licensed original benchmark contributions retain their existing MIT
grant and notice in LICENSES/MIT-original-benchmark.txt. Earlier MIT grants are not
revoked. The original MIT permission does not purport to cover third-party code.
Third-party copyrights, source-specific licenses and notices remain applicable.
The retained iSTBench GPL-3.0 terms also apply to its covered portions; those files
are not relicensed as MIT or stripped of their existing terms. GPLv3/AGPLv3 section
13 describes the special terms for a combined work. Complete source and the local
modifications are included in this source candidate.

## Zappia et al.: feature-selection benchmark infrastructure

Reference: Zappia L et al. Feature selection methods affect the performance of
scRNA-seq data integration and querying. Nature Methods 22, 834-844 (2025).
https://doi.org/10.1038/s41592-025-02624-3

Repository: https://github.com/theislab/atlas-feature-selection-benchmark
Audited version: `{atlas_head}`. Upstream license: MIT;
Copyright (c) 2022 Theis Lab. The complete upstream license is retained at
external/atlas-feature-selection-benchmark/LICENSE.

The R plotting helper is directly sourced by the benchmark figure scripts.
theme_features(), theme_features_pub() and the export helper are reused/adapted.
Local changes specify Arial-compatible typography, font/export handling and the
ggplot2 linewidth interface. The existing release checkout still contains the
unadapted helper; this candidate synchronizes it to the exact helper used by the
experiments/figures and adds an attribution comment. No figure was regenerated or
restyled to prepare this candidate.

Observed-range metric scaling and grouped aggregation follow the upstream
summarisation.R framework, with a separate Python implementation in
scripts/build_frozen_scores.py and spatial R aggregation. Ranking/heatmap,
feature-budget, and method-versus-metric dotplot designs are adapted in the
figure-specific benchmark scripts. This is distinct from claiming that every
upstream plotting script is called unchanged. The spatial task definitions and
frozen ranges/weighting are documented separately.

## iSTBench and GraphST implementation provenance

iSTBench: https://github.com/bm2-lab/iSTBench , audited commit `{ist_head}`;
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
'''
(CANDIDATE / 'THIRD_PARTY_NOTICES.md').write_text(notices)
modifications = f'''# GraphST source and benchmark modifications

Original authors: Yahui Long and colleagues. Original repository and commit:
https://github.com/JinmiaoChenLab/GraphST/tree/d62b0b7b6cd38ee285f3ac8cd67b7341a10bcc74
The benchmark obtained these files through iSTBench commit `{ist_head}`.
The pre-adaptation copies of all five Python files match the official originals.

Local sparse adaptations were recorded in the benchmark provenance on 2026-07-29:
- GraphST.py retains sparse adjacency and neighborhood matrices rather than
  converting them to dense tensors; its training input is handled accordingly.
- model.py uses sparse matrix multiplication for the neighborhood readout.
- The benchmark wrapper builds block-diagonal within-slice graphs and reduces
  reconstructed expression to 20 PCs for its common downstream evaluation.
- The original loss, optimizer and 600-epoch setting are retained.

These changes do not constitute a new GraphST algorithm. Original upstream
hashes and the local adaptation are recorded in the benchmark's
external/integration-method-references/PROVENANCE.tsv.

Revision notices were added on {DATE}. The explicit upstream AGPL-3.0 text is
retained as LICENSE.md; setup.py's inconsistent MIT metadata is recorded in the
audit rather than silently substituted for that text. This revision changes
source comments and licensing documentation, not the computational statements.
'''
(CANDIDATE / GRAPHST / 'MODIFICATIONS.md').write_text(modifications)

readme = (RELEASE / 'README.md').read_text()
readme += '''
## Attribution and source licensing (revision candidate)

The benchmark builds on Zappia et al. (Nature Methods 22, 834-844, 2025;
https://doi.org/10.1038/s41592-025-02624-3) for the observed-range scaling/grouped
summarisation framework and reused/adapted plotting infrastructure. The GraphST
implementation is obtained through iSTBench and traced to the original GraphST
source; both sources are credited in THIRD_PARTY_NOTICES.md.

This combined source release candidate is prepared under AGPL-3.0, with the
upstream GPL-3.0 and MIT component notices retained. Previously granted MIT
permissions for original benchmark contributions are preserved separately.
See LICENSE, LICENSES/MIT-original-benchmark.txt and THIRD_PARTY_NOTICES.md for
the precise scope, including conflicting upstream GraphST license metadata.
The candidate has not been pushed or deposited.
'''
(CANDIDATE / 'README.md').write_text(readme)
(CANDIDATE / 'external/README.md').write_text('''# External code and provenance

Third-party source and license records are described in ../THIRD_PARTY_NOTICES.md.
The plotting helper from atlas-feature-selection-benchmark retains its MIT
license. The iSTBench GPL-3.0 license is retained, and the vendored GraphST source
also includes the explicit upstream AGPL-3.0 LICENSE.md and a modification notice.
The combined source release candidate is not presented as MIT-only.

See svg-method-references/PROVENANCE.tsv and
integration-method-references/PROVENANCE.tsv for method versions and source hashes.
''')

audit = {'release_base_commit': head, 'candidate_path': str(CANDIDATE), 'candidate_published': False,
         'atlas_upstream_commit': atlas_head, 'istbench_upstream_commit': ist_head,
         'istbench_publication_doi': '10.1186/s13059-025-03796-z',
         'istbench_publication': 'Dong et al., Genome Biology 26, 318 (2025)',
         'distinct_benchmarkst_publication_doi': '10.1186/s13059-024-03361-0',
         'graphst_upstream_commit': 'd62b0b7b6cd38ee285f3ac8cd67b7341a10bcc74',
         'graphst_explicit_license': 'AGPL-3.0', 'graphst_setup_metadata': 'MIT',
         'graphst_license_metadata_conflict': True, 'istbench_license': 'GPL-3.0', 'atlas_license': 'MIT',
         'candidate_combined_distribution_basis': 'AGPL-3.0-only',
         'prior_original_benchmark_MIT_grant_retained': True,
         'graphst_files_checked': rows, 'atlas_original_plotting_sha256': digest(theme_original),
         'atlas_adapted_plotting_sha256': digest(theme_current),
         'atlas_release_plotting_sha256': digest(theme_released),
         'atlas_release_matches_upstream': theme_released == theme_original,
         'atlas_release_matches_runtime_before_candidate': theme_released == theme_current,
         'upstream_license_sha256': digest(license_text), 'upstream_setup_sha256': digest((UPSTREAM / 'setup.py').read_bytes())}
(CANDIDATE / 'docs/revision_major5_source_audit.json').write_text(json.dumps(audit, indent=2))
(OUT / 'source_audit.json').write_text(json.dumps(audit, indent=2))

tracked = git(RELEASE, 'ls-files').decode().splitlines()
py_count = 0
for name in tracked:
    original = RELEASE / name
    candidate = CANDIDATE / name
    if name.endswith('.py'):
        assert ast.dump(ast.parse(candidate.read_text())) == ast.dump(ast.parse(original.read_text())), name
        py_count += 1
    elif name.endswith(('.R', '.yaml', '.yml', '.tsv', '.csv', '.pdf', '.png', '.svg')):
        if name == str(ATLAS / theme):
            assert candidate.read_bytes() == theme_header + theme_current
        elif name != 'CITATION.cff':
            assert candidate.read_bytes() == original.read_bytes(), name
assert git(RELEASE, 'status', '--porcelain') == b'', 'Original release checkout changed'
for p in git(CANDIDATE, 'ls-files', '--others', '--exclude-standard').decode().splitlines():
    subprocess.run(['git', '-C', str(CANDIDATE), 'add', '--intent-to-add', p], check=True)
(OUT / 'release_candidate.patch').write_bytes(git(CANDIDATE, 'diff', '--binary', 'HEAD'))
(OUT / 'candidate_diffstat.txt').write_bytes(git(CANDIDATE, 'diff', '--stat', 'HEAD'))
verification = {'all_checks_passed': True, 'python_files_AST_unchanged': py_count,
                'original_release_checkout_clean': True, 'original_data_configs_figures_unchanged': True,
                'plotting_helper_synced_to_actual_experiment_version': True,
                'candidate_only': True, 'published': False,
                'upstream_graphst_files_match_istbench_committed_baseline': 5,
                'runtime_graphst_files_match_release': 5}
(OUT / 'candidate_verification.json').write_text(json.dumps(verification, indent=2))
print(json.dumps(verification, indent=2))
print((OUT / 'candidate_diffstat.txt').read_text())
