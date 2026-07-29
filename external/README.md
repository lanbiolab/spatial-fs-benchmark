# External Code And Provenance

The public release contains only the third-party files required by the tracked
workflow:

- `atlas-feature-selection-benchmark/` provides the plotting helper reused by
  the reference-aligned R figures. Its upstream README and license are retained.
- `iSTBench/Benchmark/RunModel/GraphST/GraphST/` contains the minimal GraphST
  source files used by the integration wrapper. The iSTBench GPL-3.0 license and
  README are retained at `iSTBench/`.
- `svg-method-references/PROVENANCE.tsv` records the official SVG software and
  versions used by the wrappers.
- `integration-method-references/PROVENANCE.tsv` records the verified upstream
  GraphST commit, original file hashes, and the sparse adaptations used here.

Third-party code remains subject to its upstream license. The project-level MIT
license applies only to original benchmark code.
