# Analysis workflows

Use the entry points below to locate each analysis. The primary ranking uses five Visium datasets; embryo and additional-platform analyses are reported separately.

| Analysis | Entry points |
| --- | --- |
| Five-dataset ranking and separate embryos | `revision_r1_major1_five_20261105/` |
| Random-2,000 | `results/revision_major2_random2000_20261105/` |
| Unranked-gene sensitivity | `results/revision_major3_missing_ranks_20261105/scripts/` |
| External markers and section-held-out evaluation | `results/revision_major4_independent_validation_20261106/scripts/` |
| Individual-held-out validation | `results/revision_r2_major2_20261106/scripts/` |
| Actual-count budget comparisons | `results/revision_r2_major6_20261107/code/` |
| Required-gene retention | `results/revision_r2_major4_20261107/code/` |
| Adult mouse-brain platform extension | `results/revision_r3_comment5_20261107/scripts/` |
| Paired spatial aggregation | `results/revision_r2_m1_scale_match_20261111/` |
| Final figure rendering | `publication/final_main/`, `publication/final_supplement/` |

See [CODE_INDEX.tsv](../CODE_INDEX.tsv) for the full index. Download the required inputs, configure local paths and adjust device/thread settings before execution. Run experiment steps before their summary and plotting steps.

Directory names retain the identifiers used by the execution scripts so that relative paths remain consistent.
