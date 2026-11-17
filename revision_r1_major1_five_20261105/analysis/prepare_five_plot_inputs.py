from pathlib import Path
import shutil
r=Path("/home/xd/y/Benchmark/revision_r1_major1_five_20261105"); root=Path("/home/xd/y/Benchmark")
p=(root/"scripts/generate_supplemental_robustness_figure.R").read_text()
p=p.replace('    "integration_weight_0.6" = "Integration\\nweight 0.6",\n    "labelled_datasets_only" = "Annotated\\ndatasets"','    "integration_weight_0.6" = "Integration\\nweight 0.6"')
p=p.replace('filter(.data$NScenarios == 6)','filter(.data$NScenarios == 5)')
assert "labelled_datasets_only" not in p and "NScenarios == 5" in p
(r/"analysis/plot_s2_five.R").write_text(p)
d=r/"analysis/robustness_inputs";d.mkdir(exist_ok=True)
for src,dst in [("five_dataset_sensitivity_ranks.tsv","scaling_weight_sensitivity_ranks.tsv"),("five_dataset_bootstrap_ranks.tsv","dataset_bootstrap_rank_intervals.tsv"),("five_dataset_sensitivity_concordance.tsv","scaling_weight_concordance.tsv")]:shutil.copy2(r/"analysis"/src,d/dst)
