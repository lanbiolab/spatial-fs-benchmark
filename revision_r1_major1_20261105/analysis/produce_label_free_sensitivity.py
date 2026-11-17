from pathlib import Path
import pandas as pd
root=Path("/home/xd/y/Benchmark"); out=root/"revision_r1_major1_20261105/analysis"
src=root/"results/spatial_svg_rebuild_v1/frozen_scores"
keep=["MouseBrainSerialSections","STOmics0224","E8p5Embryo","E9p5Embryo"]
metrics=["iLISI","CHAOS","PAS","silhouette","Ratio"]
settings=pd.read_csv(src/"representative_settings.tsv",sep="\t",dtype={"n_features":str})
settings=settings[settings.dataset.isin(keep) & (settings.fs_method!="wilcoxon")]
d=pd.read_csv(src/"setting_metric_seed_summary.tsv",sep="\t",dtype={"n_features":str})
d=d.merge(settings[["dataset","fs_method","n_features","MethodGroup"]],on=["dataset","fs_method","n_features"],how="inner",validate="many_to_one")
d=d[d.metric_name.isin(metrics)].copy()
assert set(d.dataset)==set(keep)
d.sort_values(["dataset","integration_method","metric_name","fs_method"]).to_csv(out/"four_dataset_label_free_sensitivity.tsv",sep="\t",index=False,na_rep="NA")
print("Label-free metric rows",len(d));print(d.groupby(["dataset","metric_name"]).size().to_string())
