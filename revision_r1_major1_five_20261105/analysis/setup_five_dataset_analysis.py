from pathlib import Path
import pandas as pd
r=Path("/home/xd/y/Benchmark"); rev=r/"revision_r1_major1_five_20261105"
out=rev/"frozen_five";out.mkdir(parents=True,exist_ok=True);(rev/"analysis").mkdir(exist_ok=True)
keep=["DLPFC","MouseBrainSerialSections","STOmics0212","STOmics0218","STOmics0224"]
src=r/"results/spatial_svg_rebuild_v1/frozen_scores"
for name in ["representative_settings.tsv","representative_task_scores.tsv","setting_task_scores.tsv","scaled_setting_metrics.tsv","frozen_metric_ranges.tsv","setting_metric_seed_summary.tsv"]:
 d=pd.read_csv(src/name,sep="\t",dtype={"n_features":str});d=d[d.dataset.isin(keep)];assert set(d.dataset)==set(keep),(name,set(d.dataset));d.to_csv(out/name,sep="\t",index=False,na_rep="NA")
a=pd.read_csv(out/"representative_task_scores.tsv",sep="\t",dtype={"n_features":str});g=pd.read_csv(r/"results/graphst_canonical_v1/summary/graphst_setting_task_scores.tsv",sep="\t",dtype={"n_features":str});g=g[g.dataset.isin(keep)];g.to_csv(out/"graphst_setting_task_scores.tsv",sep="\t",index=False,na_rep="NA")
combined=pd.concat([a,g],ignore_index=True);combined.to_csv(out/"combined_integrator_scores.tsv",sep="\t",index=False,na_rep="NA");combined[combined.MethodGroup.isin(["expression_driven","spatially_informed"])].to_csv(out/"competitive_integrator_scores.tsv",sep="\t",index=False,na_rep="NA")
print("5 dataset source files ready",len(a),len(g))
