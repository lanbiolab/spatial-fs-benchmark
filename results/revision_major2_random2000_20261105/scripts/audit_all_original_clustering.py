from pathlib import Path
import json,csv,math,time
from spatial_fs_benchmark.config import load_benchmark_config
from spatial_fs_benchmark.data.io import load_dataset
from spatial_fs_benchmark.data.preprocess import preprocess_dataset
from spatial_fs_benchmark.benchmark.experiment import read_feature_selection,read_integration_result
from spatial_fs_benchmark.tasks import build_task
root=Path("/home/xd/y/Benchmark");rev=root/"results/revision_major2_random2000_20261105";manifest=json.loads((rev/"config_manifest.json").read_text());rows=[]
for x in manifest:
 if "subset" in x["name"]:continue
 cfg=load_benchmark_config(root/x["source_config"]); ds=preprocess_dataset(load_dataset(cfg.datasets[0]),cfg.datasets[0].preprocess);base=Path(cfg.output_dir)/ds.name.lower();t=next(t for t in cfg.tasks if t.name=="clustering_eval");task=build_task(t.name,**t.params)
 for seed in cfg.seeds:
  sel=read_feature_selection(base/"feature_selection/random/n500"/f"seed{seed}"/"selected_features.json")
  for integ in cfg.integration_methods:
   run=base/integ.name/"random/n500"/f"seed{seed}";emb=read_integration_result(run/"embedding.npz",integ.name);orig={a["metric_name"]:a["metric_value"] for a in json.loads((run/"clustering_eval_records.json").read_text())};new=task.evaluate(ds,sel,emb,random_seed=seed).metrics;assert set(orig)==set(new)
   for name,val in new.items():
    a=orig[name];b=float(val);equal=(a is None and math.isnan(b)) or (a is not None and math.isclose(float(a),b,rel_tol=1e-7,abs_tol=1e-7));rows.append(dict(dataset=ds.name,integrator=integ.name,seed=seed,metric=name,original=a,recomputed=b,abs_difference=abs(float(a)-b) if a is not None else "",matches=equal))
   print(ds.name,integ.name,"seed",seed,"matched",sum(z["matches"] for z in rows[-len(new):]),"of",len(new),flush=True)
 del ds
out=rev/"legacy_clustering_all_seed_check.tsv"
with out.open("w",encoding="utf8",newline="") as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t");w.writeheader();w.writerows(rows)
print("TOTAL",len(rows),"matching",sum(x["matches"] for x in rows),flush=True)
if not all(x["matches"] for x in rows):
 for row in rows:
  if not row["matches"]:print("DIFF",row,flush=True)
 raise SystemExit(2)
