from pathlib import Path
import argparse,json,csv,math,time
from spatial_fs_benchmark.config import load_benchmark_config
from spatial_fs_benchmark.data.io import load_dataset
from spatial_fs_benchmark.data.preprocess import preprocess_dataset
from spatial_fs_benchmark.benchmark.experiment import read_feature_selection,read_integration_result
from spatial_fs_benchmark.tasks import build_task
p=argparse.ArgumentParser();p.add_argument("datasets",nargs="+");args=p.parse_args()
root=Path("/home/xd/y/Benchmark");rev=root/"results/revision_major2_random2000_20261105";manifest=json.loads((rev/"config_manifest.json").read_text());items=[x for x in manifest if x["name"] in args.datasets];assert len(items)==len(args.datasets)
rows=[]
for x in items:
 cfg=load_benchmark_config(root/x["source_config"]); ds=preprocess_dataset(load_dataset(cfg.datasets[0]),cfg.datasets[0].preprocess);base=Path(cfg.output_dir)/ds.name.lower();t=[t for t in cfg.tasks if t.name=="clustering_eval"][0];task=build_task(t.name,**t.params);sel=read_feature_selection(base/"feature_selection/random/n500/seed0/selected_features.json")
 for integ in cfg.integration_methods:
  before=time.monotonic();run=base/integ.name/"random/n500/seed0";emb=read_integration_result(run/"embedding.npz",integ.name);old=json.loads((run/"clustering_eval_records.json").read_text());orig={a["metric_name"]:a["metric_value"] for a in old};new=task.evaluate(ds,sel,emb,random_seed=0).metrics;assert set(orig)==set(new),(x["name"],integ.name,orig.keys(),new.keys())
  for name,val in new.items():
   a=orig[name];b=float(val);equal=(a is None and math.isnan(b)) or (a is not None and math.isclose(float(a),b,rel_tol=1e-7,abs_tol=1e-7));rows.append(dict(dataset=ds.name,integrator=integ.name,metric=name,original=a,recomputed=b,abs_difference=abs(float(a)-b) if a is not None else "",matches=equal))
  print(ds.name,integ.name,"matched",sum(z["matches"] for z in rows[-len(new):]),"of",len(new),"time",round(time.monotonic()-before,1),flush=True)
 del ds
out=rev/("original_clustering_audit_"+"_".join(args.datasets)+".tsv")
with out.open("w",encoding="utf8",newline="") as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t");w.writeheader();w.writerows(rows)
print("TOTAL",len(rows),"matching",sum(x["matches"] for x in rows),flush=True)
if not all(x["matches"] for x in rows):
 for row in rows:
  if not row["matches"]:print("DIFF",row,flush=True)
 raise SystemExit(2)
