from pathlib import Path
import argparse,json,sys,csv
from spatial_fs_benchmark.config import load_benchmark_config
root=Path("/home/xd/y/Benchmark"); rev=root/"results/revision_major2_random2000_20261105"
parser=argparse.ArgumentParser();parser.add_argument("name");args=parser.parse_args();items=json.loads((rev/"config_manifest.json").read_text());item=next(x for x in items if x["name"]==args.name)
old=load_benchmark_config(root/item["source_config"]);new=load_benchmark_config(root/item["generated_config"]);dataset=old.datasets[0].name.lower();checked=0
for seed in old.seeds:
 oldsel=Path(old.output_dir)/dataset/"feature_selection"/"random"/"n500"/f"seed{seed}"/"selected_features.json"
 newsel=Path(new.output_dir)/dataset/"feature_selection"/"random"/"n2000"/f"seed{seed}"/"selected_features.json"
 a=json.loads(oldsel.read_text());b=json.loads(newsel.read_text());assert len(b["feature_names"])==2000 and b["feature_names"][:500]==a["feature_names"],(dataset,seed,"selection")
 for integrator in old.integration_methods:
  ap=Path(old.output_dir)/dataset/integrator.name/"random"/"n500"/f"seed{seed}";bp=Path(new.output_dir)/dataset/integrator.name/"random"/"n2000"/f"seed{seed}"
  am=json.loads((ap/"embedding.meta.json").read_text());bm=json.loads((bp/"embedding.meta.json").read_text());assert (bp/"embedding.npz").is_file()
  for key in ["dataset","integrator_name","integrator_params","integrator_implementation_version","selector_name","selector_params","seed"]:assert am.get(key)==bm.get(key),(dataset,seed,integrator.name,key,am.get(key),bm.get(key))
  assert am["n_features"]==500 and bm["n_features"]==2000 and bm["selection_feature_hash"]==b["metadata"]["feature_name_hash"]
  for task in old.tasks:
   ar=json.loads((ap/(task.name+"_records.json")).read_text());br=json.loads((bp/(task.name+"_records.json")).read_text());at=json.loads((ap/(task.name+"_records.meta.json")).read_text());bt=json.loads((bp/(task.name+"_records.meta.json")).read_text())
   for key in ["dataset","integrator_name","integrator_params","task_name","task_params","task_implementation_version","selector_name","selector_params","seed"]:
    if key=="task_implementation_version" and at.get(key) is None and task.name=="clustering_eval":
     assert bt.get(key)=="v2_kmeans_slice_balanced_spatial_metrics"
     evidence=[]
     for f in rev.glob("original_clustering_audit_*.tsv"):
      with f.open(encoding="utf8",newline="") as handle:evidence.extend(csv.DictReader(handle,delimiter="\t"))
     match=[row for row in evidence if row["dataset"].lower()==dataset and row["integrator"]==integrator.name]
     assert len(match)==5 and all(row["matches"]=="True" for row in match),(dataset,integrator.name,"missing independent original clustering audit")
    else:assert at.get(key)==bt.get(key),(dataset,seed,integrator.name,task.name,key)
   assert {x["metric_name"] for x in ar}=={x["metric_name"] for x in br},(dataset,seed,task.name,"metric keys")
   assert all(x["n_features"]==2000 and x["effective_n_features"]==2000 and x["random_seed"]==seed for x in br)
   checked+=len(br)
 print(dataset,"seed",seed,"verified",flush=True)
print("PASS",args.name,"metric records",checked,flush=True)
