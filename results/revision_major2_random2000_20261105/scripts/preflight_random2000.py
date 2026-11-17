from pathlib import Path
import json,csv,time,sys
from spatial_fs_benchmark.config import load_benchmark_config
from spatial_fs_benchmark.data.io import load_dataset
from spatial_fs_benchmark.data.preprocess import preprocess_dataset
from spatial_fs_benchmark.feature_selection.random_baseline import RandomSelector
root=Path("/home/xd/y/Benchmark");rev=root/"results/revision_major2_random2000_20261105"
manifest=json.loads((rev/"config_manifest.json").read_text());rows=[]
for item in manifest:
 t=time.monotonic();cfg=load_benchmark_config(root/item["source_config"]);ds=preprocess_dataset(load_dataset(cfg.datasets[0]),cfg.datasets[0].preprocess)
 for seed in cfg.seeds:
  old=Path(cfg.output_dir)/ds.name.lower()/"feature_selection"/"random"/"n500"/("seed"+str(seed))/"selected_features.json"
  prior=json.loads(old.read_text());same=RandomSelector().select(ds,500,random_seed=seed)
  assert same.feature_names==prior["feature_names"],f"random seed reconstruction mismatch: {ds.name}, seed {seed}"
  newer=RandomSelector().select(ds,2000,random_seed=seed);assert len(newer.feature_names)==2000 and newer.feature_names[:500]==same.feature_names
  rows.append(dict(dataset=ds.name,seed=seed,old_500_matches_original=True,new_2000_starts_with_old_500=True,n_2000=len(newer.feature_names),source_feature_hash=prior["metadata"]["feature_name_hash"],new_feature_hash=newer.metadata["feature_name_hash"]))
 print(ds.name,"3 seeds verified, seconds",round(time.monotonic()-t,1),flush=True);del ds
with (rev/"random_selection_preflight.tsv").open("w",newline="",encoding="utf8") as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t");w.writeheader();w.writerows(rows)
print("PASS: all",len(rows),"random selections matched the stored original",flush=True)
