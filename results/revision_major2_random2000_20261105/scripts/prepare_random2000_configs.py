from pathlib import Path
import copy,hashlib,json,yaml
root=Path("/home/xd/y/Benchmark");revision=root/"results/revision_major2_random2000_20261105"
configdir=revision/"configs";configdir.mkdir(parents=True,exist_ok=True)
items=[("dlpfc", "configs/rebuild_v1/non_spatial_downstream/canonical/dlpfc.yaml"),("mouse_brain_serial_sections","configs/rebuild_v1/non_spatial_downstream/canonical/mouse_brain_serial_sections.yaml"),("stomics_0212","configs/rebuild_v1/non_spatial_downstream/canonical/stomics_0212.yaml"),("stomics_0218","configs/rebuild_v1/non_spatial_downstream/canonical/stomics_0218.yaml"),("stomics_0224","configs/rebuild_v1/non_spatial_downstream/canonical/stomics_0224.yaml"),("stomics_0212_immune_subset","configs/rebuild_v1/biological_subsets/non_spatial/stomics_0212_immune_subset.yaml"),("stomics_0212_epithelial_subset","configs/rebuild_v1/biological_subsets/non_spatial/stomics_0212_epithelial_subset.yaml")]
manifest=[]
for name,relative in items:
 source=root/relative;orig=yaml.safe_load(source.read_text());random=[m for m in orig["feature_selection_methods"] if m["name"]=="random"];assert len(random)==1 and random[0].get("n_features")==[500]
 new=copy.deepcopy(orig);changed=copy.deepcopy(random[0]);changed["n_features"]=[2000];new["feature_selection_methods"]=[changed];new["output_dir"]=str(revision/"runs"/name)
 for k in orig:
  if k not in ("feature_selection_methods","output_dir"):assert orig[k]==new[k],(name,k)
 assert {k:v for k,v in random[0].items() if k!="n_features"}=={k:v for k,v in changed.items() if k!="n_features"};assert orig["seeds"]==[0,1,2]
 dest=configdir/(name+".yaml");dest.write_text(yaml.safe_dump(new,sort_keys=False),encoding="utf8")
 manifest.append(dict(name=name,source_config=relative,source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),generated_config=str(dest.relative_to(root)),generated_sha256=hashlib.sha256(dest.read_bytes()).hexdigest(),dataset=orig["datasets"],seeds=orig["seeds"],integrators=[x["name"] for x in orig["integration_methods"]],tasks=[x["name"] for x in orig["tasks"]],n_features_before=500,n_features_after=2000))
(revision/"config_manifest.json").write_text(json.dumps(manifest,indent=2,ensure_ascii=False)+"\n");print(json.dumps(manifest,indent=2))
