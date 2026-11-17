import argparse
import json
import types
from pathlib import Path

from spatial_fs_benchmark.config import load_benchmark_config
from spatial_fs_benchmark.benchmark.experiment import read_feature_selection
import spatial_fs_benchmark.benchmark.runner as runner_module

ROOT=Path('/home/xd/y/Benchmark')
OUT=ROOT/'results/revision_major3_missing_ranks_20261105'

def main():
    p=argparse.ArgumentParser()
    p.add_argument('dataset');p.add_argument('method');p.add_argument('seed',type=int)
    p.add_argument('integrator',choices=['scvi','cellcharter','graphst'])
    p.add_argument('version',choices=['missing_excluded','original_zero'])
    a=p.parse_args()
    config_root='graphst_canonical' if a.integrator=='graphst' else 'canonical'
    cfg=load_benchmark_config(ROOT/f'configs/rebuild_v1/downstream/{config_root}/{a.dataset}.yaml')
    cfg.feature_selection_methods=[x for x in cfg.feature_selection_methods if x.name==a.method]
    cfg.integration_methods=[x for x in cfg.integration_methods if x.name==a.integrator]
    assert len(cfg.feature_selection_methods)==len(cfg.integration_methods)==1
    cfg.seeds=[a.seed];cfg.n_features=[2000];cfg.feature_selection_methods[0].n_features=[2000]
    selection_seed=a.seed if a.method=='spatialde' else 0
    folder=OUT/'selection'/a.dataset/a.method/f'seed{selection_seed}'
    assert (folder/'DONE').exists()
    saved=folder/f'{a.version}_n2000.json'
    selection=read_feature_selection(saved)
    orig=ROOT/f'results/spatial_svg_rebuild_v1/{a.dataset}/{cfg.datasets[0].name.lower()}/{a.integrator}/{a.method}/n2000/seed{a.seed}'
    orig_selection=ROOT/f'results/spatial_svg_rebuild_v1/{a.dataset}/{cfg.datasets[0].name.lower()}/feature_selection/{a.method}/n2000/seed{selection_seed}/selected_features.json'
    original=read_feature_selection(orig_selection)
    job=OUT/'downstream'/a.dataset/a.integrator/a.method/f'seed{a.seed}'/a.version
    job.mkdir(parents=True,exist_ok=True)
    required=[t.name for t in cfg.tasks]
    if selection.feature_names == original.feature_names:
        for task in required:
            assert (orig/f'{task}_records.json').exists(), (orig,task)
        (job/'manifest.json').write_text(json.dumps({'mode':'reuse_unchanged_original','original_dir':str(orig),'selection':str(saved),'tasks':required},indent=2))
        (job/'DONE').write_text('ok\n'); print('Reused exactly identical original selected feature order',flush=True); return
    cfg.output_dir=str(job/'result')
    factory=runner_module.build_feature_selector
    def build(name,**params):
        selector=factory(name,**params)
        selector.implementation_version=str(getattr(selector,'implementation_version',''))+'_m3_'+a.version
        def select(self,dataset,n_features,random_seed=0):
            assert n_features==2000 and random_seed==selection_seed
            assert all(dataset.adata.var_names[idx]==gene for idx,gene in zip(selection.feature_indices,selection.feature_names))
            return read_feature_selection(saved)
        selector.select=types.MethodType(select,selector)
        return selector
    runner_module.build_feature_selector=build
    runner_module.BenchmarkRunner(cfg).run()
    result=job/'result'/cfg.datasets[0].name.lower()/a.integrator/a.method/'n2000'/f'seed{a.seed}'
    for task in required:
        path=result/f'{task}_records.json'
        assert path.exists(),f'Incomplete {path}'
        data=json.loads(path.read_text());assert data, f'Empty {path}'
    (job/'manifest.json').write_text(json.dumps({'mode':'recomputed','result_dir':str(result),'selection':str(saved),'tasks':required},indent=2))
    (job/'DONE').write_text('ok\n')
    print('Downstream results verified',a.dataset,a.method,a.seed,a.integrator,a.version,flush=True)

if __name__=='__main__':
    import fcntl
    import sys
    if len(sys.argv) != 6:
        raise SystemExit('Expected dataset method seed integrator version')
    _,ds,method,seed,integ,version=sys.argv
    lock_dir=OUT/'downstream'/ds/integ/method/f'seed{seed}'/version
    lock_dir.mkdir(parents=True,exist_ok=True)
    with (lock_dir/'job.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        if (lock_dir/'DONE').exists(): print('Job already complete; no recomputation',flush=True)
        else: main()
        if (lock_dir/'DONE').exists():
            (lock_dir/'FAILED').unlink(missing_ok=True)
