import json,os,subprocess,time
from pathlib import Path
root=Path('/home/xd/y/Benchmark');out=root/'results/revision_major3_missing_ranks_20261105'
ds='mouse_brain_serial_sections';methods=['somde','spatialde','sparkx','nnsvg'];integrators=['scvi','cellcharter','graphst']
env=os.environ.copy();env.update(PYTHONPATH='src',CUDA_VISIBLE_DEVICES='3',OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2',MKL_NUM_THREADS='2',NUMBA_NUM_THREADS='2')
attempted=set();print('Second worker: mouse brain paired downstream jobs',flush=True)
while True:
    progressed=False
    for method in methods:
        for seed in range(3):
            ss=seed if method=='spatialde' else 0
            select=out/'selection'/ds/method/f'seed{ss}'
            if not (select/'DONE').exists():continue
            info=json.loads((select/'summary.json').read_text());versions=['missing_excluded']+([] if info['old_reproduction_exact_order'] else ['original_zero'])
            for integ in integrators:
                for version in versions:
                    key=(method,seed,integ,version);job=out/'downstream'/ds/integ/method/f'seed{seed}'/version
                    if (job/'DONE').exists() or key in attempted:continue
                    attempted.add(key);job.mkdir(parents=True,exist_ok=True)
                    cmd=[str(root/'.conda-env/bin/python'),'-u',str(out/'scripts/run_downstream.py'),ds,method,str(seed),integ,version]
                    print('START',key,flush=True)
                    with (job/'run.log').open('a') as log:ret=subprocess.run(cmd,cwd=root,env=env,stdout=log,stderr=subprocess.STDOUT).returncode
                    print('FINISH',key,'exit',ret,flush=True)
                    if ret:(job/'FAILED').write_text(f'exit {ret}\n')
                    progressed=True
    if all((out/'selection'/ds/m/'seed0'/'DONE').exists() for m in methods) and not progressed:break
    if not progressed:time.sleep(20)
print('Mouse downstream worker finished',flush=True)
