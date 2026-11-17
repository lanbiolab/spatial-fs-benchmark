import json
import os
import subprocess
import time
from pathlib import Path

root=Path('/home/xd/y/Benchmark');out=root/'results/revision_major3_missing_ranks_20261105'
datasets=['mouse_brain_serial_sections','dlpfc','stomics_0212','stomics_0218','stomics_0224']
methods=['morans_i','somde','spatialde','sparkx','nnsvg']
integrators=['scvi','cellcharter','graphst']
env=os.environ.copy();env.update(PYTHONPATH='src',CUDA_VISIBLE_DEVICES='3',OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2',MKL_NUM_THREADS='2',NUMBA_NUM_THREADS='2')
attempted=set()
print('Downstream queue: three original integrators, five primary datasets; one job on GPU 3 at a time',flush=True)
while True:
    progressed=False
    for dataset in datasets:
        for method in methods:
            for seed in range(3):
                selection_seed=seed if method=='spatialde' else 0
                selection=out/'selection'/dataset/method/f'seed{selection_seed}'
                if not (selection/'DONE').exists(): continue
                summary=json.loads((selection/'summary.json').read_text())
                versions=['missing_excluded']
                if not summary['old_reproduction_exact_order']:
                    versions.append('original_zero')
                for integ in integrators:
                    for version in versions:
                        key=(dataset,method,seed,integ,version)
                        job=out/'downstream'/dataset/integ/method/f'seed{seed}'/version
                        if (job/'DONE').exists() or key in attempted: continue
                        attempted.add(key);job.mkdir(parents=True,exist_ok=True)
                        cmd=[str(root/'.conda-env/bin/python'),'-u',str(out/'scripts/run_downstream.py'),dataset,method,str(seed),integ,version]
                        print('START',key,flush=True)
                        with (job/'run.log').open('a') as log:
                            ret=subprocess.run(cmd,cwd=root,env=env,stdout=log,stderr=subprocess.STDOUT).returncode
                        print('FINISH',key,'exit',ret,flush=True)
                        if ret: (job/'FAILED').write_text(f'exit {ret}\n')
                        progressed=True
    selection_log=out/'selection_queue.log'
    queue_done=selection_log.exists() and 'Queue finished' in selection_log.read_text(errors='replace')
    if queue_done and not progressed: break
    if not progressed: time.sleep(15)
print('Downstream queue finished; check FAILED and uncompleted selections',flush=True)
