import concurrent.futures as cf,json,os,subprocess,time
from pathlib import Path
ROOT=Path('/home/xd/y/Benchmark');OUT=ROOT/'results/revision_r2_major2_20261106'
runner=OUT/'scripts/run_individual.py'
protocol=json.loads((OUT/'analysis_protocol.json').read_text())
core=[m for m in protocol['methods'] if m not in ['nnsvg','spatialde','somde']]
groups=[('core',core,4),('spatialde_somde',['spatialde','somde'],2),('nnsvg',['nnsvg'],4)]
env=os.environ.copy();env.update(PYTHONPATH=str(ROOT/'src'),OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2',MKL_NUM_THREADS='2',NUMEXPR_NUM_THREADS='2',RETICULATE_PYTHON=str(ROOT/'.conda-env/bin/python'))
def job(group,methods,fold):
    logfile=OUT/'logs'/(fold['dataset']+'__'+fold['individual']+'__'+group+'.log');logfile.parent.mkdir(exist_ok=True)
    cmd=[str(ROOT/'.conda-env/bin/python'),str(runner),'--dataset',fold['dataset'],'--individual',fold['individual'],'--methods',','.join(methods)]
    with logfile.open('a') as f:
        result=subprocess.run(cmd,cwd=ROOT,env=env,stdout=f,stderr=subprocess.STDOUT)
    print(time.strftime('%F %T'),group,fold['dataset'],fold['individual'],'exit',result.returncode,flush=True)
    return {'group':group,'dataset':fold['dataset'],'individual':fold['individual'],'exit':result.returncode,'log':str(logfile)}
def group_run(group,methods,workers):
    with cf.ThreadPoolExecutor(max_workers=workers) as pool:
        return list(pool.map(lambda f:job(group,methods,f),protocol['folds']))
with cf.ThreadPoolExecutor(max_workers=3) as pool:
    results=list(pool.map(lambda x:group_run(*x),groups))
(OUT/'queue_completion.json').write_text(json.dumps(results,indent=2))
print('ALL JOB GROUPS EXITED',flush=True)
subprocess.run([str(ROOT/'.conda-env/bin/python'),str(OUT/'scripts/summarize_individual.py')],cwd=ROOT,env=env,check=True)
