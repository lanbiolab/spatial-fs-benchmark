"""Claim pending heavy folds with the unchanged experiment runner and per-method locks."""
from pathlib import Path
import concurrent.futures as cf
import json,os,subprocess,time,shlex

ROOT=Path('/home/xd/y/Benchmark');OUT=ROOT/'results/revision_r2_major2_20261106';RUNNER=OUT/'scripts/run_individual.py'
spec=json.loads((OUT/'analysis_protocol.json').read_text())
active=set()
for line in subprocess.check_output(['ps','-u','xd','-o','args='],text=True).splitlines():
    if str(RUNNER) not in line or '--dataset' not in line:continue
    args=shlex.split(line)
    ds=args[args.index('--dataset')+1];donor=args[args.index('--individual')+1]
    methods=args[args.index('--methods')+1].split(',')
    for method in methods:active.add((ds,donor,method))
env=os.environ.copy();env.update(PYTHONPATH=str(ROOT/'src'),OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2',MKL_NUM_THREADS='2',NUMEXPR_NUM_THREADS='2',RETICULATE_PYTHON=str(ROOT/'.conda-env/bin/python'))
def status(fold,method):
    p=OUT/'folds'/(fold['dataset']+'__individual_'+fold['individual'])/'status'/(method+'.json')
    return json.loads(p.read_text()).get('status') if p.exists() else 'not_started'
def candidates(methods):
    selected=[]
    for f in reversed(spec['folds']):
        if any((f['dataset'],f['individual'],m) in active for m in methods):continue
        if all(status(f,m) in ['complete','not_applicable'] for m in methods):continue
        selected.append(f)
    return selected
groups=[('extra_nnsvg',['nnsvg'],4),('extra_spatialde_somde',['spatialde','somde'],2)]
planned={name:candidates(methods) for name,methods,_ in groups}
(OUT/'additional_parallel_plan.json').write_text(json.dumps({'created_at':time.time(),'existing_active_methods_excluded':sorted(active),'planned_groups':planned,'settings_unchanged':True,'mutual_exclusion':'exclusive advisory lock per fold/method; cached complete results reused'},indent=2))
def one(name,methods,fold):
    if all(status(fold,m) in ['complete','not_applicable'] for m in methods):return
    logfile=OUT/'logs'/(fold['dataset']+'__'+fold['individual']+'__'+name+'.log')
    cmd=[str(ROOT/'.conda-env/bin/python'),str(RUNNER),'--dataset',fold['dataset'],'--individual',fold['individual'],'--methods',','.join(methods)]
    with logfile.open('a') as log:r=subprocess.run(cmd,cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT)
    print(time.strftime('%F %T'),name,fold['dataset'],fold['individual'],'exit',r.returncode,flush=True)
def group(name,methods,workers):
    with cf.ThreadPoolExecutor(max_workers=workers) as pool:list(pool.map(lambda f:one(name,methods,f),planned[name]))
with cf.ThreadPoolExecutor(max_workers=len(groups)) as pool:list(pool.map(lambda g:group(*g),groups))
print('ADDITIONAL QUEUE FINISHED',flush=True)
subprocess.run([str(ROOT/'.conda-env/bin/python'),str(OUT/'scripts/summarize_individual.py')],cwd=ROOT,env=env,check=True)
