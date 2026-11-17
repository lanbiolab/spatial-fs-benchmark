from pathlib import Path
import os,sys,json,queue,time,subprocess,concurrent.futures,threading
ROOT=Path('/home/xd/y/Benchmark');OUT=ROOT/'results/revision_r2_m1_scale_match_20261111'
pilot='--pilot' in sys.argv;accelerate='--accelerate' in sys.argv;tailwork='--tailwork' in sys.argv;references='--references' in sys.argv;mode='pilot' if pilot else ('references' if references else ('tailwork' if tailwork else ('accelerate' if accelerate else 'full')))
jobs=json.loads((OUT/'jobs.json').read_text())
if pilot:jobs=[j for j in jobs if (j['arm'],j['method']) in [('pseudo_visium','morans_i'),('pseudo_visium','random'),('pseudo_visium','somde'),('fine_footprint','morans_i')]]
else:
    gate=json.loads((OUT/'pilot_queue_finished.json').read_text());assert gate['all_completed'],gate
    # Start the slower spatial/R procedures early to avoid a long serial tail.
    heavy={'nnsvg':0,'spatialde':1,'seurat_sct':2,'scPNMF':3,'singleCellHaystack':4,'dubstepr':5}
    jobs.sort(key=lambda j:(heavy.get(j['method'],10),j['method'],j['arm']))
    if accelerate:
        running={(s['arm'],s['method']) for p in (OUT/'logs').glob('full__*.status.json') for s in [json.loads(p.read_text())] if s['state'] in ['running','completed']}
        jobs=[j for j in reversed(jobs) if (j['arm'],j['method']) not in running]
    if tailwork:
        running={(s['arm'],s['method']) for p in (OUT/'logs').glob('*.status.json') for s in [json.loads(p.read_text())] if s['state'] in ['running','completed']}
        jobs=[j for j in jobs if j['method'] in ['scPNMF','singleCellHaystack','dubstepr'] and (j['arm'],j['method']) not in running]
    if references:jobs=[j for j in jobs if j['method']=='all_features']
q=queue.Queue();[q.put(j) for j in jobs];lock=threading.Lock();states=[]
def worker(gpu):
    while True:
        try:j=q.get_nowait()
        except queue.Empty:return
        name=j['arm']+'__'+j['method'];status=OUT/'logs'/f'{mode}__{name}.status.json'
        state={**j,'gpu':gpu,'state':'running','started_at':time.strftime('%Y-%m-%dT%H:%M:%S%z')};status.write_text(json.dumps(state,indent=2))
        env=os.environ.copy();env.update(CUDA_VISIBLE_DEVICES=str(gpu),OMP_NUM_THREADS='4',OPENBLAS_NUM_THREADS='4',MKL_NUM_THREADS='4',NUMBA_NUM_THREADS='4',PYTHONPATH=str(ROOT/'src'),PYTHONUNBUFFERED='1')
        print('START',mode,name,gpu,flush=True)
        try:
            with (OUT/'logs'/f'{mode}__{name}.log').open('w') as f:
                cmd=[str(ROOT/'.conda-env/bin/python'),'-u',str(OUT/'scripts/run_job.py'),j['config']]+(['--pilot'] if pilot else [])
                proc=subprocess.Popen(cmd,cwd=ROOT,env=env,stdout=f,stderr=subprocess.STDOUT,start_new_session=True)
                state['pid']=proc.pid;status.write_text(json.dumps(state,indent=2));rc=proc.wait()
            state.update(state='completed' if rc==0 else 'failed_or_incomplete',exit_code=rc)
        except Exception as e:state.update(state='failed_or_incomplete',error=repr(e))
        state['ended_at']=time.strftime('%Y-%m-%dT%H:%M:%S%z');status.write_text(json.dumps(state,indent=2))
        with lock:states.append(state);(OUT/f'{mode}_queue_status.json').write_text(json.dumps(states,indent=2))
        print('END',mode,name,state['state'],flush=True);q.task_done()
with concurrent.futures.ThreadPoolExecutor(2 if tailwork or references else 4) as ex:list(ex.map(worker,range(2 if tailwork or references else 4)))
(OUT/f'{mode}_queue_finished.json').write_text(json.dumps({'jobs':states,'all_completed':all(x['state']=='completed' for x in states)},indent=2))
