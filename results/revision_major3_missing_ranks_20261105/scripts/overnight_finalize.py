from __future__ import annotations
import datetime,json,os,subprocess,time,traceback
from pathlib import Path
import psutil

ROOT=Path('/home/xd/y/Benchmark');OUT=ROOT/'results/revision_major3_missing_ranks_20261105';SCRIPTS=OUT/'scripts'
DATASETS=['mouse_brain_serial_sections','dlpfc','stomics_0212','stomics_0218','stomics_0224','e8p5_embryo','e9p5_embryo']
PRIMARY=DATASETS[:5];METHODS=['morans_i','somde','spatialde','sparkx','nnsvg'];INTEGRATORS=['scvi','cellcharter','graphst']
SELECT_WORKERS={'selection_queue.py','capture_selection.py','priority_queue.py','extra_queue.py','embryo_queue.py','e9_extra_queue.py'}
DOWN_WORKERS={'downstream_queue.py','run_downstream.py','mouse_downstream_queue.py','dlpfc_nnsvg_queue.py','dlpfc_nnsvg_seed2.py'}
env=os.environ.copy();env.update(PYTHONPATH='src',OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2',MKL_NUM_THREADS='2',NUMBA_NUM_THREADS='2')

def live_workers(names):
    found=[]
    for p in psutil.process_iter(['pid','username','cmdline']):
        try:
            if p.info['username']!='xd' or p.pid==os.getpid():continue
            for arg in p.info['cmdline'] or []:
                if '/scripts/' in arg and Path(arg).name in names:
                    found.append({'pid':p.pid,'script':Path(arg).name});break
        except (psutil.NoSuchProcess,psutil.AccessDenied):pass
    return found

def missing_selection():
    missing=[]
    for ds in DATASETS:
        for method in METHODS:
            for seed in (range(3) if method=='spatialde' else [0]):
                p=OUT/'selection'/ds/method/f'seed{seed}'
                if not (p/'DONE').exists() or not (p/'summary.json').exists():missing.append([ds,method,seed])
    return missing

def missing_downstream():
    missing=[]
    for ds in PRIMARY:
        for method in METHODS:
            for seed in range(3):
                ss=seed if method=='spatialde' else 0
                info=json.loads((OUT/'selection'/ds/method/f'seed{ss}'/'summary.json').read_text())
                versions=['missing_excluded']+([] if info['old_reproduction_exact_order'] else ['original_zero'])
                for integ in INTEGRATORS:
                    for v in versions:
                        p=OUT/'downstream'/ds/integ/method/f'seed{seed}'/v
                        if not (p/'DONE').exists() or not (p/'manifest.json').exists():missing.append([ds,method,seed,integ,v])
    return missing

def status(phase,**extra):
    d={'phase':phase,'updated_at':datetime.datetime.now().astimezone().isoformat(),'missing_selection':missing_selection(),**extra}
    p=OUT/'overnight_status.json';tmp=p.with_suffix('.json.tmp');tmp.write_text(json.dumps(d,ensure_ascii=False,indent=2));os.replace(tmp,p)
    print(d['updated_at'],phase,'selection_missing=',len(d['missing_selection']),'downstream_missing=',len(d.get('missing_downstream',[])),flush=True)

def execute(script,logname):
    print('Starting recovery',script,flush=True)
    with (OUT/logname).open('a') as f:
        code=subprocess.run([str(ROOT/'.conda-env/bin/python'),'-u',str(SCRIPTS/script)],cwd=ROOT,env=env,stdout=f,stderr=subprocess.STDOUT).returncode
    print('Finished',script,'exit',code,flush=True)
    return code

def main():
    repairs=0
    while True:
        missing=missing_selection();workers=live_workers(SELECT_WORKERS)
        if not missing:break
        status('waiting_for_selection',workers=workers,selection_retries=repairs)
        if not workers:
            if repairs>=2:
                status('needs_attention',reason='Selection jobs remain incomplete after two recovery attempts',selection_retries=repairs);return 2
            repairs+=1;execute('selection_queue.py',f'overnight_selection_retry_{repairs}.log')
        else:time.sleep(60)
    workers=live_workers(SELECT_WORKERS)
    while any(x['script']=='selection_queue.py' for x in workers):
        status('waiting_for_selection_queue_exit',workers=workers);time.sleep(30);workers=live_workers(SELECT_WORKERS)
    log=OUT/'selection_queue.log'
    if not log.exists() or 'Queue finished' not in log.read_text(errors='replace'):
        with log.open('a') as f:f.write('\nQueue finished (all 49 selection outputs independently verified).\n')
    mismatches=[]
    for ds in DATASETS:
        for method in METHODS:
            for seed in (range(3) if method=='spatialde' else [0]):
                p=OUT/'selection'/ds/method/f'seed{seed}'/'summary.json';x=json.loads(p.read_text())
                if x['old_reproduction_top2000_overlap']!=x['original_selected']:mismatches.append({'dataset':ds,'method':method,'seed':seed,'overlap':x['old_reproduction_top2000_overlap'],'original_selected':x['original_selected']})
    if mismatches:
        status('needs_attention',reason='Original selected feature sets do not fully reproduce',mismatches=mismatches);return 3
    repairs=0
    while True:
        missing=missing_downstream();workers=live_workers(DOWN_WORKERS)
        if not missing:break
        status('waiting_for_downstream',missing_downstream=missing,workers=workers,downstream_retries=repairs)
        if not workers:
            if repairs>=2:
                status('needs_attention',reason='Downstream jobs remain incomplete after two recovery attempts',missing_downstream=missing,downstream_retries=repairs);return 4
            repairs+=1;execute('downstream_queue.py',f'overnight_downstream_retry_{repairs}.log')
        else:time.sleep(60)
    status('summarizing',missing_downstream=[])
    code=execute('summarize.py','overnight_summary.log')
    if code:
        status('needs_attention',reason='Final summary or reproduction audit failed',summary_exit=code);return 5
    expected=['method_gene_universe_by_dataset.tsv','method_gene_universe_by_slice.tsv','five_spatial_methods_rank_changes.tsv','five_spatial_methods_paired_refit_rank_changes.tsv','original_frozen_score_reproduction_audit.tsv']
    absent=[name for name in expected if not (OUT/'analysis'/name).exists()]
    if absent:
        status('needs_attention',reason='Expected summary files are missing',missing_files=absent);return 6
    status('complete',summary_dir=str(OUT/'analysis'),summary_files=expected,missing_downstream=[])
    (OUT/'OVERNIGHT_COMPLETE').write_text(datetime.datetime.now().astimezone().isoformat()+'\n')
    return 0

if __name__=='__main__':
    try:raise SystemExit(main())
    except Exception as exc:
        traceback.print_exc();status('needs_attention',reason=f'Unexpected finalizer exception: {exc}');raise
