from pathlib import Path
import time,json,subprocess,sys
R=Path('/home/xd/y/Benchmark');O=R/'results/revision_r3_comment5_20261107';S=O/'candidate_data/st_counts.tsv.gz.status.json';t=time.monotonic()
while True:
 try:state=json.loads(S.read_text())
 except (FileNotFoundError,json.JSONDecodeError):state={}
 if state.get('state')=='downloaded':break
 if state.get('state')=='failed':raise RuntimeError('ST download failed; no ST benchmark started')
 if time.monotonic()-t>12*3600:raise RuntimeError('ST data not available after 12h; inspect download')
 time.sleep(20)
subprocess.run([str(R/'.conda-env/bin/python'),'-u',str(O/'scripts/prepare_st_data.py')],cwd=R,check=True)
# Share the four GPUs with the Stereo-seq queue by waiting for that queue to exit.
while not (O/'stereo_queue_finished.json').exists():
 if time.monotonic()-t>48*3600:raise RuntimeError('Stereo queue completion absent after 48h')
 time.sleep(30)
subprocess.run([str(R/'.conda-env/bin/python'),'-u',str(O/'scripts/run_st_queue.py')],cwd=R,check=True)
