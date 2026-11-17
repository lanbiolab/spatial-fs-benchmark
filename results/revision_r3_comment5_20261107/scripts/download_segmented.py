from pathlib import Path
import requests,time,json,hashlib,concurrent.futures,shutil
OUT=Path('/home/xd/y/Benchmark/results/revision_r3_comment5_20261107');D=OUT/'candidate_data'
jobs=[('st_counts.tsv.gz','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147747/suppl/GSE147747_expr_raw_counts_table.tsv.gz',163382763,None),('openst_adult_mouse_by_cell.h5ad.tar.gz','https://bimsbstatic.mdc-berlin.de/rajewsky/openst-public-data/adult_mouse_by_cell.h5ad.tar.gz',683465819,'c837f3184d25fb4842a0e27bde1b3941')]
def download(job):
 name,url,size,expected=job;p=D/name;status=D/(name+'.status.json');partial=D/(name+'.partial');parts=D/(name+'.parts');parts.mkdir(exist_ok=True);start=time.monotonic()
 if p.exists():return
 prefix=partial.stat().st_size if partial.exists() else 0
 assert prefix<=size
 def get(span):
  lo,hi=span;dest=parts/f'{lo:012d}_{hi:012d}';length=hi-lo+1
  if dest.exists() and dest.stat().st_size==length:return
  for attempt in range(4):
   try:
    with requests.get(url,headers={'Range':f'bytes={lo}-{hi}'},stream=True,timeout=(25,60)) as r:
     assert r.status_code==206 and r.headers.get('Content-Range')==f'bytes {lo}-{hi}/{size}',(r.status_code,r.headers.get('Content-Range'))
     with dest.with_suffix('.partial').open('wb') as f:
      for chunk in r.iter_content(256*1024):
       if chunk:f.write(chunk)
    assert dest.with_suffix('.partial').stat().st_size==length
    dest.with_suffix('.partial').replace(dest);return
   except Exception:
    if attempt==3:raise
    time.sleep(3)
 spans=[(lo,min(lo+4*1024*1024,size)-1) for lo in range(prefix,size,4*1024*1024)]
 meta={'name':name,'url':url,'state':'downloading_segmented','content_length':size,'resumed_prefix_bytes':prefix,'workers':4,'expected_md5':expected}
 try:
  with concurrent.futures.ThreadPoolExecutor(4) as pool:
   fs=[pool.submit(get,s) for s in spans]
   while not all(f.done() for f in fs):
    done=sum((hi-lo+1) for (lo,hi),f in zip(spans,fs) if f.done() and f.exception() is None)
    meta.update(bytes=prefix+done,elapsed_seconds=round(time.monotonic()-start));status.write_text(json.dumps(meta,indent=2));time.sleep(10)
   for f in fs:f.result()
  with partial.open('ab') as f:
   for lo,hi in spans:
    with (parts/f'{lo:012d}_{hi:012d}').open('rb') as src:shutil.copyfileobj(src,f)
  assert partial.stat().st_size==size
  md5=hashlib.md5();sha=hashlib.sha256()
  with partial.open('rb') as f:
   for b in iter(lambda:f.read(8*1024*1024),b''):md5.update(b);sha.update(b)
  if expected:assert md5.hexdigest()==expected,'MD5 mismatch'
  partial.replace(p);meta.update(state='downloaded',bytes=size,md5=md5.hexdigest(),sha256=sha.hexdigest(),elapsed_seconds=round(time.monotonic()-start))
 except Exception as e:meta.update(state='failed',error=repr(e))
 status.write_text(json.dumps(meta,indent=2));print(json.dumps(meta),flush=True)
with concurrent.futures.ThreadPoolExecutor(2) as ex:list(ex.map(download,jobs))
