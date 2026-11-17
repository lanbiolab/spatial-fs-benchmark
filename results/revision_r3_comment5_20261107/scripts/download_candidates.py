from pathlib import Path
import requests,time,json,hashlib,concurrent.futures
OUT=Path('/home/xd/y/Benchmark/results/revision_r3_comment5_20261107');D=OUT/'candidate_data';D.mkdir(parents=True,exist_ok=True)
jobs=[
('st_metadata.tsv.gz','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147747/suppl/GSE147747_meta_table.tsv.gz',None),
('st_counts.tsv.gz','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147747/suppl/GSE147747_expr_raw_counts_table.tsv.gz',None),
('stereo_mouse_brain.h5ad','https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000058/stomics/Mouse_brain.h5ad',None),
('openst_adult_mouse_by_cell.h5ad.tar.gz','https://bimsbstatic.mdc-berlin.de/rajewsky/openst-public-data/adult_mouse_by_cell.h5ad.tar.gz','c837f3184d25fb4842a0e27bde1b3941')]
def fetch(job):
 name,url,expected=job;p=D/name;status=D/(name+'.status.json');meta={'name':name,'url':url,'expected_md5':expected,'state':'started','started_at':time.strftime('%Y-%m-%dT%H:%M:%S%z')};status.write_text(json.dumps(meta,indent=2));start=time.monotonic()
 try:
  if p.exists():meta.update(state='already_downloaded',bytes=p.stat().st_size);return meta
  partial=D/(name+'.partial');n=0;md5=hashlib.md5();sha=hashlib.sha256();last=0
  with requests.get(url,stream=True,timeout=(25,45)) as r:
   r.raise_for_status();length=int(r.headers.get('Content-Length',0));meta.update(content_length=length,content_type=r.headers.get('Content-Type'))
   if length>2_500_000_000:raise ValueError('Candidate file exceeds 2.5GB per-file audit limit')
   with partial.open('wb') as f:
    for chunk in r.iter_content(1024*1024):
     if not chunk:continue
     n+=len(chunk)
     if n>2_500_000_000 or time.monotonic()-start>1800:raise ValueError('Candidate download size/time audit limit reached')
     f.write(chunk);md5.update(chunk);sha.update(chunk)
     if time.monotonic()-last>15:
      meta.update(state='downloading',bytes=n,elapsed_seconds=round(time.monotonic()-start));status.write_text(json.dumps(meta,indent=2));last=time.monotonic()
   if length and length!=n:raise ValueError('Incomplete byte count')
   if expected and expected!=md5.hexdigest():raise ValueError('MD5 mismatch')
   partial.replace(p);meta.update(state='downloaded',bytes=n,md5=md5.hexdigest(),sha256=sha.hexdigest(),elapsed_seconds=round(time.monotonic()-start))
 except Exception as e:meta.update(state='failed',error=str(e),elapsed_seconds=round(time.monotonic()-start))
 finally:status.write_text(json.dumps(meta,indent=2));print(json.dumps(meta),flush=True)
 return meta
with concurrent.futures.ThreadPoolExecutor(4) as ex:rows=list(ex.map(fetch,jobs))
(OUT/'download_manifest.json').write_text(json.dumps(rows,indent=2))
