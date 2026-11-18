#!/usr/bin/env python3
"""Compare mapped/LZ4 persistence in the same current library on sealed study inputs.
Usage: python3 benchmarks/mapped_storage_study.py BUILD_DIRECTORY OUTPUT_DIRECTORY
Requires build/mapped-cache-2026-09-22 from mapped_cache_study.py and its follow-up.
"""
import hashlib,json,os,platform,random,shutil,statistics,subprocess,sys
from pathlib import Path
from time import perf_counter
ROOT=Path(__file__).resolve().parents[1]
build=Path(sys.argv[1]).resolve();out=Path(sys.argv[2]).resolve();out.mkdir(parents=True,exist_ok=False)
prior=ROOT/'build/mapped-cache-2026-09-22'
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def write(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def run(cmd):
    cmd=list(map(str,cmd));start=perf_counter();p=subprocess.run(cmd,text=True,capture_output=True)
    with (out/'commands.jsonl').open('a') as f:f.write(json.dumps(dict(command=cmd,returncode=p.returncode,wall_seconds=perf_counter()-start,stdout=p.stdout,stderr=p.stderr))+'\n')
    p.check_returncode();return p.stdout
# Reuse the fully consuming, equally warmed worker instead of a new benchmark path.
worker=out/'worker.cpp';shutil.copy2(prior/'equal-warmup-worker.cpp',worker)
hts=ROOT/'build/research-inputs/samtools-1.24/htslib-1.24'
run(['c++','-std=c++20','-O3','-DNDEBUG','-include','cstdint','-pthread','-I'+str(ROOT/'include'),'-I'+str(ROOT/'build/_deps/ankerl_unordered_dense-src/include'),'-I'+str(hts),worker,build/'libtraceon_core.a',build/'_deps/zlib-ng-build/libz.a',build/'third_party/lz4/build/cmake/liblz4.a',hts/'libhts.a','-lz','-lbz2','-llzma','-ldeflate','-lcurl','-lcrypto','-lm','-o',out/'worker'])
write(out/'environment.json',dict(platform=platform.platform(),cpu=Path('/proc/cpuinfo').read_text(),affinity=sorted(os.sched_getaffinity(0)),core_sha256=sha(build/'libtraceon_core.a'),compiler=run(['c++','--version']),head=run(['git','-C',ROOT,'rev-parse','HEAD']).strip()))
shutil.copy2(__file__,out/Path(__file__).name)
for source in [* (ROOT/'src').glob('*'),*(ROOT/'include').glob('*')]:
    if source.is_file():
        dest=out/'source'/source.relative_to(ROOT);dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,dest)
datasets=[d for d in json.loads((prior/'datasets.json').read_text()) if d['name'].endswith('-100MB')]
rows=[];memory=[];builds=[]
for d in datasets:
    source=prior/d['name']/'reads.fa';ids=prior/d['name']/'queries.txt'
    assert sha(source)==d['input_sha256'] and sha(ids)==d['queries_sha256']
    for arm in ['mapped','lz4']:
        target=out/(d['name']+'-'+arm+'.traceon')
        start=perf_counter();run(['taskset','-c','8-11',build/'traceon_make_cache',source,target]+(['--mapped'] if arm=='mapped' else []))
        builds.append(dict(dataset=d['name'],arm=arm,wall_seconds=perf_counter()-start,bytes=target.stat().st_size,sha256=sha(target)))
        write(out/'builds.json',builds)
    for rep in range(6):
        arms=['mapped','lz4'];random.Random(20260922+rep).shuffle(arms)
        for arm in arms:
            target=out/(d['name']+'-'+arm+'.traceon')
            result=json.loads(run(['taskset','-c','8',out/'worker','query','cache',target,ids,d['count'],d['length']]))
            assert result['sha256']==d['query_sha256'] and result['records']==10000
            rows.append(dict(dataset=d['name'],arm=arm,rep=rep,**result));write(out/'queries.json',rows)
        for nproc in [1,4]:
            for arm in arms:
                processes=[];results=[];samples=[]
                try:
                    for i in range(nproc):
                        cmd=list(map(str,['taskset','-c',8+i,out/'worker','hold','cache',out/(d['name']+'-'+arm+'.traceon'),'-',d['count'],d['length']]))
                        p=subprocess.Popen(cmd,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True);processes.append(p)
                    for p in processes:
                        line=p.stdout.readline()
                        if not line:raise RuntimeError(p.stderr.read())
                        r=json.loads(line);assert r['sha256']==d['full_sha256'] and r['records']==d['count'];results.append(r)
                    for p in processes:
                        smaps=Path(f'/proc/{p.pid}/smaps_rollup').read_text()
                        sample={k:int(v.split()[0]) for k,v in (line.split(':',1) for line in smaps.splitlines()[1:])}
                        samples.append({k:sample[k] for k in ['Rss','Pss','Pss_Anon','Pss_File']})
                    memory.append(dict(dataset=d['name'],arm=arm,rep=rep,nproc=nproc,workers=results,memory_kib=samples,total_kib={k:sum(s[k] for s in samples) for k in samples[0]}));write(out/'memory.json',memory)
                finally:
                    for p in processes:
                        if p.poll() is None:p.stdin.write('quit\n');p.stdin.flush()
                    for p in processes:
                        _,err=p.communicate(timeout=30)
                        if p.returncode:raise RuntimeError(err)
    print(d['name']+' complete',flush=True)
from summarize_paired import interval
summary=[]
for d in datasets:
    for arm in ['mapped','lz4']:
        r=[x for x in rows if x['dataset']==d['name'] and x['arm']==arm]
        summary.append(dict(dataset=d['name'],arm=arm,setup_ms=statistics.median(x['setup_ms'] for x in r),consume_ms=statistics.median(x['consume_ms'] for x in r),pss_mib={n:statistics.median(x['total_kib']['Pss']/1024 for x in memory if x['dataset']==d['name'] and x['arm']==arm and x['nproc']==n) for n in [1,4]}))
write(out/'summary.json',summary);write(out/'inputs.json',datasets)
paired=[]
for d in datasets:
    for metric in ['setup_ms','consume_ms']:
        a={r['rep']:r[metric] for r in rows if r['dataset']==d['name'] and r['arm']=='lz4'};b={r['rep']:r[metric] for r in rows if r['dataset']==d['name'] and r['arm']=='mapped'}
        diff=[a[i]-b[i] for i in range(6)];paired.append(dict(dataset=d['name'],metric=metric,mean_saving_ms=statistics.mean(diff),exploratory_95pct_mean_saving_ms=interval(diff),wins=sum(x>0 for x in diff)))
write(out/'paired.json',paired)
write(out/'manifest.json',{str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file() and p.name!='manifest.json'})
print(json.dumps(summary,indent=2))
