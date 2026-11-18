#!/usr/bin/env python3
"""Local speed comparison using existing sealed inputs and the current mapped build.
Usage: python3 benchmarks/local_speed_study.py OUTPUT_DIRECTORY
"""
import gc,hashlib,json,os,platform,random,shutil,statistics,struct,subprocess,sys
from pathlib import Path
from summarize_paired import interval
ROOT=Path(__file__).resolve().parents[1];out=Path(sys.argv[1]).resolve();out.mkdir(parents=True,exist_ok=False)
prior=ROOT/'build/mapped-cache-2026-09-22';storage=ROOT/'build/mapped-storage-2026-09-22';build=Path('/tmp/traceon-mapped-check-20260922')
def write(path,value):path.write_text(json.dumps(value,indent=2)+'\n')
def sha(path):
    with path.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def run(cmd):
    cmd=list(map(str,cmd));p=subprocess.run(cmd,text=True,capture_output=True,timeout=120)
    with (out/'commands.jsonl').open('a') as f:f.write(json.dumps(dict(command=cmd,returncode=p.returncode,stdout=p.stdout,stderr=p.stderr))+'\n')
    p.check_returncode();return p.stdout
# Reuse the existing worker. The consuming variant adds only a start barrier and
# shared monotonic timestamps, so aggregate concurrent throughput can be measured.
worker=out/'lookup.cpp';shutil.copy2(ROOT/'benchmarks/mapped_cache_worker.cpp',worker)
source=(prior/'equal-warmup-worker.cpp').read_text()
needle='    begin=Clock::now();\n    for(size_t i=0;i<queries;++i)'
assert needle in source
source=source.replace(needle,'    std::cout<<"READY"<<std::endl; std::string go; if(!std::getline(std::cin,go)||go!="RUN")throw std::runtime_error("barrier");\n'+needle)
source=source.replace('double consume=elapsed(begin);','auto end=Clock::now(); double consume=std::chrono::duration<double,std::milli>(end-begin).count();')
needle='<<",\\\"consume_ms\\\":"<<consume'
assert needle in source
source=source.replace(needle,'<<",\\\"start_ns\\\":"<<std::chrono::duration_cast<std::chrono::nanoseconds>(begin.time_since_epoch()).count()<<",\\\"end_ns\\\":"<<std::chrono::duration_cast<std::chrono::nanoseconds>(end.time_since_epoch()).count()'+needle)
(out/'consume.cpp').write_text(source)
hts=ROOT/'build/research-inputs/samtools-1.24/htslib-1.24'
for name in ['lookup','consume']:
    run(['c++','-std=c++20','-O3','-DNDEBUG','-include','cstdint','-pthread','-I'+str(ROOT/'include'),'-I'+str(ROOT/'build/_deps/ankerl_unordered_dense-src/include'),'-I'+str(hts),out/(name+'.cpp'),build/'libtraceon_core.a',build/'_deps/zlib-ng-build/libz.a',build/'third_party/lz4/build/cmake/liblz4.a',hts/'libhts.a','-lz','-lbz2','-llzma','-ldeflate','-lcurl','-lcrypto','-lm','-o',out/name])
assert sha(build/'libtraceon_core.a')==json.loads((storage/'environment.json').read_text())['core_sha256']
write(out/'environment.json',dict(platform=platform.platform(),cpu=Path('/proc/cpuinfo').read_text(),affinity=sorted(os.sched_getaffinity(0)),core_sha256=sha(build/'libtraceon_core.a'),htslib_sha256=sha(hts/'libhts.a'),compiler=run(['c++','--version']),head=run(['git','-C',ROOT,'rev-parse','HEAD']).strip()))
shutil.copy2(__file__,out/Path(__file__).name);shutil.copy2(ROOT/'benchmarks/summarize_paired.py',out/'summarize_paired.py')
datasets=json.loads((storage/'inputs.json').read_text());write(out/'inputs.json',datasets)
sealed=json.loads((storage/'manifest.json').read_text());rows=[];lookups=[]
for d in datasets:
    name=d['name'];folder=out/name;folder.mkdir();fasta=prior/name/'reads.fa';ids=prior/name/'queries.txt'
    assert sha(fasta)==d['input_sha256'] and sha(ids)==d['queries_sha256']
    for arm in ['mapped','lz4']:
        f=storage/(name+'-'+arm+'.traceon');assert sha(f)==sealed[f.name]
    all_ids=ids.read_text().splitlines();jobs=[];records={}
    with fasta.open('rb') as f:
        while header:=f.readline():records[header[1:].strip().decode()]=f.readline().strip()
    for i in range(4):
        selected=all_ids[i*10000:(i+1)*10000];query=folder/f'ids-{i}.txt';query.write_text('\n'.join(selected)+'\n');digest=hashlib.sha256()
        for key in selected:
            for value in [key.encode(),records[key]]:digest.update(struct.pack('<Q',len(value)));digest.update(value)
        jobs.append(dict(path=str(query),sha256=digest.hexdigest(),ids_sha256=sha(query)))
    del records,all_ids;gc.collect();write(folder/'jobs.json',jobs)
    # Warm file pages consistently; no root cache dropping and no cold-cache claim.
    for f in [fasta,storage/(name+'-mapped.traceon'),storage/(name+'-lz4.traceon')]:sha(f)
    for rep in range(6):
        arms=['mapped','lz4'];random.Random(90+rep).shuffle(arms)
        for arm in arms:
            r=json.loads(run(['taskset','-c','8',out/'lookup','query','cache',storage/(name+'-'+arm+'.traceon'),ids,d['count'],d['length']]))
            assert r['length_sum']==1000000*d['length'] and r['sha256']==d['query_sha256']
            lookups.append(dict(dataset=name,rep=rep,arm=arm,**r));write(out/'lookups.json',lookups)
        for nproc in [1,4]:
            arms=['mapped','lz4','htslib'];random.Random(20260922+rep).shuffle(arms)
            for arm in arms:
                processes=[];results=[];commands=[]
                try:
                    for i in range(nproc):
                        path=fasta if arm=='htslib' else storage/(name+'-'+arm+'.traceon')
                        cmd=list(map(str,['taskset','-c',8+i,out/'consume','query','indexed' if arm=='htslib' else 'cache',path,jobs[i]['path'],d['count'],d['length']]))
                        commands.append(cmd);processes.append(subprocess.Popen(cmd,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True))
                    for p in processes:
                        line=p.stdout.readline()
                        if line!='READY\n':raise RuntimeError((line,p.stderr.read()))
                    for p in processes:p.stdin.write('RUN\n');p.stdin.flush()
                    for i,p in enumerate(processes):
                        stdout,stderr=p.communicate(timeout=60)
                        if p.returncode:raise RuntimeError(stderr)
                        r=json.loads(stdout);assert r['sha256']==jobs[i]['sha256'] and r['records']==10000;results.append(r)
                    start=min(r['start_ns'] for r in results);end=max(r['end_ns'] for r in results);ms=(end-start)/1e6
                    rows.append(dict(dataset=name,arm=arm,rep=rep,nproc=nproc,commands=commands,workers=results,cohort_ms=ms,records_per_second=nproc*10000/(ms/1000),start_skew_ms=(max(r['start_ns'] for r in results)-start)/1e6))
                    write(out/'cohorts.json',rows)
                finally:
                    for p in processes:
                        if p.poll() is None:p.kill();p.communicate()
    print(name+' complete',flush=True)
summary=[];paired=[]
for d in datasets:
    for nproc in [1,4]:
        for arm in ['mapped','lz4','htslib']:
            selected=[r for r in rows if r['dataset']==d['name'] and r['nproc']==nproc and r['arm']==arm]
            summary.append(dict(dataset=d['name'],nproc=nproc,arm=arm,median_ms=statistics.median(r['cohort_ms'] for r in selected),median_records_per_second=statistics.median(r['records_per_second'] for r in selected),median_setup_ms=statistics.median(w['setup_ms'] for r in selected for w in r['workers']),range_ms=[min(r['cohort_ms'] for r in selected),max(r['cohort_ms'] for r in selected)]))
        for baseline in ['lz4','htslib']:
            a={r['rep']:r['cohort_ms'] for r in rows if r['dataset']==d['name'] and r['nproc']==nproc and r['arm']==baseline};b={r['rep']:r['cohort_ms'] for r in rows if r['dataset']==d['name'] and r['nproc']==nproc and r['arm']=='mapped'}
            diff=[a[i]-b[i] for i in range(6)];paired.append(dict(dataset=d['name'],nproc=nproc,baseline=baseline,mean_saving_ms=statistics.mean(diff),exploratory_95pct_mean_saving_ms=interval(diff),wins=sum(x>0 for x in diff)))
write(out/'summary.json',summary);write(out/'paired.json',paired)
write(out/'manifest.json',{str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file() and p.name!='manifest.json'})
print(json.dumps(summary,indent=2))
