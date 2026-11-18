#!/usr/bin/env python3
"""Bounded local v1 mmap replication. Run from repository root; outputs are uncommitted."""
import hashlib, json, os, platform, random, shutil, statistics, struct, subprocess, time
from pathlib import Path
ROOT=Path.cwd(); OUT=ROOT/'build/mapped-cache-2026-09-22'; OUT.mkdir(parents=True,exist_ok=True)
OLD=OUT/'v1'; BUILD=Path('/tmp/traceon-release-check-20260922'); HTS=ROOT/'build/research-inputs/samtools-1.24/htslib-1.24'
def write(path,obj): path.write_text(json.dumps(obj,indent=2)+'\n')
def sha(path): return hashlib.file_digest(path.open('rb'),'sha256').hexdigest()
def run(cmd):
    p=subprocess.run(list(map(str,cmd)),capture_output=True,text=True)
    with (OUT/'commands.jsonl').open('a') as f: f.write(json.dumps(dict(command=list(map(str,cmd)),returncode=p.returncode,stdout=p.stdout,stderr=p.stderr))+'\n')
    p.check_returncode();return p.stdout
for path in run(['git','ls-tree','-r','--name-only','v1.0.0']).splitlines():
    if path.startswith('include/') or path in ['src/Cache.cpp','src/SmartStrategy.cpp','src/IEncodingStrategy.cpp','third_party/robin_hood.h']:
        target=OLD/path;target.parent.mkdir(parents=True,exist_ok=True)
        target.write_bytes(subprocess.check_output(['git','show','v1.0.0:'+path]))
assert (OLD/'third_party/robin_hood.h').exists()
worker=ROOT/'benchmarks/mapped_cache_worker.cpp'
common=['c++','-std=c++20','-O3','-DNDEBUG','-include','cstdint','-pthread','-I'+str(HTS)]
libs=[str(HTS/'libhts.a'),'-lz','-lbz2','-llzma','-ldeflate','-lcurl','-lcrypto','-lm']
run(common+['-I'+str(OLD),'-I'+str(OLD/'include'),worker]+[OLD/'src'/s for s in ['Cache.cpp','SmartStrategy.cpp','IEncodingStrategy.cpp']]+libs+['-o',OUT/'v1-worker'])
run(common+['-Iinclude','-Ibuild/_deps/ankerl_unordered_dense-src/include',worker,BUILD/'libtraceon_core.a',BUILD/'_deps/zlib-ng-build/libz.a',BUILD/'third_party/lz4/build/cmake/liblz4.a']+libs+['-o',OUT/'current-worker'])
write(OUT/'environment.json',dict(platform=platform.platform(),python=platform.python_version(),cpu=Path('/proc/cpuinfo').read_text(),affinity=sorted(os.sched_getaffinity(0)),compiler=run(['c++','--version']),v1=run(['git','rev-parse','v1.0.0^{}']).strip(),head=run(['git','rev-parse','HEAD']).strip(),current_diff_sha256=hashlib.sha256(subprocess.check_output(['git','diff'])).hexdigest(),libraries={str(p):sha(p) for p in [BUILD/'libtraceon_core.a',HTS/'libhts.a']}))
shutil.copy2(worker,OUT/worker.name);shutil.copy2(__file__,OUT/Path(__file__).name)
def field(h,s): h.update(struct.pack('<Q',len(s)));h.update(s)
def memory(pid):
    data={k:int(v.split()[0]) for k,v in (line.split(':',1) for line in Path(f'/proc/{pid}/smaps_rollup').read_text().splitlines()[1:])}
    return {k:data[k] for k in ['Rss','Pss','Pss_Anon','Pss_File','Private_Clean','Private_Dirty','Shared_Clean','Shared_Dirty']}
rows=[]; memrows=[]; datasets=[]
for length,mb in [(150,10),(150,100),(30000,10),(30000,100)]:
    name=f'len{length}-{mb}MB'; folder=OUT/name;folder.mkdir(exist_ok=True)
    count=mb*1000000//length; fasta=folder/'reads.fa'; idsfile=folder/'queries.txt'
    rng=random.Random(20260922); sequences=[]; full=hashlib.sha256()
    translation=bytes.maketrans(bytes(range(256)),bytes(b'ACGT'[i%4] for i in range(256)))
    with fasta.open('wb') as f:
        for i in range(count):
            seq=rng.randbytes(length).translate(translation);id=f'read{i}'.encode();sequences.append(seq)
            f.write(b'>'+id+b'\n'+seq+b'\n');field(full,id);field(full,seq)
    queries=[rng.randrange(count) for _ in range(200000)];idsfile.write_text(''.join(f'read{i}\n' for i in queries))
    expected=hashlib.sha256()
    for i in queries[:10000]:field(expected,f'read{i}'.encode());field(expected,sequences[i])
    del sequences,queries
    meta=dict(name=name,count=count,length=length,input_bytes=fasta.stat().st_size,input_sha256=sha(fasta),queries_sha256=sha(idsfile),full_sha256=full.hexdigest(),query_sha256=expected.hexdigest(),builds={})
    for arm in ['v1','current','indexed']:
        exe=OUT/('v1-worker' if arm=='v1' else 'current-worker'); kind='indexed' if arm=='indexed' else 'cache'
        cache=folder/(arm+'.cache'); cmd=['taskset','-c','8-11',exe,'build',kind,fasta,cache,count,length]
        meta['builds'][arm]=json.loads(run(cmd));meta['builds'][arm]['bytes']=(Path(str(fasta)+'.fai') if arm=='indexed' else cache).stat().st_size
    datasets.append(meta);write(OUT/'datasets.json',datasets)
    for rep in range(6):
        arms=['v1','current','indexed'];random.Random(500+rep).shuffle(arms)
        for arm in arms:
            exe=OUT/('v1-worker' if arm=='v1' else 'current-worker');kind='indexed' if arm=='indexed' else 'cache'; inp=fasta if arm=='indexed' else folder/(arm+'.cache')
            result=json.loads(run(['taskset','-c','8',exe,'query',kind,inp,idsfile,count,length]))
            assert result['sha256']==meta['query_sha256'] and result['records']==10000
            if arm!='indexed':assert result['length_sum']==1000000*length
            rows.append(dict(dataset=name,rep=rep,arm=arm,**result));write(OUT/'query-results.json',rows)
        # Fresh execs; no inherited cache. All payload bytes read before PSS sampling.
        for nproc in [1,4]:
            for arm in arms:
                processes=[];results=[]
                try:
                    for i in range(nproc):
                        exe=OUT/('v1-worker' if arm=='v1' else 'current-worker');kind='indexed' if arm=='indexed' else 'cache';inp=fasta if arm=='indexed' else folder/(arm+'.cache')
                        p=subprocess.Popen(list(map(str,['taskset','-c',8+i,exe,'hold',kind,inp,'-',count,length])),stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True);processes.append(p)
                    for p in processes:
                        line=p.stdout.readline()
                        if not line:raise RuntimeError(p.stderr.read())
                        result=json.loads(line);assert result['sha256']==meta['full_sha256'] and result['records']==count;results.append(result)
                    samples=[memory(p.pid) for p in processes]
                    memrows.append(dict(dataset=name,rep=rep,arm=arm,nproc=nproc,workers=results,memory_kib=samples,total_kib={k:sum(s[k] for s in samples) for k in samples[0]}));write(OUT/'memory-results.json',memrows)
                finally:
                    for p in processes:
                        if p.poll() is None:p.stdin.write('quit\n');p.stdin.flush()
                    for p in processes:
                        _,err=p.communicate(timeout=30)
                        if p.returncode:raise RuntimeError(err)
    print(name+' complete',flush=True)
summary=[]
for d in datasets:
    for arm in ['v1','current','indexed']:
        r=[r for r in rows if r['dataset']==d['name'] and r['arm']==arm]
        summary.append(dict(dataset=d['name'],arm=arm,setup_ms=statistics.median(x['setup_ms'] for x in r),consume_ms=statistics.median(x['consume_ms'] for x in r),lookup_mops=statistics.median(1000/x['lookup_ms'] for x in r) if arm!='indexed' else None,pss_mib={n:statistics.median(x['total_kib']['Pss']/1024 for x in memrows if x['dataset']==d['name'] and x['arm']==arm and x['nproc']==n) for n in [1,4]}))
write(OUT/'summary.json',summary)
write(OUT/'manifest.json',{str(p.relative_to(OUT)):sha(p) for p in OUT.rglob('*') if p.is_file() and p.name!='manifest.json'})
print(json.dumps(summary,indent=2))
