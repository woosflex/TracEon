#!/usr/bin/env python3
"""Local pyfastx/compression/cgroup controls; see docs/broader-record-protocol.md."""
import gzip
import hashlib
import json
import os
from pathlib import Path
import random
import shutil
import statistics
import subprocess
import sys
import time
from resident_record_study import (ROOT, BUILD, CPUS, PREFIX, SAMTOOLS, Worker,
                                  cohort, expected, save, validate)
from indexed_mapping_study import sequence_records
from mapping_baseline import sha256, run_measured
from shared_mapping_pilot import start_server
from summarize_paired import interval

OLD = ROOT / 'build/resident-record-control-2026-09-22'
PYTHON = ROOT / 'build/broader-validation-venv/bin/python'
SITE = ROOT / 'build/broader-validation-venv/lib/python3.14/site-packages'
os.environ['PYTHONPATH'] = str(SITE)
ARMS = ['indexed', 'pyfastx-plain', 'pyfastx-gzip', 'remote']
CASES = [(1, 64, 'random'), (4, 1, 'random'), (4, 64, 'random'), (4, 64, 'ordered')]
# Three permutations and their reverses guarantee 3/3 pairwise precedence.
ORDERS = [ARMS, list(reversed(ARMS)), ARMS[1:]+ARMS[:1],
          list(reversed(ARMS[1:]+ARMS[:1])), ARMS[2:]+ARMS[:2], list(reversed(ARMS[2:]+ARMS[:2]))]


def index_command(reads, fxi, fastq):
    code = "import pyfastx,sys; k={'index_file':sys.argv[2]}; k.update({} if sys.argv[3]=='1' else {'uppercase':False}); x=(pyfastx.Fastq if sys.argv[3]=='1' else pyfastx.Fasta)(sys.argv[1],**k); print(len(x))"
    return PREFIX + [str(PYTHON), '-c', code, str(reads), str(fxi), str(int(fastq))]


def compressed(source, target):
    start = time.perf_counter()
    with source.open('rb') as i, target.open('wb') as o, gzip.GzipFile(fileobj=o, mode='wb', compresslevel=6, mtime=0) as z:
        shutil.copyfileobj(i, z, 1024*1024)
    elapsed = time.perf_counter()-start
    h = hashlib.sha256()
    with gzip.open(target, 'rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''): h.update(block)
    if h.hexdigest() != sha256(source): raise ValueError('gzip round trip')
    return {'seconds': elapsed, 'bytes': target.stat().st_size, 'sha256': sha256(target)}


def prepare(out):
    original = json.loads((OLD/'manifest.json').read_text())
    manifest = {'datasets': {}, 'python': sys.version, 'cases': CASES, 'orders': ORDERS,
                'prior_manifest_sha256': sha256(OLD/'manifest.json')}
    for name in ('drosophila', 'synthetic1m', 'long'):
        d = dict(original['datasets'][name]); dest=out/name; dest.mkdir()
        src=Path(d['reads']); z=dest/'source.gz'
        d['gzip_generation']=compressed(src,z); d['gzip']=str(z)
        d['fxi']={}; d['new_builds']=[]
        for form, source in [('plain',src), ('gzip',z)]:
            for trial in range(3):
                index=dest/f'{form}-{trial}.fxi'
                m=run_measured(index_command(source,index,d['fastq']),dest,f'{form}-index-{trial}')
                m.update(arm=f'pyfastx-{form}', trial=trial, path=str(index), sha256=sha256(index), bytes=index.stat().st_size)
                d['new_builds'].append(m)
                if trial==0:d['fxi'][form]=str(index)
        for trial in range(3):
            cache=dest/f'gzip-cache-{trial}'
            m=run_measured(PREFIX+[str(BUILD/'traceon'),'build',str(z),str(cache)],dest,f'gzip-cache-build-{trial}')
            m.update(arm='gzip-cache',trial=trial,path=str(cache),sha256=sha256(cache),bytes=cache.stat().st_size)
            d['new_builds'].append(m)
            if sha256(cache)!=sha256(Path(d['cache'])):raise ValueError('plain/gzip cache artifact differs')
        jobs={}; needed=set(); n=512 if name=='long' else 1024
        for session in range(6):
            for pattern in ('random','ordered'):
                group=[]
                for worker,oldjob in enumerate(d['jobs'][f'{session}-{pattern}'][:4]):
                    ids=Path(oldjob['ids']).read_text().splitlines()[:n]; needed.update(ids)
                    path=dest/f's{session}-{pattern}-w{worker}.ids';path.write_text('\n'.join(ids)+'\n')
                    group.append({'ids':str(path), '_ids':ids})
                jobs[f'{session}-{pattern}']=group
        records={k:(s,q) for k,s,q in sequence_records(src,d['fastq']) if k in needed}
        for group in jobs.values():
            for j in group:j['expected']=expected(records,j.pop('_ids'))
        d['jobs']=jobs;manifest['datasets'][name]=d
        print('prepared',name,flush=True)
    save(out/'manifest.json',manifest)


def command(out, arm, d, port=1, snapshot='0'*64):
    py=arm.startswith('pyfastx')
    form=arm.split('-')[-1] if py else 'plain'
    reads=d['gzip'] if form=='gzip' else d['reads']
    index=d['fxi'][form] if py else d['fai']
    return PREFIX+[str(out/'worker'),'pyfastx' if py else arm,reads,index,str(int(d['fastq'])),str(port),snapshot]


def measure(out,d,dest,arm,case,jobs,server=None,port=1,setup=None):
    n,b,_=case;workers=[]
    try:
        for w in range(n):workers.append(Worker(command(out,arm,d,port,(setup or {}).get('snapshot','0'*64)),dest/f'{arm}-w{w}.stderr'))
        passes={}
        for phase in ('warmup','timing','memory'):
            passes[phase]=cohort(workers,jobs[:n],b,dest,f'{arm}-{phase}',server.pid if arm=='remote' else None,phase=='memory')
        return {'arm':arm,'case':case,'worker_setup':[w.setup for w in workers],'passes':passes}
    finally:
        for w in workers:w.close()


def smoke(out):
    dest=out/'smoke';dest.mkdir();checks=[]
    for fastq in (False,True):
        sub=dest/str(int(fastq));sub.mkdir()
        records={'id:1/2':('aCGTN','!~ABC' if fastq else ''),'last':('tT','I#' if fastq else '')}
        reads=sub/'source';reads.write_text(''.join(f'@{k} description\n{s}\n+\n{q}\n' if fastq else f'>{k} description\n{s}\n' for k,(s,q) in records.items()))
        compressed(reads,sub/'source.gz')
        d={'reads':str(reads),'gzip':str(sub/'source.gz'),'fastq':fastq,'fai':str(sub/'fai'),'cache':str(sub/'cache'),'fxi':{f:str(sub/f'{f}.fxi') for f in ('plain','gzip')}}
        run_measured([str(BUILD/'traceon'),'build',d['reads'],d['cache']],sub,'cache')
        run_measured([str(SAMTOOLS),'fqidx' if fastq else 'faidx','--fai-idx',d['fai'],d['reads']],sub,'fai')
        for form in ('plain','gzip'):run_measured(index_command(Path(d['gzip'] if form=='gzip' else d['reads']),Path(d['fxi'][form]),fastq),sub,form)
        ids=['id:1/2','last','id:1/2']*87+['last'];p=sub/'ids';p.write_text('\n'.join(ids)+'\n');job={'ids':str(p),'expected':expected(records,ids)}
        server,port,setup=start_server(BUILD/'traceon',Path(d['cache']),sub,supported_cli=True,cpus=CPUS)
        try:
            for arm in ARMS:
                for b in (1,16,64,256):
                    target=sub/f'{arm}-{b}';target.mkdir()
                    checks.append(measure(out,d,target,arm,(1,b,'fixture'),[job],server,port,setup))
        finally:server.terminate();server.wait(timeout=10)
    good=checks[-1]['passes']['timing']['workers'][0]
    validate(good,job['expected'])
    bad=dict(good,sha256='0'*64)
    try:validate(bad,job['expected'])
    except ValueError:pass
    else:raise ValueError('corruption accepted')
    save(dest/'passed.json',{'checks':checks,'corruption_rejected':True})


def seal(out):
    if (out/'seal.json').exists():raise ValueError('already sealed')
    if not json.loads((out/'smoke/passed.json').read_text())['corruption_rejected']:raise ValueError('smoke')
    # Match ordered queries to exactly the random selection, before sealing.
    manifest=json.loads((out/'manifest.json').read_text())
    shutil.copy2(out/'manifest.json',out/'preparation-manifest.json')
    for d in manifest['datasets'].values():
        needed={key for session in range(6) for j in d['jobs'][f'{session}-random'] for key in Path(j['ids']).read_text().splitlines()}
        records={};ranks={}
        for rank,(key,seq,qual) in enumerate(sequence_records(Path(d['reads']),d['fastq'])):
            if key in needed:records[key]=(seq,qual);ranks[key]=rank
        for session in range(6):
            for old,new in zip(d['jobs'][f'{session}-random'],d['jobs'][f'{session}-ordered']):
                ids=sorted(Path(old['ids']).read_text().splitlines(),key=ranks.__getitem__)
                Path(new['ids']).write_text('\n'.join(ids)+'\n');new['expected']=expected(records,ids)
    save(out/'manifest.json',manifest)
    paths={out/'worker',out/'manifest.json',BUILD/'traceon',out/'smoke/passed.json'}
    for name in ('broader_record_study.py','broader_record_worker.cpp','pyfastx_missing_probe.py','summarize_broader_records.py','resident_record_study.py','indexed_mapping_study.py','mapping_baseline.py','shared_mapping_pilot.py','subset_mapping_pilot.py','summarize_paired.py'):
        p=ROOT/'benchmarks'/name;paths.add(p);target=out/'source'/name;target.parent.mkdir(exist_ok=True);shutil.copy2(p,target)
    paths.add(ROOT/'docs/broader-record-protocol.md')
    shutil.copy2(ROOT/'docs/broader-record-protocol.md',out/'source/broader-record-protocol.md')
    for p in SITE.glob('pyfastx*.so'):paths.add(p)
    paths.add(out/'missing-probe.json')
    paths.update((ROOT/'build/research-inputs/pyfastx-wheels').glob('*.whl'))
    for d in manifest['datasets'].values():
        if sha256(Path(d['reads'])) != d['source_sha256']:raise ValueError('original input changed')
        paths.update(Path(d[k]) for k in ('reads','gzip','cache','fai'));paths.update(map(Path,d['fxi'].values()))
        paths.update(Path(j['ids']) for group in d['jobs'].values() for j in group)
    save(out/'seal.json',{str(p):sha256(p) for p in sorted(paths)})


def verify(out):
    for p,h in json.loads((out/'seal.json').read_text()).items():
        if sha256(Path(p))!=h:raise ValueError(f'changed sealed artifact: {p}')


def run(out):
    verify(out);os.sched_setaffinity(0,{16});m=json.loads((out/'manifest.json').read_text())
    with (out/'runs.jsonl').open('x') as stream:
        for name,d in m['datasets'].items():
            for session in range(6):
                dest=out/name/f'session-{session}';dest.mkdir()
                save(dest/'host-before.json',{k:Path('/proc'/Path(k)).read_text() for k in ('meminfo','vmstat','loadavg')})
                server,port,setup=start_server(BUILD/'traceon',Path(d['cache']),dest,supported_cli=True,cpus=CPUS);save(dest/'server.json',setup)
                cases=list(enumerate(CASES));random.Random(92026+session).shuffle(cases)
                try:
                    for ci,case in cases:
                        target=dest/f'case-{ci}';target.mkdir()
                        for arm in ORDERS[session]:
                            try:
                                result=measure(out,d,target,arm,case,d['jobs'][f'{session}-{case[2]}'],server,port,setup)
                            except Exception as e:
                                stream.write(json.dumps({'dataset':name,'session':session,'case':case,'arm':arm,'status':'failed','error':repr(e)})+'\n');stream.flush();raise
                            result.update(dataset=name,session=session);stream.write(json.dumps(result)+'\n');stream.flush()
                    print('completed',name,session,flush=True)
                finally:
                    server.terminate();server.wait(timeout=10)
                    save(dest/'host-after.json',{k:Path('/proc'/Path(k)).read_text() for k in ('meminfo','vmstat','loadavg')})
    verify(out)


def cg_snapshot():
    c=Path('/sys/fs/cgroup')/Path('/proc/self/cgroup').read_text().strip().split('::')[1].lstrip('/')
    return {n:(c/n).read_text() for n in ('memory.max','memory.swap.max','memory.current','memory.swap.current','memory.peak','memory.events','memory.stat','memory.pressure')}


def pressure_one(out,dest,arm,session):
    os.sched_setaffinity(0,{16});d=dict(json.loads((out/'manifest.json').read_text())['datasets']['synthetic1m']);d['fxi']=dict(d['fxi'])
    server=None;local=dest/'local';local.mkdir();save(dest/'before.json',cg_snapshot())
    try:
        keys=['cache'] if arm=='remote' else ['reads','fai'] if arm=='indexed' else ['gzip' if arm.endswith('gzip') else 'reads']
        for key in keys:
            target=local/key
            with Path(d[key]).open('rb') as src,target.open('wb') as dst:shutil.copyfileobj(src,dst,1024*1024)
            d[key]=str(target)
        if arm.startswith('pyfastx'):
            form=arm.split('-')[1];target=local/'index.fxi'
            with Path(d['fxi'][form]).open('rb') as src,target.open('wb') as dst:shutil.copyfileobj(src,dst,1024*1024)
            d['fxi'][form]=str(target)
        save(dest/'after-copy.json',cg_snapshot());port=1;setup=None
        if arm=='remote':server,port,setup=start_server(BUILD/'traceon',Path(d['cache']),dest,supported_cli=True,cpus=CPUS)
        r=measure(out,d,dest,arm,(4,64,'random'),d['jobs'][f'{session}-random'],server,port,setup)
        r.update(status='passed',server_setup=setup);save(dest/'result.json',r)
    except Exception as e:save(dest/'result.json',{'status':'failed','error':repr(e),'arm':arm})
    finally:
        save(dest/'after.json',cg_snapshot())
        if server and server.poll() is None:server.terminate();server.wait(timeout=10)
        shutil.rmtree(local)


def pressure(out):
    verify(out);base=out/'pressure';base.mkdir()
    for limit in (256,768,1536):
        for session in range(6):
            for arm in ORDERS[session]:
                dest=base/f'm{limit}-s{session}-{arm}';dest.mkdir()
                unit=f'traceon-pressure-{os.getpid()}-{limit}-{session}-{arm}'
                cmd=['systemd-run','--user','--quiet','--wait','--pipe',f'--unit={unit}',f'--property=MemoryMax={limit}M','--property=MemorySwapMax=0','--property=OOMPolicy=continue','--property=RuntimeMaxSec=180',sys.executable,str(Path(__file__).resolve()),'pressure-one',str(out),str(dest),arm,str(session)]
                start=time.perf_counter()
                try:r=subprocess.run(cmd,capture_output=True,text=True,timeout=200)
                except subprocess.TimeoutExpired as e:
                    subprocess.run(['systemctl','--user','stop',unit],capture_output=True,timeout=20)
                    r=subprocess.CompletedProcess(cmd,124,'',repr(e))
                diagnostic=subprocess.run(['systemctl','--user','show',unit,'-p','Result','-p','ExecMainStatus','-p','MemoryPeak','-p','OOMPolicy'],capture_output=True,text=True)
                (dest/'unit-diagnostics.txt').write_text(diagnostic.stdout+diagnostic.stderr)
                (dest/'unit.stdout').write_text(r.stdout);(dest/'unit.stderr').write_text(r.stderr)
                save(dest/'unit.json',{'command':cmd,'returncode':r.returncode,'seconds':time.perf_counter()-start})
                # A killed controller is also an observed failure, never a dropped arm.
                if not (dest/'result.json').exists():save(dest/'result.json',{'status':'unit-failed','arm':arm})
            print('pressure completed',limit,session,flush=True)
    verify(out)


if __name__=='__main__':
    action=sys.argv[1];out=Path(sys.argv[2]).resolve()
    if action=='pressure-one':pressure_one(out,Path(sys.argv[3]),sys.argv[4],int(sys.argv[5]))
    else:globals()[action](out)
