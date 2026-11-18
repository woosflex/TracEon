#!/usr/bin/env python3
"""Reconstruct the first reported benchmark; never label adapted timings exact history."""
import hashlib,json,os,platform,random,re,struct,subprocess,sys,time
from pathlib import Path
import statistics as st
ROOT=Path(__file__).resolve().parents[1]
HTS=ROOT/'build/research-inputs/samtools-1.24/htslib-1.24'
SAM=HTS.parent/'samtools'
CPU='8'
def save(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(cmd,out,label):
    r=subprocess.run(cmd,capture_output=True,text=True)
    (out/(label+'.stdout')).write_text(r.stdout);(out/(label+'.stderr')).write_text(r.stderr)
    if r.returncode:raise RuntimeError(f'{label}: exit {r.returncode}; see stderr')
    return r.stdout

def expected(records,ids):
    h=hashlib.sha256()
    for k in ids:
        for v in (k,records[k]):
            b=v.encode();h.update(struct.pack('<Q',len(b)));h.update(b)
    return {'records':len(ids),'bytes':sum(len(records[k]) for k in ids),'sha256':h.hexdigest()}
def check(row,wanted):
    for k,v in wanted.items():
        if row.get(k)!=v:raise ValueError(f'wrong {k}')
def compile_all(out):
    src=out/'fa95f1d';commands=[]
    for opt in ('O0','O3'):
        d=out/opt;d.mkdir()
        flags=['c++','-std=c++20','-'+opt,'-include','cstdint','-I'+str(src)]+(['-DNDEBUG'] if opt=='O3' else [])
        core=[str(src/n) for n in ('Cache.cpp','PlainTextStrategy.cpp','TwoBitDnaStrategy.cpp','IEncodingStrategy.cpp')]
        for name in ('benchmark_disk','benchmark_traceon','verified'):
            if name=='verified':
                cmd=flags+['-I'+str(HTS),str(ROOT/'benchmarks/earliest_claim_worker.cpp')]+core+[str(HTS/'libhts.a'),'-lz','-lbz2','-llzma','-ldeflate','-lcurl','-lcrypto','-lpthread','-lm']
            else:cmd=flags+[str(src/(name+'.cpp'))]+(core if name=='benchmark_traceon' else [])
            cmd+=['-o',str(d/name)];commands.append(cmd);run(cmd,d,'compile-'+name)
    save(out/'compile-commands.json',commands)
def prepare(out):
    # Preserve first commit's missing-implementation failure after resolving header compatibility.
    src=out/'499d7f8';cmd=['c++','-std=c++20','-O0','-include','cstdint','-I'+str(src)]+[str(src/n) for n in ('benchmark_traceon.cpp','Cache.cpp','PlainTextStrategy.cpp','TwoBitDnaStrategy.cpp')]+['-o',str(out/'original-compat-attempt')]
    r=subprocess.run(cmd,capture_output=True,text=True);save(out/'original-compat-build.json',{'command':cmd,'returncode':r.returncode,'stderr':r.stderr})
    for n in ('Cache.cpp','benchmark_disk.cpp','benchmark_traceon.cpp'):
        if (src/n).read_bytes()!=(out/'fa95f1d'/n).read_bytes():raise ValueError('historical algorithm changed')
    compile_all(out)
    rng=random.Random(20260922);translation=bytes('ACGT'[i%4].encode()[0] for i in range(256));records={f'seq_{i:05}':rng.randbytes(1000).translate(translation).decode() for i in range(50000)}
    for layout,width in [('single',1000),('wrapped',80)]:
        p=out/f'{layout}.fa'
        with p.open('w') as f:
            for k,s in records.items():f.write('>'+k+'\n'+'\n'.join(s[i:i+width] for i in range(0,len(s),width))+'\n')
    ids=out/'all.ids';ids.write_text('\n'.join(records)+'\n');save(out/'all-expected.json',expected(records,list(records)))
    for session in range(6):
        keys=random.Random(1300+session).sample(list(records),64);p=out/f's{session}.ids';p.write_text('\n'.join(keys)+'\n');save(out/f's{session}-expected.json',expected(records,keys))
        for n in (1,2,3,4,8,16):
            p=out/f's{session}-n{n}.ids';p.write_text('\n'.join(keys[:n])+'\n');save(out/f's{session}-n{n}-expected.json',expected(records,keys[:n]))
    # Small smoke with the same algorithms, independently predicted payloads.
    small={'a':'aCGTN','b':'ACGT'*23};p=out/'smoke.fa';p.write_text(''.join('>'+k+'\n'+'\n'.join(s[i:i+7] for i in range(0,len(s),7))+'\n' for k,s in small.items()));q=out/'smoke.ids';keys=['a','b','a'];q.write_text('\n'.join(keys)+'\n')
    run([str(SAM),'faidx','--fai-idx',str(out/'smoke.fai'),str(p)],out,'smoke-index');smoke=[]
    for opt in ('O0','O3'):
        for arm in ('scan','cache','indexed','roundtrip'):
            result=json.loads(run(['taskset','-c',CPU,str(out/opt/'verified'),arm,str(p),str(q),str(out/('smoke.fai' if arm=='indexed' else 'smoke.cache'))],out,f'smoke-{opt}-{arm}'));check(result,expected(small,keys));smoke.append(result)
    save(out/'smoke.json',smoke)
    started=time.perf_counter();run([str(SAM),'faidx','--fai-idx',str(out/'single.fai'),str(out/'single.fa')],out,'index-build');save(out/'index-build.json',{'seconds':time.perf_counter()-started})
    paths=[p for p in out.rglob('*') if p.is_file() and p.suffix not in ('.stdout','.stderr')]
    paths += [ROOT/'benchmarks/earliest_claim_study.py',ROOT/'benchmarks/earliest_claim_worker.cpp',ROOT/'docs/earliest-claims-protocol.md',HTS/'libhts.a',SAM]
    save(out/'seal.json',{str(p):digest(p) for p in paths})
    save(out/'environment.json',{'platform':platform.platform(),'python':sys.version,'compiler':subprocess.check_output(['c++','--version'],text=True),'cpu':CPU,'cpuinfo':Path('/proc/cpuinfo').read_text(),'memory':Path('/proc/meminfo').read_text(),'generator':'Python Random20260922 randbytes translated modulo4; 50k unique IDs×1000bases; equivalent single/80-column layout'})
def verify(out):
    for p,h in json.loads((out/'seal.json').read_text()).items():
        if digest(Path(p))!=h:raise ValueError('changed sealed artifact '+p)
def execute(out):
    verify(out);raw=out/'raw';raw.mkdir()
    # Full-content validation precedes benchmark results.
    correctness=[]
    for opt in ('O0','O3'):
        for layout in ('single','wrapped'):
            for arm in ('cache','roundtrip'):
                row=json.loads(run(['taskset','-c',CPU,str(out/opt/'verified'),arm,str(out/f'{layout}.fa'),str(out/'all.ids'),str(raw/'roundtrip.cache')],raw,f'full-{opt}-{layout}-{arm}'))
                check(row,json.loads((out/'all-expected.json').read_text()))
                if row['cache_count']!=50000:raise ValueError('cache count')
                correctness.append(dict(row,opt=opt,layout=layout,arm=arm))
    save(out/'correctness.json',correctness)
    legacy=[]
    for opt in ('O0','O3'):
        for layout in ('single','wrapped'):
            for key in ('seq_00000','seq_25000','seq_49999'):
                for rep in range(6):
                    row={'opt':opt,'layout':layout,'key':key,'session':rep}
                    for arm in (('disk','traceon') if rep%2==0 else ('traceon','disk')):
                        text=run(['taskset','-c',CPU,str(out/opt/f'benchmark_{arm}'),str(out/f'{layout}.fa'),key],raw,f'legacy-{opt}-{layout}-{key}-{rep}-{arm}')
                        values=[float(v) for v in re.findall(r'(?:took:|load time:)\s*([0-9.eE+-]+)\s*ms',text)]
                        if len(values)!=(1 if arm=='disk' else 2):raise ValueError('legacy parse')
                        if arm=='disk':row['scan_ms']=values[0]
                        else:row.update(load_ms=values[0],lookup_ms=values[1])
                    legacy.append(row)
    save(out/'legacy.json',legacy)
    comparisons=[]
    for session in range(6):
        arms=['scan','cache','indexed'];arms=arms if session%2==0 else list(reversed(arms))
        for arm in arms:
            row=json.loads(run(['taskset','-c',CPU,str(out/'O3/verified'),arm,str(out/'single.fa'),str(out/f's{session}.ids'),str(out/'single.fai')],raw,f'control-{session}-{arm}'));check(row,json.loads((out/f's{session}-expected.json').read_text()));comparisons.append(dict(row,session=session,arm=arm))
    save(out/'controls.json',comparisons)
    amortization=[]
    for opt in ('O0','O3'):
        for n in (1,2,3,4,8,16):
            for session in range(6):
                for arm in (('scan','cache') if session%2==0 else ('cache','scan')):
                    row=json.loads(run(['taskset','-c',CPU,str(out/opt/'verified'),arm,str(out/'single.fa'),str(out/f's{session}-n{n}.ids'),str(out/'single.fai')],raw,f'amort-{opt}-{n}-{session}-{arm}'));check(row,json.loads((out/f's{session}-n{n}-expected.json').read_text()));amortization.append(dict(row,opt=opt,n=n,session=session,arm=arm))
    save(out/'amortization.json',amortization);verify(out)
    print('Completed:72 legacy pairs,18 comparator cohorts,144 amortization arms; outputs verified')
def summarize(out):
    legacy=json.loads((out/'legacy.json').read_text());result={'legacy':[],'controls':[],'amortization':[]}
    for opt in ('O0','O3'):
        for layout in ('single','wrapped'):
            g=[r for r in legacy if r['opt']==opt and r['layout']==layout];a={'opt':opt,'layout':layout}
            for k in ('scan_ms','load_ms','lookup_ms'):a[k]={'median':st.median(r[k] for r in g),'range':[min(r[k] for r in g),max(r[k] for r in g)]}
            a['median_paired_scan_lookup_ratio']=st.median(r['scan_ms']/r['lookup_ms'] for r in g);result['legacy'].append(a)
    controls=json.loads((out/'controls.json').read_text())
    for arm in ('scan','cache','indexed'):
        g=[r for r in controls if r['arm']==arm];result['controls'].append({'arm':arm,'median_setup_ms':st.median(r['setup_ms'] for r in g),'median_query_ms':st.median(r['query_ms'] for r in g),'query_range_ms':[min(r['query_ms'] for r in g),max(r['query_ms'] for r in g)]})
    amort=json.loads((out/'amortization.json').read_text())
    for opt in ('O0','O3'):
        for n in (1,2,3,4,8,16):
            g={arm:sorted([r for r in amort if r['opt']==opt and r['n']==n and r['arm']==arm],key=lambda r:r['session']) for arm in ('scan','cache')}
            result['amortization'].append({'opt':opt,'n':n,'scan_median_ms':st.median(r['total_ms'] for r in g['scan']),'cache_median_ms':st.median(r['total_ms'] for r in g['cache']),'cache_wins':sum(a['total_ms']>b['total_ms'] for a,b in zip(g['scan'],g['cache'])),'pairs':6})
    save(out/'summary.json',result);print(json.dumps(result,indent=2))
if __name__=='__main__':globals()[sys.argv[1]](Path(sys.argv[2]).resolve())
