#!/usr/bin/env python3
"""Exploratory indexed break-even follow-up, chosen after the first rerun."""
import json,random,sys
from pathlib import Path
import statistics as st
from earliest_claim_study import run,save,digest,check,expected,verify,CPU

def main(base):
    verify(base);out=base/'indexed-followup';out.mkdir()
    records={};key=None
    for line in (base/'single.fa').read_text().splitlines():
        if line.startswith('>'):key=line[1:];records[key]=''
        else:records[key]+=line
    if len(records)!=50000 or any(len(s)!=1000 for s in records.values()):raise ValueError('input')
    plan={}
    for count in (1000,5000,10000,50000):
        for session in range(6):
            ids=random.Random(52026+session).choices(list(records),k=count);p=out/f'n{count}-s{session}.ids';p.write_text('\n'.join(ids)+'\n')
            plan[f'{count}-{session}']={'ids':str(p),'expected':expected(records,ids)}
    save(out/'plan.json',{'classification':'Exploratory follow-up chosen after initial64-query control; not an independent confirmatory test','counts':[1000,5000,10000,50000],'sessions':6,'jobs':plan})
    paths=[Path(__file__).resolve(),base/'single.fa',base/'single.fai',base/'O3/verified',base/'seal.json',out/'plan.json']+[Path(j['ids']) for j in plan.values()]
    hashes={str(p):digest(p) for p in paths};save(out/'seal.json',hashes)
    rows=[]
    for count in (1000,5000,10000,50000):
        for session in range(6):
            job=plan[f'{count}-{session}']
            for arm in (('cache','indexed') if session%2==0 else ('indexed','cache')):
                r=json.loads(run(['taskset','-c',CPU,str(base/'O3/verified'),arm,str(base/'single.fa'),job['ids'],str(base/'single.fai')],out,f'n{count}-s{session}-{arm}'))
                check(r,job['expected']);rows.append(dict(r,n=count,session=session,arm=arm))
        save(out/'runs.json',rows)
    for p,h in hashes.items():
        if digest(Path(p))!=h:raise ValueError('changed artifact')
    summary=[]
    for count in (1000,5000,10000,50000):
        groups={arm:sorted([r for r in rows if r['n']==count and r['arm']==arm],key=lambda r:r['session']) for arm in ('cache','indexed')}
        summary.append({'requests':count,'cache_median_total_ms':st.median(r['total_ms'] for r in groups['cache']),'indexed_median_total_ms':st.median(r['total_ms'] for r in groups['indexed']),'cache_median_query_ms':st.median(r['query_ms'] for r in groups['cache']),'indexed_median_query_ms':st.median(r['query_ms'] for r in groups['indexed']),'cache_wins':sum(c['total_ms']<i['total_ms'] for c,i in zip(groups['cache'],groups['indexed'])),'pairs':6})
    save(out/'summary.json',summary);print(json.dumps(summary,indent=2));verify(base)
if __name__=='__main__':main(Path(sys.argv[1]).resolve())
