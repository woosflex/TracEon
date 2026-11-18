#!/usr/bin/env python3
"""Revalidate broader-study evidence and export all comparisons, including failures."""
import csv
import json
from pathlib import Path
import statistics as st
import sys
from broader_record_study import ARMS, CASES, validate, save
from summarize_paired import interval


def main(out):
    manifest=json.loads((out/'manifest.json').read_text())
    runs=[json.loads(x) for x in (out/'runs.jsonl').read_text().splitlines()]
    if len(runs)!=288:raise ValueError(f'expected 288 rows, got {len(runs)}')
    totals={'arm_observations':len(runs),'worker_processes':0,'record_instances_all_passes':0,'record_instances_timing':0}
    indexed={}
    def check(r,d,session):
        jobs=d['jobs'][f'{session}-{r["case"][2]}'][:r['case'][0]]
        for phase,p in r['passes'].items():
            if len(p['workers'])!=len(jobs):raise ValueError('worker count')
            for row,job in zip(p['workers'],jobs):validate(row,job['expected'])
            span=max(w['end'] for w in p['workers'])-min(w['start'] for w in p['workers'])
            if span!=p['active_seconds']:raise ValueError('span mismatch')
    for r in runs:
        key=(r['dataset'],tuple(r['case']),r['session'],r['arm'])
        if key in indexed:raise ValueError('duplicate cell')
        indexed[key]=r;check(r,manifest['datasets'][r['dataset']],r['session'])
        totals['worker_processes']+=len(r['worker_setup'])
        for phase,p in r['passes'].items():
            n=sum(w['records'] for w in p['workers']);totals['record_instances_all_passes']+=n
            if phase=='timing':totals['record_instances_timing']+=n
    results=[]
    for name in manifest['datasets']:
        for case in CASES:
            groups={arm:[indexed[name,case,s,arm] for s in range(6)] for arm in ARMS}
            row={'dataset':name,'case':case,'arms':{},'comparisons':{}}
            for arm,rows in groups.items():
                row['arms'][arm]={'median_ms':st.median(r['passes']['timing']['active_seconds'] for r in rows)*1000,
                    'median_pss_mib':st.median(r['passes']['memory']['memory_peak_pss_kib'] for r in rows)/1024}
            for arm in ARMS[:-1]:
                differences=[(a['passes']['timing']['active_seconds']-b['passes']['timing']['active_seconds'])*1000 for a,b in zip(groups[arm],groups['remote'])]
                row['comparisons'][arm]={'paired_mean_saving_ms':st.mean(differences),'exploratory_95pct_interval_ms':interval(differences),'differences_ms':differences}
            results.append(row)
    pressure=[];pressure_records=0
    def counters(text):return {a:int(b) for a,b in (line.split() for line in text.splitlines())}
    for limit in (256,768,1536):
        for arm in ARMS:
            trials=[]
            for s in range(6):
                dest=out/'pressure'/f'm{limit}-s{s}-{arm}'
                r=json.loads((dest/'result.json').read_text());u=json.loads((dest/'unit.json').read_text())
                t={'session':s,'status':r['status'],'unit_returncode':u['returncode']}
                if r['status']=='passed':
                    check(r,manifest['datasets']['synthetic1m'],s)
                    pressure_records+=sum(w['records'] for p in r['passes'].values() for w in p['workers'])
                    t.update(ms=r['passes']['timing']['active_seconds']*1000,pss_mib=r['passes']['memory']['memory_peak_pss_kib']/1024)
                else:t['error']=r.get('error')
                if (dest/'after.json').exists():
                    a=json.loads((dest/'after.json').read_text());b=json.loads((dest/'before.json').read_text());c=json.loads((dest/'after-copy.json').read_text()) if (dest/'after-copy.json').exists() else None
                    if int(a['memory.max'])!=limit*1024*1024 or int(a['memory.swap.max'])!=0:raise ValueError('wrong memory limit')
                    t.update(peak_mib=int(a['memory.peak'])/1024**2,events=counters(a['memory.events']),stat=counters(a['memory.stat']),swap_bytes=int(a['memory.swap.current']),pressure=a['memory.pressure'])
                    if c:t['after_copy_stat']=counters(c['memory.stat'])
                trials.append(t)
            good=[t for t in trials if t['status']=='passed']
            pressure.append({'limit_mib':limit,'arm':arm,'completed':len(good),'attempted':6,
                             'median_ms':st.median(t['ms'] for t in good) if good else None,'trials':trials})
    totals['pressure_verified_record_instances']=pressure_records
    setup={}
    for name,d in manifest['datasets'].items():
        setup[name]={'gzip_generation':d['gzip_generation'],'builds':{}}
        for arm in ('pyfastx-plain','pyfastx-gzip','gzip-cache'):
            builds=[b for b in d['new_builds'] if b['arm']==arm]
            setup[name]['builds'][arm]={'median_seconds':st.median(b['wall_seconds'] for b in builds),'range_seconds':[min(b['wall_seconds'] for b in builds),max(b['wall_seconds'] for b in builds)],'artifact_bytes':builds[0]['bytes']}
    save(out/'summary.json',{'totals':totals,'main':results,'pressure':pressure,'setup':setup})
    with (out/'all-conditions.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['dataset','workers','batch','pattern','arm','median_ms','median_pss_mib','mean_saving_vs_traceon_ms','interval_low','interval_high'])
        for r in results:
            for arm,m in r['arms'].items():
                c=r['comparisons'].get(arm,{})
                w.writerow([r['dataset'],*r['case'],arm,m['median_ms'],m['median_pss_mib'],c.get('paired_mean_saving_ms'),*c.get('exploratory_95pct_interval_ms',[None,None])])
    print(json.dumps(totals))

if __name__=='__main__':main(Path(sys.argv[1]))
