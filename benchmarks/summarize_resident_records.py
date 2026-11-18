#!/usr/bin/env python3
"""Validate completed evidence and export setup/primary summaries plus full CSV."""
import csv
import json
from pathlib import Path
import statistics
import sys

from resident_record_study import CASES, validate
from mapping_baseline import sha256


def main(out):
    manifest = json.loads((out / 'manifest.json').read_text())
    analysis = json.loads((out / 'analysis.json').read_text())
    rows = [json.loads(line) for line in (out / 'runs.jsonl').read_text().splitlines()]
    if len(rows) != len(manifest['datasets']) * 6 * len(CASES) * 2:
        raise ValueError('incomplete study')
    totals = {'process_instances': 0, 'passes': 0, 'records_all_passes': 0,
              'records_timing_passes': 0, 'bytes_timing_passes': 0}
    for row in rows:
        data = manifest['datasets'][row['dataset']]
        jobs = data['jobs'][f'{row["session"]}-{row["case"]["pattern"]}'][:row['case']['processes']]
        totals['process_instances'] += len(jobs)
        for phase, result in row['passes'].items():
            if len(result['workers']) != len(jobs):
                raise ValueError('worker count mismatch')
            totals['passes'] += 1
            for worker, job in zip(result['workers'], jobs):
                validate(worker, job['expected'])
                totals['records_all_passes'] += worker['records']
                if phase == 'timing':
                    totals['records_timing_passes'] += worker['records']
                    totals['bytes_timing_passes'] += worker['bytes']
    setup = {}
    for name, data in manifest['datasets'].items():
        builds = {arm: [b['wall_seconds'] for b in data['builds'] if b['arm'] == arm] for arm in ('cache', 'fai')}
        restarts = [json.loads((out / name / f'session-{s}/server-setup.json').read_text()) for s in range(6)]
        if any(r['records'] != data['count'] for r in restarts):
            raise ValueError('server dataset count mismatch')
        primaries = {arm: sorted([r for r in rows if r['dataset'] == name and r['case'] ==
                     {'processes': 4, 'batch': 64, 'pattern': 'random'} and r['arm'] == arm], key=lambda r: r['session']) for arm in ('indexed', 'remote')}
        # Initialization-only projection: identical worker binary/process launch
        # overhead is not measured here; it is excluded from both worker terms.
        extra = [statistics.median(builds['cache']) - statistics.median(builds['fai']) + restarts[s]['startup_seconds']
                 + sum(w['ready_seconds'] for w in primaries['remote'][s]['worker_setup'])
                 - sum(w['ready_seconds'] for w in primaries['indexed'][s]['worker_setup']) for s in range(6)]
        primary = next(a for a in analysis if a['dataset'] == name and a['case'] == {'processes': 4, 'batch': 64, 'pattern': 'random'})
        saving = primary['mean_indexed_minus_remote_seconds']
        low = primary['exploratory_95pct_mean_interval_seconds'][0]
        setup[name] = {'source_records': data['count'], 'source_bytes': data['source_bytes'], 'length': data['length'],
            'cache_build_median_seconds': statistics.median(builds['cache']), 'cache_build_range_seconds': [min(builds['cache']), max(builds['cache'])],
            'fai_build_median_seconds': statistics.median(builds['fai']), 'fai_build_range_seconds': [min(builds['fai']), max(builds['fai'])],
            'server_restart_median_seconds': statistics.median(r['startup_seconds'] for r in restarts),
            'server_restart_range_seconds': [min(r['startup_seconds'] for r in restarts), max(r['startup_seconds'] for r in restarts)],
            'initialization_only_primary_extra_seconds': extra,
            'projected_primary_break_even_cohorts': max(0, statistics.mean(extra)) / saving if low > 0 and saving > 0 else None,
            'projection_limit': 'Not an observed break-even. Median construction costs plus measured library initialization/readiness; common worker process-launch remainder omitted. Three builds and six sessions only.'}
    with (out / 'all-conditions.csv').open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['dataset', 'processes', 'batch', 'pattern', 'indexed_ms', 'traceon_ms',
                         'indexed_pss_mib', 'traceon_pss_mib', 'mean_saving_ms', 'interval_low_ms', 'interval_high_ms'])
        for a in analysis:
            writer.writerow([a['dataset'], a['case']['processes'], a['case']['batch'], a['case']['pattern'],
                a['median_active_seconds']['indexed'] * 1000, a['median_active_seconds']['remote'] * 1000,
                a['median_memory_pss_kib']['indexed'] / 1024, a['median_memory_pss_kib']['remote'] / 1024,
                a['mean_indexed_minus_remote_seconds'] * 1000, *[v * 1000 for v in a['exploratory_95pct_mean_interval_seconds']]])
    report = {'totals': totals, 'setup': setup,
              'interpretation': 'Session-level exploratory evidence; primary metric includes common SHA-256 consumption; not universal tool ranking.',
              'analysis_sha256': sha256(out / 'analysis.json'), 'raw_sha256': sha256(out / 'runs.jsonl'),
              'summary_source_sha256': sha256(Path(__file__))}
    (out / 'evidence-summary.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main(Path(sys.argv[1]).resolve())
