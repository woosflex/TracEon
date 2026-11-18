#!/usr/bin/env python3
"""Describe paired cohort timings; bootstrap intervals are exploratory with six runs."""
import argparse
import json
import random
import statistics
from pathlib import Path


def interval(values):
    """Percentile bootstrap of paired mean differences, in original time units."""
    if len(values) < 2:
        raise ValueError("at least two paired observations are required")
    rng = random.Random(42)
    means = sorted(statistics.mean(rng.choices(values, k=len(values))) for _ in range(10000))
    return [means[249], means[9749]]


def summarize(report):
    samples = report['samples']
    rows = []
    for processes, count in sorted({(s['processes'], s['count']) for s in samples}):
        selected = [s for s in samples if (s['processes'], s['count']) == (processes, count)]
        groups = {}
        for s in selected:
            key = (s['arm'], s['repetition'])
            if key in groups:
                raise ValueError('duplicate arm/repetition')
            groups[key] = s
        repetitions = sorted({s['repetition'] for s in selected})
        for baseline in ('indexed', 'shared'):
            for stage in ('retrieval', 'workflow'):
                def seconds(s):
                    return s['retrieval']['wall_seconds'] if stage == 'retrieval' else s['workflow_seconds']
                differences = [seconds(groups[baseline, r]) - seconds(groups['batched', r]) for r in repetitions]
                row = dict(processes=processes, count=count, baseline=baseline, stage=stage,
                           pairs=len(differences), mean_saving_seconds=statistics.mean(differences),
                           saving_range_seconds=[min(differences), max(differences)],
                           exploratory_mean_saving_95pct_bootstrap_seconds=interval(differences))
                if baseline == 'indexed' and stage == 'retrieval':
                    extra = max(0, report['cache_build']['wall_seconds'] + report['server_setup']['startup_seconds']
                                - report['fai_build']['wall_seconds'])
                    saving = row['mean_saving_seconds']
                    row['extra_setup_seconds'] = extra
                    row['projected_break_even_cohorts'] = extra / saving if saving > 0 else None
                rows.append(row)
    return rows


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('summary', type=Path)
    args = parser.parse_args()
    print(json.dumps(summarize(json.loads(args.summary.read_text())), indent=2))
