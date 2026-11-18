#!/usr/bin/env python3
"""Record-service claim check: indexed versus scalar/batched resident retrieval.

Controlled affinity, balanced order, warm cache, exact payload/order checks.
No mapper is involved. Both one-time preparation and per-cohort costs are kept.
"""
import argparse
import itertools
import json
import os
from pathlib import Path
import platform
import random
import statistics
import subprocess

from indexed_mapping_study import sequence_records
from mapping_baseline import run_measured, sha256
from shared_mapping_pilot import run_group, start_server


def choose_ids(names, count, pattern, seed):
    rng = random.Random(42 if pattern == 'repeated' else seed)
    if pattern == 'duplicates':
        return rng.choices(names[:32], k=count)
    ids = rng.sample(names, count)
    if pattern == 'ordered':
        chosen = set(ids)
        return [name for name in names if name in chosen]
    return ids


def validate_records(path, fastq, expected):
    """Reject payload, quality, identifier, order or multiplicity mismatches."""
    if list(sequence_records(path, fastq)) != expected:
        raise ValueError(f'retrieved records differ from source: {path}')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('traceon', 'samtools', 'reads', 'output'):
        p.add_argument('--' + name, type=Path, required=True)
    p.add_argument('--fastq', action='store_true')
    p.add_argument('--cpus', type=int, nargs='+', required=True)
    p.add_argument('--controller-cpu', type=int, required=True)
    a = p.parse_args()
    allowed = os.sched_getaffinity(0)
    if not set(a.cpus) <= allowed or len(set(a.cpus)) != len(a.cpus) or a.controller_cpu not in allowed or a.controller_cpu in a.cpus:
        p.error('worker CPUs must be distinct and available; controller needs a separate available CPU')
    os.sched_setaffinity(0, {a.controller_cpu})
    paths = {k: v.resolve(strict=True) for k, v in vars(a).items() if isinstance(v, Path) and k != 'output'}
    out = a.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    records = {}
    for key, sequence, quality in sequence_records(paths['reads'], a.fastq):
        if key in records:
            raise ValueError('source IDs must be unique for indexed comparison')
        records[key] = (sequence, quality)
    names = list(records)
    if len(names) < 1000:
        raise ValueError('source must contain at least 1000 unique records')
    prefix = ['taskset', '-c', ','.join(map(str, a.cpus))]
    operation = 'fqidx' if a.fastq else 'faidx'
    cache, fai = out / 'reads.traceon', out / 'reads.fai'
    manifest = {'scope': __doc__, 'platform': platform.platform(), 'cpus': a.cpus,
                'controller_cpu': a.controller_cpu, 'fastq': a.fastq, 'records': len(names),
                'inputs': {k: {'path': str(v), 'sha256': sha256(v)} for k, v in paths.items()},
                'harness_sha256': {n: sha256(Path(__file__).with_name(n)) for n in
                    ('record_reuse_study.py', 'indexed_mapping_study.py', 'mapping_baseline.py', 'shared_mapping_pilot.py', 'subset_mapping_pilot.py')},
                'policy': 'Nonexclusive affinity pool shared by server and workers; controller separate. Warm OS cache. Indexed PSS excludes idle server; shared PSS includes server. Six balanced repetitions and one warmup per cell. CLI startup and output writing included.',
                'samtools_version': subprocess.check_output([str(paths['samtools']), '--version'], text=True),
                'cpu_topology': subprocess.check_output(['lscpu', '-e=CPU,CORE,SOCKET,ONLINE'], text=True)}
    (out / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    cache_build = run_measured(prefix + [str(paths['traceon']), 'build', str(paths['reads']), str(cache)], out, 'cache-build')
    fai_build = run_measured(prefix + [str(paths['samtools']), operation, '--fai-idx', str(fai), str(paths['reads'])], out, 'fai-build')
    server, port, setup = start_server(paths['traceon'], cache, out, supported_cli=True, cpus=a.cpus)
    samples = []
    orders = list(itertools.permutations(('indexed', 'scalar', 'batched')))
    try:
        if setup['records'] != len(records) or setup['cpu_affinity'] != sorted(a.cpus):
            raise ValueError('server record count or CPU affinity differs from requested setup')
        for pattern, processes, count in itertools.product(('changing', 'ordered', 'repeated', 'duplicates'), (1, 4), (10, 1000)):
            for repetition in range(7):
                jobs = []
                for worker in range(processes):
                    ids = choose_ids(names, count, pattern, 42 + count * 100 + repetition * 10 + worker)
                    path = out / f'{pattern}-p{processes}-n{count}-r{repetition}-w{worker}.ids'
                    path.write_text(''.join(key + '\n' for key in ids))
                    jobs.append((path, [(key, *records[key]) for key in ids]))
                for arm in orders[(repetition - 1) % 6]:
                    label = f'{pattern}-p{processes}-n{count}-r{repetition}-{arm}'
                    commands = []
                    for ids, _ in jobs:
                        if arm == 'indexed':
                            command = [str(paths['samtools']), operation, '--fai-idx', str(fai), '-n', '0', '-r', str(ids), str(paths['reads'])]
                        else:
                            command = [str(paths['traceon']), 'query', '--port', str(port), '--expect-snapshot', setup['snapshot'], '--ids', str(ids)]
                            if a.fastq:
                                command.append('--fastq')
                            if arm == 'batched':
                                command.extend(['--batch-size', '64'])
                        commands.append(prefix + command)
                    result = run_group(commands, out, label, server.pid if arm != 'indexed' else None)
                    output_hashes = []
                    for worker, (_, expected) in enumerate(jobs):
                        path = out / f'{label}-{worker}.stdout'
                        validate_records(path, a.fastq, expected)
                        output_hashes.append(sha256(path))
                    if server.poll() is not None:
                        raise RuntimeError('resident server exited during study')
                    result.update(pattern=pattern, processes=processes, count=count, repetition=repetition, arm=arm, output_sha256=output_hashes)
                    samples.append(result)
                    with (out / 'validated.jsonl').open('a') as stream:
                        stream.write(json.dumps(result) + '\n')
            print(pattern, processes, count, 'passed', flush=True)
    finally:
        server.terminate()
        try:
            server.wait(timeout=10)
        except subprocess.TimeoutExpired:
            server.kill()
            server.wait()
        server.stdout.close()
    if server.returncode != 0:
        raise RuntimeError(f'resident server did not shut down cleanly: {server.returncode}')
    aggregates = []
    extra = cache_build['wall_seconds'] + setup['startup_seconds'] - fai_build['wall_seconds']
    for pattern, processes, count in itertools.product(('changing', 'ordered', 'repeated', 'duplicates'), (1, 4), (10, 1000)):
        groups = {arm: [s for s in samples if (s['pattern'], s['processes'], s['count'], s['arm']) == (pattern, processes, count, arm) and s['repetition']]
                  for arm in ('indexed', 'scalar', 'batched')}
        differences = [x['wall_seconds'] - y['wall_seconds'] for x, y in zip(groups['indexed'], groups['batched'])]
        saving = statistics.mean(differences)
        aggregates.append({'pattern': pattern, 'processes': processes, 'count': count,
            'median_seconds': {k: statistics.median(s['wall_seconds'] for s in v) for k, v in groups.items()},
            'median_peak_pss_kib': {k: statistics.median(s['sampled_peak_aggregate_pss_kib'] for s in v) for k, v in groups.items()},
            'mean_indexed_minus_batched_seconds': saving, 'paired_saving_range_seconds': [min(differences), max(differences)],
            'extra_setup_seconds': extra, 'projected_break_even_cohorts': max(0, extra) / saving if saving > 0 else None,
            'projected_batched_minus_indexed_total_seconds': {str(n): extra - n * saving for n in (1, 10, 100, 1000)}})
    report = {'cache_build': cache_build, 'fai_build': fai_build, 'server_setup': setup,
              'cache_sha256': sha256(cache), 'fai_sha256': sha256(fai), 'cache_bytes': cache.stat().st_size,
              'fai_bytes': fai.stat().st_size, 'samples': samples, 'aggregates': aggregates}
    (out / 'summary.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(aggregates, indent=2))


if __name__ == '__main__':
    main()
