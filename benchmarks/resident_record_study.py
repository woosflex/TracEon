#!/usr/bin/env python3
"""Prospective persistent-reader control. See docs/resident-record-protocol.md.

Usage: resident_record_study.py prepare|smoke|seal|run|analyze OUTPUT
All execution is local. Run prepare/smoke/seal before run; retain every artifact.
"""
import argparse
from contextlib import ExitStack
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import platform
import random
import selectors
import statistics
import struct
import subprocess
import sys
import threading
import time

from indexed_mapping_study import sequence_records
from mapping_baseline import run_measured, sha256
from shared_mapping_pilot import pss_kib, process_cpu_seconds, start_server
from summarize_paired import interval

ROOT = Path(__file__).resolve().parents[1]
BUILD = Path('/tmp/traceon-release-check-20260922')
HTS = ROOT / 'build/research-inputs/samtools-1.24/htslib-1.24'
SAMTOOLS = HTS.parent / 'samtools'
CPUS = [8, 9, 10, 11]
PREFIX = ['taskset', '-c', ','.join(map(str, CPUS))]
CASES = ([{'processes': n, 'batch': b, 'pattern': 'random'} for n, b in itertools.product((1, 4, 8), (1, 64))]
         + [{'processes': 4, 'batch': b, 'pattern': 'random'} for b in (16, 256)]
         + [{'processes': 4, 'batch': 64, 'pattern': p} for p in ('ordered', 'hotspot')])


def save(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n')


def expected(records, ids):
    h = hashlib.sha256()
    total = found = 0
    def field(value):
        b = value.encode()
        h.update(struct.pack('<Q', len(b))); h.update(b)
    for key in ids:
        field(key)
        present = key in records
        h.update(bytes([present]))
        if present:
            seq, qual = records[key]
            field(seq); field(qual)
            total += len(seq) + len(qual)
            found += 1
    return {'records': len(ids), 'found': found, 'bytes': total, 'sha256': h.hexdigest()}


def validate(result, wanted):
    for key in ('records', 'found', 'bytes', 'sha256'):
        if result.get(key) != wanted[key]:
            raise ValueError(f'output validation failed: {key}')


def check_prepared(out, manifest):
    hashes = {out / 'resident_record_worker': manifest['worker_sha256'],
              BUILD / 'traceon': manifest['traceon_sha256'],
              HTS / 'libhts.a': manifest['htslib_sha256'],
              SAMTOOLS: manifest['samtools_sha256'],
              ROOT / 'docs/resident-record-protocol.md': manifest['protocol_sha256']}
    for data in manifest['datasets'].values():
        hashes[Path(data['reads'])] = data['source_sha256']
        for arm in ('cache', 'fai'):
            build = next(b for b in data['builds'] if b['arm'] == arm and b['trial'] == 0)
            hashes[Path(data[arm])] = build['sha256']
    for path, digest in hashes.items():
        if sha256(path) != digest:
            raise ValueError(f'prepared artifact changed: {path}')


def seal(out):
    """Freeze execution provenance separately from the original preparation."""
    if (out / 'execution-manifest.json').exists() or (out / 'runs.jsonl').exists():
        raise ValueError('execution already sealed or measurements exist')
    manifest = json.loads((out / 'manifest.json').read_text())
    smoke_result = json.loads((out / 'smoke/passed.json').read_text())
    if not smoke_result.get('corruption_rejected') or not smoke_result.get('optimized_corruption_rejected'):
        raise ValueError('normal and optimized smoke validation must pass before seal')
    check_prepared(out, manifest)
    names = ('resident_record_study.py', 'resident_record_worker.cpp',
             'indexed_mapping_study.py', 'mapping_baseline.py', 'shared_mapping_pilot.py',
             'subset_mapping_pilot.py', 'summarize_paired.py')
    paths = [ROOT / 'benchmarks' / name for name in names]
    paths.append(ROOT / 'docs/resident-record-protocol.md')
    hashes = {}
    for path in paths:
        relative = path.relative_to(ROOT)
        target = out / 'execution-source' / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        content = path.read_bytes()
        with target.open('xb') as stream:
            stream.write(content)
        hashes[str(relative)] = hashlib.sha256(content).hexdigest()
    ids = {job['ids']: sha256(Path(job['ids'])) for data in manifest['datasets'].values()
           for group in data['jobs'].values() for job in group}
    save(out / 'execution-manifest.json', {
        'preparation_manifest_sha256': sha256(out / 'manifest.json'),
        'smoke_sha256': sha256(out / 'smoke/passed.json'),
        'python': sys.version, 'source_sha256': hashes, 'request_sha256': ids})


def prepare(out):
    if (out / 'manifest.json').exists():
        raise ValueError('prepared study already exists; use a new directory')
    out.mkdir(parents=True, exist_ok=True)
    compile_command = ['c++', '-std=c++20', '-O3', '-DNDEBUG', '-I' + str(ROOT / 'include'), '-I' + str(HTS),
        str(ROOT / 'benchmarks/resident_record_worker.cpp'), str(BUILD / 'libtraceon_remote.a'),
        str(BUILD / 'libtraceon_core.a'), str(HTS / 'libhts.a'), str(BUILD / '_deps/zlib-ng-build/libz.a'),
        str(BUILD / 'third_party/lz4/build/cmake/liblz4.a'), '-lcrypto', '-lbz2', '-llzma', '-ldeflate',
        '-lcurl', '-lm', '-lpthread', '-o', str(out / 'resident_record_worker')]
    with (out / 'compile.log').open('w') as log:
        subprocess.run(compile_command, check=True, stdout=log, stderr=subprocess.STDOUT)
    real = ROOT / 'build/research-inputs/GSM461176_untreat_single_subset.fastq'
    if hashlib.md5(real.read_bytes()).hexdigest() != 'bba714969365c904694aa9d7c1199ffa':
        raise ValueError('Zenodo dataset MD5 mismatch')
    synthetic = out / 'synthetic-1000000.fq'
    small = out / 'synthetic-100000.fq'
    rng = random.Random(20260922)
    bases = bytes(b'ACGT'[i % 4] for i in range(256))
    qualities = bytes(33 + i % 41 for i in range(256))
    with synthetic.open('wb') as big, small.open('wb') as little:
        for i in range(1000000):
            record = b'@synthetic_%07d\n' % i + rng.randbytes(150).translate(bases) + b'\n+\n' + rng.randbytes(150).translate(qualities) + b'\n'
            big.write(record)
            if i < 100000:
                little.write(record)
    cookbook = ROOT / 'build/research-inputs/cookbook'
    datasets = [('short', cookbook / 'ecoli_mason_1.fq', True),
                ('long', cookbook / 'ecoli_p6_25x_canu.fa', False),
                ('drosophila', real, True), ('synthetic100k', small, True), ('synthetic1m', synthetic, True)]
    manifest = {'protocol_sha256': sha256(ROOT / 'docs/resident-record-protocol.md'),
                'platform': platform.platform(), 'cpus': CPUS, 'controller_cpu': 16,
                'compile_command': compile_command, 'compiler': subprocess.check_output(['c++', '--version'], text=True),
                'worker_sha256': sha256(out / 'resident_record_worker'), 'traceon_sha256': sha256(BUILD / 'traceon'),
                'htslib_sha256': sha256(HTS / 'libhts.a'), 'samtools_sha256': sha256(SAMTOOLS),
                'harness_sha256': {str(p.relative_to(ROOT)): sha256(p) for p in
                    [Path(__file__), ROOT / 'benchmarks/resident_record_worker.cpp', ROOT / 'benchmarks/shared_mapping_pilot.py', ROOT / 'benchmarks/mapping_baseline.py', ROOT / 'benchmarks/indexed_mapping_study.py']},
                'cpu_topology': subprocess.check_output(['lscpu', '-e=CPU,CORE,SOCKET,ONLINE'], text=True),
                'initial_meminfo': Path('/proc/meminfo').read_text(), 'cases': CASES, 'sessions': 6, 'datasets': {},
                'biological_source': {'doi': '10.5281/zenodo.61771', 'url': 'https://zenodo.org/records/61771/files/GSM461176_untreat_single_subset.fastq?download=1', 'md5': 'bba714969365c904694aa9d7c1199ffa'},
                'synthetic': {'seed': 20260922, 'length': 150, 'algorithm': 'Python Random.randbytes translated modulo 4 bases and modulo 41 quality values; 100k is prefix of 1m; not biological replication'}}
    for name, path, fastq in datasets:
        dest = out / name
        dest.mkdir()
        builds = []
        for trial in range(3):
            for arm in (('cache', 'fai') if trial % 2 == 0 else ('fai', 'cache')):
                target = dest / f'{arm}-{trial}'
                command = ([str(BUILD / 'traceon'), 'build', str(path), str(target)] if arm == 'cache' else
                           [str(SAMTOOLS), 'fqidx' if fastq else 'faidx', '--fai-idx', str(target), str(path)])
                result = run_measured(PREFIX + command, dest, f'{arm}-{trial}')
                builds.append(result | {'arm': arm, 'trial': trial, 'sha256': sha256(target), 'bytes': target.stat().st_size})
        records = {}
        for key, seq, qual in sequence_records(path, fastq):
            if key in records:
                raise ValueError(f'duplicate source ID in {name}: {key}')
            records[key] = (seq, qual)
        names = list(records)
        lengths = [len(v[0]) for v in records.values()]
        jobs = {}
        count = 4096 if name == 'long' else 8192
        hotspot = random.Random(20260922).sample(names, 32)
        for session, pattern in itertools.product(range(6), ('random', 'ordered', 'hotspot')):
            group = []
            for worker in range(8 if pattern == 'random' else 4):
                local = random.Random(8100 + session * 100 + worker)
                if pattern == 'hotspot':
                    ids = local.choices(hotspot, k=count)
                else:
                    positions = local.sample(range(len(names)), count)
                    if pattern == 'ordered':
                        positions.sort()
                    ids = [names[i] for i in positions]
                ids_file = dest / f's{session}-{pattern}-w{worker}.ids'
                ids_file.write_text(''.join(x + '\n' for x in ids))
                group.append({'ids': str(ids_file), 'expected': expected(records, ids)})
            jobs[f'{session}-{pattern}'] = group
        manifest['datasets'][name] = {'reads': str(path), 'fastq': fastq, 'source_sha256': sha256(path),
            'source_bytes': path.stat().st_size, 'count': len(records), 'length': {'min': min(lengths), 'median': statistics.median(lengths), 'max': max(lengths)},
            'cache': str(dest / 'cache-0'), 'fai': str(dest / 'fai-0'), 'builds': builds, 'jobs': jobs}
        save(out / 'preparation.json', manifest)
        print(name, 'prepared', len(records), flush=True)
        del records, names, lengths
    save(out / 'manifest.json', manifest)


class Worker:
    def __init__(self, command, log):
        self.log = log.open('wb')
        self.child = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=self.log, bufsize=0)
        self.buffer = b''
        try:
            self.setup = json.loads(self.line())
        except BaseException:
            self.close()
            raise

    def line(self):
        deadline = time.monotonic() + 120
        while b'\n' not in self.buffer:
            with selectors.DefaultSelector() as selector:
                selector.register(self.child.stdout, selectors.EVENT_READ)
                if not selector.select(max(0, deadline - time.monotonic())):
                    raise TimeoutError('worker timed out')
            block = os.read(self.child.stdout.fileno(), 65536)
            if not block:
                raise RuntimeError('worker exited: inspect stderr')
            self.buffer += block
        line, self.buffer = self.buffer.split(b'\n', 1)
        return line.decode()

    def send(self, text):
        self.child.stdin.write((text + '\n').encode())

    def close(self):
        if self.child.poll() is None:
            self.child.stdin.close()
            try:
                self.child.wait(timeout=5)
            except subprocess.TimeoutExpired:
                self.child.kill(); self.child.wait()
        self.child.stdout.close(); self.log.close()


def cohort(workers, jobs, batch, dest, label, server_pid, memory):
    for i, (worker, job) in enumerate(zip(workers, jobs)):
        worker.send(f'PREP\t{job["ids"]}\t{batch}\t{dest / (label + "-w" + str(i) + ".us")}')
    for worker in workers:
        if worker.line() != 'PREPARED':
            raise RuntimeError('worker barrier failed')
    pids = [w.child.pid for w in workers] + ([server_pid] if server_pid else [])
    samples, errors = [], []
    def sample_once():
        samples.append({'time': time.monotonic(), 'pss_kib': sum(pss_kib(p) for p in pids)})
    def sample_loop():
        try:
            while not stopped.wait(.005):
                sample_once()
        except BaseException as error:
            errors.append(error)
    stopped = threading.Event()
    monitor = None
    if memory:
        sample_once()
        monitor = threading.Thread(target=sample_loop); monitor.start()
    server_cpu = process_cpu_seconds(server_pid) if server_pid else 0
    try:
        for worker in workers:
            worker.send('RUN')
        rows = [json.loads(w.line()) for w in workers]
        server_cpu = process_cpu_seconds(server_pid) - server_cpu if server_pid else 0
    finally:
        stopped.set()
        if monitor:
            monitor.join()
    if memory:
        sample_once()
    if errors:
        raise errors[0]
    for row, job in zip(rows, jobs):
        validate(row, job['expected'])
    span = max(r['end'] for r in rows) - min(r['start'] for r in rows)
    latencies = []
    for i, row in enumerate(rows):
        values = [float(v) for v in (dest / f'{label}-w{i}.us').read_text().splitlines()]
        if len(values) != row['requests'] or any(not math.isfinite(v) or v < 0 for v in values):
            raise ValueError('invalid request latency observations')
        latencies.extend(values)
    if samples:
        save(dest / f'{label}-memory.json', samples)
    return {'label': label, 'workers': rows, 'active_seconds': span, 'records_per_second': sum(r['records'] for r in rows) / span,
            'bytes_per_second': sum(r['bytes'] for r in rows) / span, 'server_cpu_seconds': server_cpu,
            'total_cpu_seconds': server_cpu + sum(r['cpu_seconds'] for r in rows),
            'request_median_us': statistics.median(latencies),
            'memory_peak_pss_kib': max((s['pss_kib'] for s in samples), default=None),
            'memory_samples': len(samples)}


def smoke(out):
    dest = out / 'smoke'
    dest.mkdir()
    reads = dest / 'fixture.fq'
    records = {'id:1/2': ('aCGTN', '!~ABC'), 'other': ('TT', '#$'), 'last': ('ACGT'*32, 'I'*128)}
    reads.write_text(''.join(f'@{k} description\n{s}\n+\n{q}\n' for k, (s, q) in records.items()))
    cache, fai = dest / 'cache', dest / 'fai'
    run_measured([str(BUILD / 'traceon'), 'build', str(reads), str(cache)], dest, 'build')
    run_measured([str(SAMTOOLS), 'fqidx', '--fai-idx', str(fai), str(reads)], dest, 'index')
    ids = ['id:1/2', 'missing', 'other', 'id:1/2', 'last'] * 53
    path = dest / 'ids'; path.write_text('\n'.join(ids) + '\n')
    job = {'ids': str(path), 'expected': expected(records, ids)}
    server, port, setup = start_server(BUILD / 'traceon', cache, dest, supported_cli=True, cpus=CPUS)
    checks = []
    try:
        for arm in ('indexed', 'remote'):
            w = Worker(PREFIX + [str(out / 'resident_record_worker'), arm, str(reads), str(fai), '1', str(port), setup['snapshot']], dest / f'{arm}.stderr')
            try:
                for batch in (1, 16, 64, 256):
                    checks.append(cohort([w], [job], batch, dest, f'{arm}-{batch}', server.pid if arm == 'remote' else None, True))
            finally:
                w.close()
            if w.child.returncode:
                raise RuntimeError('worker shutdown failed')
        broken = dict(checks[0]['workers'][0], sha256='0'*64)
        try:
            validate(broken, job['expected'])
        except ValueError:
            pass
        else:
            raise RuntimeError('corruption was accepted')
        optimized_code = '''import json, sys
sys.path.insert(0, sys.argv[1])
from resident_record_study import validate
broken, wanted = json.loads(sys.stdin.read())
try:
    validate(broken, wanted)
except ValueError:
    print("optimized corruption rejected")
else:
    raise RuntimeError("optimized validation accepted corruption")
'''
        optimized = subprocess.run([sys.executable, '-O', '-c', optimized_code, str(ROOT / 'benchmarks')],
                                   input=json.dumps([broken, job['expected']]), text=True,
                                   capture_output=True, check=True)
        (dest / 'optimized-check.stdout').write_text(optimized.stdout)
        (dest / 'optimized-check.stderr').write_text(optimized.stderr)
        save(dest / 'passed.json', {'checks': checks, 'corruption_rejected': True,
                                  'optimized_corruption_rejected': True,
                                  'optimized_check': {'command': optimized.args, 'returncode': optimized.returncode}})
    finally:
        server.terminate(); server.wait(timeout=10); server.stdout.close()
    if server.returncode:
        raise RuntimeError('server shutdown failed')


def run(out):
    manifest = json.loads((out / 'manifest.json').read_text())
    if not (out / 'execution-manifest.json').exists() or (out / 'runs.jsonl').exists():
        raise ValueError('require sealed execution and no existing measurements')
    execution = json.loads((out / 'execution-manifest.json').read_text())
    if sha256(out / 'manifest.json') != execution['preparation_manifest_sha256']:
        raise ValueError('original preparation manifest changed')
    if sha256(out / 'smoke/passed.json') != execution['smoke_sha256']:
        raise ValueError('smoke evidence changed')
    for relative, digest in execution['source_sha256'].items():
        if sha256(ROOT / relative) != digest:
            raise ValueError(f'source changed since execution seal: {relative}')
        if sha256(out / 'execution-source' / relative) != digest:
            raise ValueError(f'sealed source snapshot changed: {relative}')
    for path, digest in execution['request_sha256'].items():
        if sha256(Path(path)) != digest:
            raise ValueError(f'request IDs changed since execution seal: {path}')
    check_prepared(out, manifest)
    for name, data in manifest['datasets'].items():
        for session in range(6):
            dest = out / name / f'session-{session}'; dest.mkdir()
            save(dest / 'host-before.json', {'loadavg': Path('/proc/loadavg').read_text(), 'meminfo': Path('/proc/meminfo').read_text(), 'vmstat': Path('/proc/vmstat').read_text()})
            server, port, setup = start_server(BUILD / 'traceon', Path(data['cache']), dest, supported_cli=True, cpus=CPUS)
            save(dest / 'server-setup.json', setup)
            cases = list(enumerate(CASES)); random.Random(20260922 + session).shuffle(cases)
            save(dest / 'schedule.json', cases)
            try:
                for case_index, case in cases:
                    jobs = data['jobs'][f'{session}-{case["pattern"]}'][:case['processes']]
                    for arm in (('indexed', 'remote') if (session + case_index) % 2 == 0 else ('remote', 'indexed')):
                        label = f'c{case_index}-{arm}'
                        workers = []
                        try:
                            for i in range(case['processes']):
                                workers.append(Worker(PREFIX + [str(out / 'resident_record_worker'), arm, data['reads'], data['fai'],
                                    str(int(data['fastq'])), str(port), setup['snapshot']], dest / f'{label}-w{i}.stderr'))
                            passes = {}
                            for mode in ('warmup', 'timing', 'memory'):
                                passes[mode] = cohort(workers, jobs, case['batch'], dest, label + '-' + mode,
                                                      server.pid if arm == 'remote' else None, mode == 'memory')
                            row = {'dataset': name, 'session': session, 'case_index': case_index, 'case': case, 'arm': arm,
                                   'worker_setup': [w.setup for w in workers], 'passes': passes}
                            with (out / 'runs.jsonl').open('a') as stream:
                                stream.write(json.dumps(row) + '\n')
                        finally:
                            for w in workers:
                                w.close()
                        if any(w.child.returncode for w in workers):
                            raise RuntimeError('worker failed during shutdown')
                    print(name, session, case_index, 'passed', flush=True)
            finally:
                server.terminate(); server.wait(timeout=10); server.stdout.close()
            if server.returncode:
                raise RuntimeError('server shutdown failed')
            save(dest / 'host-after.json', {'loadavg': Path('/proc/loadavg').read_text(), 'meminfo': Path('/proc/meminfo').read_text(), 'vmstat': Path('/proc/vmstat').read_text()})


def analyze(out):
    rows = [json.loads(line) for line in (out / 'runs.jsonl').read_text().splitlines()]
    manifest = json.loads((out / 'manifest.json').read_text())
    if len(rows) != len(manifest['datasets']) * 6 * len(CASES) * 2:
        raise ValueError('incomplete experiment')
    result = []
    for name in manifest['datasets']:
        for i, case in enumerate(CASES):
            groups = {arm: sorted([r for r in rows if r['dataset'] == name and r['case_index'] == i and r['arm'] == arm], key=lambda r: r['session']) for arm in ('indexed', 'remote')}
            if any([r['session'] for r in g] != list(range(6)) for g in groups.values()):
                raise ValueError('missing or duplicate sessions')
            diffs = [x['passes']['timing']['active_seconds'] - y['passes']['timing']['active_seconds'] for x, y in zip(groups['indexed'], groups['remote'])]
            result.append({'dataset': name, 'case': case, 'sessions': 6,
                'median_active_seconds': {a: statistics.median(r['passes']['timing']['active_seconds'] for r in g) for a, g in groups.items()},
                'median_records_per_second': {a: statistics.median(r['passes']['timing']['records_per_second'] for r in g) for a, g in groups.items()},
                'median_memory_pss_kib': {a: statistics.median(r['passes']['memory']['memory_peak_pss_kib'] for r in g) for a, g in groups.items()},
                'median_cpu_seconds': {a: statistics.median(r['passes']['timing']['total_cpu_seconds'] for r in g) for a, g in groups.items()},
                'mean_indexed_minus_remote_seconds': statistics.mean(diffs), 'range_seconds': [min(diffs), max(diffs)],
                'exploratory_95pct_mean_interval_seconds': interval(diffs),
                'memory_pass_over_timing_ratio': {a: statistics.median(r['passes']['memory']['active_seconds'] / r['passes']['timing']['active_seconds'] for r in g) for a, g in groups.items()}})
    save(out / 'analysis.json', result)
    print(json.dumps([r for r in result if r['case'] == {'processes': 4, 'batch': 64, 'pattern': 'random'}], indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('action', choices=('prepare', 'smoke', 'seal', 'run', 'analyze'))
    p.add_argument('output', type=Path)
    args = p.parse_args()
    allowed = os.sched_getaffinity(0)
    if not set(CPUS + [16]) <= allowed:
        p.error('CPU pool unavailable; explicitly adapt and record protocol for this host')
    os.sched_setaffinity(0, {16})
    globals()[args.action](args.output.resolve())
