#!/usr/bin/env python3
"""Linux resident/fork index reuse: native fresh vs native resident vs TracEon resident.

Workers inherit a parent's index; arbitrary exec-launched processes cannot attach.
The native-resident control distinguishes sharing gains from lookup-backend gains.
"""
import argparse
import difflib
import itertools
import json
import os
from pathlib import Path
import platform
import selectors
import shutil
import subprocess
import threading
import time

from mapping_baseline import sha256, sam_summary
from minimap_lookup_pilot import patch as patch_lookup
from shared_mapping_pilot import pss_kib, run_group


def patch_main(text, traceon):
    def replace(old, new):
        nonlocal text
        assert text.count(old) == 1, old
        text = text.replace(old, new)
    replace('int main(int argc, char *argv[])',
            '#include "resident_minimap_pilot.h"\n#include <malloc.h>\n' +
            ('void traceon_pilot_discard_khash(mm_idx_t *mi);\n' if traceon else '') +
            'int main(int argc, char *argv[])')
    replace('if ((opt.flag & MM_F_OUT_SAM) && idx_rdr->n_parts == 1)',
            'if (!getenv("TRACEON_RESIDENT_PILOT") && (opt.flag & MM_F_OUT_SAM) && idx_rdr->n_parts == 1)')
    hook = '''
        if (getenv("TRACEON_RESIDENT_PILOT")) {
            int result;
            if (!mm_idx_reader_eof(idx_rdr) || idx_rdr->n_parts != 1 ||
                !(opt.flag & MM_F_OUT_SAM) || opt.split_prefix || argc - o.ind != 2) {
                fprintf(stderr, "resident pilot requires single-part SAM mapping and one query argument\\n");
                return 1;
            }
            DISCARD
            malloc_trim(0);
            result = resident_pilot(mi, &opt, n_threads, rg, argc, argv);
            mm_idx_destroy(mi);
            mm_idx_reader_close(idx_rdr);
            return result;
        }
'''.replace('DISCARD', 'traceon_pilot_discard_khash(mi);' if traceon else '')
    replace('\t\tif (alt_list) mm_idx_alt_read(mi, alt_list);',
            '\t\tif (alt_list) mm_idx_alt_read(mi, alt_list);' + hook)
    return text


def build(src, lib, out):
    root = Path(__file__).resolve().parents[1]
    binaries = {}
    for arm in ('native', 'traceon'):
        dest = out / arm
        shutil.copytree(src, dest, ignore=shutil.ignore_patterns('.git', '*.o', '*.a', 'minimap2'))
        for filename in ('main.c', 'index.c'):
            original = (src / filename).read_text()
            changed = original
            if filename == 'main.c':
                changed = patch_main(original, arm == 'traceon')
            elif arm == 'traceon':
                changed = patch_lookup(original) + '''
/* Called only after map-option updates/statistics and exhaustive copy checks. */
void traceon_pilot_discard_khash(mm_idx_t *mi)
{
    uint32_t i;
    for (i = 0; i < 1U<<mi->b; ++i) {
        kh_destroy(idx, (idxhash_t*)mi->B[i].h);
        mi->B[i].h = 0;
    }
}
'''
            (dest / filename).write_text(changed)
            (out / f'{arm}-{filename}.patch').write_text(''.join(difflib.unified_diff(original.splitlines(True), changed.splitlines(True), fromfile=f'a/{filename}', tofile=f'b/{filename}')))
        command = ['make', '-j4', f'INCLUDES=-I{root / "include"} -I{root / "benchmarks"}']
        if arm == 'traceon':
            command.append(f'LIBS={lib} -lstdc++ -lm -lz -lpthread')
        with (out / f'{arm}-build.log').open('w') as stream:
            subprocess.run(command, cwd=dest, stdout=stream, stderr=subprocess.STDOUT, check=True)
        binaries[arm] = dest / 'minimap2'
    return binaries


class Resident:
    def __init__(self, command, out, label):
        started = time.perf_counter()
        self.log = (out / f'{label}.stderr').open('wb')
        self.child = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                      stderr=self.log, env={**os.environ, 'TRACEON_RESIDENT_PILOT': '1'},
                                      start_new_session=True, bufsize=0)
        self.buffer = b''
        try:
            assert self.line() == f'READY {self.child.pid}'
            assert len(list(Path(f'/proc/{self.child.pid}/task').iterdir())) == 1
            self.setup = {'seconds': time.perf_counter() - started, 'pid': self.child.pid,
                          'ready_pss_kib': pss_kib(self.child.pid), 'command': command}
        except BaseException:
            self.close()
            raise

    def line(self):
        deadline = time.monotonic() + 60
        while b'\n' not in self.buffer:
            with selectors.DefaultSelector() as selector:
                selector.register(self.child.stdout, selectors.EVENT_READ)
                if not selector.select(max(0, deadline - time.monotonic())):
                    raise TimeoutError('resident protocol timed out')
            block = os.read(self.child.stdout.fileno(), 65536)
            if not block:
                raise RuntimeError('resident terminated before completing protocol')
            self.buffer += block
        result, self.buffer = self.buffer.split(b'\n', 1)
        return result.decode()

    def cohort(self, reads, outputs, out, label):
        stopped, samples, errors = threading.Event(), [], []
        parent = self.child.pid
        def sample():
            try:
                while not stopped.is_set():
                    children = [int(p) for p in Path(f'/proc/{parent}/task/{parent}/children').read_text().split()]
                    samples.append({'seconds': time.perf_counter() - started,
                                    'worker_pids': children,
                                    'aggregate_pss_kib': pss_kib(parent) + sum(pss_kib(p) for p in children)})
                    stopped.wait(0.005)
            except BaseException as error:
                errors.append(error)
        for path in [*reads, *outputs]:
            assert not any(c in str(path) for c in '\t\r\n')
        request = str(len(reads)) + '\n' + ''.join(f'{r}\t{o}\n' for r, o in zip(reads, outputs))
        started = time.perf_counter()
        monitor = threading.Thread(target=sample)
        monitor.start()
        children = []
        try:
            self.child.stdin.write(request.encode())
            while True:
                line = self.line()
                if line.startswith('DONE '):
                    assert line == 'DONE 0', line
                    break
                fields = line.split()
                assert fields[0] == 'CHILD' and len(fields) == 6, line
                assert int(fields[2]) == 0, line
                children.append({'pid': int(fields[1]), 'user_seconds': float(fields[3]),
                                 'system_seconds': float(fields[4]), 'peak_rss_kib': int(fields[5])})
            elapsed = time.perf_counter() - started
        finally:
            stopped.set()
            monitor.join()
        if errors:
            raise errors[0]
        assert len(children) == len(reads)
        (out / f'{label}-memory.json').write_text(json.dumps(samples) + '\n')
        return {'label': label, 'wall_seconds': elapsed, 'children': children,
                'sampled_peak_aggregate_pss_kib': max(s['aggregate_pss_kib'] for s in samples),
                'memory_samples': len(samples), 'settled_parent_pss_kib': pss_kib(parent)}

    def close(self):
        import signal
        if self.child.poll() is None:
            self.child.stdin.close()
            try:
                self.child.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(self.child.pid, signal.SIGKILL)
                self.child.wait()
        self.child.stdout.close()
        self.log.close()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for arg in ('source', 'library', 'previous-pilot', 'output'):
        p.add_argument('--' + arg, type=Path, required=True)
    a = p.parse_args()
    src, lib, previous, out = (v.resolve() for v in (a.source, a.library, a.previous_pilot, a.output))
    out.mkdir(parents=True, exist_ok=False)
    binaries = build(src, lib, out)
    prior = json.loads((previous / 'results.json').read_text())
    manifest = {'scope': __doc__, 'platform': platform.platform(),
                'affinity': sorted(os.sched_getaffinity(0)), 'library_sha256': sha256(lib),
                'binaries': {k: sha256(v) for k, v in binaries.items()},
                'policy': 'Warm OS cache, nonexclusive host; parent included in resident PSS; both parents alive during all arms; no exec in workers.',
                'datasets': {}}
    results = []
    orders = list(itertools.permutations(('fresh', 'native-resident', 'traceon-resident')))
    for name in ('short', 'long', 'sirv'):
        prior_run = next(r for r in prior['runs'] if r['dataset'] == name)
        index, reads = map(Path, prior_run['command'][-2:])
        preset = prior_run['command'][3]
        # One mapping thread per process isolates process-count scaling.
        def command(arm):
            return [str(binaries[arm]), '-a', '-x', preset, '-t', '1', str(index), str(reads)]
        parents = {}
        try:
            for arm in ('native', 'traceon'):
                parents[arm] = Resident(command(arm), out, f'{name}-{arm}-parent')
            manifest['datasets'][name] = {'index_sha256': sha256(index), 'reads_sha256': sha256(reads),
                                          'parents': {k: v.setup for k, v in parents.items()}}
            # Exercise repeated cohorts in each same parent, including a warmup.
            for count in (1, 2, 4):
                for repetition in range(-1, 6):
                    for arm in orders[max(0, repetition)]:
                        label = f'{name}-p{count}-r{repetition}-{arm}'
                        if arm == 'fresh':
                            result = run_group([command('native')] * count, out, label)
                            paths = [out / f'{label}-{i}.stdout' for i in range(count)]
                        else:
                            paths = [out / f'{label}-{i}.sam' for i in range(count)]
                            result = parents[arm.split('-')[0]].cohort([reads] * count, paths, out, label)
                        for path in paths:
                            assert sam_summary(path) == prior_run['sam'], path
                        result.update(dataset=name, processes=count, repetition=repetition, arm=arm,
                                      canonical_sam=prior_run['sam'])
                        results.append(result)
                        with (out / 'cohorts.jsonl').open('a') as stream:
                            stream.write(json.dumps(result) + '\n')
                print(name, count, 'passed', flush=True)
        finally:
            for parent in parents.values():
                parent.close()
        (out / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    (out / 'results.json').write_text(json.dumps(results, indent=2) + '\n')


if __name__ == '__main__':
    main()
