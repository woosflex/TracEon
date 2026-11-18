#!/usr/bin/env python3
"""Bounded internal lookup substitution; retains khash for index I/O/statistics.

This tests the local kmer C API, NOT the resident sequence service. Retaining
both tables adds setup/memory cost; timings are whole-process feasibility data.
Run with a pristine minimap2 source directory, portable kmer archive and cookbook.
"""
import argparse
import difflib
import json
import platform
from pathlib import Path
import shutil
import statistics
import subprocess

from mapping_baseline import run_measured, sam_summary, sha256


def patch(source):
    def replace(old, new):
        nonlocal source
        assert source.count(old) == 1, old
        source = source.replace(old, new)

    replace('#include "khash.h"', '#include "khash.h"\n#include "kmerindex_c_api.h"')
    replace('mm128_v a;   //', 'kmerindex_t *traceon;\n\tmm128_v a;   //')
    helper = '''
/* ponytail: retain khash for statistics/serialization; this pilot measures
 * lookup compatibility, not a memory-efficient or shared index. */
static void traceon_copy(mm_idx_bucket_t *b)
{
    idxhash_t *h = (idxhash_t*)b->h;
    khint_t k;
    b->traceon = kmerindex_create();
    if (!b->traceon || !kmerindex_reserve(b->traceon, kh_size(h))) abort();
    for (k = 0; k < kh_end(h); ++k)
        if (kh_exist(h, k) && kmerindex_insert(b->traceon, kh_key(h,k), kh_val(h,k)) != 1) abort();
    /* Exhaustively validate singleton flag and value preservation before use. */
    for (k = 0; k < kh_end(h); ++k) if (kh_exist(h,k)) {
        uint64_t key;
        const uint64_t *v = kmerindex_get(b->traceon, kh_key(h,k)&~1ULL, &key);
        if (!v || key != kh_key(h,k) || *v != kh_val(h,k)) abort();
    }
}
'''
    replace('mm_idx_t *mm_idx_init(', helper + '\nmm_idx_t *mm_idx_init(')
    replace('free(mi->B[i].p);', 'kmerindex_destroy(mi->B[i].traceon);\n\t\t\tfree(mi->B[i].p);')
    start = source.index('const uint64_t *mm_idx_get(')
    end = source.index('\nvoid mm_idx_stat', start)
    source = source[:start] + '''const uint64_t *mm_idx_get(const mm_idx_t *mi, uint64_t minier, int *n)
{
    mm_idx_bucket_t *b = &mi->B[minier & ((1<<mi->b)-1)];
    uint64_t key;
    const uint64_t *v;
    *n = 0;
    if (!b->traceon) return 0;
    v = kmerindex_get(b->traceon, minier>>mi->b<<1, &key);
    if (!v) return 0;
    if (key & 1) { *n = 1; return v; }
    *n = (uint32_t)*v;
    return &b->p[*v>>32];
}
''' + source[end:]
    replace('\tb->h = h;\n', '\tb->h = h;\n\ttraceon_copy(b);\n')
    replace('\t\t\tkh_val(h, k) = x[1];\n\t\t}', '\t\t\tkh_val(h, k) = x[1];\n\t\t}\n\t\ttraceon_copy(b);')
    return source


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--library', type=Path, required=True)
    p.add_argument('--cookbook', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    src, lib, data, out = (v.resolve() for v in (a.source, a.library, a.cookbook, a.output))
    out.mkdir(parents=True, exist_ok=False)
    root = Path(__file__).resolve().parents[1]
    original = (src / 'index.c').read_text()
    changed = patch(original)
    (out / 'index.patch').write_text(''.join(difflib.unified_diff(original.splitlines(True), changed.splitlines(True), fromfile='a/index.c', tofile='b/index.c')))
    binaries = {}
    for arm in ('baseline', 'traceon'):
        dest = out / arm
        shutil.copytree(src, dest, ignore=shutil.ignore_patterns('.git', '*.o', '*.a', 'minimap2'))
        if arm == 'traceon':
            (dest / 'index.c').write_text(changed)
        command = ['make', '-j4', f'INCLUDES=-I{root / "include"}']
        if arm == 'traceon':
            command.append(f'LIBS={lib} -lstdc++ -lm -lz -lpthread')
        with (out / f'{arm}-build.log').open('w') as log:
            subprocess.run(command, cwd=dest, stdout=log, stderr=subprocess.STDOUT, check=True)
        binaries[arm] = dest / 'minimap2'
    results, inputs = [], {}
    for name, reference, reads, preset in (
        ('short', 'ecoli_ref.fa', 'ecoli_mason_1.fq', 'sr'),
        ('long', 'ecoli_ref.fa', 'ecoli_p6_25x_canu.fa', 'map-pb'),
        ('sirv', 'SIRV_E2.fa', 'SIRV_iso-seq.fq', 'splice:hq'),
    ):
        # Bound long-read execution while preserving complete FASTA records.
        read_path = data / reads
        if name == 'long':
            from subset_mapping_pilot import fasta_records
            read_path = out / 'long-1000.fa'
            import itertools
            with read_path.open('w') as stream:
                for identifier, sequence in itertools.islice(fasta_records(data / reads), 1000):
                    stream.write(f'>{identifier}\n{sequence}\n')
        index = out / f'{name}.mmi'
        inputs[name] = {'reference_sha256': sha256(data / reference), 'reads_sha256': sha256(read_path), 'preset': preset}
        run_measured([str(binaries['baseline']), '-x', preset, '-d', str(index), str(data / reference)], out, f'{name}-index')
        expected = None
        for threads in (1, 4):
            for repetition in range(6):
                for arm in (('baseline', 'traceon') if repetition % 2 == 0 else ('traceon', 'baseline')):
                    label = f'{name}-{threads}-{repetition}-{arm}'
                    result = run_measured([str(binaries[arm]), '-a', '-x', preset, '-t', str(threads), str(index), str(read_path)], out, label)
                    summary = sam_summary(out / f'{label}.stdout')
                    if expected is None:
                        expected = summary
                    assert summary == expected, (label, summary, expected)
                    result.update(dataset=name, threads=threads, repetition=repetition, arm=arm, sam=summary)
                    results.append(result)
        # Exercise construction and serialization in the patched mapper too.
        label = f'{name}-traceon-fasta'
        run_measured([str(binaries['traceon']), '-a', '-x', preset, '-t', '4', '-d', str(out / f'{name}-traceon.mmi'), str(data / reference), str(read_path)], out, label)
        assert sam_summary(out / f'{label}.stdout') == expected
        assert sha256(index) == sha256(out / f'{name}-traceon.mmi')
    report = {'scope': __doc__, 'platform': platform.platform(), 'inputs': inputs, 'source_index_sha256': sha256(src / 'index.c'), 'library_sha256': sha256(lib), 'binary_sha256': {k: sha256(v) for k, v in binaries.items()}, 'runs': results}
    (out / 'results.json').write_text(json.dumps(report, indent=2) + '\n')
    for name in ('short', 'long', 'sirv'):
        for threads in (1, 4):
            print(name, threads, {arm: round(statistics.median(r['wall_seconds'] for r in results if r['dataset'] == name and r['threads'] == threads and r['arm'] == arm), 4) for arm in binaries}, flush=True)


if __name__ == '__main__':
    main()
