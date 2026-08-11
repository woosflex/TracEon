#!/usr/bin/env python3
"""Sequential A/B: new-code vs old-code (20483dc) on 1000MB GZ + controls.
Writes results to a file so nothing is lost. Run from benchmarks/."""
import benchmark_runner as br
import os, sys

SCEN = [('WGS_GZ', 'fastq.gz', 150),
        ('PacBio_GZ', 'fastq.gz', 15000),
        ('Nanopore_GZ', 'fastq.gz', 30000),
        ('RefGenome_GZ', 'fasta.gz', 1000000),
        ('WGS', 'fastq', 150)]

def measure(label, fmt, avg, runs=3):
    times = []
    r = None
    for _ in range(runs):
        r = br.run_benchmark(label, fmt, 1000, avg)
        times.append(r['TracEon'][1])
    return sorted(times)[len(times)//2], times, (r['Restore'] if r else 0)

def main():
    out = open('/tmp/ab_results.txt', 'w')
    def log(*a):
        s = ' '.join(str(x) for x in a)
        print(s, flush=True)
        out.write(s + '\n'); out.flush()

    log('=== A/B: 1000MB, median of 3 ===')
    log(f'{"SCENARIO":<14} {"NEW(s)":<8} {"NEW runs":<30} {"OLD(s)":<8} {"OLD runs":<30} {"RESTORE_NEW":<12} {"RESTORE_OLD":<12}')
    results = {}

    # Leg 1: NEW code (current working tree binary)
    br.TRACEON_BINARY = '../build/traceon_driver'
    log('\n[leg 1] NEW code:', br.TRACEON_BINARY)
    for label, fmt, avg in SCEN:
        med, runs, rst = measure(label, fmt, avg)
        results[label] = {'new': med, 'new_runs': runs, 'new_restore': rst}
        log(f'{label:<14} {med:<8.3f} {str([round(t,3) for t in runs]):<30} {0:<8} {"-":<30} {rst:<12.3f} {-1:<12}')

    # Leg 2: OLD code (pre-fix 20483dc binary)
    br.TRACEON_BINARY = '/tmp/traceon_head/build/traceon_driver'
    log('\n[leg 2] OLD code:', br.TRACEON_BINARY)
    for label, fmt, avg in SCEN:
        med, runs, rst = measure(label, fmt, avg)
        r = results[label]
        r['old'] = med; r['old_runs'] = runs; r['old_restore'] = rst
        delta = ((med - r['new']) / r['new'] * 100) if r['new'] else 0
        log(f'{label:<14} {r["new"]:<8.3f} {str([round(t,3) for t in r["new_runs"]]):<30} {med:<8.3f} {str([round(t,3) for t in runs]):<30} {r["new_restore"]:<12.3f} {rst:<12.3f}  delta={delta:+.1f}%')

    log('\n=== DONE ===')
    out.close()

if __name__ == '__main__':
    main()
