#!/usr/bin/env python3
"""Check changed queries, worker failure, output protection and protocol bounds.

Usage: python benchmarks/check_resident_minimap_pilot.py PILOT_DIR PREVIOUS_DIR
"""
import json
from pathlib import Path
import subprocess
import sys

from mapping_baseline import sam_summary
from shared_minimap_index_pilot import Resident


def main():
    out, previous = (Path(p).resolve() for p in sys.argv[1:])
    check = out / 'lifecycle-checks'
    check.mkdir(exist_ok=False)
    prior = json.loads((previous / 'results.json').read_text())
    run = next(r for r in prior['runs'] if r['dataset'] == 'sirv')
    index, reads = map(Path, run['command'][-2:])
    subset = check / 'different-query.fq'
    # Cookbook fixture is four-line FASTQ; use a different, complete ten-read job.
    lines = reads.read_text().splitlines(True)
    assert lines[0].startswith('@') and lines[2].startswith('+')
    subset.write_text(''.join(lines[:40]))
    base = [str(out / 'native/minimap2'), '-a', '-x', 'splice:hq', '-t', '1', str(index)]
    expected = check / 'expected.sam'
    with expected.open('wb') as stdout, (check / 'expected.stderr').open('wb') as stderr:
        subprocess.run(base + [str(subset)], stdout=stdout, stderr=stderr, check=True)
    for arm in ('native', 'traceon'):
        command = [str(out / arm / 'minimap2'), *base[1:], str(reads)]
        parent = Resident(command, check, f'{arm}-changed-query')
        try:
            for number, query, wanted in ((0, reads, run['sam']), (1, subset, sam_summary(expected)), (2, reads, run['sam'])):
                output = check / f'{arm}-{number}.sam'
                parent.cohort([query], [output], check, f'{arm}-{number}')
                assert sam_summary(output) == wanted
        finally:
            parent.close()
        assert parent.child.returncode == 0
        for case in ('missing-query', 'existing-output', 'too-many-workers'):
            parent = Resident(command, check, f'{arm}-{case}')
            output = check / f'{arm}-{case}.sam'
            if case == 'existing-output':
                output.write_text('preserve this\n')
            query = check / 'absent.fq' if case == 'missing-query' else reads
            request = '5\n' if case == 'too-many-workers' else f'1\n{query}\t{output}\n'
            try:
                parent.child.stdin.write(request.encode())
                parent.child.stdin.close()
                parent.child.wait(timeout=10)
                assert parent.child.returncode != 0, case
                if case == 'existing-output':
                    assert output.read_text() == 'preserve this\n'
            finally:
                parent.close()
    (check / 'passed.json').write_text(json.dumps({'changed_query_cohorts': 6, 'failure_cases': 6, 'passed': True}) + '\n')
    print('Six changing-query cohorts and six failure checks passed.')


if __name__ == '__main__':
    main()
