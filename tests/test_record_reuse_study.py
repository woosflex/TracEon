"""Record validation must remain active under optimized Python execution."""
from pathlib import Path
import subprocess
import sys
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'benchmarks'))
from record_reuse_study import choose_ids, validate_records


def check():
    expected = [('a', 'AC', 'II'), ('b', 'GT', 'JJ'), ('a', 'AC', 'II')]
    cases = [expected, [('a', 'AA', 'II'), *expected[1:]],
             [('a', 'AC', 'IJ'), *expected[1:]],
             [('other', 'AC', 'II'), *expected[1:]],
             list(reversed(expected[:2])) + expected[2:], expected[:2]]
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / 'records.fq'
        for index, records in enumerate(cases):
            path.write_text(''.join(f'@{name}\n{seq}\n+\n{qual}\n'
                                    for name, seq, qual in records))
            try:
                validate_records(path, True, expected)
            except ValueError:
                if index == 0:
                    raise
            else:
                if index != 0:
                    raise RuntimeError(f'corrupted record output accepted: case {index}')
    names = [str(i) for i in range(100)]
    if not set(choose_ids(names, 1000, 'duplicates', 42)) <= set(names[:32]):
        raise RuntimeError('historical duplicates pattern changed')


if __name__ == '__main__':
    check()
    if '--optimized-child' not in sys.argv:
        subprocess.run([sys.executable, '-O', str(Path(__file__).resolve()),
                        '--optimized-child'], check=True)
        print('record validation passed, including corrupted output under python -O')
