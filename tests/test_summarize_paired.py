from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'benchmarks'))
from summarize_paired import interval, summarize
assert interval([2, 2, 2]) == [2, 2]
assert interval([-1, 1])[0] < 0 < interval([-1, 1])[1]
samples = [dict(arm=arm, repetition=r, processes=1, count=10,
                retrieval={'wall_seconds': seconds}, workflow_seconds=seconds + 4)
           for r in (1, 2) for arm, seconds in [('indexed', 3), ('shared', 2), ('batched', 1)]]
report = dict(samples=samples, cache_build={'wall_seconds': 8},
              server_setup={'startup_seconds': 1}, fai_build={'wall_seconds': 1})
rows = summarize(report)
assert rows[0]['mean_saving_seconds'] == 2
assert rows[0]['projected_break_even_cohorts'] == 4
assert rows[2]['mean_saving_seconds'] == 1
print('paired analysis checks passed')
