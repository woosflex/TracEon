# Bounded batching: development measurements

Date: 2026-09-22. These are local working-tree measurements, not release results.

The existing indexed mapping harness now optionally compares samtools indexed
retrieval, scalar TracEon queries, and TracEon batches of up to 64 IDs. Each
dataset uses one resident server across all independent querying processes.
Six measured repetitions cover all six arm-order permutations, following one
unmeasured warmup. Inputs, executable hashes, CPU topology, commands, query IDs,
raw measurements, and harness hashes are recorded in the artifacts.

The SIRV FASTQ example uses 10 and 50 IDs per worker; the E. coli simulated
short-read example uses 100 and 1,000. Both use 1, 2, and 4 workers. Workers and
server share CPU affinity 8–11; the controller uses CPU 16. CPUs are not reserved,
frequency is uncontrolled, and the OS page cache is warm. Short reads use one
mate as unpaired input. Minimap2 is unchanged. Every arm's retrieved IDs,
sequences, qualities, and order match the source, and canonical SAM results match
across arms. Both servers shut down successfully.

At four workers, median cohort wall times across six measured repetitions:

| Dataset / IDs per worker | Stage | Indexed | Scalar | Batched |
|---|---|---:|---:|---:|
| SIRV / 50 | Retrieval | 14.27 ms | 9.61 ms | 9.38 ms |
| SIRV / 50 | Retrieval + mapping | 363.83 ms | 356.50 ms | 355.79 ms |
| Short reads / 1,000 | Retrieval | 27.46 ms | 28.93 ms | 10.40 ms |
| Short reads / 1,000 | Retrieval + mapping | 111.65 ms | 114.10 ms | 94.52 ms |

The short-read result supports retaining opt-in batching. The small SIRV
difference does not establish a material benefit. Cohort times include process
startup and output; they are not individual request latencies. Six samples do
not establish robust uncertainty bounds or a general mapping speedup. Larger
datasets, long reads, batch latency distributions, and release-revision repeats
remain necessary. Setup construction and server startup are measured separately;
they must be included when estimating deployment cost and reuse break-even.

Sampled aggregate retrieval PSS for the short-read row peaked at 8,200 KiB
indexed, 9,927 KiB scalar, and 9,286 KiB batched. Shared arms include the server;
indexed excludes the idle server. The controller and OS page cache are excluded.
Sampling every 5 ms may miss short peaks. This is not evidence of lower memory
use than indexed retrieval.

[Archived manifests, summaries, raw measurements, IDs, and harness sources](results/batch-study-2026-09-22.tar.gz)
have SHA-256 `c535dd8c3f688a173f6b310bb020d34179f3c41bc55e9901f13827e1a0b57a64`.
Full generated read/SAM outputs remain in `build/batch-sirv-2026-09-22/` and
`build/batch-short-2026-09-22/`; they are not bundled in this compact archive.
The archive is development evidence, not a complete release reproduction bundle.

## Paired variation and setup projection

`python3 benchmarks/summarize_paired.py STUDY/summary.json` compares arms within
each repetition. Positive differences mean time saved by batching. It reports
the observed range, mean, and an exploratory percentile bootstrap interval for
the mean (10,000 resamples of paired differences, seed 42). With only six pairs,
these intervals cannot characterize machine-to-machine variability, temporal
drift, or uncertainty in setup costs; they are not publication-grade confidence
in a general speedup. The script's arithmetic has a runnable synthetic check in
`tests/test_summarize_paired.py`.

For four workers requesting 1,000 short reads each, batching saved 14.05–20.55 ms
per retrieval cohort relative to indexed access across all six pairs; mean
saving was 17.12 ms. The exploratory mean interval is 15.30–19.02 ms. Relative
to scalar access, the mean saving was 18.31 ms, with observed range
17.32–19.44 ms. These paired means differ from differences between the medians
reported above.

The measured extra cache-build-plus-server-startup cost over FAI construction
was 36.81 ms on this short-read dataset. Dividing that cost by the paired mean
retrieval saving projects break-even at 2.15 four-worker cohorts (three complete
cohorts). This is a projection from one setup measurement and warm retrievals,
not a directly observed deployment break-even or an uncertainty bound.

For SIRV's four-worker, 50-ID workflow, batching versus indexed access had an
exploratory mean-saving interval of −0.43 to 15.28 ms; versus scalar access it
was −6.57 to 12.53 ms. Both span zero. This evidence does not resolve a full
workflow benefit on that dataset. The corresponding retrieval-only setup
projection is 3.96 cohorts, with the same limitations as above.

These analyses are generated in each study's `paired-analysis.json`; they were
added after the compact archive above and are reproducible from its summaries.

## Long-read follow-up

The same balanced design completed for the E. coli corrected PacBio FASTA
example: 100 and 1,000 IDs per worker, 1/2/4 workers, six measured repetitions,
and one warmup. All read and canonical-SAM comparisons passed; server exit was
zero. At four workers requesting 1,000 reads each:

| Stage | Indexed | Scalar | Batched |
|---|---:|---:|---:|
| Median retrieval | 79.45 ms | 44.98 ms | 40.89 ms |
| Median retrieval + mapping | 4.427 s | 4.401 s | 4.421 s |
| Peak sampled retrieval PSS | 11,770 KiB | 131,664 KiB | 138,161 KiB |

Batched retrieval saved a paired mean 39.44 ms versus indexed access and 3.43 ms
versus scalar access. Full-workflow exploratory intervals span zero against
both baselines. Extra setup cost was 52.85 s; retrieval savings project break-even
at approximately 1,341 complete four-worker cohorts, or 5.36 million requested
record instances. This projection inherits the single-setup and warm-cache
limitations above. Resident memory and setup costs remain substantial trade-offs.

[Long-read evidence archive](results/batch-long-study-2026-09-22.tar.gz), including
paired analysis and harness sources, has SHA-256
`f7befe06ed00d7d5bbc47e86dfcdfacf3a033527433347892af34e28cacdbe2f`.
Full outputs remain in `build/batch-long-2026-09-22/`.
