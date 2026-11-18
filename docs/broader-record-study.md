# Broader record-access claim validation — 2026-09-22

This extension tests competitive performance beyond the earlier native HTSlib
control, ordinary gzip access, and constrained-memory feasibility. The user
confirmed that no second machine is available: cross-host reproduction remains
unvalidated. No evidence here establishes a new mapping paradigm or exhaustive
novelty across sequence-access systems.

## Design

The [prospective protocol](broader-record-protocol.md) specifies three inputs:
Drosophila 100k FASTQ, synthetic 1m FASTQ, and corrected E. coli long FASTA.
Four arms use persistent native HTSlib 1.24, pyfastx 2.3.1 plain, pyfastx 2.3.1 gzip,
or TracEon. One native driver embeds CPython to call pyfastx's public API; all
arms use the same C++ SHA-256 consumption, timing and output-validation loop.
Python initialization belongs to pyfastx startup and PSS. This is an API/backend
integration comparison, not a comparison of language-independent algorithms or
ordinary Python application throughput. Public indexed API behavior is described
in the [official pyfastx documentation](https://pyfastx.readthedocs.io/en/latest/api_reference.html).

Four conditions span one/four workers, scalar/batch 64 and random/source-ordered
requests. Each input has six restarted server sessions. Three arm permutations
and their reverses balance each pair's precedence 3/3; arm position counts are
not equal. Case order is shuffled deterministically. Each worker requests 1,024
records (512 long), with ordered cases using the identical random selection
sorted by source position. Preparation initially truncated previously sorted
lists; execution sealing corrected them before measurement and preserved the
original preparation manifest. Timing includes digest consumption, excludes
initialization and output-file writes. Separate memory passes measure aggregate
process PSS; filesystem page-cache memory is not included. Warmup/timing/memory
passes reuse the same IDs. The replication unit is the session; six-pair bootstrap
intervals are exploratory and not multiplicity-adjusted significance tests.

The gzip stream is ordinary gzip level 6, not BGZF. Decompression must reproduce
source bytes; three gzip-derived TracEon builds must match plain cache artifacts.
Three pyfastx index builds per representation record setup and storage costs.
Gzip indexed retrieval performs decompression that TracEon's resident decoded
snapshot has already paid for. A speedup here does not mean equivalent memory,
storage, or setup costs. The synthetic million-read HTSlib arm remains a stress
test beyond its documented recommended FASTQ scope; pyfastx is the added
appropriate indexed alternative.

Memory-limit experiments use one systemd user service per arm/session, with
controller, server and four workers under 256/768/1536 MiB limits and zero swap.
Fresh input/index/cache copies are made by ordinary read/write inside each
service so new file pages are charged there, addressing the kernel's
[file-page ownership rules](https://docs.kernel.org/admin-guide/cgroup-v2.html#memory-ownership).
Copies occur before retrieval timing. The controller is included in cgroup
memory, while sampled process PSS includes workers and the server only. Peak
cgroup memory includes copying and startup; it is not a query-only metric.
Successful limited runs have warmup/timing/memory passes. This tests limited
working-set feasibility; a limit does not guarantee sustained reclaim during
every measured query window. Failed/OOM/time-limited
runs remain in the denominator. Cgroup events, pressure, file/anonymous memory,
swap and unit diagnostics distinguish observed reclaim and failures; boundary
counters cannot reconstruct every transient event. No host-wide eviction or
pressure allocation is performed. Limits are tested only on synthetic 1m.

## Correctness boundary discovered before timing

The first smoke test failed on a missing FASTQ identifier. A separate minimal
Python reproduction on the installed pyfastx 2.3.1 wheel showed that after
`Fastq['missing']` raises `KeyError`, subsequent known-present lookups on that
handle also raise `KeyError`, for both plain and gzip. FASTA handles recovered.
This is a version-specific observed limitation, not a claim about all pyfastx
versions or inputs. The benchmark does not patch the competitor or count these
failed deliveries as speed gains. All performance requests are present IDs;
a second smoke suite checks mixed case, punctuation, qualities, duplicates,
partial batches and all batch sizes. A post-run probe reproduced the FASTQ failure on both complete Drosophila and
synthetic-million inputs, while the long FASTA reader recovered. Source/index
hashes remained unchanged. The failed fixture and reproduction scripts
are retained alongside successful checks. TracEon missing-ID recovery remains
covered by the earlier service tests.

## Results

All 288 warm-study arm observations completed: 936 persistent workers, 864
validated cohorts and **2,396,160 exact record deliveries** across three passes.
The 72 limited-memory trials produced 60 completed runs and 737,280 additionally
validated deliveries; 12 runs failed with recorded OOM kills. These are workload
observations and correctness counts, not independent biological replicates.

Primary condition: four workers, batch 64, random IDs. Median verified-retrieval
milliseconds; initialization excluded. Positive paired savings favor TracEon.

| Input | HTSlib ms | pyfastx plain ms | pyfastx gzip ms | TracEon ms | Mean saving vs plain pyfastx ms [exploratory 95% interval] |
|---|---:|---:|---:|---:|---:|
| drosophila | 19.03 | 15.55 | 1,054.47 | 2.55 | 13.01 [12.58, 13.46] |
| synthetic1m | 27.18 | 18.64 | 3,407.78 | 4.31 | 14.57 [13.97, 15.18] |
| long | 15.47 | 15.41 | 584.78 | 21.66 | -6.48 [-7.80, -5.34] |

TracEon median time is lower in **6/12 conditions versus plain pyfastx**,
**7/12 versus HTSlib**, and **12/12 versus gzip pyfastx**. This is a descriptive
count, not a corrected family of significance tests. All primary paired
short-record differences favor TracEon; all primary long-record differences
favor both plain indexed readers. Scalar synthetic1m access illustrates why
comparator choice matters: TracEon is faster than HTSlib in this condition but
slower than plain pyfastx. The earlier study's scalar losses apply to its own
larger-cohort conditions, not every possible workload. The full CSV retains
every condition and paired interval.

| Primary input | HTSlib PSS MiB | pyfastx plain PSS MiB | pyfastx gzip PSS MiB | TracEon PSS MiB |
|---|---:|---:|---:|---:|
| drosophila | 55.60 | 40.23 | 55.30 | 36.07 |
| synthetic1m | 554.96 | 40.51 | 90.63 | 439.52 |
| long | 27.47 | 39.79 | 65.26 | 141.73 |

The million-record service uses approximately **10.8× plain pyfastx process
PSS**. Sharing a decoded cache therefore does not establish superior memory
efficiency against a suitable indexed alternative.

### Limited-memory outcomes

Each entry is successful/attempted trials, with median retrieval ms among
successful trials in parentheses. No timing is assigned to failed trials.

| Cgroup limit | HTSlib | pyfastx plain | pyfastx gzip | TracEon |
|---|---:|---:|---:|---:|
| 256 MiB | 0/6 (OOM) | 6/6 (373.07) | 6/6 (4,915.28) | 0/6 (OOM) |
| 768 MiB | 6/6 (403.97) | 6/6 (18.87) | 6/6 (3,417.84) | 6/6 (3.25) |
| 1536 MiB | 6/6 (23.84) | 6/6 (17.46) | 6/6 (3,453.27) | 6/6 (3.25) |

At 256 MiB, each HTSlib failure records two OOM kills and each TracEon failure
one; pyfastx completes with limit/reclaim events but no OOM kills. All arms
complete at 768/1536 MiB. These limits bound cgroup-charged memory, including
new file pages and controller costs. They are not universal minimum-memory
requirements. Timings across limits are descriptive: limit groups run in fixed
order and page-cache/reclaim behavior differs. They show that feasibility and
speed under a memory cap must be evaluated separately from unrestricted PSS.

### Setup and storage

| Input | pyfastx plain index median s | pyfastx gzip index median s | Gzip-derived TracEon cache median s | Source / gzip MiB |
|---|---:|---:|---:|---:|
| drosophila | 0.201 | 0.338 | 1.144 | 20.77 / 3.00 |
| synthetic1m | 1.157 | 3.907 | 43.975 | 308.04 / 162.02 |
| long | 0.197 | 1.097 | 50.083 | 111.25 / 32.51 |

These are three build trials per representation; build ranges, artifact sizes,
gzip generation time and worker/server initialization are retained in raw data.
Gzip generation is not charged to retrieval. The primary retrieval gains do
not by themselves prove a deployment break-even; this extension does not run
thousands of cohorts to observe amortization.

**Claim disposition:** batched random short-record gains extend to the tested
pyfastx alternatives and gzip inputs on this host. General speed superiority,
general memory superiority, and long-record superiority are rejected by
counterexamples. Cross-host reproducibility, broad architectural novelty and
general mapper-engine replacement remain unestablished.

The [evidence archive](results/broader-record-validation-2026-09-22.tar.gz) and
[SHA-256](results/broader-record-validation-2026-09-22.tar.gz.sha256) retain the
complete CSV/JSON matrix, raw passes, memory counters, failures, source seals,
preparation history, dependency wheel and correctness probes. Large datasets,
cache/index files and executables are omitted; hashes and commands identify
their local counterparts. A separate assistant independently rechecked main
digests, counts, medians and 256 MiB outcomes; this is internal review, not
external replication.


## Reproduction and limitations

See the [environment record](../outputs/broader-record-environment.md). Tools and
paths are local and pinned, not a portable release package. A clean output
requires compiling `broader_record_worker.cpp` with the recorded command and
running `broader_record_study.py prepare`, `smoke`, `seal`, `run`, and `pressure`.
The known missing-ID probe is run separately; its fixture is retained. The
archive's exact preparation source and sealed execution source preserve changes
made between preparation and measurement. `summarize_broader_records.py`
revalidates every successful recorded delivery and exports all comparisons.

This developer-led, single-host study is not external replication, external
research adoption, a complete competitor survey, a cold-cache benchmark or a
WAN evaluation. Memory pressure is controlled only within the tested cgroups;
CPUs 8–11 are shared and not exclusively reserved. The practical conclusion must
be conditioned on request locality, batching, payload length, reuse, storage
format and available memory.
