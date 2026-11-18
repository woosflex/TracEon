# Local speed tests — 22 September 2026

**Mapped TracEon retains fast local queries while reducing restore time.** Its
sequence-consuming queries beat persistent HTSlib 1.24 on both tested inputs,
with one and four workers. Mapped versus compressed TracEon query performance
is broadly comparable; this experiment does not establish a general query-speed
improvement from mapping alone.

## Complete sequence queries

Each worker requests 10,000 random, present identifiers, reads every sequence
byte, and computes the same length-framed ID/sequence SHA-256. HTSlib buffers
are hashed directly, with no additional string copy. Each process opens its
reader once, warms the same IDs, and waits at a barrier before measurement.
Four-process cohorts request 40,000 records total, using separate portions of
the deterministic query list. Timings exclude build, process launch, restore,
and warmup. There is no network service in these tests.

Median elapsed query time across six repetitions (ms):

| Dataset (~100 MB payload) | Processes | Mapped TracEon | LZ4 TracEon | HTSlib | HTSlib / mapped time |
|---|---:|---:|---:|---:|---:|
| 150-base records | 1 | 8.81 | 8.25 | 124.82 | 14.17× |
| 150-base records | 4 | 7.96 | 10.89 | 156.50 | 19.67× |
| 30-kb records | 1 | 136.40 | 136.21 | 227.89 | 1.67× |
| 30-kb records | 4 | 164.75 | 161.25 | 272.39 | 1.65× |

Mapped was faster than HTSlib in every paired run in all four conditions.
Exploratory paired-mean bootstrap intervals also stayed positive. Against LZ4,
three of four intervals included zero. The four-worker short-record interval
was narrowly positive (0.03–2.96 ms saved, mean 1.19 ms), with five wins in six
pairs; treat that as a local exploratory result, not a general speed claim.

## Borrowed-view lookup speed

A separate test precomputes 200,000 random IDs, warms them, then performs five
passes (one million lookups per run). It validates returned lengths and
accumulates them so the work is observable. It does **not read sequence bytes**;
these rates must not be presented as full-record throughput or compared directly
with the consuming-query timings above. HTSlib has no equivalent borrowed-view
operation in this test and is therefore omitted.

Median million lookups/second, six runs:

| Dataset | Mapped TracEon | LZ4 TracEon |
|---|---:|---:|
| 150-base records | 3.26 | 3.37 |
| 30-kb records | 33.27 | 36.50 |

## Restore / index-open costs

Median single-process setup time from the consuming-query workers (ms):

| Dataset | Mapped restore | LZ4 restore | HTSlib FAI open |
|---|---:|---:|---:|
| 150-base records | 151.02 | 219.55 | 349.83 |
| 30-kb records | 19.51 | 64.13 | 1.33 |

HTSlib opens the long-record index substantially faster than either full cache.
Cache preparation and validation must be amortized over enough useful queries;
fast hot-path lookup alone is not end-to-end workflow evidence. Build times and
file-size trade-offs are reported separately in the
[mapped-storage implementation validation](mapped-storage-validation.md).

## Method and checks

- Same current static library for both cache modes, checked against the preceding
  study's library hash. No production code changes for this experiment.
- Two sealed synthetic single-line ACGT FASTA inputs: 666,666 × 150 bases and
  3,333 × 30,000 bases. Neither measures FASTQ quality retrieval or real analysis.
- Six repetitions with reproducibly shuffled arm order. Warm OS file caches;
  workloads restricted to CPU 8 for one process and CPUs 8–11 for four.
  Intel Core Ultra 5 125H, one nonexclusive Linux host; CPU frequency uncontrolled.
- Four-worker elapsed time spans the earliest worker start through the latest
  finish on the shared monotonic clock. Maximum observed start skew: 0.171 ms.
  Startups and reader setup finish before the barrier.
- 72 consuming-query cohorts, 180 independent workers, 24 million length-only
  lookups, and **2,040,000 independently checksum-verified sequence deliveries**
  including the separate lookup workers' content-validation passes. All passed.
- The checksum-consuming loop includes hashing overhead and is not a mapper or
  researcher-defined workflow. No new pyfastx or gzip comparison, cold-cache
  result, independent-host replication, or universal speedup is claimed.
- Medians of rates and medians of times are separate statistics; ratios in the
  primary table are explicitly ratios of median times. Individual observations,
  timing ranges, and paired intervals are preserved in the raw evidence.

## Reproduce

```sh
python3 benchmarks/local_speed_study.py build/NEW-local-speed-run
```

Requires the existing mapped build at `/tmp/traceon-mapped-check-20260922`,
HTSlib 1.24, and sealed inputs/caches from the preceding two mapped studies.
The output directory must be new. The runner reuses the existing benchmark
worker; the consuming variant adds only a start barrier and monotonic timestamps.
Commands, source, input/reference digests, environment, and raw results are saved.

Full local evidence: `build/local-speed-2026-09-22/`.
Compact [evidence archive](results/local-speed-2026-09-22.tar.gz) and
[SHA-256](results/local-speed-2026-09-22.tar.gz.sha256). Binaries and the already
archived/regenerable datasets are excluded; source and provenance hashes remain.

Everything remains uncommitted.
