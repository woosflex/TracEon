# v1 mapped-cache validation — 22 September 2026

**The original local-cache direction has a measurable advantage in these tests.** Both v1 and the current embedded library reduced repeated sequence retrieval time against a persistent HTSlib reader. v1 additionally shared mapped payload pages across independent processes and restored faster than the current compressed format. This supports investigating safe mapped persistence; it does not establish a novel unmet need, mapper replacement, or universal superiority.

Four deterministic FASTA datasets, six repetitions per condition, one Intel Core Ultra 5 125H host. Dataset sizes below mean decimal MB of sequence payload; actual files include headers. All complete-record digest checks passed: **67,769,820 record deliveries** across 144 query observations and 144 memory cohorts (360 memory workers). Separately, the original driver completed 24 million length-only lookups with zero misses; the common worker completed another 48 million.

## Closest rerun of the historical lookup claim

The original v1 driver was compiled unchanged against tag `v1.0.0` (`600734f94bf657a20d74195527f5dec5c3a5eea3`), including its Robin Hood hash-map header. It selects actual stored IDs, checks membership, then times borrowed views and lengths. These rates **do not measure reading sequence bytes**. Original input files and historical host conditions are unavailable.

| Dataset | Historical README, million lookups/s | Original driver rerun, median million/s |
|---|---:|---:|
| 150 bp, 10 MB | 40–55 | 12.79 |
| 150 bp, 100 MB | 12–18 | 3.23 |
| 30 kb, 10 MB | 50–60 | 30.22 |
| 30 kb, 100 MB | 25–35 | 27.65 |

One of four medians falls within its historical range. The other three are lower on these inputs and this host. That is partial reproduction, not evidence that the original timings were fabricated. The previously identified broken pyfastx comparison remains unsuitable evidence for a 54× claim; this experiment does not rerun or validate that ratio.

## Complete sequence consumption

The controlled follow-up first issues the same 10,000 untimed IDs to every arm, then times the same 10,000 sequence lookups with length-framed ID/sequence SHA-256. Native HTSlib buffers are hashed directly. Medians exclude construction, startup, restore/open, and warmup. Each arm uses its natural API: borrowed cache views or allocated indexed-file fetches.

| Dataset | v1, ms | Current local library, ms | HTSlib 1.24, ms |
|---|---:|---:|---:|
| 150 bp, 10 MB | 7.73 | 7.08 | 76.78 |
| 150 bp, 100 MB | 9.66 | 10.01 | 132.64 |
| 30 kb, 10 MB | 133.25 | 132.76 | 208.50 |
| 30 kb, 100 MB | 145.14 | 144.19 | 268.87 |

Both cache arms won all six paired comparisons against HTSlib on every dataset. Exploratory bootstrap intervals of paired mean savings also stayed positive (`paired-analysis.json`). These observations support repeated local record access; a hash-consuming microbenchmark is not an alignment workflow. v1/current differences in complete consumption are comparatively small and inconsistent in direction. The initial unequal-warmup observations remain in the archive but are not the primary payload comparison.

## Restore and independent-process memory

Warm-host restore/index-open medians from the primary experiment:

| Dataset | v1 restore, ms | Current restore, ms | HTSlib FAI open, ms |
|---|---:|---:|---:|
| 150 bp, 10 MB | 9.20 | 20.23 | 30.10 |
| 150 bp, 100 MB | 79.47 | 214.45 | 376.61 |
| 30 kb, 10 MB | 0.85 | 13.34 | 0.31 |
| 30 kb, 100 MB | 7.38 | 77.06 | 1.73 |

v1 rebuilds its private hash table, so restore is not literally instantaneous or O(1). The many-record input is substantially more expensive than the long-record input despite similar payload size. No 1 GB restore claim was tested.

Every memory worker was independently exec’d, read and verified the entire dataset, then remained resident while PSS was sampled. Entries below are **one process → four-process aggregate**, MiB, median of six cohorts:

| Dataset | v1 PSS | Current PSS | HTSlib PSS |
|---|---:|---:|---:|
| 150 bp, 10 MB | 23.51 → 59.21 | 18.81 → 71.90 | 11.32 → 42.63 |
| 150 bp, 100 MB | 192.76 → 441.11 | 160.82 → 639.97 | 77.18 → 306.53 |
| 30 kb, 10 MB | 12.09 → 17.50 | 12.33 → 46.41 | 2.84 → 8.66 |
| 30 kb, 100 MB | 98.56 → 105.75 | 98.52 → 390.94 | 3.28 → 10.59 |

For 100 MB long records, v1 uses **72.95% less aggregate PSS** than the current library across four workers (105.75 versus 390.94 MiB). Its aggregate file PSS is about 97.77 MiB and anonymous PSS 7.98 MiB; the current implementation instead has about 388.16 MiB anonymous PSS. This directly supports the mapped-payload sharing mechanism. With many short records, duplicated indexes reduce the benefit: v1 still uses 441.11 MiB versus 639.97 MiB for current, but is larger than current in a single process.

HTSlib has the smallest process PSS in every condition. It need not keep every decoded sequence resident. PSS excludes unmapped filesystem page-cache pages and the controller, so these numbers are not total machine or cgroup memory. We have not demonstrated universal memory superiority over indexed readers.

## Costs and limits

Cache creation was measured once per input; these values are descriptive, not repeated performance estimates:

| Dataset | v1 load + save, s | Current load + save, s | FAI build, s | v1 / current cache, MB |
|---|---:|---:|---:|---:|
| 150 bp, 10 MB | 0.096 | 4.537 | 0.115 | 11.39 / 4.32 |
| 150 bp, 100 MB | 0.697 | 42.675 | 1.088 | 114.56 / 43.39 |
| 30 kb, 10 MB | 0.035 | 0.069 | 0.037 | 10.00 / 5.67 |
| 30 kb, 100 MB | 0.157 | 47.201 | 0.343 | 100.06 / 39.26 |

Compression trades smaller artifacts for construction/decode cost and private resident payloads. These query medians alone do not prove end-to-end savings for short-lived jobs; build amortization needs an explicit reuse workload. The inputs are synthetic single-line sequence-only FASTA, with no qualities, ambiguous bases, missing-ID performance test, gzip access, pyfastx/BioPython control, or second-host replication. Six repeats quantify local variation, not generality. CPU frequency and unrelated host activity were not controlled. The different map implementations and safety work prevent attributing every v1/current difference solely to compression.

The evidence supports a bounded next direction: evaluate a **validated uncompressed mapped-cache option for the current local API**, preserving integrity checks and lifetime guarantees. Do not roll back to v1’s unchecked binary reader. General bindings and an actual repeated-record research workflow remain needed to demonstrate user value; mapper-engine replacement remains unsupported.

## Reproduction and evidence

- [Protocol](mapped-cache-protocol.md), including the two exploratory follow-ups.
- [Main runner](../benchmarks/mapped_cache_study.py) and [common worker](../benchmarks/mapped_cache_worker.cpp).
- [Direction audit](project-direction-audit.md) and [earliest prototype rerun](earliest-claims-rerun.md).
- Full local inputs, binaries, exact commands, source snapshots, digests, raw observations, and analyses: `build/mapped-cache-2026-09-22/`.
- Compact evidence: [archive](results/mapped-cache-2026-09-22.tar.gz) and [SHA-256](results/mapped-cache-2026-09-22.tar.gz.sha256). Regenerable large input/cache files and binaries are excluded from the compact archive; their hashes remain in its full manifest.

No production code was changed. All work remains uncommitted.
