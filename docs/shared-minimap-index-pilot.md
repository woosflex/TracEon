# Cache once, map in many processes: internal-index experiment

On 2026-09-22, a Linux resident-parent prototype loaded a minimap2 reference index
once and reused it across 49 distinct mapping workers per parent. All 441 mapping
jobs matched the original mapper's canonical SAM output. Sharing reduced repeated
startup and index-memory costs in some workloads. The native minimap2 resident
control achieved comparable or better results, so this does **not** establish a
general advantage for TracEon's lookup implementation.

## What was tested

Three arms used the same prebuilt minimap2 index and read files:

1. **Fresh native:** launch ordinary minimap2 processes, each loading its index.
2. **Resident native:** initialize minimap2 once, then fork mapping workers.
3. **Resident TracEon:** initialize once, replace internal `mm_idx_get` lookups
   with TracEon's existing kmer C API, then fork mapping workers.

The resident TracEon parent exhaustively checks every copied key/value against
the native table, freezes the TracEon maps, computes mapping options/statistics,
then destroys the duplicate native minimizer tables before announcing readiness.
Both resident variants call `malloc_trim(0)` after setup. Posting lists and
reference sequence storage remain minimap2's implementation. Children inherit
the index through Linux fork/copy-on-write, perform real mapping, flush SAM, and
exit without destroying the inherited index. The parent remains single-threaded;
mapping threads run only in children.

**Boundary:** this demonstrates cache-once reuse by forked worker processes.
Unrelated programs launched with `exec` cannot attach to this index. It is not a
portable shared-memory file format, network lookup service, or shipped TracEon
CLI feature. The resident sequence-record RPC service is not used.

## Results

Each dataset used 1, 2 and 4 simultaneous workers, one mapping thread per worker,
six measured repetitions, and a warmup cohort per arm/process count. All six
execution-order permutations were used. Each worker mapped the entire indicated
read set; the table shows the time until all four workers finished.

| Four-worker workload | Fresh native seconds / MiB | Resident native seconds / MiB | Resident TracEon seconds / MiB |
|---|---:|---:|---:|
| 10,000 E. coli simulated short reads per worker | 0.1553 / 130.22 | 0.1016 / 50.16 | 0.1060 / 72.46 |
| 1,000 corrected E. coli PacBio reads per worker | 3.0146 / 1004.45 | 2.8820 / 931.42 | 2.8463 / 955.80 |
| 53 SIRV Iso-Seq reads per worker | 0.2822 / 245.32 | 0.2725 / 230.40 | 0.2868 / 234.65 |

Entries are medians. Memory is the median sampled peak **aggregate PSS**, including
the relevant resident parent and all its workers. PSS apportions shared pages;
summing worker RSS would double-count them. Sampling every 5 ms can miss brief
peaks, and process reads are sequential rather than perfectly simultaneous.

For four short-read workers, resident TracEon used about 44% less aggregate memory
and finished about 32% sooner than fresh native processes. Native residency saved
more memory and had a lower median time. For long reads, per-worker alignment
memory dominated, so sharing saved a much smaller fraction. SIRV showed no
consistent time advantage. One short-read worker used **more** memory with
resident TracEon: 56.08 MiB versus 33.98 MiB for fresh native.

The paired four-worker TracEon-minus-native-resident comparison does not support
a consistent backend speed advantage: the time saved by TracEon ranged from
-14.8 to +4.9 ms for short reads, -68.7 to +74.8 ms for long reads, and -45.8 to
+10.5 ms for SIRV. These are observed ranges, not confidence intervals.

Resident startup is excluded from cohort times and recorded separately:

| Dataset | Native startup seconds / ready PSS MiB | TracEon startup seconds / ready PSS MiB |
|---|---:|---:|
| Short | 0.0343 / 28.61 | 0.0791 / 50.45 |
| Long | 0.0323 / 21.83 | 0.0687 / 43.64 |
| SIRV | 0.0068 / 6.68 | 0.0115 / 10.31 |

These are single startup observations, not repeated estimates. Index construction
was common setup using the earlier pilot's `.mmi` files. Each resident parent
served 21 cohorts / 49 workers, with exactly one index-load log entry. Settled
parent PSS varied by at most 4 KiB across those cohorts.

## Correctness and limitations

- All 189 cohorts / 441 jobs passed SAM equivalence checks: 378 measured jobs
  and 63 warmup jobs. All mandatory alignment fields, qualities, optional tag
  values and reference headers were compared; only `@PG`, record order and tag
  order were ignored. Equivalence is to the original mapper, not independently
  established biological truth.
- Every resident cohort's memory samples observed the requested simultaneous
  worker count. The 294 resident jobs used distinct worker PIDs within each
  parent, with no index reload between jobs.
- Six additional cohorts changed query files and switched back in the same
  parent. Six failure checks covered missing reads, existing output protection,
  and rejection of more than four workers; all passed.
- This is a small local feasibility study on a nonexclusive host with warm,
  uncontrolled OS caches and no CPU affinity. Both resident parents stayed alive
  during comparisons; reported PSS includes only the arm's required processes.
  Small lifecycle checks overlapped part of the two-worker long-read runs; that
  cell should not be used for performance conclusions. The four-worker results
  above did not overlap those checks.
- Only single-part indexes and single-end jobs were tested. The controlled stdin
  driver is a benchmark mechanism, not a production service. Arbitrary process
  attachment, crash recovery, reference replacement, larger genomes and paired
  reads remain untested.

The supported conclusion is that resident index reuse works and can benefit
repeated mapping. A claim that TracEon's engine is superior needs further evidence
against the resident-native control, especially given TracEon's larger index
footprint here. A future independently attachable index should retain that
control and measure total memory, setup amortization and mapping equivalence.

## Reproduce

```sh
python benchmarks/shared_minimap_index_pilot.py \
  --source build/research-inputs/minimap2-79c9cc186b95f50bd899f69b48eba995ced810c6 \
  --library /tmp/traceon-release-check-20260922/libtraceon_kmer.a \
  --previous-pilot build/minimap-lookup-pilot-2026-09-22 \
  --output build/shared-minimap-index-new
python benchmarks/check_resident_minimap_pilot.py \
  build/shared-minimap-index-new build/minimap-lookup-pilot-2026-09-22
```

Requires Linux, a C/C++ toolchain and the inputs from the
[internal lookup pilot](minimap-lookup-pilot.md). Each output directory must be
new. Full raw SAM, patched source trees and binaries remain under
`build/shared-minimap-index-2026-09-22/`. The
[evidence archive](results/shared-minimap-index-2026-09-22.tar.gz) preserves source
patches, harness/header snapshots, build logs, commands, hashes, cohort metrics,
memory samples and lifecycle-check evidence; large SAM files are regenerable.
