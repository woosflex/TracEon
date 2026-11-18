# Installed shared-cache workflow versus indexed access

Date: 2026-09-21. Status: preliminary development-tree measurements; not a release
benchmark or a general mapping-speedup claim.

Three public cookbook workloads completed using the installed, snapshot-pinned
CLI. Each study kept one server resident across 84 independent querying processes.
All 189,840 returned record instances matched source IDs, sequences, qualities,
and request order. The 504 mapper executions passed 252 paired canonical-SAM
comparisons, including all warmups. FASTQ qualities were retained throughout.

## Inputs and execution

| Workload | Source records | Payload | Mapper preset | IDs per worker |
|---|---:|---|---|---|
| SIRV Iso-Seq example | 53 | FASTQ, 184–2,492 bases | `splice:hq` | 10, 50 |
| E. coli corrected PacBio reads | 13,124 | FASTA, approximately 116.7 MB source | `map-pb` | 100, 1,000 |
| E. coli Mason mate-1 simulation | 10,000 | FASTQ, 100 bases | `sr` | 100, 1,000 |

The simulated mate-1 input was deliberately mapped as unpaired reads; these
results do not evaluate paired-end behavior or biological mapping accuracy.
All inputs come from the [minimap2 cookbook archive](https://github.com/lh3/minimap2/releases/download/v2.10/cookbook-data.tgz).
Source-file SHA-256 values, versions, binary hashes, and every command are in the
archived manifests. Build inputs were:

| Archive | SHA-256 |
|---|---|
| `cookbook-data.tgz` | `2cb3626a3c1cf1d26102aa3a98d0f73c7eac10f40428e4f6075ab9e2f179ac31` |
| minimap2 commit `79c9cc186b95f50bd899f69b48eba995ced810c6` source tar.gz | `6219d3ac9f0c667481075ff122108e7e587eb1b71c14a5713b12ff7bd601b233` |
| samtools 1.24 release tar.bz2 | `89b2a440123eeaa400392ce1736e7d60ce9041843027d76819753c5a8246bfdd` |

minimap2 2.30-r1287 used `make -j4`; samtools/HTSlib 1.24 used
`./configure --without-curses --disable-libcurl` followed by `make -j4` (GCC 16.2.1).
TracEon was the installed Release candidate, binary SHA-256
`dae5ca6ea5f89d5f6ea61c32009eadb82eafe2da9f80db8a21e65720dcb76b4a`.

Every cohort used 1, 2, or 4 processes, five measured repetitions, and one warmup.
Arm order alternated; five repetitions leave a 3/2 order imbalance, so these
small-sample results need balanced repeats before publication. Worker commands
and the server shared CPUs 8–11 (four distinct efficiency cores on this host).
The controller/sampler used CPU 16. Affinity was verified for the server and
included explicitly in every child command. Cores were not reserved and CPU
frequency was not fixed. Inputs were already read for validation, so OS page
cache was warm; no cache was dropped.

The comparator is [samtools indexed retrieval](https://www.htslib.org/doc/samtools-faidx.html),
using faidx for FASTA and fqidx for FASTQ. Its FAI was built once per study.
Each mapper used the same prebuilt minimap2 reference index. The shared arm
restored/hashed its cache once, and each worker pinned the expected snapshot.
The legacy per-process restore adapter was not part of this comparison.

## Observed results

The table shows the largest queried subset at four workers. Time is the median
for the whole concurrent retrieval group, including process startup. PSS is the
largest sampled simultaneous process total across the five repetitions; the
shared total includes its server.

| Workload | Indexed retrieval (ms) | Shared retrieval (ms) | Indexed PSS (MiB) | Shared PSS (MiB) | Indexed full workflow (ms) | Shared full workflow (ms) |
|---|---:|---:|---:|---:|---:|---:|
| SIRV, 50 IDs/worker | 13.03 | 9.58 | 2.26 | 4.48 | 320.52 | 338.92 |
| E. coli long, 1,000 IDs/worker | 77.94 | 44.89 | 11.52 | 119.74 | 3,750.80 | 3,796.05 |
| E. coli simulated short, 1,000 IDs/worker | 29.43 | 28.69 | 8.01 | 9.54 | 108.98 | 106.71 |

Shared retrieval was faster at these reported medians, but process memory was
higher than indexed access. The two long-read examples had slower median full
workflows despite faster extraction; mapping variation dominated the small
retrieval saving. The short-read difference was small. These observations do
not establish a dependable end-to-end mapping improvement.

Process wall timing blocks in `wait4`; the 5 ms PSS sampler runs independently
and does not impose a 5 ms polling floor on time. A separate 20-repeat observer
check on the small indexed command measured 3.75 ms median with sampling versus
3.56 ms without it; short measurements remain sensitive to instrumentation and
scheduling. Raw child CPU counters and server CPU deltas are retained. Server
CPU uses `/proc` clock ticks and is coarse for brief operations. Per-request tail
latency and detailed CPU/copy/checksum profiles remain outstanding.

PSS excludes the controller, filesystem page-cache pages not mapped into a
process, and the idle TracEon server during indexed-arm accounting. The server
physically stayed resident during both arms. Sampling can miss short peaks.
`wait4` RSS also includes process creation effects; the PSS totals above are the
reported memory comparison. None of these totals is total machine memory.

## Setup and amortization

| Workload | TracEon cache build (s) | Server restore + digest/start (s) | FAI build (s) |
|---|---:|---:|---:|
| SIRV | 0.0088 | 0.0061 | 0.0034 |
| E. coli long | 53.4336 | 0.1614 | 0.1845 |
| E. coli simulated short | 0.0228 | 0.0125 | 0.0158 |

Cache construction is once per artifact; restore/digest is once per server
lifetime. Neither is charged to every query process. Both belong in deployment
costs. The shared reference-index construction is recorded separately.

Using only the four-worker, 1,000-ID long-read median retrieval difference,
recovering the approximately 53.41 seconds of extra setup would take roughly
1,617 such groups (about 6.47 million requested record instances). Across the
long-read cohorts the estimates ranged from approximately 1,617 to 5,726 groups.
These are illustrative projections from median differences, not measured reuse
runs or confidence bounds, and they ignore the observed mapping variation.
Do not promise a break-even point for an arbitrary workflow from these values.

## Artifacts and remaining checks

[Compressed evidence bundle](results/development-study-2026-09-21.tar.gz) contains
raw resource/timing samples, manifests, result summaries, query ID lists, and
mapping harness sources verified against their recorded hashes. Full read/SAM
outputs remain under the corresponding `build/indexed-*-2026-09-21/` directories;
the public data, ID lists, and exact commands permit regeneration.

See [workflow instructions](research-workflow.md) and
[sustained resource checks](service-resource-study.md).
Before a release/publication claim: collect balanced, adequately repeated timing
samples and per-request tails; profile round trips/copies/CRC/server CPU; decide
on batching from that evidence; rerun at an exact released revision. The example
datasets do not establish performance for millions of short reads or large
human-reference workloads. In particular, samtools itself cautions against fqidx
for millions of short reads; a scale study should include an appropriate
alternative indexed reader rather than exploiting that limitation.

Evidence bundle SHA-256:
`a6d436cc8004880b2a4d72f8d0e89aa811798b9ee9090d5758e860dbd7e55a46`.
