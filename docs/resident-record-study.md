# Persistent indexed-reader control and validity audit

This study tests whether earlier CLI-level retrieval gains survive keeping an
established indexed reader open. It is a developer-led local experiment, not
external replication or a claim of universal superiority. The
[prospective protocol](resident-record-protocol.md) was written before these
measurements; [literature notes](../outputs/record-access-literature.md) document
the comparator choice and relevant prior work.

## Results

All 600 arm observations completed: five datasets × ten conditions × six paired
sessions × two arms. Thirty server lifetimes and 2,520 persistent worker
processes delivered 18,579,456 record instances during timing. Across warmup,
timing and memory passes, all 55,738,368 record instances matched independently
computed digests, counts and presence flags. These deliveries are correctness
coverage, not millions of independent statistical replicates.

The prespecified primary comparison uses four workers, batch 64, random IDs.
Times are median verified-retrieval cohort milliseconds; PSS is median sampled
aggregate process MiB, including the service in its arm.

| Input | HTSlib ms | TracEon ms | Mean paired saving ms [exploratory 95% interval] | HTSlib PSS MiB | TracEon PSS MiB |
|---|---:|---:|---:|---:|---:|
| E. coli short, 10k | 55.26 | 6.63 | 49.18 [48.50, 50.28] | 17.15 | 17.68 |
| E. coli long, 13,124 | 73.18 | 88.83 | −17.28 [−25.74, −11.17] | 20.91 | 163.07 |
| Drosophila, 100k | 98.81 | 7.33 | 92.01 [90.14, 94.23] | 49.42 | 30.26 |
| Synthetic, 100k | 93.74 | 11.88 | 82.16 [80.66, 83.73] | 49.95 | 55.13 |
| Synthetic, 1m | 133.85 | 14.58 | 119.88 [118.20, 121.10] | 549.65 | 434.86 |

Batched random short-record gains survive removing repeated CLI initialization.
The full matrix has 28 positive mean savings (two intervals include zero) and
22 negative savings. **Scalar access loses all 15 conditions.** Long records
have eight losses and two unresolved comparisons, with no interval-supported
win. Ordered batch-64 access favors HTSlib on short E. coli (4.49 versus 6.56 ms),
long E. coli (36.13 versus 87.87 ms), and synthetic 100k (8.40 versus 11.68 ms).
These are descriptive comparisons, not multiplicity-adjusted hypothesis tests.
Memory is also workload-dependent: the primary service uses less process PSS
on Drosophila and synthetic 1m, but more on the other three inputs.

| Input | Cache build median s | FAI build median s | Server restart median s |
|---|---:|---:|---:|
| E. coli short | 0.027 | 0.017 | 0.013 |
| E. coli long | 49.534 | 0.196 | 0.178 |
| Drosophila | 0.773 | 0.099 | 0.067 |
| Synthetic 100k | 4.144 | 0.109 | 0.094 |
| Synthetic 1m | 41.634 | 1.602 | 1.182 |

Construction has three trials per artifact; restart has six. An
initialization-only projection includes construction, server readiness and
library initialization, but omits the common worker process-launch remainder.
Under stable primary-cohort costs it projects roughly 0.05, 6.3, 48.2 and 323
cohorts for short, Drosophila, synthetic 100k and synthetic 1m respectively.
These are not observed deployment break-even points. No finite long-read
break-even is supported because primary retrieval is slower as well as setup
being more expensive. Raw ranges and the signed cost calculation are preserved
in `evidence-summary.json`.

Host snapshots show approximately 7,015–7,217 MiB available memory and no new
swap-out, normal allocation-stall or OOM-kill increments during recorded session
windows. Swap was already occupied; 204 swap-in pages and 214 major faults were
recorded system-wide. These counters include other activity and cannot exclude
transient pressure or attribute a result to swapping. This was not a controlled
memory-pressure experiment.

The evidence supports **batched random short-record access on this host**;
scalar, long-record and some ordered workloads favor persistent indexed access.
It does not establish a universally better query engine.


## Why this control matters

The [previous matrix](record-reuse-study.md) compared a resident TracEon server
with newly launched samtools commands. That is a valid measurement of those
CLI workflows, but includes repeated index initialization in the indexed arm.
More repetitions of that comparison cannot isolate a retrieval advantage.
HTSlib exposes reusable native indexed readers, so the new experiment keeps
both HTSlib handles and TracEon client connections open.

The literature also limits novelty claims. Pyfastx offers reusable disk-backed
sequence access, and MetaProFi's FAQIndex uses LMDB for sequence-file indexing
shared across processes. This experiment does not benchmark those systems or
establish that cache-once access is new. See the
[HTSlib publication](https://doi.org/10.1093/gigascience/giab007),
[pyfastx publication](https://doi.org/10.1093/bib/bbaa368), and
[MetaProFi publication](https://doi.org/10.1093/bioinformatics/btad101).

## Methods and interpretation

Both arms use the same C++20 worker executable. The indexed arm loads one HTSlib
1.24 handle per process and uses exact-name full-sequence and quality access.
Returned HTSlib buffers are consumed without an unnecessary string copy.
The service arm maintains one snapshot-pinned TracEon connection per process.
Batch 1 uses scalar access; larger groups use protocol batching. This compares
native indexed access with a network service, not two equivalent network servers.
IPC overhead and differences in data representation are part of the comparison.

The five inputs are 10,000 simulated E. coli short reads, 13,124 corrected E. coli
long reads, a separately sourced 100,000-read Drosophila RNA-seq teaching subset,
and synthetic 100,000/1,000,000-record FASTQ collections. The biological subset's
published MD5 was verified against the
[Zenodo record](https://doi.org/10.5281/zenodo.61771). Its selection procedure is
not established by the deposit description. Synthetic records have unique IDs,
150 bases and varied qualities; the 100,000-record collection is a prefix of the
million-record one. These synthetic inputs test scale, not biological diversity.
The million-record HTSlib arm is a stress test beyond samtools' recommended FASTQ
scope and cannot establish superiority over indexing tools designed for that
scale. [Samtools explicitly documents this limitation](https://www.htslib.org/doc/samtools-fqidx.html).

Each dataset has six independently restarted server sessions. Ten conditions
cover 1/4/8 processes, batch sizes 1/16/64/256, random, source-ordered and
distributed-hotspot requests. Case order is randomized; each condition has three
indexed-first and three service-first sessions. The primary condition, selected
before measurement, is four processes and batch 64 with random IDs. Each worker
requests 8,192 records, or 4,096 for long reads. The same request lists are used
for warmup, timing and memory passes; random lists change between sessions and
workers, not between those three passes.

Initialization and ID-list loading occur before a READY/RUN barrier. Primary
time spans the earliest native worker start to the latest finish and includes
SHA-256 consumption of all returned fields. Output-file writing occurs after
the native finish timestamp. The secondary lookup timer excludes digest work,
includes result destruction, and is not interchangeable with primary verified
retrieval time. Every record contributes its requested ID, presence status and
length-framed sequence/quality bytes to a digest independently computed in
Python. Errors stop the experiment rather than excluding the observation.

Primary timing passes have no PSS sampler. Separate memory passes sample before,
during and after work while readers remain alive; the TracEon server is included
in its arm. This avoids missing the footprint of a worker that exits before a
sampler starts. Memory is **sampled aggregate process PSS**, not total physical
memory: filesystem page-cache memory is omitted. Server allocations can be
retained across cases. The later memory-pass timing is a sensitivity check, not
an isolated estimate of sampler overhead, because pass order is fixed.

Workload CPUs 8–11 are four distinct reported cores on this host; the controller
uses CPU16. Cores are not exclusively reserved. Host load, memory and swap
counters are recorded per session. Warm OS caches are intended. Server CPU is
read from `/proc` ticks, so small deltas are quantized. No cold-cache, controlled
memory-pressure, WAN, cross-machine or downstream mapping result is implied.

The replication unit is a restarted session. Confidence intervals are percentile
bootstrap intervals over six paired session-level differences (10,000 resamples,
seed 42). They are exploratory with this small sample, not confirmatory tests.
Requests and workers are not counted as independent statistical replicates.
The entire 50-condition matrix is retained rather than selecting only wins.
These controls follow the transparency and comparator-selection principles in
[Weber et al.](https://doi.org/10.1186/s13059-019-1738-8); following a guideline
does not itself certify scientific validity.

## Correctness and reproducibility

The fixture smoke checks exact names containing punctuation, mixed sequence
case, qualities, missing IDs, duplicate queries and all four batch sizes,
including partial final batches. Normal and optimized Python both reject a
deliberately wrong digest. The earlier CLI harness was also corrected to use
explicit validation errors instead of removable assertions; its regression
test changes IDs, sequence, qualities, order and multiplicity. Old benchmark
artifacts were preserved rather than retroactively rewritten.

Preparation builds each cache and FAI three times. Artifact/tool/input hashes
and compile commands are recorded. A separate execution seal captures every
benchmark source/helper, the protocol and all request-file hashes. The original
preparation source remains archived because provenance checks were tightened
after preparation started and before steady-state measurements. Runtime
preflight rejects changed binaries, inputs, selected cache/index files, request
lists or sealed sources.

```sh
python benchmarks/resident_record_study.py prepare build/resident-control-new
python benchmarks/resident_record_study.py smoke build/resident-control-new
python benchmarks/resident_record_study.py seal build/resident-control-new
python benchmarks/resident_record_study.py run build/resident-control-new
python benchmarks/resident_record_study.py analyze build/resident-control-new
python benchmarks/summarize_resident_records.py build/resident-control-new
python tests/test_record_reuse_study.py
```

The script uses the pinned local input/tool locations recorded in the
[environment note](../outputs/resident-record-environment.md). Adapt and record
those paths and CPU choices when reproducing on another host. `smoke` and `run`
require loopback sockets. Outputs must be new. Local full artifacts are under
`build/resident-record-control-2026-09-22/`.

The [compact evidence archive](results/resident-record-control-2026-09-22.tar.gz)
and [SHA-256 checksum](results/resident-record-control-2026-09-22.tar.gz.sha256)
include the full 50-condition CSV, analysis JSON, raw timing and memory samples,
request lists, manifests, smoke checks, build logs, and sealed benchmark sources.
Large datasets, cache/index files and executables remain local; their hashes and
regeneration commands are recorded. The post-analysis exporter is archived
separately from the prospectively sealed execution sources. An independent
assistant review rechecked recorded correctness, timing spans and arm balance;
this is internal review, not external replication. Updated manuscripts pass
Pandoc citation processing; previous PDF renders predate these revisions.

## Subsequent broader validation

The [competitor/compression/memory extension](broader-record-study.md) adds
pyfastx plain/gzip and isolated memory limits. Its separate experiment confirms
batched random short-record gains over pyfastx on this host, while exposing a
substantial million-record memory disadvantage and long-record timing losses.
Do not merge its smaller query cohorts numerically with this experiment.

## What remains before broader claims

This strengthens a scoped local comparison; it does not make the project
submission-ready. The subsequent study adds pyfastx, compressed-input and constrained-memory
controls. Independent workloads on another host and other appropriate
alternatives remain necessary for stronger generalization. A researcher-developed
workflow and final release reproduction are distinct from correctness testing
on public examples. The manuscripts must keep those evidence gaps explicit.
