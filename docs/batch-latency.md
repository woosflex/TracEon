# Batch request latency

Development profiles, 2026-09-22. Each request retrieves 64 IDs from one resident
snapshot. Independent clients have no local cache copy. Each configuration uses
six repetitions of 2,000 measured requests per process, preceded by 1,000 warmup
requests per process. Server and clients share CPUs 8–11; the controller uses 16.
CPU frequency is uncontrolled and CPUs are not exclusive. The data and binaries
are identified by hashes in the manifests.

The values below are medians of the six per-run latency statistics, in
microseconds **per batch request**. They include batch-key construction, protocol
checks, owned payload creation and destruction. Handshake and warmup are excluded.
Within each process, requests are sequential; different processes start their
native timing windows independently. This is a closed-loop workload, not a
measurement of queueing delay under a fixed external arrival rate.

| Dataset | Processes | Median | p95 | p99 |
|---|---:|---:|---:|---:|
| Short FASTQ | 1 | 43.44 | 49.38 | 59.73 |
| Short FASTQ | 2 | 44.96 | 51.88 | 61.68 |
| Short FASTQ | 4 | 42.06 | 48.23 | 55.25 |
| Long FASTA | 1 | 427.53 | 741.26 | 1,132.92 |
| Long FASTA | 2 | 646.88 | 917.87 | 1,575.02 |
| Long FASTA | 4 | 1,169.22 | 1,706.65 | 2,498.97 |

Short batches contain 12,800 sequence-plus-quality bytes. Long batches average
approximately 560–570 kB, depending on sampled IDs. Their higher latency and
degradation with concurrency show why a fixed batch size is not a universal
optimum. Per-batch percentiles divided by 64 are not individual-record latency
percentiles. Keep batching opt-in and choose it for the demonstrated workload.

At four processes, median server CPU per cohort, including warmup, was 0.23 s
for short reads and 8.58 s for long reads. Each such cohort comprises 12,000
requests including warmup, or 768,000 returned record instances. Linux CPU
accounting has a 10 ms tick. Peak sampled aggregate process PSS was 18,073 KiB
and 149,855 KiB respectively; this includes the server, clients' ID pools and
timing vectors, but excludes the controller and OS page cache. Five-millisecond
sampling can miss brief peaks.

Both profiles completed with valid request counts, CRC-checked responses, no
missing IDs, and clean server shutdown. Exact source-payload and mapping parity
are covered separately by the [mapping study](batch-study.md); the timed profiler
checks counts and byte totals rather than comparing every payload to the source.

[Raw profiles and harness sources](results/batch-profiles-2026-09-22.tar.gz):
SHA-256 `ca81c2fa51520006680fc47b6bb69dd9c04059edd2de2f88016fb0c36016c692`.
This is development evidence; release-revision reproduction remains open.
