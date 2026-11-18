# Prospective persistent-reader experiment

Written 2026-09-22 before the new persistent-reader measurements. This is a
prospective local protocol, not a public preregistration. Earlier CLI results
informed the questions. Preserve this file's hash in the experiment manifest.

**Question:** do observed record-service benefits survive keeping an established
indexed reader open? Correct output is required; a TracEon speedup is not.

**Primary contrast:** persistent TracEon batch-64 versus persistent HTSlib 1.24
handles, four independent worker processes, changing random record IDs. Use the
same compiled C++ driver, ordered complete sequence/quality results, identical
query lists, and the same four-CPU workload budget including the TracEon server.
Each indexed worker owns a reusable handle. TracEon workers keep TCP connections
open. This is a native indexed-reader control, not an HTSlib network server.

**Inputs:** prior E. coli short/long examples; independently obtained biological
Drosophila RNA-seq teaching subset (Zenodo 61771, GSM461176); deterministic
synthetic 100,000- and 1,000,000-record FASTQ fixtures (150 bases, varied qualities,
unique IDs). Synthetic fixtures test identifier-count/byte stress, not biological
replication. The million-read HTSlib case is explicitly a stress case beyond
samtools' recommended FASTQ scope, not the sole competitive baseline.

**Matrix:** 1/4/8 processes × batch 1/64 with random queries; at four processes,
also batch 16/256 and source-ordered/distributed-hotspot batch-64 queries.
Each worker requests 8,192 records (4,096 for long reads). Six independently
restarted sessions per dataset, randomized cell order with a recorded seed,
balanced arm order across sessions. Independent sessions are the replication
unit; workers and requests are not treated as independent experiments.

**Measurement:** readers initialize and prepare ID lists before a READY/RUN
barrier. Primary active-window time spans the first worker's native start through
the last native finish and includes SHA-256 consumption of every returned field.
Secondary per-request timings exclude digest computation and include ownership
release. Digest framing includes ID, presence and length-prefixed sequence/quality,
validated independently in Python. Keep raw request times; avoid strong p99
claims when a batch size yields few requests. Run one unmeasured warmup and one
timing pass without a memory sampler, followed by a separate memory pass using
the same IDs. Sample process PSS before, during and after memory passes; include
the server, retain process CPU and server CPU. Active-window timing excludes
writing latency files. Compare memory-pass timing as an observer sensitivity
check. Page-cache memory is outside process PSS.

**Setup:** build each cache and file index three times, alternating order, then
retain one artifact. Record size/hashes and construction time/CPU/RSS. Measure
server/reader initialization each session/cell, not in steady-state retrieval.
Report signed setup differences. Any projected break-even uses session-level
paired savings and stays explicitly a projection; if intervals include zero,
do not present a finite break-even as established.

**Analysis:** report all conditions, primary paired effect and bootstrap intervals
over six session-level differences (seed fixed), raw ranges and ratios. Intervals
are descriptive/exploratory with six sessions, not confirmatory significance
tests. Dataset-specific results and failures remain visible. No post-hoc winning
subset is the main result. CPU placement is controlled but not exclusive, warm
OS cache is intended, and host load/swap state are recorded; one host is not
cross-machine validation. Avoid cold-cache, memory-pressure, distributed-service,
novelty or universal performance claims.

**Correctness smoke:** exact IDs including colon/slash, case, varied qualities,
missing IDs, repeated requests, batch boundaries, and deliberate wrong-digest
rejection under normal and optimized Python. Preserve errors and do not silently
drop failed observations. Out-of-scope source formats unsupported by an indexed
baseline are reported, not forced into an invalid equivalence comparison.
