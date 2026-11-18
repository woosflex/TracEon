# Historical TracEon performance notes

> **Status: quarantined historical evidence.** The removed tables were point
> measurements from older development/release lines and were not produced by
> the corrected publication-grade harness. Raw repeated samples, complete
> input manifests, equivalent-output validation, and clean-checkout
> reproduction are incomplete. They must not be used as current speedup,
> throughput, latency, memory, or production-readiness claims.

The latest published release is v2.3.0. The active candidate is 3.0.0-dev and
has no released performance guarantee. Current single-host studies are recorded
in the [historical rerun](docs/earliest-claims-rerun.md),
[mapped-cache investigation](docs/mapped-cache-study.md), and
[current mapped-storage validation](docs/mapped-storage-validation.md).
Local-library and TCP-service results concern different access paths.

Retained historical reports are available under `benchmarks/`, including
`benchmarks/v2.3.0_benchmark_report.md`. Their numerical values are provenance
only and require rerunning with the corrected methodology before citation.

## Required future methodology

A current benchmark must:

1. pin exact source/tool revisions and input accession URLs plus SHA-256 values;
2. fail on every non-zero subprocess result;
3. use the same valid identifiers and verify equal hit counts and payload hashes;
4. separate ingest/build, restore/open, sequential iteration, and random lookup;
5. distinguish borrowed-view results from equivalent owned-output work;
6. retain raw warmup and repeated samples and report a defined statistic; and
7. record CPU, RAM, storage, kernel, compiler, affinity/governor, and page-cache policy.

Until that run exists, use the benchmark executables only as diagnostics.
