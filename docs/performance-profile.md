# Performance profile — validation required

> **Current status:** no publication-grade performance profile is established
> for 3.0.0-dev. Earlier throughput thresholds were based on
> historical point measurements without sufficient raw provenance and have
> been removed from current guidance.

Use `read_scaling`, `remote_bench`, and the benchmark programs only for local
diagnostics. Do not interpret their output as a regression gate or comparative
claim until the harness has pinned inputs, exact revisions, equivalent output
validation, raw repeated samples, and a declared environment/page-cache policy.

See [`PERFORMANCE.md`](../PERFORMANCE.md) for the required future methodology and
[`benchmarks/v2.3.0_benchmark_report.md`](../benchmarks/v2.3.0_benchmark_report.md)
for explicitly historical provenance.
