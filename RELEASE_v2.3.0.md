# TracEon v2.3.0 "Harpe" — Release Notes

**Date:** 2026-08-20  
**Release type:** Minor — remote read access, read-scaling validation, and API hardening

## Highlights

### Remote read access

TracEon can now serve an immutable v4 cache over a bounded, length-prefixed
TCP protocol. The POSIX-only `traceon_remote` library provides `TraceonServer`
and `TraceonClient`; responses are CRC32C-verified and the protocol supports
handshake, lookup, presence, statistics, and clean disconnect messages.

The protocol has no authentication or encryption. Bind it to loopback or use
it only on a protected network; it is not suitable for direct internet
exposure.

Deployment support includes:

- `traceon_make_cache` for creating v4 caches from FASTA/FASTQ input
- `remote_bench` for local, server, and remote benchmarking
- A Docker Compose testbed under `docker/`
- [ADR-006](docs/architecture/ADR-006-traceon-remote-access.md), which defines
  the wire protocol, lifecycle, and trust boundary

### Multi-core reads

The new `read_scaling` benchmark measures concurrent `getView()` throughput
against one immutable cache. On the documented Intel Core Ultra 5 125H host,
reads scaled nearly linearly through eight threads (about 7×) and reached
about 9.5× at 14 threads. Hot read-path atomics are now isolated on cache
lines, preventing false sharing in lifecycle diagnostics.

### API and integrity hardening

- `Cache` now exposes `loadGzipFile()`, `clearCache()`, `getAllKeys()`,
  `getQuality()`, and `getDetectedFormat()`.
- Explicit `loadGzipFile()` calls reject non-GZIP input and retain documented
  failure-atomicity behavior.
- `kmerindex_reserve()` rejects implausible sizes before allocation.
- Remote shutdown, hostname resolution, and `remote_bench` option parsing were
  corrected before release.

## Compatibility

The release is additive. The v4 `.traceon` format and existing core API remain
compatible with v2.0–v2.2. Remote access is optional and is disabled on
Windows.

## Verification

- Core, k-mer, C linkage, and remote tests are registered with CTest.
- Seven libFuzzer targets cover binary-cache, GZIP, FASTA/FASTQ, k-mer, and
  TRKI input surfaces.
- The standalone read-scaling benchmark verifies hits and checksums while
  measuring aggregate throughput.

Reproduced Release-build measurements and exact commands are recorded in the
[v2.3.0 benchmark report](benchmarks/v2.3.0_benchmark_report.md).

See [CHANGELOG.md](CHANGELOG.md) for the complete change list.
