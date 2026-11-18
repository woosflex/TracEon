# TracEon roadmap

Updated 2026-09-22 after historical and current local-cache validation.

## Direction

**Prepare once, open locally, query repeatedly from many processes.**

TracEon is an embedded FASTA/FASTQ record cache. Its core operation is retrieving
whole sequence payloads by identifier, with borrowed views for local C++ callers
and owned sequence/quality access when needed. An immutable mapped cache lets
independent processes share file-backed payload pages without a resident daemon.
Each process still owns its index and pays validation/index construction costs.

Optimize correct record delivery, repeated-query latency, preparation and restore
costs, and aggregate process memory against persistent indexed readers. Preserve
integrity, atomic replacement, and explicit borrowed-view lifetimes. Compression
remains a storage trade-off; the optional TCP service serves deployments needing
remote access. It is not the primary architecture or the default evaluation path.

The [direction audit](docs/project-direction-audit.md),
[earliest prototype rerun](docs/earliest-claims-rerun.md), and
[v1 mapped-cache study](docs/mapped-cache-study.md) explain this decision.
Historical results support a local reuse benefit and mapped payload sharing;
they do not establish universal superiority, a new alignment paradigm, or
replacement of mapper seeding/chaining/alignment engines.

## Current implementation slice

- Add explicit portable v5 mapped persistence to the current local API and cache
  builders, reusing the checked parser, CRC32C, and atomic writer.
- Retain compressed writes as the compatibility default; restore auto-detects
  storage. Do not revive unsafe historical readers or change local query semantics.
- Verify exact FASTA/FASTQ payloads, empty/missing records, both index modes,
  corruption/truncation rejection, conversions, and mapped-reader lifetime across
  atomic replacement. Exercise installation and existing compressed/remote paths.
- Compare mapped and compressed storage with the **same current implementation**,
  independent processes, full-content checks, preparation costs, and PSS.
  Historical v1/current differences are not a substitute for this control.

This slice is implemented and locally validated: [results and limitations](docs/mapped-storage-validation.md).
All work remains an uncommitted 3.0.0-dev candidate. See [unreleased changes](CHANGELOG.md).

## Next priorities

| Order | Work | Completion condition |
|---|---|---|
| Done locally | Validate the mapped local path | 23 regression checks pass; the same-library comparison verifies 40,439,940 delivered records and lower aggregate PSS across four mapped readers. Single-host limits remain. |
| 2 | Prove usefulness in a repeated-record workflow | A realistic workload repeatedly selects records from one stable dataset. Compare persistent HTSlib and pyfastx, include preparation/restore and useful computation, and retain losses. No mapper internals need replacing. |
| 3 | Make the local API accessible | Build a narrow binding for the workflow's actual language if C++ integration is insufficient. Specify ownership, missing IDs, and sequence/quality semantics before expanding the API. |
| 4 | Release a reproducible candidate | Review the uncommitted changes, settle compatibility, run CI/downstream checks at an exact revision, and package source, tools, and evidence. |
| 5 | Open-source release & community | Release reproducible candidate, update public documentation and benchmarks, and support community integrations. |

A second machine is unavailable. Mark all current measurements as single-host;
independent-host replication remains pending rather than silently assumed.

## Existing capabilities and evidence to preserve

- Immutable snapshots, concurrent local reads, complete FASTQ payload access,
  structural input validation, portable v5 persistence and v4 migration reads.
- Atomic POSIX replacement and fault-injection checks; relocatable CMake exports,
  core-only builds, lifecycle contracts, and downstream consumer verification.
- Optional TCP service with snapshot identities, bounded batches, deadlines,
  connection limits, and owned client records. Maintain its correctness and
  compatibility without making service expansion the main roadmap.
- [Persistent-reader service control](docs/resident-record-study.md),
  [broader competitor/compression/memory controls](docs/broader-record-study.md),
  [record reuse](docs/record-reuse-study.md), [batching](docs/batch-study.md),
  [service resources](docs/service-resource-study.md), and
  [native mapper control](docs/shared-minimap-index-pilot.md).
  These remain evidence about their measured access paths, including negative
  results. They do not settle local borrowed-view performance.

## Deferred

Mapper-engine replacement, experimental k-mer expansion, distributed sharding,
eviction, live mutable snapshots, new network transports, a broad binding suite,
and speculative SIMD tuning. Revisit only when the repeated-record workflow
identifies a measured need. A mapped artifact is immutable: update it through
atomic replacement, never by modifying or truncating a file held by readers.
