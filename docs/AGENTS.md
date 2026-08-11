# docs/AGENTS.md — Architecture & Design Documentation

This domain owns architecture decision records (ADRs), performance profiles, and design documentation.

## Purpose

Document and justify major architectural decisions, performance characteristics, and design trade-offs. Serve as reference for understanding why the codebase is structured as it is.

## Ownership

Maintainer creates/updates ADRs when architectural decisions are made. ADRs are immutable records (once closed, they don't change).

## Local Contracts

### Documentation Standards
- **ADR Format**: Architecture Decision Record (title, status, decision, alternatives considered, consequences, risks)
- **Performance Profile**: Expected throughput ranges, regression thresholds, hardware notes
- **Audience**: Future contributors, reviewers, and maintainers

### Current ADRs
- **ADR-001: Lock-Free Reads** — C++20 atomics, memory ordering, acquire/release semantics
- **ADR-002: GZIP Integration** — zlib-ng, pre-allocation heuristics, streaming decompression
- **ADR-003: SIMD Parsing & Hash Map Optimization** — SIMD boundary scanning, ankerl::unordered_dense, pre-reserved maps
- **ADR-004: Parallel GZIP & LZ4 Binary Cache** — concatenated-stream parallel decompression, LZ4 binary cache compression (v1.3.0)
- **ADR-005: `.traceon` v4 Binary Format & CRC32C** — v4 header layout, whole-payload CRC32C coverage and streaming order, hardware dispatch, reader hardening + failure atomicity (v2.0.0)

### Performance Profile
- **Location**: `performance-profile.md`
- **Specifies**: Expected OPS/s ranges by dataset size/type, regression thresholds (green 10%, yellow 10–20%, red >20%)
- **Hardware**: Always include processor, RAM, OS, compiler in results
- **Updates**: Revised when baseline changes significantly (new architecture version, compiler upgrade, hardware platform)

## Work Guidance

### Creating a New ADR
1. Use numbered sequence (ADR-004, ADR-005, etc.)
2. File: `docs/architecture/ADR-NNN-short-title.md`
3. Structure:
   ```
   # ADR-NNN: Title
   **Status**: Proposed/Accepted/Superseded
   **Decision**: One paragraph summary
   **Alternatives Considered**: List with brief rationale for rejection
   **Consequences**: Impact on codebase, performance, maintenance
   **Risks**: Known limitations or edge cases
   ```
4. Commit message: Link ADR in PR description
5. Once merged: Status = Accepted (immutable thereafter)

### Updating Performance Profile
1. Run full benchmark matrix (benchmarks/benchmark_runner.py)
2. If results differ >10% from expected: update expected ranges in `performance-profile.md`
3. Document hardware used (CPU, RAM, OS, compiler, -O3 flags)
4. Include date of measurement
5. Note if architectural change caused shift (e.g., "ankerl::unordered_dense reduced load time 86%")

### Referencing ADRs in Code
- Use comment with ADR link in performance-critical sections:
  ```cpp
  // See ADR-001: Lock-free reads with acquire/release semantics
  if (data_ready_.load(std::memory_order_acquire)) { ... }
  ```
- Reference in CLAUDE.md and child AGENTS.md where architectural invariants apply

## Verification

- All open ADRs are readable and internally consistent
- Performance profile matches observed benchmarks ±10%
- No contradictions between ADRs (e.g., ADR-001 and ADR-003 must align on concurrency model)

## Child DOX Index

- **[architecture/AGENTS.md](architecture/AGENTS.md)** (if needed) — ADR governance and naming conventions

For now, no child AGENTS.md; flat structure is sufficient.
