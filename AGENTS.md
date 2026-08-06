# AGENTS.md — TracEon Project Root

This is the DOX rail for TracEon: High-Performance Genomic Data Cache. It establishes project-wide contracts, build standards, and workflow rules.

## Purpose

TracEon is a zero-copy, lock-free C++20 genomic data caching library. This DOX tree ensures all work maintains:
- Architectural integrity (immutable-after-load, lock-free reads)
- Performance invariants (zero-copy, pre-reserved maps, SIMD dispatch)
- Build and test contracts (Release mode mandatory, Catch2 test structure)
- Documentation standards (ADRs for major decisions, performance profiles)

## Ownership

**Maintainer**: Adnan Raza (@woosflex)  
**Primary Domains**: C++20 core implementation, performance optimization, architecture decisions

## Global Contracts

### Build Standard
- **CMake 3.20+**, **C++20 required**
- **Release mode is MANDATORY** for any performance measurement or testing. Debug builds give misleading results.
- Dependencies: zlib-ng v2.2.2, ankerl::unordered_dense v4.4.0, Catch2 v3.6.0 (all via FetchContent or vendored)
- Runtime dispatch macros: `TRACEON_HAS_AVX2` (x86_64), `TRACEON_HAS_NEON` (ARM64) set by CMake

### Core Invariants (Non-Negotiable)
1. **Zero-copy**: All `getView()` calls return `std::string_view` directly into arena memory
2. **Immutable after load**: Data never changes post-parse; this enables lock-free reads
3. **Lock-free reads**: Single `std::atomic<bool> data_ready_` (acquire/release) guards load completion; no mutexes on read path
4. **Pre-reserved maps**: Multithreaded parsers pre-reserve `ankerl::unordered_dense` maps to avoid mid-parse rehashing
5. **GZIP: pre-allocate + direct-write**: zlib-ng decompresses directly into `text_arena_` (pre-sized at `compressed_size × 3` with OOM guard)
6. **Binary cache format versioning**: `.traceon` v3 format (`"TRO\x03"`, LZ4 Frame, streamed) is what `saveBinary()` writes; v2 (`"TRO\x02"`, LZ4 block) and v1 (`"TRO\x01"`, uncompressed) still load for backward compat
7. **Parallel GZIP**: Concatenated GZIP streams are independently decompressed using zlib-ng `inflate()` API; single-stream files use unmodified `loadGzipSingleStream()` path

### Workflow Standards
- All unit tests via Catch2 in `tests/` directory
- Tests must run from `build/` directory (relative path assumption)
- Performance validation via `benchmarks/benchmark_runner.py` before release
- Regression check via `benchmarks/check_regression.py` gates merges
- Architecture decisions documented as ADRs in `docs/architecture/`

## Work Guidance

### Before Editing Code
1. Identify which domain(s) you're touching (use Child DOX Index below)
2. Read the nearest child AGENTS.md for that domain
3. Review relevant ADRs in `docs/architecture/` if touching core logic
4. Verify CLAUDE.md build/test commands are still accurate after changes

### After Editing Code
1. Verify the change doesn't break any core invariant
2. Update relevant child AGENTS.md if scope/responsibility changes
3. Update CLAUDE.md if build/test commands are affected
4. If architectural change: create or update an ADR
5. Run full test suite from `build/` directory before committing

### Performance-Sensitive Areas
- `src/SmartStrategy.cpp` — Core arena, parsing, lock-free logic
- `include/SimdUtils.h` — SIMD dispatch (AVX2/NEON)
- `include/MapDefs.h` — Hash map selection and tuning
- GZIP decompression buffer management in `src/SmartStrategy.cpp`

## Child DOX Index

- **[tests/AGENTS.md](tests/AGENTS.md)** — Unit test framework (Catch2), test execution, verification contracts
- **[benchmarks/AGENTS.md](benchmarks/AGENTS.md)** — Performance benchmarking, regression checking, Python environment
- **[src/AGENTS.md](src/AGENTS.md)** — Core implementation (SmartStrategy, GZIP, SIMD, lock-free logic) and headers (Cache, MapDefs, SimdUtils)
- **[docs/AGENTS.md](docs/AGENTS.md)** — Architecture documentation (ADRs, performance profiles, decision records)
- **[examples/AGENTS.md](examples/AGENTS.md)** — Usage examples and demonstrations
