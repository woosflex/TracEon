# src/AGENTS.md — Core Implementation & Headers

This domain owns the core engine, zero-copy architecture, lock-free mechanisms, GZIP integration, SIMD dispatch, and public API definitions.

## Purpose

Implement and maintain:
- **Arena allocation and zero-copy access** (SmartStrategy, text_arena_)
- **Lock-free concurrent reads** (std::atomic, memory_order_acquire/release)
- **GZIP decompression** (zlib-ng integration, pre-sized buffers)
- **Binary cache persistence** (LZ4 compression, format versioning, mmap)
- **SIMD-accelerated boundary scanning** (AVX2, NEON, scalar fallback)
- **Efficient hashing** (ankerl::unordered_dense, pre-reserved maps)
- **Public API surface** (Cache, IEncodingStrategy)

## Ownership

Core implementation is owned by the project maintainer. Performance-critical sections require benchmarking validation before merge.

## Local Contracts

### File Responsibilities

**Implementation (src/):**
- `SmartStrategy.cpp` — Arena, lock-free data_ready_ flag, multithreaded parsing, GZIP decompression, binary cache (streaming LZ4 Frame, v1/v2/v3 format detection, mmap), SIMD dispatch
- `Cache.cpp` — High-level API (loadFile, getView, save, restore, getFastqRecord), routes to SmartStrategy
- `IEncodingStrategy.cpp` — Strategy interface/base class (minimal, mostly headers)

**Headers (include/):**
- `Cache.h` — Public API surface
- `SmartStrategy.h` — Core engine definition, arena/map access (for advanced users)
- `SimdUtils.h` — SIMD character-search function with runtime dispatch
- `MapDefs.h` — ankerl::unordered_dense typedefs and pre-reservation logic
- `RecordTypes.h`, `DecodedRecordTypes.h` — Type definitions
- `IEncodingStrategy.h` — Strategy interface
- `TracEon.h` — Single-include convenience header

### Core Invariants (Binding)
1. **Zero-copy**: `getView()` must return `std::string_view` pointing into `text_arena_` (never allocate)
2. **Immutable**: No changes to `text_arena_` or map after `data_ready_ = true`
3. **Lock-free**: `data_ready_` is the only inter-thread synchronization; reads use acquire semantics, load completion uses release
4. **Pre-reserved maps**: Parsers call `ankerl::unordered_dense::map::reserve(capacity)` before inserting to avoid rehashing
5. **GZIP: allocate-then-decompress**: Pre-allocate `text_arena_` at `compressed_size × 3`, decompress directly in, shrink_to_fit when done
6. **Binary cache format versioning**: `loadBinary()` detects magic byte version (`0x01` v1 uncompressed, `0x02` v2 LZ4 block, `0x03` v3 streaming LZ4 Frame) and routes to the appropriate decompressor; `saveBinary()` always writes v3

### Performance-Sensitive Code Paths
- `SmartStrategy::loadFile()` — Parse initiation, GZIP detection, chunking for multithreaded parsing
- `SmartStrategy::loadGzipFile()` — zlib-ng decompression loop, buffer pre-sizing, OOM guarding
- `SmartStrategy::loadGzipInternal()` — Dispatcher: detects stream count, routes to parallel or single-stream path
- `SmartStrategy::scanGzipStreams()` — O(n) scan for GZIP stream boundaries; drives parallel dispatch decision
- `SmartStrategy::loadGzipParallel()` — Per-stream inflate() in parallel threads; merges into text_arena_
- `SmartStrategy::loadBinary()` — Format version detection, LZ4 decompression (v3 streaming Frame, v2 block, v1 uncompressed), mmap + parsing, record parsing
- `SmartStrategy::saveBinary()` — Payload serialization, LZ4 compression, header + data write
- `SmartStrategy::getView()` — Must be O(1) map lookup; no allocations
- `simd_find_char()` in SimdUtils.h — Boundary scanning (newline, `>`, `@`, `+`); runtime dispatch critical
- `SmartStrategy::parseFasta[Threaded]()` — Map pre-reservation and insertion; SIMD boundary scanning integration
- `SmartStrategy::parseFastq[Threaded]()` — Same as FASTA plus quality score handling

## Work Guidance

### Before Modifying SmartStrategy
1. Review ADR-001 (lock-free reads) and ADR-003 (SIMD + hashing)
2. Verify changes don't introduce heap allocations in `getView()` path
3. Check that `data_ready_` ordering assumptions still hold (load = release, read = acquire)
4. For GZIP changes: review pre-allocation logic and OOM guard in `loadGzipFile()`

### Before Modifying SIMD (SimdUtils.h)
1. Verify runtime dispatch via `__builtin_cpu_supports()` is still sound
2. Test on both AVX2 and non-AVX2 hardware (or simulator)
3. Ensure scalar fallback is bit-identical to SIMD variant
4. Benchmark: SIMD should be ≥2x scalar on 100MB files

### Before Modifying Hash Map (MapDefs.h)
1. Verify `ankerl::unordered_dense` pre-reservation heuristic remains valid
2. If changing capacity estimates (`chunk_size / 500` for FASTA, `/600` for FASTQ): re-run benchmarks
3. Test multi-threaded insertion with pre-reserved maps (should not rehash mid-parse)

### Before Modifying Binary Cache (loadBinary / saveBinary)
1. Understand v1 (uncompressed), v2 (LZ4 block), and v3 (streaming LZ4 Frame) formats (see ADR-002 "Binary Cache Compression" and ADR-004)
2. Binary cache format versioning is critical: bumping magic version byte is backward-incompatible for v1 readers
3. When deserializing: always check magic version byte first, route to appropriate decompressor
4. When serializing: save to v3 format (streaming LZ4 Frame — the only format `saveBinary()` writes today)
5. Test round-trip: save → restore → verify all sequences match
6. Benchmark: ensure compression ratio is ~3x for typical data and decompression overhead is <0.03s

### Adding New Features
- Keep arena allocation model (no new malloc in hot paths)
- Maintain immutable-after-load contract
- Extend API via Cache.h only (SmartStrategy is internal)
- Include unit tests in `tests/` and benchmark integration tests in `benchmarks/`

## Verification

- All unit tests pass (Catch2, run from build/ directory)
- No memory leaks (valgrind optional)
- Performance benchmarks within expected ranges (docs/performance-profile.md)
- Code follows existing C++20 patterns (const-correctness, RAII, no exceptions)

## Child DOX Index

No child domains in src/. `include/` is logically paired with `src/` for API/implementation coherence.
