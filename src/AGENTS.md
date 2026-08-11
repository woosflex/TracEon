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
- `SmartStrategy.cpp` — Arena, lock-free data_ready_ flag, multithreaded parsing, GZIP decompression + integrity validation (gzerror/gzoffset), binary cache (v4 streaming LZ4 Frame + whole-payload CRC32C, no legacy readers, size/count hardening + failure atomicity, mmap), SIMD dispatch
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
2. **Immutable**: No changes to `text_arena_` or map after `data_ready_ = true`; enforced via `data_loaded_` — `addEntry()`/`set()` after a load throws `std::logic_error`, and `clearCache()` reopens the build phase
3. **Lock-free**: `data_ready_` is the only inter-thread synchronization; reads use acquire semantics, load completion uses release
4. **Pre-reserved maps**: Parsers call `ankerl::unordered_dense::map::reserve(capacity)` before inserting to avoid rehashing
5. **GZIP: allocate-then-decompress**: Pre-allocate `text_arena_` at `compressed_size × 3`, decompress directly in, shrink_to_fit when done
6. **Binary cache format versioning**: v2.0.0 ships **v4 as the only binary format** (`"TRO\x04"`: magic, codec flags, index mode, logical length u64-LE, frame length u64-LE, CRC32C u32-LE, then a streamed LZ4 Frame). `loadBinary()` rejects v1/v2/v3 with 'unsupported cache version; regenerate with v2.0.0' — there are no legacy readers. `saveBinary()` always writes v4. CRC32C covers the ENTIRE uncompressed payload + canonical header fields, updated incrementally as chunks stream through LZ4F (no second pass). Hardware dispatch: x86 SSE4.2 / AArch64 crc32 / table fallback
7. **Strict FASTQ framing (one terminator per boundary)**: Both FASTQ parsers consume a strict 4-line positional cycle (header / seq / '+' / qual) and advance exactly one line terminator per boundary — collapsing newline runs would skip empty sequence/quality lines. Empty quality lines parse as `qual=''`; empty sequence lines no longer break the cycle. The multithreaded chunk-boundary classifier uses a FORWARD check (genuine header iff line+2 begins with '+') with the same one-terminator rule
8. **GZIP integrity (gzerror/gzoffset)**: `loadGzipFile()` throws `std::runtime_error` when `gzread()` returns 0 with a non-`Z_OK` `gzerror()` (truncated stream), when `gzoffset()` < file size (trailing garbage after the last member; `scanGzipStreams()` gives exact coverage for files ≥ 1 MiB), and when the parallel path hits `Z_BUF_ERROR` with `avail_in == 0` (truncated tail member). Partial data is never served as complete
9. **Binary-cache hardening (v4 reader)**: implausible sizes are rejected before allocation (undersized logical length < 8, 32-bit `size_t` guard, OOM guard + `bad_alloc` catch); all pointer bounds checks use the overflow-safe subtraction form (`n > end - ptr`); record counts are bounded by remaining payload bytes (≥12 B/record GENOME, ≥20 B NGS); decompression requires EXACT logical length and EXACT frame termination; whole-payload CRC32C must match — a modified/truncated/crafted cache fails loudly. **Failure atomicity**: an invalid load never publishes `data_ready_` or partial map state (cache resets to pristine empty)

### Performance-Sensitive Code Paths
- `SmartStrategy::loadFile()` — Parse initiation, GZIP detection, chunking for multithreaded parsing
- `SmartStrategy::loadGzipFile()` — zlib-ng decompression loop, buffer pre-sizing, OOM guarding
- `SmartStrategy::loadGzipInternal()` — Dispatcher: detects stream count, routes to parallel or single-stream path
- `SmartStrategy::scanGzipStreams()` — O(n) scan for GZIP stream boundaries; drives parallel dispatch decision
- `SmartStrategy::loadGzipParallel()` — Per-stream inflate() in parallel threads; merges into text_arena_
- `SmartStrategy::loadBinary()` — v4-only: header validation (magic/version/codec/mode), streaming LZ4 Frame decompression with exact-length + exact-frame-termination checks, streaming CRC32C verification, mmap + parsing, record parsing with count-bounded reserve
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
5. FASTQ/FASTA parser changes must preserve the one-terminator-per-boundary rule in both parse paths AND the chunk-boundary classifier; the classifier is FORWARD-looking (candidate `@` is a genuine header iff line+2 begins with `+`) by contract
6. GZIP changes must preserve integrity rejection: consult `gzerror()` when `gzread()` returns 0, verify `gzoffset()` equals file size at EOF, and treat `Z_BUF_ERROR` with `avail_in == 0` as truncation in `loadGzipParallel()`

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
1. Understand the v4 format (header layout, CRC32C coverage, LZ4 Frame streaming) — v1/v2/v3 are REJECTED, never read
2. Binary cache format versioning is critical: bumping magic version byte is backward-incompatible for v1 readers
3. When deserializing: always check magic version byte first, route to appropriate decompressor
4. When serializing: save to v4 format (streaming LZ4 Frame + CRC32C — the only format `saveBinary()` writes)
5. Test round-trip: save → restore → verify all sequences match
6. Benchmark: ensure compression ratio is ~3x for typical data and decompression overhead is <0.03s
7. Sizes/counts read from the file are attacker-controlled: reject implausible logical/frame lengths before allocation, bound record count by remaining payload bytes (≥12 B GENOME / ≥20 B NGS per record), keep all bounds checks in subtraction form (`n > end - ptr`), and verify the whole-payload CRC32C before publishing data

### Adding New Features
- Keep arena allocation model (no new malloc in hot paths)
- Maintain immutable-after-load contract
- Extend API via Cache.h only (SmartStrategy is internal)
- Include unit tests in `tests/` and benchmark integration tests in `benchmarks/`

## Verification

- All unit tests pass (Catch2, run from build/ directory)
- No memory leaks (valgrind optional)
- Performance benchmarks within expected ranges (docs/performance-profile.md)
- Code follows existing C++20 patterns (const-correctness, RAII); exceptions are used only for error reporting on invalid/corrupt input (per the integrity contracts above)

## Child DOX Index

No child domains in src/. `include/` is logically paired with `src/` for API/implementation coherence.
