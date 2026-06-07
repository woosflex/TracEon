# Changelog

All notable changes to TracEon will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased] — v1.1.0 "Bakuya" (in progress)

### ✨ Added

- **zlib-ng Integration**: Replaced system zlib with zlib-ng v2.2.2 via CMake FetchContent in zlib-compat mode. Provides ~10-15% faster inflate through SIMD-optimized CRC and decompression routines.
- **Pre-Size with OOM Guard**: Pre-allocate `text_arena_` to `min(compressed_size × 3, available_memory × 0.25)` before decompression. Uses `std::filesystem::file_size()` for accurate compressed size and platform-specific memory queries (`/proc/meminfo`, `sysctl`, `GlobalMemoryStatusEx`). On OOM, throws `std::bad_alloc` with diagnostic.
- **Direct-Write Decompression**: Eliminated `temp_buffer` entirely. Decompressed data written directly into `text_arena_` with geometric growth fallback. Retains 1MB chunk buffer for `gzread` source reads.
- **shrink_to_fit**: Release excess arena capacity after decompression to minimize steady-state memory.

### 🔧 Changed

- **GZIP Load Time**: Reduced from 0.28s to 0.251s for 100MB GZIP (-29ms, 10% improvement)
  - zlib-ng SIMD inflate: -25ms
  - Pre-size (no realloc): -3ms
  - Direct-write (no move): -1ms
- **GZIP Memory Peak**: Reduced from ~366MB to ~266MB for 100MB GZIP (27% reduction) by eliminating `temp_buffer` co-allocation
- **Removed `third_party/zlib/`**: zlib-ng is now fetched at build time via FetchContent
- **Update `third_party/` notes**: zlib-ng no longer vendored in-tree

### 🐛 Fixed

- No bug fixes in this phase — pure performance optimization

### 🧪 Testing

- Validated GZIP load times match v1.1.0 targets (0.251s ± 0.01s)
- Memory profiling confirmed 266MB peak RSS for 132MB compressed file
- All existing GZIP test cases pass (auto-detect, explicit load, magic bytes)

---

## [1.0.1] — Patch Release

### 🐛 Fixed (Critical)

- **`addEntry()` dangling string_views (UB)**: `SequenceView` members were views into
  function-parameter strings that went out of scope immediately. Replaced with
  `std::deque<std::string> manual_store_` whose elements have reference stability
  (push_back never invalidates existing references per the C++ standard).

- **`loadBinary()` mmap with no error handling**: `open()` / `fstat()` / `mmap()` failures
  on Linux, and `CreateFileA()` / `MapViewOfFile()` failures on Windows, all went
  unchecked — leading to reads from invalid/NULL pointers (crash/UB). Full error
  propagation added via `std::runtime_error`.

- **`loadBinary()` no bounds checking**: Reading from the mmap buffer advanced `ptr`
  without verifying it stayed within `[data, data+size)`. A corrupt or truncated binary
  cache file would cause out-of-bounds reads. A `safe_advance(n)` lambda now throws
  on any overrun, and a record-count sanity cap (1B) prevents runaway loops.

- **Unaligned reads in `loadBinary()`**: `*reinterpret_cast<const uint32_t*>(ptr)` is
  undefined behaviour on strict-alignment architectures (ARM). All integer reads now use
  `std::memcpy` which the compiler optimises away on x86/x64 while remaining portable.

- **Multi-line FASTA sequences contained embedded newlines**: NCBI-style wrapped FASTA
  (60 chars/line) produced sequence `string_view`s with internal `\n` chars. A new
  `normalizeFastaArena()` pass compacts the arena in-place before any parser runs, so
  every sequence is guaranteed to be a single contiguous string. Works for both
  single-threaded and multithreaded parse paths.

### ✨ Added

- **Heterogeneous lookup on hot read path**: `GenomeIndex` now uses `StringHashMap`
  (from `MapDefs.h`) with a `TransparentStringHash` and `std::equal_to<>`. `getView()`
  can now call `map.find(sequence_id)` with a `string_view` directly — no `std::string`
  heap allocation per lookup. Requires C++20; falls back to `std::unordered_map` without
  robin_hood.

- **Binary cache format versioning**: Magic changed from `"MMAP"` to `"TRO\x01"`.
  Opening an old v1.0 binary cache now throws a descriptive `std::runtime_error` instead
  of silently producing garbage data.

### 🔧 Changed

- **GZIP buffer growth**: Removed the single upfront `3×compressed_size` reservation
  (inaccurate for highly compressible or incompressible data). The decompression loop
  now doubles `temp_buffer` capacity geometrically when needed, matching standard
  container growth behaviour.

- **Thread count guard**: Both multithreaded parser templates now use
  `std::max(1u, std::thread::hardware_concurrency())` to avoid division-by-zero on
  platforms where `hardware_concurrency()` returns 0.

- **`clearInternal()`**: Now also clears `manual_store_`, freeing strings added via
  `addEntry()` when the cache is reloaded or destroyed.

### 🧪 Testing

- 15 new assertions across 7 new test cases: multi-line FASTA normalization,
  `addEntry` stable-view correctness, `loadFile` / `loadBinary` on missing files,
  truncated binary, old-format rejection, initial state validation.
- Total: **43 assertions in 12 test cases** (was 28 in 6).

---

## [1.0.0] - 2025-12-16 "Avalon"

### 🎯 Codename
**Avalon** - Like EMIYA's absolute defense from Fate/stay night, this release provides performance isolation through lock-free concurrency and zero-copy architecture.

### ✨ Added

#### Core Features
- **Native GZIP Support** ([ADR-002](docs/architecture/ADR-002-gzip-integration.md))
  - Automatic format detection via `.gz` extension or magic bytes (`0x1f 0x8b`)
  - Streaming decompression with zlib (1MB chunks)
  - Minimal performance impact: < 7% overhead on lookups after decompression
  - Public APIs: `loadFile()` (auto-detect), `loadGzipFile()` (explicit)

- **Lock-Free Concurrent Reads** ([ADR-001](docs/architecture/ADR-001-lock-free-reads.md))
  - C++20 atomics with `memory_order_acquire`/`release` semantics
  - Immutable-after-load design eliminates all read locks
  - `std::atomic<bool> data_ready_` replaces `std::shared_mutex` on hot path
  - 2x potential throughput improvement on concurrent workloads

- **Robin Hood Flat Hash Maps**
  - Switched from `std::unordered_map` to `robin_hood::unordered_flat_map`
  - Open addressing provides ~15% better cache locality
  - Pre-allocation strategy: Reserve 125% of estimated capacity to prevent rehashing

- **Memory-Mapped Binary Cache**
  - Custom `.traceon` format with magic bytes verification (`MMAP`)
  - Instant restoration via `mmap()` / `MapViewOfFile()` (sub-second for GB files)
  - Cross-platform support (Linux, macOS, Windows)

#### Architecture
- **Zero-Copy Data Access**
  - Single arena allocation (`std::vector<char>` text_arena_)
  - All sequences returned as `std::string_view` (pointer + length)
  - No heap fragmentation

- **Hybrid Key Storage**
  - `std::string` keys (leverages SSO for short IDs < 15 chars)
  - `std::string_view` values (zero-copy)
  - Trade-off: 24 bytes overhead per key for improved cache locality

- **Forced GenomeIndex Mode**
  - Disabled NGSIndex (hash-only) despite being implemented
  - Rationale: SSO cache locality outweighs hash speed for typical genomic IDs

#### Performance Optimizations
- **Multithreaded Parsing**
  - Automatic parallelization for files > 10MB
  - Scales linearly up to 8 cores
  - Thread-local caches merged without lock contention

- **Capacity Pre-Allocation**
  - Parsers reserve 125% of estimated record count
  - Prevents rehashing during insertion
  - Heuristics: ~100 bytes/record (FASTA), ~150 bytes/record (FASTQ)

### 📊 Performance

#### Benchmarks (Intel Core Ultra 5 125H, 16GB RAM)

**Small Datasets (10MB, L3 Cache Resident):**
- WGS Short Reads: 40-55M OPS/s
- PacBio Long Reads: 50-60M OPS/s
- Reference Genomes: 40-50M OPS/s

**Medium Datasets (100MB, RAM Latency):**
- WGS Short Reads: 12-18M OPS/s (plain), 11-17M OPS/s (GZIP)
- PacBio Long Reads: 25-35M OPS/s
- Reference Genomes: 15-25M OPS/s

**Large Datasets (500MB-1GB):**
- WGS Short Reads: 10-15M OPS/s
- PacBio Long Reads: 20-30M OPS/s
- Reference Genomes: 12-20M OPS/s

**Memory Efficiency:**
- 10MB dataset: ~25 MB RSS
- 100MB dataset: ~180 MB RSS (5x more efficient than BioPython)
- 500MB dataset: ~900 MB RSS

**Binary Cache:**
- Save: 0.08s (100MB), 0.25s (1GB)
- Restore: ~0.00s (instant via mmap)

### 🔧 Changed

#### API
- `loadFile()` now auto-detects GZIP format (extension + magic bytes)
- Renamed internal `clear()` → `clearInternal()` to fix deadlock
- Made `getView()` lock-free (removed `std::shared_lock`)

#### Internal
- Replaced `std::shared_mutex` with `std::atomic<bool>` on read path
- GZIP decompression uses streaming with geometric buffer growth
- Optimized chunk size: 1MB for GZIP decompression

### 🐛 Fixed

#### Critical
- **Lock Overhead**: Eliminated 25-35% of execution time on read-heavy workloads
- **Performance Regression**: Fixed debug build masking as release (6M → 40M OPS/s)
- **HashMap Load Factor**: Pre-allocation prevents rehashing (was causing 20% degradation)
- **Deadlock in clearCache()**: Refactored to use internal helper

#### Minor
- Missing error handling for `tellg()` failures in GZIP size estimation
- Race condition in format auto-detection (moved detection outside lock)
- Test path dependencies (tests require running from build directory)

### 📚 Documentation

#### Added
- **ADR-001**: Lock-Free Reads (memory ordering justification)
- **ADR-002**: GZIP Integration (design decisions, rejected alternatives)
- **Performance Profile** with regression thresholds and troubleshooting
- **Comprehensive Benchmark Guide** (`benchmarks/README.md`)
- Single-include header (`TracEon.h`) with Doxygen comments
- Working example code (`examples/simple_usage.cpp`)

#### Updated
- README with GZIP feature highlights and accurate system specs
- Architecture diagrams showing lock-free data flow
- API documentation for `loadFile()` / `loadGzipFile()`

### 🧪 Testing

#### Added
- Lock-free concurrency tests (4-thread validation)
- GZIP test suite (8 assertions):
  - Explicit load (`loadGzipFile()`)
  - Auto-detect via extension (`.gz`)
  - Auto-detect via magic bytes (`0x1f 0x8b`)
- Zero-copy memory semantics verification
- Architecture validation (GenomeIndex usage)

#### Coverage
- 100% of public API methods
- Edge cases: empty files, malformed records, truncated streams
- Concurrency: Multi-threaded read validation

### 🔄 Migrations

**No breaking changes** - v1.0.0 is the initial public release.

### 🏗️ Build System

#### Dependencies
- **Required**: CMake 3.20+, C++20 compiler (GCC 10+, Clang 12+, MSVC 2019+)
- **Bundled**: 
  - robin_hood.h v3.11.5 (MIT License)
  - zlib 1.2.11 (zlib License)
  - lz4 1.9.4 (BSD License)
- **Test**: Catch2 v3.6.0 (fetched via CMake FetchContent)

#### Platforms
- ✅ Linux x86_64 (primary development platform)
- ✅ macOS ARM64 (Apple Silicon M1/M2/M3)
- ⚠️ Windows x64 (experimental, MSVC 2019+ required)
- 🔄 ARM64 Linux (untested but should work)

### ⚠️ Known Limitations

1. **RAM Latency Wall**: 3-4x degradation when dataset exceeds L3 cache (~16-32MB)
   - This is **hardware physics** (100ns RAM vs 10ns cache), not a software bug
   - Mitigation: Binary cache reduces load time, not lookup time

2. **Immutable After Load**: Cannot incrementally update cache
   - Required for lock-free reads
   - Alternative: Reload entire cache (fast with mmap)

3. **Single-Node Only**: No distributed caching
   - Target: Local workstation/server workloads
   - Future: v2.1.0 distributed architecture

4. **Test Path Dependency**: Tests must run from build directory
   - Workaround: `cd build && ./unit_tests`

5. **RecordTypes.h Structure**: `SequenceView` defined in `SmartStrategy.h`
   - Architectural decision for internal implementation details
   - Public API abstracts this via `Cache.h`

---

## Version History

| Version | Date | Codename | Key Features |
|---------|------|----------|--------------|
| 1.0.0 | 2025-12-16 | Avalon | Zero-copy, Lock-free, GZIP support, Memory mapping |
| 1.1.0 | 2026-Q1 | Bakuya | zlib-ng integration, Pre-size + Direct-write, Optimized GZIP |
| 1.2.0 | 2026-Q2 | Caladbolg | Binary cache compression |
| 2.0.0 | 2026-Q3 | Durandal | C API, Streaming |
| 2.1.0 | 2026-Q4 | Excalibur | Distributed caching |

---

## Comparison with State-of-the-Art

### vs BioPython
- **Performance**: 26-54x faster on long reads
- **Memory**: 5x more efficient (180MB vs >900MB for 100MB dataset)
- **Concurrency**: Lock-free vs single-threaded

### vs PyFastX
- **Architecture**: Memory-resident vs disk-cached
- **Long Reads**: 26-54x faster (TracEon excels with sparse indexing)
- **Short Reads**: Competitive (disk caching can be faster for specific workloads)
- **Concurrency**: Native multi-threading vs limited parallelism

---

## References

### External Resources
- [C++ Memory Model](https://en.cppreference.com/w/cpp/atomic/memory_order)
- Preshing, *Acquire and Release Semantics* (2012)
- Boehm & Adve, *Foundations of the C++ Concurrency Memory Model* (2008)
- [Robin Hood Hashing](https://github.com/martinus/robin-hood-hashing)
- [What Every Programmer Should Know About Memory](https://people.freebsd.org/~lstewart/articles/cpumemory.pdf) (Drepper, 2007)

### Comparison Tools
- BioPython: https://biopython.org/
- PyFastX: https://github.com/lmdu/pyfastx
- SeqKit: https://bioinf.shenwei.me/seqkit/

---

**Status:** ✅ Production Release  
**Next Version:** v1.2.0 "Caladbolg" (Q2 2026)

*"Trace On" - Projecting legendary performance from genomic data across eons.* ⚔️