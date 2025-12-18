# Changelog

All notable changes to TracEon will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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

## [Unreleased]

### Planned for v1.1.0 "Bakuya" (Q1 2026)
*Twin swords - Dual optimization paths*

- SIMD prefetching for hash probes (target: 20M OPS/s on 100MB RAM datasets)
- Optimized GZIP buffering with better chunk management
- Adaptive capacity estimation for different file types
- Benchmark against production datasets from collaborators

### Planned for v1.2.0 "Caladbolg" (Q2 2026)
*Rainbow sword - Spectrum of compression*

- Binary cache compression using LZ4 (target: 3x size reduction)
- Transparent decompression for reference genomes
- Smart compression algorithm selection based on data characteristics

### Planned for v2.0.0 "Durandal" (Q3 2026)
*Peerless - Language-agnostic access*

- C API for Python/R/Julia bindings
- Streaming API for datasets larger than RAM
- Foreign Function Interface (FFI) with zero-copy guarantees

### Planned for v2.1.0 "Excalibur" (Q4 2026)
*Holy sword - Distributed power*

- Distributed caching across multiple nodes
- NUMA-aware architecture for >8 core systems
- Network protocol for shared genomic databases

---

## Version History

| Version | Date | Codename | Key Features |
|---------|------|----------|--------------|
| 1.0.0 | 2025-12-16 | Avalon | Zero-copy, Lock-free, GZIP support, Memory mapping |
| 1.1.0 | 2026-Q1 | Bakuya | SIMD prefetching, Optimized GZIP |
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
**Next Version:** v1.1.0 "Bakuya" (Q1 2026)

*"Trace On" - Projecting legendary performance from genomic data across eons.* ⚔️