# TracEon: High-Performance Genomic Data Cache

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![C++20](https://img.shields.io/badge/C%2B%2B-20-blue.svg)](https://isocpp.org/)
[![Version](https://img.shields.io/badge/Release-v2.2.0-blueviolet.svg)]()
[![Performance](https://img.shields.io/badge/Performance-81M%20OPS%2Fs-brightgreen.svg)]()
[![GZIP](https://img.shields.io/badge/GZIP-Native%20Support-orange.svg)]()
[![BinaryCache](https://img.shields.io/badge/BinaryCache-v4%20CRC32C%20Integrity-green.svg)]()

**TracEon** (v2.2.0) is a **zero-copy, lock-free** genomic data caching library written in modern C++20. It accelerates bioinformatics pipelines by providing microsecond-latency random access to FASTA/FASTQ datasets, including native `.gz` support, CRC32C-protected LZ4 binary caches, and **workload-optimized index modes** for genomic and NGS applications.

*"Trace On" - A nod to Fate/stay night, re-contextualized for tracing biological data across eons.*

---

## 🎯 What's New in v2.2.0

- **Aligner validation matrix COMPLETE** — 5 real tools, all PAF byte-identical: minimap2, Winnowmap, mm2-fast, minigraph (first khashl adapter), BLEND. Every fork: opt-in `TRACEON=1` table backend + `tcache` mmap'd prebuilt index (zero-rebuild load). Recurring wins −13.6% to −80%; full numbers in `RELEASE_v2.2.0.md`

## 🎯 What's New in v2.1.0

- **Fuzzing gate**: 7 libFuzzer targets (v4 loader, gzip, FASTQ/FASTA, kmer encode, kmer C API, TRKI) + CI workflow — the last verification layer is now closed
- **kmer C API hardened**: atomic freeze flag (thread-safe concurrent lookup) + allocation-free error path (total no-throw boundary)
- **Real-world aligner validation**: `woosflex/minimap2` + `woosflex/Winnowmap` forks — `tcache v2` mmap'd prebuilt index. Recurring load+map: minimap2 −13.6% (full chr1), Winnowmap −80% (chr1-50Mb); PAF byte-identical to stock everywhere

## 🎯 What's New in v2.0.0

> **Breaking release.** v2.0.0 makes a clean break in the binary format and
> lifecycle contract (no public adopters were identified; see the
> [design review](outputs/traceon-v2-design-review.md)). Read the
> [Migration / Breaking Changes](#migration--breaking-changes) section before
> upgrading. The experimental k-mer C API ships in v2.0.0+ and is used by the
> aligner forks (minimap2/Winnowmap `traceon-backend`).

### **`.traceon` v4 Binary Format** 📦
- `saveBinary()` / `Cache::save()` now write **v4** (`"TRO\x04"`): a 26-byte
  header (magic, codec flags, index mode, logical payload length u64-LE,
  compressed frame length u64-LE, CRC32C u32-LE) followed by a streamed LZ4
  Frame.
- **Whole-payload CRC32C**: the checksum covers the entire uncompressed
  logical payload plus the canonical header fields, computed **incrementally**
  as chunks pass through the LZ4F compressor/decompressor — no second
  full-payload pass, no extra allocation.
- **Hardware-accelerated CRC32C** with runtime dispatch: x86 SSE4.2 `crc32`
  instruction, AArch64 `crc32cx/crc32cw`, and a portable table fallback
  (equivalence-tested).
- **v4 is the ONLY readable format.** v1/v2/v3 caches are rejected with
  "unsupported cache version; regenerate with v2.0.0". Corrupt, truncated,
  or modified caches are rejected: wrong magic/version, unsupported
  codec/mode, wrong logical length, checksum mismatch, implausible sizes,
  and count-bombs all fail loudly — never silently loading garbage.
- All v2/v3 hardening semantics are preserved in the v4 reader: implausible
  sizes rejected before allocation, subtraction-form bounds checks, count-
  bounded reserve, OOM guard + `bad_alloc` catch, and **failure atomicity**
  (an invalid load never publishes partial state).

### **Lifecycle Contract (Reader Quiescence)** 🔒
- `getView()` returns a **non-owning** `std::string_view` valid only while the
  same loaded snapshot is installed. It is **dangling** once `clearCache()`,
  reload, or destruction begins — stop using views before clearing/reloading.
- Concurrent lock-free reads against a live snapshot remain fully supported;
  the write-side mutex serializes writers only (it is not reader reclamation).
- `TRACEON_DEBUG_LIFECYCLE` (opt-in) adds a reader counter diagnostic.

### **Integrity Hardening** 🛡️ (carried from fix/p0-integrity)
- Strict FASTQ 4-line framing (empty seq/qual lines, `'@'`/`'+'` in quality).
- GZIP truncation / trailing garbage / truncated tail members throw instead
  of serving partial data.
- `set()`/`addEntry()` after a load throws `std::logic_error`.

### **Zero-copy + Performance** ⚡
- `std::string_view` index keys, transparent wyhash-style hashing, and
  pre-reserved thread-local maps (no allocation on insert/read hot paths).

---

## 🚨 Migration / Breaking Changes (v1.x → v2.0.0)

1. **Binary cache:** `.traceon` files are **v4 only**. v1/v2/v3 caches are
   rejected at load time with "unsupported cache version; regenerate with
   v2.0.0". Regenerate caches after upgrading:
   ```cpp
   TracEon::Cache old;              // old library
   // ... load + rebuild from source text/GZIP ...
   old.save("genome.traceon");      // now writes v4
   ```
2. **Input integrity:** truncated/corrupt GZIP and trailing garbage after the
   last GZIP member now throw `std::runtime_error`; partial data is never
   published.
3. **Mutation/lifecycle:** `set()`/`addEntry()` after a load throws
   `std::logic_error`; call `clearCache()` only after all readers stop using
   views.
4. **View/storage:** `getView()` remains zero-copy and non-owning; its result
   must not be used after `clearCache()`, reload, or destruction (see the
   [Lifecycle contract](#-lifecycle-contract-reader-quiescence-) below).
5. **Source/ABI:** `std::string_view`-backed index keys require recompilation;
   **no ABI stability is promised** across 1.x → 2.x.
6. **Deferred surface:** the experimental k-mer C API and `TRKI` are **not**
   shipped in v2.0.0 (planned separately).

---

## 🔒 Lifecycle Contract (Reader Quiescence)

`getView()` returns a `std::string_view` pointing directly into the immutable
snapshot's arena (`text_arena_` / `manual_store_`). Because the view is
**non-owning**, it is valid only while the same loaded snapshot stays installed:

- ✅ Concurrent lock-free reads against a **live** snapshot are safe.
- ❌ Using a view after `clearCache()`, a reload (`loadFile`/`loadGzipFile`/
  `restore`), or destruction is **undefined behavior** — the backing arena has
  been destroyed.
- The write-side mutex serializes writers; it does **not** reclaim memory
  from readers and is deliberately not taken on the read path.
- Sequential clear/reload cycles are fully supported and tested.
- For debugging coordinated misuse, build with `-DTRACEON_DEBUG_LIFECYCLE`:
  `clearInternal()` warns if a `getView()` lookup overlaps a teardown
  (diagnostic only — not synchronization).

See [ADR-001](docs/architecture/ADR-001-lock-free-reads.md) for the full
memory-ordering and quiescence contract.

### **v1.3.0 "Hrunting" Highlights** (Previous Release)

### **LZ4 Binary Cache Compression** 📦
- **Format v2** with LZ4 compression for `.traceon` binary cache files
- **3x size reduction** (~105MB → ~35MB for 100MB FASTA)
- **Backward compatible**: Auto-detection of v1 (uncompressed) and v2 (LZ4) formats
- **Minimal overhead**: +0.02s decompression (LZ4 @ 1-2 GB/s)
- Binary cache sizes now practical for large datasets and network transfer

### **v1.2.0 "Caladbolg" Highlights** (Previous Release)

### **SIMD-Accelerated Parsing** ⚡
- New `simd_find_char()` function with **AVX2** (32 bytes/iter), **NEON** (16 bytes/iter), and scalar fallback
- Runtime dispatch via `__builtin_cpu_supports("avx2")` — zero overhead on non-AVX2 CPUs
- Integrated into FASTA and FASTQ parsers for newline and record-boundary (`>`, `@`, `+`) scanning
- Biggest win on large-reference multi-line FASTA where the sequence scan dominates

### **ankerl::unordered_dense Hash Map** 📊
- Replaced `robin_hood::unordered_flat_map` with **ankerl::unordered_dense::map** (Swiss-table design)
- Better cache locality, **~0.8s faster** insert for 100MB datasets
- Pre-reserved thread-local maps in multithreaded parsers eliminate mid-parse rehashing
- Conservative heuristic: `chunk_size / 500` (FASTA), `chunk_size / 600` (FASTQ) with 1.25× safety margin

### **Reduced Memory Footprint** 📦
- **Peak RSS**: 185MB (30% reduction from 263MB in v1.1.0)
- Memory savings come from ankerl's more compact table representation vs Robin Hood

### **Fixed: normalizeFastaArena() Trailing Newline** 🐛
- `resize()` zero-initialized trailing `\n` when input lacked trailing newline, causing false trailing-data detection
- Fixed by growing vector by 1, writing `\n` past original end, then trimming

### **v1.1.0 "Bakuya" Highlights** (Previous Release)

### **zlib-ng Integration** ⚡
- Replaced system zlib with **zlib-ng v2.2.2** (SIMD-optimized inflate)
- Drop-in replacement via zlib-compat mode — API unchanged
- **10% faster GZIP decompression** (0.251s vs 0.28s for 100MB)
- Vendored via CMake FetchContent — no system dependency required

### **Optimized Buffer Pre-Allocation** 📦
- **Pre-size + Direct-Write**: Decompress directly into `text_arena_` — no intermediate `temp_buffer`
- **3x heuristic**: Pre-allocate `compressed_size × 3` to cover typical GZIP ratios
- **OOM Guard**: Caps pre-allocation at 25% of available memory
- **27% lower memory peak** during load (266MB vs 366MB for 100MB file)
- Geometric growth as fallback for under-estimation

### **Native GZIP Support** ✨ (v1.0.0, improved in v1.1.0)
- Automatic detection via file extension (`.gz`) or magic bytes (`0x1f 0x8b`)
- Streaming decompression preserves zero-copy architecture
- No performance degradation on queries: GZIP lookups as fast as plain text

### **Lock-Free Concurrent Access** 🔓 (v1.0.0)
- C++20 atomics with acquire/release semantics
- 2x throughput improvement on multi-threaded workloads
- Zero mutex contention on read-heavy operations

### **Production-Ready Performance** ⚡ (v1.0.0)
- **54x faster** than PyFastX on compressed long reads
- **60-80% memory reduction** vs BioPython
- Sub-second binary cache restoration via memory mapping

---

## 🚀 Performance Highlights

### System Specifications
**Test Platform:** Intel Core Ultra 5 125H (14 cores, 16GB LPDDR5 @ 7467MT/s)  
**Operating System:** Ubuntu 24.04 LTS  
**Compiler:** GCC 13.3.0  
**Build Configuration:** Release (`-O3`, C++20)

### Real-World Throughput (Synthetic Benchmark Data)

#### Short Reads (WGS Illumina)
| Metric | 10MB Dataset | 100MB Dataset | Notes |
|--------|--------------|---------------|-------|
| **Random Lookups** | 40-55M OPS/s | 12-18M OPS/s | L3 cache vs RAM latency |
| **GZIP Lookups** | 35-50M OPS/s | 11-17M OPS/s | Minimal overhead |
| **Memory Usage** | ~25 MB | ~185 MB | 5x efficient vs BioPython |
| **Load Time (plain)** | <0.05s | ~0.15s | Includes parsing |
| **Load Time (GZIP v1.2.0)** | <0.05s | **~0.245s** | SIMD + ankerl optimizations |

#### Long Reads (PacBio/Nanopore)
| Metric | 10MB Dataset | 100MB Dataset | Speedup vs PyFastX |
|--------|--------------|---------------|-------------------|
| **Random Lookups** | 50-60M OPS/s | 28-81M OPS/s | **26-54x faster** |
| **GZIP Lookups** | 45-55M OPS/s | 25-35M OPS/s | Maintains advantage |
| **Memory Usage** | ~20 MB | ~120 MB | Highly efficient |

**Key Insight:** Performance scales with dataset size relative to CPU cache hierarchy. 3-4x degradation when exceeding L3 cache (~16-32MB) is **hardware physics** (100ns RAM vs 10ns cache), not a software limitation.

### Binary Cache Performance
| Operation | 100MB Dataset | 1GB Dataset | Compression |
|-----------|---------------|-------------|------------|
| **Save Time (v2 LZ4)** | 0.10s | 0.30s | With LZ4 compression |
| **Restore Time (v2 LZ4)** | ~0.021s | ~0.15s | Decompress + parse |
| **Restore Time (mmap baseline)** | **~0.001s** | **~0.01s** | v1 uncompressed |
| **File Size (v1)** | ~105 MB | ~1.05 GB | Uncompressed |
| **File Size (v2)** | **~35 MB** | **~350 MB** | **3x smaller (LZ4)** |

**Note:** v1 (uncompressed) caches restore instantly via mmap. v2 (LZ4) adds ~0.02s decompression overhead but saves ~3x disk space. Binary cache format automatically detected at load time.

---

## 🏗️ Architecture Overview

```
┌─────────────┐
│ User Code   │
└──────┬──────┘
       │
       v
┌─────────────────────────────────────┐
│  Cache (Public API)                 │
│  - loadFile()   ← Auto-detect GZIP  │
│  - getView()    ← Zero-copy access  │
└──────┬──────────────────────────────┘
       │
       v
┌─────────────────────────────────────┐
│  SmartStrategy (Engine)             │
│  ┌─────────────────────────────┐   │
│  │ Arena Allocator             │   │  ← Single malloc
│  │ ┌─────────────────────────┐ │   │
│  │ │ ACGTACGT... (text_arena)│ │   │
│  │ └─────────────────────────┘ │   │
│  └─────────────────────────────┘   │
│                                     │
 │  ┌───────────────────────────────────┐   │
 │  │ ankerl::unordered_dense Map      │   │  ← Lock-free reads
 │  │ [string → SequenceView]          │   │     (C++20 atomics)
 │  └───────────────────────────────────┘   │
│                                     │
│  ┌─────────────────────────────┐   │
│  │ GZIP Decompression          │   │  ← zlib-ng (SIMD)
│  │ (On-demand, cached)         │   │
│  └─────────────────────────────┘   │
└─────────────────────────────────────┘
```

### Core Design Principles

#### 1. **Zero-Copy Data Access**
- **Arena Allocation**: Entire file content in single `std::vector<char>` buffer
- **String Views**: All sequences returned as `std::string_view` (pointer + length)
- **No Heap Fragmentation**: Single allocation eliminates memory overhead

```cpp
// Zero-copy access pattern
std::string_view seq = cache.getView("read_001");
// seq points directly into text_arena_ - no allocation!
```

#### 2. **Lock-Free Concurrent Reads** (ADR-001)
```cpp
std::atomic<bool> data_ready_; // C++20 acquire/release semantics

// Reader threads (no locks!)
if (data_ready_.load(std::memory_order_acquire)) {
    std::string_view seq = getView("read_001"); 
}
```
- **Immutable After Load**: Once parsing completes, data never changes
- **Memory Ordering**: `memory_order_acquire` on reads ensures visibility
- **2x Throughput**: Eliminates 25-35% overhead from mutex contention

#### 3. **Native GZIP Support** (ADR-002)
```cpp
// Automatic detection
cache.loadFile("genome.fastq.gz"); // Detects .gz extension

// Or explicit
cache.loadGzipFile("data.fq.gz");

// Magic bytes fallback (handles mis-named files)
cache.loadFile("genome.fasta"); // Checks 0x1f 0x8b if ext doesn't match
```

**Implementation (v1.1.0):**
- **zlib-ng v2.2.2**: SIMD-optimized inflate via FetchContent, zlib-compat mode
- **Streaming Decompression**: zlib-ng processes file in 1MB chunks
- **Pre-Size + Direct-Write**: Pre-allocate `text_arena_` using `compressed_size × 3` heuristic with OOM guard (25% of available memory), then write decompressed data directly — no intermediate `temp_buffer`
- **Geometric growth**: Falls back to 2× capacity doubling if pre-size estimate is too low
- **shrink_to_fit**: Releases excess capacity after decompression (~183MB steady-state for 100MB file)

> **⚠️ Breaking change (stricter validation):** `loadFile()` / `loadGzipFile()`
> and `restore()` now **reject** corrupt or truncated inputs instead of loading
> them silently. Truncated GZIP streams, trailing garbage after the last GZIP
> member, and implausible binary-cache sizes/counts throw `std::runtime_error`.
> `.gz` / `.traceon` files that previously loaded partial data as complete must
> be regenerated.

#### 4. **ankerl::unordered_dense Hashing** (v1.2.0)
- **Swiss-Table Design**: Better cache locality and ~0.8s faster inserts vs Robin Hood
- **Pre-allocation**: Reserve 125% of estimated capacity to prevent rehashing
- **Pre-reserved Thread-Local Maps**: Multithreaded parsers pre-reserve maps to avoid mid-parse rehashing
- **Hybrid Keys**: `std::string` keys (SSO for short IDs) with `std::string_view` values
- **Replaced**: `robin_hood::unordered_flat_map` → `ankerl::unordered_dense::map`

---

## 📦 Project Structure

```
TracEon/
├── include/               # Public API headers
│   ├── TracEon.h         # Single-include convenience header
│   ├── Cache.h           # High-level interface
│   ├── SmartStrategy.h   # Core engine (advanced users)
│   ├── SimdUtils.h       # SIMD character-search (AVX2/NEON/scalar)
│   ├── RecordTypes.h     # Type definitions
│   └── MapDefs.h         # ankerl::unordered_dense typedefs
├── src/                  # Implementation
│   ├── Cache.cpp
│   ├── SmartStrategy.cpp # Lock-free logic, GZIP, SIMD parsers
│   └── IEncodingStrategy.cpp
├── tests/                # Unit & integration tests (Catch2)
│   ├── TestHelpers.h          # Shared fixtures/helpers
│   ├── CacheTests.cpp
│   ├── MapDefsTests.cpp
│   ├── FastqTests.cpp
│   ├── test_parser_fastq.cpp  # FASTQ parser domain
│   ├── test_parser_fasta.cpp  # FASTA parser domain
│   ├── test_gzip.cpp          # GZIP/parallel-decode domain
│   ├── test_binary_cache.cpp  # `.traceon` v4 binary cache domain
│   ├── test_lifecycle.cpp     # lifecycle/immutability domain
│   └── test_api_misc.cpp      # misc API/SIMD/NGS domain
├── benchmarks/           # Performance validation
│   ├── benchmark_runner.py    # Matrix benchmark (sizes × scenarios)
│   ├── validate_real_data.py  # Real-world datasets (NCBI/SRA)
│   ├── check_regression.py    # CI regression checker
│   └── README.md              # Benchmarking guide
├── examples/             # Usage demonstrations
│   └── simple_usage.cpp
├── docs/                 # Architecture documentation
│   ├── architecture/
│   │   ├── ADR-001-lock-free-reads.md
│   │   └── ADR-002-gzip-integration.md
│   └── performance-profile.md
├── test_data/            # Sample files
│   ├── simple.fasta
│   └── simple.fastq
└── third_party/          # Vendored dependencies
    ├── lz4/              # BSD license (future: binary cache compression)
    # ankerl::unordered_dense and zlib-ng are fetched via CMake FetchContent at build time
    # robin_hood.h was removed in v1.2.0
```

---

## 🛠️ Quick Start

### Prerequisites
- **C++20 Compiler**: GCC 10+, Clang 12+, MSVC 2019+
- **CMake**: 3.20+
- **Python** (optional): For benchmarks

### Installation

#### From Source
```bash
git clone https://github.com/woosflex/TracEon.git
cd TracEon

# IMPORTANT: Release build is MANDATORY for performance
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release
cmake --build build -j

# Run tests to verify
./build/unit_tests
```

**Expected output:**
```
All tests passed (330+ assertions in 52 test cases)
```

#### CMake Integration
```cmake
# In your CMakeLists.txt
add_subdirectory(third_party/TracEon)
target_link_libraries(your_target PRIVATE traceon_core)
```

### Basic Usage

#### Example 1: Simple Lookup
```cpp
#include <TracEon.h>
#include <iostream>

int main() {
    TracEon::Cache cache;
    
    // Load FASTA/FASTQ (plain or GZIP, auto-detected)
    cache.loadFile("genome.fastq.gz");
    
    // Zero-copy access (returns string_view)
    std::string_view seq = cache.getView("read_001");
    std::cout << "Sequence: " << seq << "\n";
    
    return 0;
}
```

#### Example 2: Binary Cache Workflow
```cpp
#include <TracEon.h>

int main() {
    TracEon::Cache cache;
    
    // First run: Parse and save
    cache.loadFile("large_genome.fasta.gz");  // ~0.245s for 100MB (v1.2.0)
    cache.save("genome.bin");                  // One-time cost
    
    // Subsequent runs: Instant restore via mmap
    cache.restore("genome.bin");               // ~0.00s for 100MB
    
    // Same zero-copy access
    auto seq = cache.getView("chr1");
    
    return 0;
}
```

#### Example 3: Concurrent Access
```cpp
#include <TracEon.h>
#include <thread>
#include <vector>

void worker(const TracEon::Cache& cache, const std::vector<std::string>& ids) {
    for (const auto& id : ids) {
        auto seq = cache.getView(id);  // Lock-free!
        // Process sequence...
    }
}

int main() {
    TracEon::Cache cache;
    cache.restore("genome.bin");
    
    // Spawn multiple reader threads (no locks needed)
    std::vector<std::thread> threads;
    for (int i = 0; i < 8; ++i) {
        threads.emplace_back(worker, std::cref(cache), get_id_batch(i));
    }
    
    for (auto& t : threads) t.join();
    
    return 0;
}
```

> **⚠️ Immutability contract (ADR-001):** a loaded cache is **immutable**.
> `set()` is only legal **before** any `loadFile()`/`restore()`, or after
> `clearCache()` — calling it on a loaded cache throws `std::logic_error`.
> Concurrent reads of a *loaded* cache are lock-free and safe; build-phase
> `set()` calls must be single-threaded (they are not meant to race with
> readers).

See `examples/simple_usage.cpp` for a complete working example.

---

## 📊 Benchmarking

### Quick Performance Check
```bash
# Build benchmark driver
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release
cmake --build build -j

# Run matrix benchmark (sizes: 10MB, 100MB, 500MB, 1GB)
cd benchmarks
conda env create -f traceon_benchmark.yml
conda activate traceon_bench
python benchmark_runner.py 2>&1 | tee benchmark_results_$(date +%Y%m%d).txt
```

### Interpreting Results

#### Expected Throughput Ranges (Release Build)
| Dataset Type | File Size | L3 Cache (Hot) | RAM (Warm) | RAM (Cold) |
|--------------|-----------|----------------|------------|------------|
| WGS (short reads) | 10MB | 40-55M OPS/s | 40-50M OPS/s | 35-45M OPS/s |
| WGS (short reads) | 100MB | 20-25M OPS/s | 12-18M OPS/s | 10-15M OPS/s |
| PacBio (long reads) | 10MB | 50-60M OPS/s | 45-55M OPS/s | 40-50M OPS/s |
| PacBio (long reads) | 100MB | 35-45M OPS/s | 25-35M OPS/s | 20-30M OPS/s |
| Reference genome | 100MB | N/A | 15-25M OPS/s | 12-20M OPS/s |

**Regression Thresholds** (See `docs/performance-profile.md`):
- 🟢 Green: Within 10% of expected range
- 🟡 Yellow: 10-20% degradation (investigate)
- 🔴 Red: >20% degradation (block release)

### Continuous Integration
```bash
# Automated regression check (exits non-zero on failure)
python check_regression.py benchmark_results_{date_stamp}.txt
```

---

## 🧪 Testing

### Run All Tests
```bash
cd build
ctest --output-on-failure

# Or with verbose output
./unit_tests -d yes
```

### Run Specific Test Suites
```bash
# GZIP tests only
./unit_tests "[gzip]"

# Lock-free concurrency tests
./unit_tests "[strategy]"

# All cache API tests
./unit_tests "[Cache]"
```

### Test Coverage
- ✅ **Lock-free concurrent reads** (4-thread validation)
- ✅ **Zero-copy memory semantics** (pointer equality checks)
- ✅ **GZIP support** (explicit load, auto-detect extension, magic bytes)
- ✅ **Architecture validation** (GenomeIndex vs NGSIndex selection)
- ✅ **Edge cases** (empty files, malformed records, truncated streams)

### Memory Leak Validation (Optional)
```bash
# Requires valgrind
cd build
valgrind --leak-check=full --show-leak-kinds=all ./unit_tests
```

Expected output: `All heap blocks were freed -- no leaks are possible`

---

## 🌐 Remote Access (v2.3.0)

TracEon can serve a **loaded, immutable cache snapshot** over TCP so many
processes/containers read the same data without each re-loading it. The remote
slice is POSIX-only, dependency-free (plain sockets), and **trusted-network-only**
(no auth — do not expose to the internet).

### Build with the server enabled
```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release -DTRACEON_BUILD_SERVER=ON
cmake --build build -j --target remote_bench traceon_make_cache
```
`TRACEON_BUILD_SERVER` is **ON by default** on non-Windows; the server/client live
in the separate `traceon_remote` library (core stays socket-free).

### Serve a v4 cache
```bash
# 1. Turn a FASTA/FASTQ into a v4 .traceon cache (once)
./build/traceon_make_cache input.fasta cache.traceon

# 2. Serve it (blocking; Ctrl-C to stop)
./build/remote_bench --mode=serve --host=127.0.0.1 --port=9876 --file=cache.traceon
```

### Connect with the benchmark (loopback example)
```bash
# In-process zero-copy baseline
./build/remote_bench --mode=local --file=cache.traceon --count=20000

# Remote client, one connection per thread (8 threads x 5000)
./build/remote_bench --mode=remote --host=127.0.0.1 --port=9876 \
    --file=cache.traceon --threads=8 --count=5000 --warmup=1000
```
Every response is **CRC32C-verified** on the client; missing keys return
NOT_FOUND (never garbage); a corrupt/crafted payload fails loudly.
See `docs/architecture/ADR-006-traceon-remote-access.md` for the full protocol.

### In code
```cpp
#include "TraceonServer.h"
#include "TraceonClient.h"

TracEon::TraceonServer server(cache, 0, "127.0.0.1"); // port 0 = ephemeral
server.start();
uint16_t port = server.port();                       // actual bound port

TracEon::TraceonClient client("127.0.0.1", port);
auto seq = client.getView("seq1");   // std::optional<std::string>
bool has = client.has("seq1");
```

### Docker testbed
```bash
docker/run.sh   # builds images, generates a cache, serves + benchmarks it
```

---

## 🎓 Architecture Documentation

### Architecture Decision Records (ADRs)
- **[ADR-001: Lock-Free Reads](docs/architecture/ADR-001-lock-free-reads.md)**  
  Memory ordering guarantees, C++20 atomics justification
  
- **[ADR-002: GZIP Integration](docs/architecture/ADR-002-gzip-integration.md)**  
  Design decisions, rejected alternatives, performance trade-offs

- **[ADR-003: SIMD Parsing & Hash Map Optimization](docs/architecture/ADR-003-simd-parsing-hash-map.md)**  
  SIMD-accelerated boundary scanning, ankerl::unordered_dense, pre-reserved thread-local maps

- **[ADR-005: v4 Binary Format & CRC32C](docs/architecture/ADR-005-traceon-v4-binary-format.md)**
  Wire layout, checksum coverage/streaming order, hardware dispatch, hardening

- **[ADR-006: Remote Read Access over TCP](docs/architecture/ADR-006-traceon-remote-access.md)**
  Length-prefixed RPC, CRC32C integrity, thread-per-connection, trusted-network scope

### Performance Characteristics
- **[Performance Profile](docs/performance-profile.md)**  
  Expected throughput ranges, regression thresholds, hardware scaling

### Benchmarking Guide
- **[Benchmark README](benchmarks/README.md)**  
  How to interpret results, comparison methodology, reproducibility

---

## 🗺️ Roadmap

### ✅ v1.0.0 "Avalon" (December 2025)
*Codename: Like EMIYA's absolute defense, provides performance isolation*

- Zero-copy arena allocation
- Lock-free concurrent reads (C++20 atomics)
- Robin Hood flat hash maps
- Native GZIP support with auto-detection
- Memory-mapped binary cache

### ✅ v1.1.0 "Bakuya" (Q1 2026) — Completed
*Twin swords - Dual optimization paths*

- **zlib-ng Integration**: Replaced system zlib with zlib-ng v2.2.2 (SIMD-optimized inflate) via CMake FetchContent
- **Optimized GZIP Buffering**: Pre-size with 3x heuristic + OOM guard (25% available memory), direct-write to text_arena_ eliminates temp_buffer
- **Load Time**: 0.251s for 100MB GZIP (10% improvement over v1.0.0)
- **Memory Peak**: 266MB during load (27% reduction from 366MB)
- **Performance Goal**: Partially met — 10% improvement achieved; 20% target requires parallel decompression

### ✅ v1.2.0 "Caladbolg" (Q2 2026) — Completed
*Rainbow sword - Spectrum of optimization*

- **SIMD-Accelerated Parsing**: `simd_find_char()` with AVX2 (32 bytes/iter), NEON (16 bytes/iter), runtime dispatch
- **ankerl::unordered_dense**: Replaced Robin Hood with Swiss-table hash map — 86% load time reduction
- **Pre-reserved Thread-Local Maps**: Eliminated mid-parse rehashing in multithreaded parsers
- **Load Time**: 0.245s for 100MB GZIP (86% reduction from 1.843s v1.1.0 baseline)
- **Memory Peak**: 185MB during load (30% reduction from 263MB)
- **normalizeFastaArena() fix**: Trailing newline bug resolved

### ✅ v1.3.0 "Hrunting" (Q3 2026) — Completed
*Hound of the red plains - Relentless pursuit of speed*

- **Binary Cache Compression**: LZ4 integration for `.traceon` files (~3–5× size reduction)
- **Parallel GZIP Decompression**: Concatenated-stream detection + parallel inflate (target: 4x speedup)
- **Transparent Decompression**: On-the-fly decompression for reference genomes
- **Smart Compression**: `saveBinary()` auto-selects LZ4_HC for large DNA/RNA payloads (> 10 MiB), LZ4_default otherwise — no new dependencies, no format bump

### ✅ v1.4.0 "Caliburn" (Q3 2026) — Completed
*The sword of selection — choosing the right index mode*

- **IndexMode Selection**: `SmartStrategy(IndexMode::NGS)` / `Cache(IndexMode::NGS)` activates the hash-keyed `NGSIndex` for short-read workloads; `IndexMode::GENOME` (default) retains string-keyed `GenomeIndex`
- **`Cache::set()` data persistence**: Manually-added entries via `set()` are now serialized by `save()` and restored by `restore()`
- **Correctness hardening**: Closed test gaps for parallel GZIP, RNA/protein format detection, v1 binary format, `clearCache()` + reload, and zero-copy `Cache::getView()`

### ✅ v2.1.0 "Excalibur" (Q4 2026) — Fuzz Gate & Aligner Validation (Released)
*The golden sword — verification and real-world proof*

- **Fuzzing gate**: 7 libFuzzer targets + 83-file seed corpus + CI workflow (per-push + weekly)
- **kmer C API hardened**: atomic freeze flag (thread-safe) + allocation-free error path (total no-throw)
- **Aligner validation**: `woosflex/minimap2` + `woosflex/Winnowmap` forks with `TRACEON=1` table backend + `tcache v2` mmap'd prebuilt index — PAF byte-identical, recurring load+map −13.6% (minimap2, full chr1) / −80% (Winnowmap, chr1-50Mb)

### ✅ v2.0.0 (Q4 2026) — Integrity & v4 Format (Released)
*Codename: Durandal*

- **`.traceon` v4 binary format**: CRC32C whole-payload integrity (SSE4.2 /
  AArch64 / table dispatch), exact-length + exact-frame-termination
  validation, malformed-cache rejection, failure atomicity
- **No legacy binary readers**: v1/v2/v3 caches must be regenerated
- **Lifecycle contract**: documented reader quiescence + debug-only
  `TRACEON_DEBUG_LIFECYCLE` diagnostic
- **Integrity hardening**: strict FASTQ framing, GZIP truncation rejection,
  immutable-after-load enforcement, string_view keys + wyhash reserve perf
- **C API deferred**: the experimental k-mer C API moved out of the v2.0.0
  stable surface (exception-boundary, iterator-ownership, and `TRKI` format
  work pending — see the design review)

### ✅ v2.2.0 "Gáe Bolg" (Q4 2026) — Aligner Validation Matrix Complete (Current Release)
*The cursed spear — piercing through every target*

- **5-tool validation matrix complete**: minimap2 (−13.6% full chr1), Winnowmap (−80%), mm2-fast (−39%/−33% RSS), minigraph (first khashl adapter), BLEND (7-config identity) — all PAF byte-identical, all forks live under `woosflex/*`
- **Upstream bugs found + documented**: Winnowmap applyWeight crash + uninitialized rep_len; minigraph GCC-16 crash on own test data + gi->g_own leak; mm2-fast safestringlib GCC-16 patch

### 🎯 v2.3.0 "Harpe" (Planned)
*The divine sword that strikes regardless of distance*

- **Distributed Caching**: Shard datasets across nodes
- **NUMA-Aware Architecture**: Optimize for >8 core systems
- **Remote Access**: Network protocol for shared genomic databases
- **Performance Goal**: Linear scaling to 64+ cores

### 🔬 Research Directions (v3.0+)

- **GPU Acceleration**: CUDA kernels for sequence alignment pre-filtering
- **Persistent Memory**: Intel Optane integration for instant cold starts
- **Columnar Storage**: Separate sequence/quality storage for better compression
- **Machine Learning Integration**: Embedding cache for feature extraction pipelines

---

## 🤝 Contributing

We welcome contributions! Whether it's bug reports, feature requests, or pull requests.

### Guidelines

1. **Code Style**: Follow existing patterns (modern C++20, RAII, const-correctness)
2. **Tests Required**: All new features must include unit tests
3. **Performance Validation**: Run regression benchmark before submitting
   ```bash
   python benchmarks/check_regression.py --baseline main
   ```
4. **Documentation**: Update relevant ADRs/docs for architectural changes

### Development Setup
```bash
# Clone the repository
git clone https://github.com/woosflex/TracEon.git

# Create environment for benchmarking and 
conda env create -f ./benchmarks/traceon_benchmark.yml
conda activate traceon_bench

# Run pre-commit checks
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
cd benchmarks
python benchmark_runner.py 2>&1 | tee benchmark_results_$(date +%Y%m%d).txt
python check_regression.py benchmark_results_{date_stamp}.txt
```

### Reporting Issues

Please include:
- System details (CPU, RAM, OS, compiler version)
- Build configuration (`CMAKE_BUILD_TYPE`, GCC/Clang version)
- Minimal reproducible example
- Expected vs actual behavior

---

## 📄 License

**MIT License** - See [LICENSE](LICENSE) file for details.

This allows commercial and non-commercial use, modification, and distribution with proper attribution.

---

## 🙏 Acknowledgments

### Core Dependencies
- **[ankerl::unordered_dense](https://github.com/martinus/unordered_dense)** - Martinus (MIT License) — Swiss-table hash map, v4.4.0
- **[zlib-ng](https://github.com/zlib-ng/zlib-ng)** - zlib-ng team (zlib License) — SIMD-optimized zlib replacement, v2.2.2
- **[LZ4](https://github.com/lz4/lz4)** - Yann Collet (BSD License)
- **[Catch2](https://github.com/catchorg/Catch2)** - Testing framework (BSL-1.0)
- **[robin_hood.h](https://github.com/martinus/robin-hood-hashing)** - Martinus (MIT License, *retired in v1.2.0*)

### Inspiration & Motivation
This project emerged from frustrations during MSc Bioinformatics research:
- Existing tools either fast but memory-hungry (BioPython) or memory-efficient but unpredictable (PyFastX disk caching)
- No good middle ground for interactive genomic data exploration
- Modern C++20 features underutilized in bioinformatics

**Result:** A library that's both **fast AND efficient**, proving you don't need to compromise.

### Special Thanks
- Bioinformatics community for feedback on early prototypes
- Fate/stay night fans who appreciate the Noble Phantasm naming scheme 😄

---

## 📞 Contact & Support

**Maintainer**: Adnan Raza  
**GitHub**: [github.com/woosflex](https://github.com/woosflex)  
**Project Repository**: [github.com/woosflex/TracEon](https://github.com/woosflex/TracEon)

### Getting Help

- 🐛 **Bug Reports**: [GitHub Issues](https://github.com/woosflex/TracEon/issues)
- 💡 **Feature Requests**: [GitHub Discussions](https://github.com/woosflex/TracEon/discussions)
- 📚 **Documentation**: [GitHub Wiki](https://github.com/woosflex/TracEon/wiki)

---

## ⚖️ Known Limitations

### v1.0.0 Constraints

1. **RAM Latency Wall**: 3-4x degradation when dataset exceeds L3 cache (~16-32MB)
   - This is **hardware physics** (100ns RAM vs 10ns cache), not a software bug
   - Mitigation: Use binary cache to reduce load time (doesn't affect lookup time)

2. **Immutable After Load**: Cannot incrementally update cache (by design for lock-free reads)
   - Alternative: Reload entire cache (fast with mmap)

3. **Single-Node Only**: No distributed caching (planned for v2.1.0)
   - Target: Local workstation/server workloads

4. **Test Path Dependency**: Tests must be run from build directory
   - Workaround: `cd build && ./unit_tests`

### Platform Support

| Platform | Status | Notes |
|----------|--------|-------|
| Linux x86_64 | ✅ Fully Supported | Primary development platform |
| macOS ARM64 | ✅ Supported | Tested on Apple Silicon |
| Windows x64 | ⚠️ Experimental | MSVC 2019+ required, mmap uses Win32 API |
| ARM64 Linux | 🔄 Untested | Should work, needs validation |

---

## 🌟 Why TracEon?

### The Problem
Bioinformatics tools force a choice:
- **Fast but memory-hungry**: BioPython loads entire datasets (OOM on 100MB)
- **Memory-efficient but unpredictable**: PyFastX's disk caching has variable latency
- **Scalable but complex**: Distributed databases (Cassandra, etc.) overkill for local analysis

### The Solution: TracEon
- **Zero-copy** eliminates allocation overhead
- **Lock-free** enables true concurrent access
- **Memory-mapped cache** provides instant "load" times
- **Modern C++20** leverages decades of language evolution

### Who Should Use TracEon?

✅ **Perfect For:**
- Interactive genomic data exploration (Jupyter notebooks, R Shiny apps)
- High-throughput pipelines processing thousands of reads/second
- Multi-threaded applications needing concurrent access
- Memory-constrained environments (laptops, cloud instances)

❌ **Not Ideal For:**
- Write-heavy workloads (TracEon is read-optimized)
- Distributed systems (single-node only in v1.0.0)
- Datasets larger than available RAM (future: streaming API in a later release)

---

*"I am the bone of my code,  
C++ is my body, and atomics are my blood.  
I have optimized over a thousand benchmarks,  
Unknown to memory leaks, nor known to locks.  
Have withstood cache misses to create many fast paths,  
Yet these hands will never hold anything but pointers.  
So as I trace—Unlimited Performance Works."*

— Adnan Raza

*"Trace On" - Projecting legendary performance from genomic data across eons.* ⚔️