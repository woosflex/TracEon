# TracEon: High-Performance Genomic Data Cache

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![C++20](https://img.shields.io/badge/C%2B%2B-20-blue.svg)](https://isocpp.org/)
[![Version](https://img.shields.io/badge/Release-v1.0.0%20%22Avalon%22-blueviolet.svg)]()
[![Performance](https://img.shields.io/badge/Performance-60M%20OPS%2Fs-brightgreen.svg)]()
[![GZIP](https://img.shields.io/badge/GZIP-Native%20Support-orange.svg)]()

**TracEon** (v1.0.0 "Avalon") is a **zero-copy, lock-free** genomic data caching library written in modern C++20. It accelerates bioinformatics pipelines by providing microsecond-latency random access to FASTA/FASTQ datasets, including native `.gz` support.

*"Trace On" - A nod to Fate/stay night, re-contextualized for tracing biological data across eons.*

---

## 🎯 What's New in v1.0.0 "Avalon"

### **Native GZIP Support** ✨
- Automatic detection via file extension (`.gz`) or magic bytes (`0x1f 0x8b`)
- Streaming decompression preserves zero-copy architecture
- No performance degradation on queries: GZIP lookups as fast as plain text

### **Lock-Free Concurrent Access** 🔓
- C++20 atomics with acquire/release semantics
- 2x throughput improvement on multi-threaded workloads
- Zero mutex contention on read-heavy operations

### **Production-Ready Performance** ⚡
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
| **Memory Usage** | ~25 MB | ~180 MB | 5x efficient vs BioPython |
| **Load Time** | <0.05s | ~0.15s | Includes parsing |

#### Long Reads (PacBio/Nanopore)
| Metric | 10MB Dataset | 100MB Dataset | Speedup vs PyFastX |
|--------|--------------|---------------|-------------------|
| **Random Lookups** | 50-60M OPS/s | 25-35M OPS/s | **26-54x faster** |
| **GZIP Lookups** | 45-55M OPS/s | 20-30M OPS/s | Maintains advantage |
| **Memory Usage** | ~20 MB | ~120 MB | Highly efficient |

**Key Insight:** Performance scales with dataset size relative to CPU cache hierarchy. 3-4x degradation when exceeding L3 cache (~16-32MB) is **hardware physics** (100ns RAM vs 10ns cache), not a software limitation.

### Binary Cache Performance
| Operation | 100MB Dataset | 1GB Dataset |
|-----------|---------------|-------------|
| **Save Time** | 0.08s | 0.25s |
| **Restore Time** | **~0.00s** | **~0.00s** |
| **File Size** | ~105 MB | ~1.05 GB |

**Note:** Restore is instant via memory mapping (`mmap`/`MapViewOfFile`)

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
│  ┌─────────────────────────────┐   │
│  │ Robin Hood HashMap          │   │  ← Lock-free reads
│  │ [string → SequenceView]     │   │     (C++20 atomics)
│  └─────────────────────────────┘   │
│                                     │
│  ┌─────────────────────────────┐   │
│  │ GZIP Decompression          │   │  ← zlib streaming
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

**Implementation:**
- **Streaming Decompression**: zlib processes file in 1MB chunks
- **Efficient Buffering**: Geometric growth minimizes reallocations
- **Zero-Copy Move**: Decompressed data moved into arena (no double-copy)

#### 4. **Robin Hood Hashing**
- **Open Addressing**: ~15% better cache locality vs chaining
- **Pre-allocation**: Reserve 125% of estimated capacity to prevent rehashing
- **Hybrid Keys**: `std::string` keys (SSO for short IDs) with `std::string_view` values

---

## 📦 Project Structure

```
TracEon/
├── include/               # Public API headers
│   ├── TracEon.h         # Single-include convenience header
│   ├── Cache.h           # High-level interface
│   ├── SmartStrategy.h   # Core engine (advanced users)
│   ├── RecordTypes.h     # Type definitions
│   └── MapDefs.h         # Robin Hood typedefs
├── src/                  # Implementation
│   ├── Cache.cpp
│   ├── SmartStrategy.cpp # Lock-free logic, GZIP, parsers
│   └── IEncodingStrategy.cpp
├── tests/                # Unit & integration tests (Catch2)
│   ├── CacheTests.cpp
│   ├── SmartStrategyTests.cpp
│   └── FastqTests.cpp
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
    ├── robin_hood.h      # MIT licensed (martinus)
    ├── zlib/             # zlib license
    └── lz4/              # BSD license (future: binary cache compression)
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
All tests passed (28 assertions in 6 test cases)
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
    cache.loadFile("large_genome.fasta.gz");  // ~0.15s for 100MB
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

## 🎓 Architecture Documentation

### Architecture Decision Records (ADRs)
- **[ADR-001: Lock-Free Reads](docs/architecture/ADR-001-lock-free-reads.md)**  
  Memory ordering guarantees, C++20 atomics justification
  
- **[ADR-002: GZIP Integration](docs/architecture/ADR-002-gzip-integration.md)**  
  Design decisions, rejected alternatives, performance trade-offs

### Performance Characteristics
- **[Performance Profile](docs/performance-profile.md)**  
  Expected throughput ranges, regression thresholds, hardware scaling

### Benchmarking Guide
- **[Benchmark README](benchmarks/README.md)**  
  How to interpret results, comparison methodology, reproducibility

---

## 🗺️ Roadmap

### ✅ Completed: v1.0.0 "Avalon" (December 2025)
*Codename: Like EMIYA's absolute defense, provides performance isolation*

- Zero-copy arena allocation
- Lock-free concurrent reads (C++20 atomics)
- Robin Hood flat hash maps
- Native GZIP support with auto-detection
- Memory-mapped binary cache

### 🎯 Planned Releases

#### v1.1.0 "Bakuya" (Q1 2026)
*Twin swords - Dual optimization paths*

- **SIMD Prefetching**: Vectorized hash probing (target: 20M OPS/s on 100MB RAM datasets)
- **Optimized GZIP Buffering**: Improved chunk management and pre-allocation
- **Adaptive Capacity Estimation**: Better heuristics for different file types
- **Performance Goal**: Close gap with competitors on small datasets

#### v1.2.0 "Caladbolg" (Q2 2026)
*Rainbow sword - Spectrum of compression*

- **Binary Cache Compression**: LZ4 integration for `.traceon` files (3x size reduction)
- **Transparent Decompression**: On-the-fly decompression for reference genomes
- **Smart Compression**: Auto-select algorithm based on data characteristics
- **Performance Goal**: <100MB binary cache for 1GB reference genome

#### v2.0.0 "Durandal" (Q3 2026)
*Peerless - Language-agnostic access*

- **C API**: Foreign Function Interface for Python/R/Julia
- **Python Bindings**: Zero-copy access from NumPy/Pandas
- **R Bindings**: Integration with Bioconductor workflows
- **Streaming API**: Process datasets larger than RAM

#### v2.1.0 "Excalibur" (Q4 2026)
*Holy sword - Distributed power*

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
- **[Robin Hood Hashing](https://github.com/martinus/robin-hood-hashing)** - Martinus (MIT License)
- **[zlib](https://zlib.net/)** - Jean-loup Gailly & Mark Adler (zlib License)
- **[LZ4](https://github.com/lz4/lz4)** - Yann Collet (BSD License)
- **[Catch2](https://github.com/catchorg/Catch2)** - Testing framework (BSL-1.0)

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
- Datasets larger than available RAM (future: streaming API in v2.0.0)

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