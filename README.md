# TracEon v1.0.0 "Avalon"
## High-Performance Genomic Data Cache

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![C++20](https://img.shields.io/badge/C%2B%2B-20-blue.svg)](https://isocpp.org/)
[![Version](https://img.shields.io/badge/version-1.0.0--Avalon-blue.svg)]()
[![Performance](https://img.shields.io/badge/Performance-57M%20OPS%2Fs-brightgreen.svg)]()

TracEon is a **zero-copy, lock-free** genomic data caching library written in modern C++20. 

> *"I am the bone of my sword" - Just as EMIYA projects legendary weapons, TracEon projects ultra-fast data access from genomic files.*

---

## 🚀 Performance Highlights

### Version 1.0.0 "Avalon" (December 2025)

---

## 🚀 Performance Highlights

###  Avalon v1.0.0 (December 2025)

| Workload | TracEon 1.0.0 | PyFastX | BioPython | Speedup |
|----------|---------------|---------|-----------|---------|
| **WGS 100MB Lookup** | **13.2M OPS/s** | 37M OPS/s† | OOM | Competitive‡ |
| **PacBio 100MB Lookup** | **57.3M OPS/s** | 1.4M OPS/s | OOM | **40x faster** |
| **Memory (100MB file)** | **180 MB** | Unknown | > 900 MB | **5x more efficient** |

† *PyFastX uses disk caching; comparison is apples-to-oranges*  
‡ *TracEon's memory-resident architecture optimizes for repeated access*

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
│  - loadFile()                       │
│  - getView() ← Zero-copy            │
└──────┬──────────────────────────────┘
       │
       v
┌─────────────────────────────────────┐
│  SmartStrategy (Engine)             │
│  ┌─────────────────────────────┐   │
│  │ Arena Allocator             │   │  ← Single malloc
│  │ ┌─────────────────────────┐ │   │
│  │ │ ACGTACGT...             │ │   │
│  │ └─────────────────────────┘ │   │
│  └─────────────────────────────┘   │
│                                     │
│  ┌─────────────────────────────┐   │
│  │ Robin Hood HashMap          │   │  ← Lock-free reads
│  │ [string → string_view]      │   │
│  └─────────────────────────────┘   │
└─────────────────────────────────────┘
```

### Key Innovations

1. **Arena Allocation**: Entire file loaded into single `std::vector<char>`, eliminating heap fragmentation
2. **String Views**: All sequences accessed via `std::string_view` (pointer + length), zero allocations
3. **Lock-Free Reads**: C++20 atomics with acquire/release semantics, 2x concurrent throughput
4. **Robin Hood Hashing**: Open-addressing hash map with ~15% better cache locality vs chaining
5. **Memory Mapping**: Binary cache files mapped directly into virtual memory (instant "load")

---

## 📦 Project Structure
```
TracEon/
├── include/               # Public API headers
│   ├── TracEon.h         # Single-include header
│   ├── Cache.h           # High-level interface
│   └── SmartStrategy.h   # Core engine (advanced users)
├── src/                  # Implementation
│   ├── Cache.cpp
│   └── SmartStrategy.cpp # Lock-free logic, parsers
├── tests/                # Unit & integration tests
│   ├── CacheTests.cpp
│   └── SmartStrategyTests.cpp
├── benchmarks/           # Performance testing
│   ├── benchmark_runner.py    # Matrix benchmark
│   ├── validate_real_data.py  # Real-world datasets
│   ├── check_regression.py    # CI regression checker
│   └── README.md              # Benchmark guide
├── examples/             # Usage demonstrations
│   └── simple_usage.cpp
├── docs/                 # Architecture documentation
│   ├── architecture/
│   │   └── ADR-001-lock-free-reads.md
│   └── performance-profile.md
└── third_party/          # Vendored dependencies
    ├── robin_hood.h      # MIT licensed
    └── zlib/             # zlib license
```

---

## 🛠️ Quick Start

### Prerequisites
- C++20 compiler (GCC 10+, Clang 12+, MSVC 2019+)
- CMake 3.20+

### Build
```bash
git clone https://github.com/woosflex/TracEon.git
cd TracEon

# Release build (CRITICAL for performance)
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release
cmake --build build -j

# Run example
./build/example_usage
```

### Basic Usage
```cpp
#include <TracEon.h>
#include <iostream>

int main() {
    TracEon::Cache cache;
    
    // Load and parse (multithreaded for files > 10MB)
    cache.loadFile("genome.fastq");
    
    // Save binary cache (one-time cost)
    cache.save("genome.bin");
    
    // Future sessions: instant restore via mmap
    cache.restore("genome.bin");
    
    // Zero-copy access (recommended)
    std::string_view seq = cache.getView("read_001");
    std::cout << "Sequence: " << seq << "\n";
    
    return 0;
}
```

---

## 📊 Benchmarking

### Run Full Matrix Benchmark
```bash
pip install psutil tqdm requests
python benchmarks/benchmark_runner.py
```

See `benchmarks/README.md` for detailed interpretation guide.

### Expected Performance (Release Build)

| Dataset Type | File Size | Throughput |
|--------------|-----------|------------|
| WGS (short reads) | 10MB | 40-55M OPS/s |
| WGS (short reads) | 100MB | 12-18M OPS/s |
| PacBio (long reads) | 100MB | 45-60M OPS/s |
| Reference genome | 100MB | 15-25M OPS/s |

**Note:** Throughput degrades 3-4x when dataset exceeds L3 cache size (~16-32MB). This is **hardware physics**, not a software bug.

---

## 🧪 Testing
```bash
# Build and run all tests
cd build
ctest --output-on-failure

# Run specific test suite
./unit_tests "[strategy]"
```

Tests cover:
- Lock-free concurrent reads (4-thread validation)
- Zero-copy memory semantics
- Architecture verification (GenomeIndex vs NGSIndex)
- Edge cases (empty files, malformed records)

---

## 🎓 Architecture Deep Dives

- **[ADR-001: Lock-Free Reads](docs/architecture/ADR-001-lock-free-reads.md)** - Memory ordering guarantees
- **[Performance Profile](docs/performance-profile.md)** - Expected characteristics & regression thresholds

---

## 📈 Roadmap

- [x] **v1.0.0 "Avalon"**: Zero-copy architecture, lock-free reads, Robin Hood hashing
- [ ] **v1.1.0 "Bakuya"**: SIMD prefetching for hash probes
- [ ] **v1.2.0 "Caladbolg"**: Transparent compression for reference genomes
- [ ] **v2.0.0 "Durandal"**: C API for Python/R bindings
- [ ] **v2.1.0 "Excalibur"**: Distributed caching across nodes

---

## 🤝 Contributing

We welcome contributions! Please see:
1. Code follows existing style (C++20, modern idioms)
2. Tests pass: `ctest --output-on-failure`
3. Performance maintained: `python benchmarks/check_regression.py`
4. Documentation updated for architectural changes

---

## 📄 License

MIT License - See [LICENSE](LICENSE) file

---

## 🙏 Acknowledgments

- **Robin Hood Hashing**: Martinus' excellent implementation ([MIT License](https://github.com/martinus/robin-hood-hashing))
- **zlib**: Jean-loup Gailly & Mark Adler
- **Inspiration**: This project bridges my MSc Bioinformatics studies with systems programming passion

---

## 📧 Contact

**Developer/Maintainer**: Adnan Raza  
**Project**: [github.com/woosflex/TracEon](https://github.com/woosflex/TracEon)

*"Trace On" - A nod to Fate/stay night, re-contextualized for tracing biological data across eons.*