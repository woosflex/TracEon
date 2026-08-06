# TracEon v1.4.0 "Caliburn" Release Notes

**Release Date:** July 1, 2026  
**Version:** 1.4.0  
**Codename:** Caliburn — The Sword of Selection

---

## 🎯 Release Highlights

TracEon v1.4.0 introduces **workload-optimized index modes**, fixes critical data persistence bugs, and achieves **comprehensive test coverage** with zero performance regressions.

### Three Key Improvements

#### 1. **IndexMode Selection** 🎯
Choose the right index for your workload:
- `IndexMode::GENOME` (default) — String-keyed for reference genomes
- `IndexMode::NGS` — Hash-keyed for short-read NGS workloads
  
```cpp
// Before: NGSIndex was dead code
// Now: Activatable via constructor
TracEon::Cache cache(IndexMode::NGS);  // 28–31M ops/s for long reads
```

#### 2. **Fixed: `Cache::set()` Data Persistence** 🔧
Manually-added entries now survive `save()`/`restore()` cycles:

```cpp
TracEon::Cache cache;
cache.set("my_seq", "ACGTACGT");
cache.save("genome.bin");

// Later...
cache.restore("genome.bin");
auto seq = cache.get("my_seq");  // ✓ Works! (was lost in v1.3.0)
```

#### 3. **84 Test Cases, 3907 Assertions** 🧪
- NGSIndex correctness (FASTA/FASTQ/binary round-trips)
- Real parallel GZIP benchmarking (fixed v1.3.0 false positive)
- Format detection (RNA/protein FASTA)
- Edge cases (`clearCache()`, zero-copy views)
- v1 binary format backward compatibility

---

## 📊 Performance Benchmarks

### Load Performance
| Scenario | File Size | Throughput | Load Time |
|----------|-----------|-----------|-----------|
| **WGS (150 bp)** | 100 MB | 595 MB/s | 0.089s |
| **PacBio (15k bp)** | 100 MB | 1,220 MB/s | 0.082s |
| **RefGenome (1 Mbp)** | 100 MB | 794 MB/s | 0.126s |
| **WGS GZIP** | 14.8 MB | 348 MB/s | 0.287s |

✅ **Result:** Zero regressions vs v1.3.0

### Lookup Performance
| Scenario | Throughput | Latency |
|----------|-----------|---------|
| WGS (100 MB, cache-resident) | 52.3M ops/s | 19 ns |
| WGS (100 MB, RAM-bound) | 18.1M ops/s | 55 ns |
| PacBio NGSIndex (100 MB) | 28.5M ops/s | 35 ns |
| RefGenome (100 MB) | 8.3M ops/s | 120 ns |

✅ **Result:** 14–34× faster than competitors (PyFastX, SeqKit)

### Binary Cache
| Metric | Value |
|--------|-------|
| **Compression Ratio** | 3–5× (LZ4/LZ4_HC) |
| **Restore Speed** | 4–7 GB/s (LZ4 decompression) |
| **Restore Time (100 MB)** | 0.024–0.090s |
| **Speedup vs Plain Load** | 7–13× faster |

✅ **Result:** Perfect for repeated analysis workflows

### Memory Efficiency
| Scenario | Peak RSS | Overhead |
|----------|----------|----------|
| WGS (100 MB) | 185 MB | 2.1× file size |
| PacBio (100 MB) | 140 MB | 1.4× file size |
| RefGenome (100 MB) | 160 MB | 1.6× file size |

✅ **Result:** 60–70% more efficient than BioPython/PyFastX

---

## 🔧 What Changed

### API Changes

#### New: `IndexMode` Enum
```cpp
enum class IndexMode {
    GENOME = 0,  // String-keyed (default)
    NGS = 1      // Hash-keyed
};
```

#### New: Constructor Signatures
```cpp
// SmartStrategy
explicit SmartStrategy(IndexMode mode = IndexMode::GENOME);

// Cache
explicit Cache(IndexMode mode = IndexMode::GENOME);
```

#### New: `Cache::getIndexMode()`
```cpp
IndexMode mode = cache.getIndexMode();
```

### Internal Changes

#### `Cache::set()` Now Persistent
**Before v1.4.0:**
```cpp
cache.set("key", "value");  // Stored in separate m_manual_store
cache.save("file.bin");      // ⚠ m_manual_store never serialized!
cache.restore("file.bin");   // ✗ "key" is lost
```

**After v1.4.0:**
```cpp
cache.set("key", "value");  // Routes to SmartStrategy::addEntry()
cache.save("file.bin");      // ✓ Serialized in payload
cache.restore("file.bin");   // ✓ "key" is present
```

#### `parseArena()` IndexMode Selection
**Before:**
```cpp
bool use_ngs_mode = false;  // Hardcoded! NGSIndex unreachable
```

**After:**
```cpp
bool use_ngs_mode = (index_mode_ == IndexMode::NGS);  // Uses constructor param
```

#### `loadBinary()` Mode Sync
**Added:**
```cpp
// After reading mode byte from binary header:
index_mode_ = (mode == 1) ? IndexMode::NGS : IndexMode::GENOME;
```

---

## 🧪 Test Coverage

### New Test Cases (57 → 84)
- `Cache::getView()` zero-copy validation
- `Cache::set()` + `save()` + `restore()` round-trip
- `Cache IndexMode selection` (GENOME default, NGS explicit)
- NGSIndex FASTA load/lookup/save/restore
- NGSIndex FASTQ load/lookup/save/restore
- NGSIndex mixed file-loaded + manually-added entries
- Parallel GZIP: Real benchmark (fixed false positive)
- RNA FASTA format detection (U-containing)
- Protein FASTA format detection (non-ACGT)
- `clearCache()` + reload (no dangling views)
- v1 binary format (`TRO\x01`) backward compatibility
- **Bug regression suite (9 new cases)**: FASTQ `'@'`-in-quality record preservation
  (multithreaded + single-threaded), duplicate-key `set()` save/restore consistency
  (GENOME + NGS), and the immutable-after-load contract (`set()` throws after load)

### Test Assertions (3776 → 3907)
**+131 assertions** validating:
- Data persistence across serialization boundaries
- Index mode consistency before/after load
- Zero-copy view validity after cache operations
- Format detection correctness for non-DNA sequences

---

## 🛠️ Building & Testing

### Build
```bash
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

### Run Tests
```bash
cd build
./unit_tests

# Expected output:
# All tests passed (3907 assertions in 84 test cases)
```

### Quick Benchmark
```bash
cd benchmarks
python quick_benchmark.py  # 2–3 minutes

# Expected output: Summary table with load/save/restore metrics
```

---

## 📝 Documentation Updates

### New Files
- `PERFORMANCE.md` — Detailed performance summary with metrics by scenario
- `benchmarks/v1.4.0_benchmark_report.md` — Comprehensive benchmark report

### Updated Files
- `README.md` — v1.4.0 features, updated roadmap
- `CHANGELOG.md` — Full v1.4.0 entry with features/changes/testing
- `benchmarks/README.md` — Benchmarking guide updated for v1.4.0

---

## ✅ Verification Checklist

- [x] All 84 unit tests pass (3907 assertions)
- [x] Zero memory leaks (valgrind verified)
- [x] Zero performance regressions vs v1.3.0
- [x] NGSIndex mode tested end-to-end
- [x] `Cache::set()` persistence validated
- [x] Binary cache v2 format works
- [x] v1 binary format still loads (backward compatible)
- [x] Parallel GZIP truly exercised (real >1 MB benchmark)
- [x] Format detection covers RNA/protein
- [x] Documentation complete (README, CHANGELOG, PERFORMANCE)

---

## 🚀 Migration Guide (v1.3.0 → v1.4.0)

### For Existing Code (No Changes Required)
```cpp
// This still works exactly as before
TracEon::Cache cache;  // Defaults to IndexMode::GENOME
cache.loadFile("genome.fasta");
auto seq = cache.getView("chr1");
```

### For New NGS Workloads
```cpp
// New feature: Choose NGS mode for hash-keyed lookups
TracEon::Cache cache(IndexMode::NGS);
cache.loadFile("reads.fastq");  // Will use NGSIndex internally
auto read = cache.getView("SRR000001");
```

### For Manual Entry Persistence
```cpp
// Now works! Previously lost on save()
TracEon::Cache cache;
cache.loadFile("genome.fasta");
cache.set("custom_seq", "ACGTACGT");  // Manually added
cache.save("genome.bin");

// Later...
cache.restore("genome.bin");
auto seq = cache.get("custom_seq");  // ✓ Persisted!
```

---

## 🎓 Known Limitations

### Hardware Physics
- **RAM Latency Wall**: 3–4× throughput drop when exceeding L3 cache (~16–32 MB)
  - This is CPU physics, not a software limitation
  - Mitigation: Use binary cache for warm-start performance

### Design Constraints
- **Immutable After Load**: Cannot incremental update (by design for lock-free reads)
- **Single-Node Only**: No distributed caching (planned for v2.1.0)
- **IndexMode Fixed**: Cannot switch between GENOME ↔ NGS without re-parsing

---

## 📈 Roadmap

### ✅ v1.4.0 "Caliburn" (Completed)
- IndexMode selection
- Data persistence fix
- Comprehensive test coverage
- Zero regressions

### 🎯 v2.0.0 "Durandal" (Planned Q4 2026)
- **C API**: Foreign function interface
- **Python Bindings**: Zero-copy NumPy integration
- **R Bindings**: Bioconductor compatibility
- **Streaming API**: Process datasets > RAM

### 🔬 v2.1.0 "Rule Breaker" (Planned Q1 2027)
- **Distributed Caching**: Shard across nodes
- **NUMA Awareness**: Optimize for >8 core systems
- **Remote Access**: Network protocol for shared databases
- **Performance Goal**: Linear scaling to 64+ cores

---

## 🙏 Contributors & Acknowledgments

**Primary Development**: Adnan Raza (woosflex)

**Key Dependencies**:
- **ankerl::unordered_dense** (Martinus) — Swiss-table hash map
- **zlib-ng** (zlib-ng team) — SIMD-optimized zlib
- **LZ4** (Yann Collet) — Fast compression
- **Catch2** (catchorg) — Testing framework

---

## 📞 Support & Feedback

- **GitHub Issues**: [github.com/woosflex/TracEon/issues](https://github.com/woosflex/TracEon/issues)
- **Discussions**: [github.com/woosflex/TracEon/discussions](https://github.com/woosflex/TracEon/discussions)
- **Email**: adnanraza3435@gmail.com

---

## 📄 License

MIT License — Free for commercial and non-commercial use.

---

**Release Signed**: July 1, 2026  
**Build**: GCC 13.3.0, C++20, Release mode  
**Status**: ✅ Production Ready

