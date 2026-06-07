# TracEon Performance Profile

**Version:** 1.1.0 "Bakuya"  
**Last Updated:** June 7, 2026  
**Test Platform:** Intel Core Ultra 5 125H (14 cores, 16GB LPDDR5 @ 7467MT/s)  
**OS:** Ubuntu 24.04 LTS  
**Compiler:** GCC 13.3.0

This document defines **expected performance characteristics** for regression testing and provides guidance on interpreting benchmark results.

---

## 📊 Benchmark Targets (Release Build)

### Parsing Performance
| Dataset | Size | Expected Time | Regression Threshold |
|---------|------|---------------|---------------------|
| WGS FASTQ | 10MB | < 0.05s | > 0.10s |
| WGS FASTQ | 100MB | < 0.15s | > 0.30s |
| WGS FASTQ.gz | 100MB | < 0.25s | > 0.40s |
| PacBio FASTQ | 100MB | < 0.08s | > 0.15s |
| Reference FASTA | 100MB | < 0.08s | > 0.15s |

**Note:** GZIP decompression adds ~0.10s overhead for 100MB files (v1.1.0 with zlib-ng, down from ~0.13s in v1.0.0).

### Lookup Throughput (Random Access, 500K ops)

#### Small Datasets (L3 Cache Resident)
| Dataset | Expected OPS/s | Warning Threshold | Fail Threshold |
|---------|----------------|-------------------|----------------|
| WGS 10MB | 40-55M | < 35M | < 20M |
| PacBio 10MB | 50-60M | < 45M | < 30M |
| RefGenome 10MB | 40-50M | < 35M | < 25M |

#### Medium Datasets (RAM Latency)
| Dataset | Expected OPS/s | Warning Threshold | Fail Threshold |
|---------|----------------|-------------------|----------------|
| WGS 100MB | 12-18M | < 11M | < 8M |
| WGS 100MB.gz | 11-17M | < 10M | < 7M |
| PacBio 100MB | 25-35M | < 20M | < 15M |
| RefGenome 100MB | 15-25M | < 12M | < 8M |

#### Large Datasets (Memory Intensive)
| Dataset | Expected OPS/s | Warning Threshold | Fail Threshold |
|---------|----------------|-------------------|----------------|
| WGS 500MB | 10-15M | < 9M | < 6M |
| PacBio 500MB | 20-30M | < 18M | < 12M |
| RefGenome 1GB | 12-20M | < 10M | < 7M |

### Memory Efficiency
| Dataset | Expected RSS | Warning Threshold | Fail Threshold |
|---------|--------------|-------------------|----------------|
| WGS 10MB | < 30 MB | > 40 MB | > 60 MB |
| WGS 100MB | < 200 MB | > 250 MB | > 300 MB |
| WGS 500MB | < 900 MB | > 1.1 GB | > 1.3 GB |
| PacBio 100MB | < 130 MB | > 160 MB | > 180 MB |
| WGS 100MB.gz (load peak) | < 280 MB | > 320 MB | > 370 MB |

---

## 🔬 Known Performance Characteristics

### 1. The L3 Cache Boundary (~16-32MB)

**Observation:** Throughput drops 3-4x when dataset exceeds L3 cache size.

| Cache State | Dataset Size | Typical Throughput | RAM Latency |
|-------------|--------------|-------------------|-------------|
| **Hot (L3 Cache)** | < 16MB | 40-60M OPS/s | ~10ns per access |
| **Warm (RAM)** | 16-100MB | 12-25M OPS/s | ~100ns per access |
| **Cold (RAM)** | > 100MB | 10-20M OPS/s | ~100ns + TLB misses |

**This is not a bug.** It's fundamental RAM physics:
- L3 cache: ~10ns latency
- Main RAM: ~100ns latency (10x slower)
- TLB misses: Additional ~50-100ns penalty

**Visual Example:**
```
10MB dataset (fits in L3):
  CPU → L3 Cache (10ns) → Data ✅ FAST (50M OPS/s)

100MB dataset (RAM access):
  CPU → L3 Miss → RAM (100ns) → Data ⚠️ SLOWER (15M OPS/s)
```

**Mitigation:** Binary cache reduces load time, but **cannot** eliminate RAM latency during lookups.

---

### 2. Load Factor Sensitivity (Robin Hood Hashing)

Robin Hood hash maps degrade rapidly above 80% load factor.

| Load Factor | Probe Length | Lookup Performance |
|-------------|--------------|-------------------|
| < 50% | 1-2 probes | Optimal (100% baseline) |
| 50-70% | 2-3 probes | Good (90-100% baseline) |
| 70-80% | 3-5 probes | Acceptable (80-90% baseline) |
| > 80% | 5-10+ probes | Poor (< 70% baseline) |

**Our Strategy:** Parsers pre-allocate **125% of estimated capacity** to maintain < 70% load factor.

**Capacity Estimation Heuristics:**
```cpp
// FASTA: ~100 bytes per record (ID + sequence length)
size_t est_records = file_size / 100;

// FASTQ: ~150 bytes per record (ID + sequence + quality)
size_t est_records = file_size / 150;

// Reserve with safety margin
map.reserve(est_records * 1.25);
```

---

### 3. Multithreading Scaling

Parsing scales linearly up to 8 cores, then plateaus.

| Core Count | Speedup | Bottleneck |
|------------|---------|------------|
| 1 core | 1.0x (baseline) | Single-threaded |
| 2 cores | 1.9x | Good scaling |
| 4 cores | 3.7x | Good scaling |
| 8 cores | 7.2x | Excellent scaling |
| 16 cores | 8.5x | **Memory bandwidth** |
| 32 cores | 9.0x | **Diminishing returns** |

**Limiting Factors:**
1. **Memory Bandwidth Saturation**: DDR5 @ 7467MT/s ≈ 60GB/s shared across all cores
2. **Arena Allocation**: Single-threaded `std::vector` growth
3. **Merge Phase**: Thread-local caches merged sequentially

**Recommendation:** Don't expect gains beyond 8-core systems without NUMA-aware architecture (planned for v2.1.0).

---

### 4. GZIP Decompression Overhead

| Operation | Plain Text | GZIP (v1.0.0) | GZIP (v1.1.0) | Improvement |
|-----------|-----------|---------------|---------------|-------------|
| Load 100MB | 0.15s | 0.28s | **0.251s** | -10% |
| Load peak memory | 180 MB | 366 MB | **266 MB** | -27% |
| Lookup 500K ops | 15M OPS/s | 14M OPS/s | 14M OPS/s | 0% |
| Steady-state memory | 180 MB | 183 MB | 183 MB | 0% |

**Key Insights:**
- **Load Time (v1.1.0)**: 0.251s for 100MB GZIP (10% improvement over v1.0.0 via zlib-ng SIMD inflate + pre-sizing)
- **Load Memory (v1.1.0)**: 266MB peak (27% reduction via direct-write to `text_arena_`, eliminating `temp_buffer`)
- **Lookup Performance**: Unchanged — once decompressed, GZIP and plain text data are identical in memory
- **Steady-State Memory**: Identical across versions (`shrink_to_fit` releases excess capacity to ~183MB)

**Improvement Breakdown (v1.0.0 → v1.1.0):**
```
Total: 0.28s → 0.251s (-29ms, 10%)
  ├─ zlib-ng SIMD inflate:    -25ms  (CPU-optimized CRC + inflate)
  ├─ Pre-size (no realloc):    -3ms  (eliminated 7 reallocations)
  └─ Direct-write (no move):   -1ms  (eliminated std::move)
```

**Why Lookup Speed Is Nearly Identical:**
Once decompressed, GZIP and plain text data occupy the same `text_arena_` buffer. All subsequent accesses are zero-copy via `std::string_view`.

---

### 5. String Key Storage (SSO Optimization)

**Small String Optimization (SSO)** keeps short keys (< 15 chars) in the map key itself, avoiding heap allocations.

| Key Length | Storage | Cache Performance |
|------------|---------|-------------------|
| < 15 chars | Stack (SSO) | Excellent (hot cache) |
| 15-23 chars | Heap (single alloc) | Good (cold cache line) |
| > 23 chars | Heap (multiple allocs) | Degraded (multiple cache misses) |

**Typical Genomic IDs:**
- **Short:** `read_001` (8 chars) → SSO ✅
- **Medium:** `chr1:1000000-2000000` (21 chars) → Single heap allocation ⚠️
- **Long:** `@HWI:1:FLOWCELL:1:1101:1234:2345 1:N:0:ACGT` (50+ chars) → Multiple allocations ❌

**Recommendation:** Keep sequence IDs short (<15 chars) when possible for optimal performance.

---

## 🎯 Regression Test Integration

### Automated Regression Checking

To prevent performance regressions in CI/CD:

```bash
# 1. Run benchmark and capture results
./build/traceon_driver lookup test.bin 1000000 prefix 100000 > result.txt

# 2. Check against baseline
python benchmarks/check_regression.py result.txt --baseline v1.1.0
```

**Exit Codes:**
- `0`: Performance within acceptable range (< 10% degradation)
- `1`: Warning threshold exceeded (10-20% degradation) - investigate
- `2`: Failure threshold exceeded (> 20% degradation) - block merge

### Example CI Configuration

**GitHub Actions:**
```yaml
- name: Performance Regression Check
  run: |
    cmake -B build -DCMAKE_BUILD_TYPE=Release
    cmake --build build --target traceon_driver
    cd benchmarks
    python benchmark_runner.py --quick
    python check_regression.py --baseline ${{ github.base_ref }}
```

**GitLab CI:**
```yaml
performance:
  script:
    - cmake -B build -DCMAKE_BUILD_TYPE=Release
    - cmake --build build -j
    - cd benchmarks && python check_regression.py
  rules:
    - if: '$CI_MERGE_REQUEST_ID'
```

---

## 🐛 Troubleshooting Common Performance Issues

### Issue 1: Getting 6M OPS/s instead of 40M+

**Symptoms:**
- Benchmark shows 5-10M OPS/s on 10MB dataset
- Expected: 40-55M OPS/s

**Causes:**
1. **Debug Build**: Check `CMAKE_BUILD_TYPE`
   ```bash
   grep CMAKE_BUILD_TYPE build/CMakeCache.txt
   # Expected: CMAKE_BUILD_TYPE:STRING=Release
   ```

2. **Missing Robin Hood**: Check CMake output
   ```bash
   grep "robin_hood.h" build/CMakeFiles/*.log
   # Expected: Found robin_hood.h at third_party/
   ```

3. **CPU Throttling**: Check CPU governor
   ```bash
   cat /sys/devices/system/cpu/cpu*/cpufreq/scaling_governor
   # Expected: performance (not powersave)
   ```

**Fix:**
```bash
# 1. Clean rebuild in Release mode
rm -rf build
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j

# 2. Set performance governor (Linux)
sudo cpupower frequency-set -g performance

# 3. Re-run benchmark
./build/unit_tests
```

---

### Issue 2: Inconsistent Performance (High Variance)

**Symptoms:**
- Performance varies 2-3x between runs
- Example: 20M OPS/s, then 8M OPS/s, then 25M OPS/s

**Causes:**
1. **Background Processes**: Other apps consuming RAM/CPU
2. **Thermal Throttling**: CPU overheating
3. **TLB Thrashing**: OS swapping TLB entries

**Fix:**
```bash
# 1. Close background apps
# 2. Monitor CPU temperature
sensors | grep Core

# 3. Run multiple iterations and take median
for i in {1..5}; do
  ./build/unit_tests | grep "OPS/s"
done | sort -n | sed -n '3p'  # Median of 5 runs
```

---

### Issue 3: Out of Memory (OOM) on Large Files

**Symptoms:**
- Process killed with `exit code 137` (SIGKILL)
- `dmesg` shows "Out of memory: Killed process"

**Causes:**
1. **Insufficient RAM**: Dataset exceeds available memory
2. **Memory Fragmentation**: Arena allocation fails
3. **Swap Disabled**: No virtual memory fallback

**Fix:**
```bash
# 1. Check available memory
free -h

# 2. Enable swap (temporary)
sudo dd if=/dev/zero of=/swapfile bs=1G count=8
sudo mkswap /swapfile
sudo swapon /swapfile

# 3. Use binary cache to reduce memory footprint
./build/traceon_driver save large.fasta large.bin
# Future loads use mmap (doesn't count toward RSS)
```

---

### Issue 4: Slow GZIP Decompression

**Symptoms:**
- GZIP load takes > 0.40s for 100MB file
- Expected: 0.20-0.25s (v1.1.0)

**Causes:**
1. **Disk I/O Bottleneck**: Slow storage (HDD, network drive)
2. **zlib-ng Build Issue**: FetchContent may have fallen back to system zlib
3. **File Corruption**: Damaged GZIP archive

**Diagnosis:**
```bash
# 1. Test raw read speed
time cat large.fastq.gz > /dev/null
# Should be < 0.05s for 100MB on SSD

# 2. Verify GZIP integrity
gzip -t large.fastq.gz

# 3. Check zlib-ng version (linkage)
ldd build/traceon_driver | grep zlib
# Expected: references zlib-ng (not system libz.so)

# 4. Verify build log for zlib-ng FetchContent
grep "zlib-ng" build/CMakeCache.txt
# Expected: zlib-ng source directory populated
```

---

## 📈 Performance Optimization Checklist

Before reporting a performance issue, verify:

- [ ] Build mode is `Release` (`-O3` optimization)
- [ ] Robin Hood hashing is enabled (check CMake output)
- [ ] CPU governor set to `performance` (not `powersave`)
- [ ] No background processes consuming RAM/CPU
- [ ] Dataset size appropriate for available RAM
- [ ] File storage is fast (SSD, not HDD or network)
- [ ] Multiple benchmark runs show consistent results (< 20% variance)
- [ ] System temperature within normal range (< 80°C)

---

## 🔍 Profiling Tools

### Recommended Tools for Performance Analysis

**Linux:**
```bash
# 1. perf (CPU profiling)
perf record -g ./build/unit_tests
perf report

# 2. Valgrind Callgrind (function-level analysis)
valgrind --tool=callgrind ./build/unit_tests
kcachegrind callgrind.out.*

# 3. Massif (memory profiling)
valgrind --tool=massif ./build/unit_tests
ms_print massif.out.*
```

**macOS:**
```bash
# Instruments (Xcode)
instruments -t "Time Profiler" ./build/unit_tests
```

**Windows:**
```powershell
# Visual Studio Profiler
VSPerfCmd /start:sample /output:profile.vsp
./build/Release/unit_tests.exe
VSPerfCmd /shutdown
```

---

## 📚 References

### Performance Analysis Papers
- Preshing, *Acquire and Release Semantics* (2012)
- Boehm & Adve, *Foundations of the C++ Concurrency Memory Model* (2008)
- [What Every Programmer Should Know About Memory](https://people.freebsd.org/~lstewart/articles/cpumemory.pdf) (Drepper, 2007)

### Profiling Guides
- [Linux perf Examples](https://www.brendangregg.com/perf.html)
- [Valgrind User Manual](https://valgrind.org/docs/manual/)
- [Intel VTune Profiler](https://www.intel.com/content/www/us/en/developer/tools/oneapi/vtune-profiler.html)

---

**Document Version:** 1.1.0  
**Last Validated:** June 7, 2026  
**Next Review:** July 2026 (post-v1.1.0 release)

*"Trace On" - Profiling legendary performance from genomic data.* 📊