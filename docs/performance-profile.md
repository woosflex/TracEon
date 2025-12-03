# TracEon Performance Profile

**Version:** 1.0.0 "Avalon"  
**Last Updated:** December 2025

This document defines **expected performance characteristics** for regression testing.

## Benchmark Targets (Release Build, AMD EPYC 7B12)

### Parsing Performance
| Dataset | Size | Expected Time | Threshold (Fail) |
|---------|------|---------------|------------------|
| WGS FASTQ | 100MB | < 0.15s | > 0.3s |
| PacBio FASTQ | 100MB | < 0.08s | > 0.15s |
| Reference FASTA | 100MB | < 0.08s | > 0.15s |

### Lookup Throughput (Random Access, 500K ops)
| Dataset | Expected OPS/s | Threshold (Fail) |
|---------|----------------|------------------|
| WGS 10MB | > 35M | < 20M |
| WGS 100MB | > 11M | < 8M |
| PacBio 100MB | > 40M | < 25M |

### Memory Efficiency
| Dataset | Expected RSS | Threshold (Fail) |
|---------|--------------|------------------|
| WGS 100MB | < 200 MB | > 300 MB |
| PacBio 100MB | < 120 MB | > 180 MB |

## Known Performance Characteristics

### The L3 Cache Boundary (~16-32MB)
**Observation:** Throughput drops 3-4x when dataset exceeds L3 cache size.
- **10MB datasets:** 40-55M OPS/s (hot cache)
- **100MB datasets:** 12-18M OPS/s (RAM latency)

**This is not a bug.** It's fundamental RAM physics (100ns latency vs 10ns cache).

### Load Factor Sensitivity
Robin Hood hash maps degrade rapidly above 80% load. Our parsers pre-allocate 125% of estimated capacity to maintain < 70% load factor.

### Multithreading Scaling
Parsing scales linearly up to 8 cores, then plateaus due to:
1. Memory bandwidth saturation
2. Arena allocation being single-threaded
3. Merge phase overhead

**Recommendation:** Don't expect gains beyond 8-core systems without NUMA-aware architecture.

## Regression Test Integration

To prevent performance regressions, run:
```bash
./build/traceon_driver lookup test.bin 1000000 prefix 100000 > result.txt
python benchmarks/check_regression.py result.txt
```


Expected: Script exits 0 if throughput within acceptable range, exits 1 if regression detected.