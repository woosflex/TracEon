# TracEon Benchmarking Suite

**Version:** 1.0.0 "Avalon"  
**Purpose:** Performance validation and regression testing

This directory contains tools for benchmarking TracEon against synthetic and real-world genomic datasets, comparing against state-of-the-art tools (BioPython, PyFastX, SeqKit).

---

## 📋 Quick Start

### Prerequisites

```bash
# Python dependencies
conda env create -f traceon_benchmark.yml
```

**System Requirements:**
- Python 3.8+
- At least 4GB RAM for 100MB benchmarks
- At least 16GB RAM for 1GB benchmarks

---

## 🚀 Running Benchmarks

### 1. Matrix Benchmark (Synthetic Data)

Tests TracEon across multiple file sizes and biological scenarios.

```bash
# From project root
cd TracEon

# Build traceon_driver (required)
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release
cmake --build build --target traceon_driver -j

# Run matrix benchmark
cd benchmarks
python benchmark_runner.py
```

**What it tests:**
- **Sizes:** 10MB, 100MB, 500MB, 1GB
- **Scenarios:** WGS, WGS_GZ, PacBio, PacBio_GZ, Nanopore, Nanopore_GZ, RefGenome, RefGenome_GZ
- **Metrics:** Parse time, restore time, lookup throughput, memory usage

**Expected runtime:**
- 10MB tests: ~2 minutes per scenario
- 100MB tests: ~5 minutes per scenario  
- Full matrix: 1-2 hours (32 combinations)

### 2. Real-Data Validation (NCBI Datasets)

Downloads and tests against actual genomic data from public repositories.

```bash
python validate_real_data.py
```

**Datasets tested:**
- *E. coli* K-12 genome (4.6 MB compressed)
- *C. elegans* chromosome I (15 MB compressed)
- Human chromosome 21 (38 MB compressed)

**Expected runtime:** 10-15 minutes (including downloads)

### 3. Quick Smoke Test

For rapid validation during development:

```bash
# Run only 10MB tests (fast)
python benchmark_runner.py --sizes 10

# Run specific scenario
python benchmark_runner.py --scenario WGS --sizes 10,100
```

---

## 📊 Interpreting Results

### Throughput Metrics

**Performance Tiers:**
| Throughput | Classification | Typical Cause |
|------------|----------------|---------------|
| **40M+ OPS/s** | Excellent | L3 cache-resident (< 16MB) |
| **15-30M OPS/s** | Good | RAM latency, long reads |
| **10-15M OPS/s** | Acceptable | RAM latency, short reads |
| **< 10M OPS/s** | Investigate | Possible DEBUG build or issue |

**Context matters:**
- **10MB dataset @ 50M OPS/s**: Normal (fits in L3 cache)
- **100MB dataset @ 15M OPS/s**: Normal (RAM latency wall)
- **100MB dataset @ 6M OPS/s**: Problem (check build mode)

### Memory Efficiency

**Comparison benchmarks:**
```
TracEon vs BioPython (100MB WGS):
  TracEon:   180 MB  (baseline)
  BioPython: 900 MB+ (5x worse)
```

**Interpretation:**
- **< 2x file size**: Excellent (zero-copy working)
- **2-3x file size**: Good (some overhead acceptable)
- **> 3x file size**: Investigate (possible memory leak)

### Load Time Analysis

**Components:**
```
Total Load Time = I/O + Decompression + Parsing + Index Building

For 100MB FASTQ.gz:
  I/O:            ~0.05s (SSD)
  Decompression:  ~0.10s (zlib)
  Parsing:        ~0.05s (multithreaded)
  Index Building: ~0.05s (Robin Hood)
  Total:          ~0.25s
```

**Regression thresholds:**
- **Expected:** 0.20-0.30s for 100MB GZIP
- **Warning:** > 0.35s (investigate)
- **Failure:** > 0.50s (block merge)

---

## 🔍 Benchmark Output Explained

### Example Output
```
  FILE SIZE: 10 MB
SCENARIO     | RECS      | SEQTK  | SEQKIT | BIOPY  | TRACEON | RESTORE | TRACEON OPS  | PYFASTX OPS  | MEM TRA  | MEM PFX  | vs FILE 
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  Generating bench_WGS_10MB.fastq (29,959 recs, ~150bp)...
WGS          | 29,959    | 0.02   | 0.05   | 0.38   | 0.02    | 0.01    | 47,267,637   | 23,611,258  | 22       | 11       | +148%   
  Generating bench_WGS_GZ_10MB.fastq.gz (29,959 recs, ~150bp)...
WGS_GZ       | 29,959    | 0.03   | 0.04   | 0.33   | 0.04    | 0.01    | 48,562,392   | 23,725,049  | 23       | 11       | +1434%  
  Generating bench_PacBio_10MB.fastq (348 recs, ~15000bp)...
PacBio       | 348       | 0.07   | 0.03   | 0.19   | 0.01    | 0.00    | 44,291,758   | 1,108,627    | 1        | 12       | -92%
```

**Column meanings:**
- **RECS:** Number of sequences in generated test file
- **SEQTK:** SeqTK processing time (FASTQ only, baseline)
- **SEQKIT:** SeqKit stats time (all formats, baseline)
- **BIOPY:** BioPython parse time
- **TRACEON:** TracEon load time (parse + index)
- **RESTORE:** TracEon binary cache restore time (mmap)
- **TRACEON OPS:** TracEon Lookups per second of Random Queries
- **PYFASTX OPS:** PyFastX lookups per second
- **MEM TRA:** TracEon peak memory usage (RSS)
- **MEM PFX:** PyFastX peak memory usage (RSS)
- **vs FILE:** Memory utilised versus the file size

**Special markers:**
- **`-`**: Tool doesn't support this format
- **`(BioPy ERR)`**: BioPython crashed (OOM or exception)

---

## 🎯 Performance Expectations (Release Build)

### Reference Values (Intel Core Ultra 5 125H)

#### Small Datasets (L3 Cache Resident)
| Dataset Type | File Size | Expected OPS/s | Memory | Load Time |
|--------------|-----------|----------------|--------|-----------|
| WGS (short reads) | 10MB | 40-55M | ~25 MB | < 0.05s |
| PacBio (long reads) | 10MB | 50-60M | ~20 MB | < 0.03s |
| RefGenome | 10MB | 40-50M | ~22 MB | < 0.03s |

#### Medium Datasets (RAM Latency)
| Dataset Type | File Size | Expected OPS/s | Memory | Load Time |
|--------------|-----------|----------------|--------|-----------|
| WGS (short reads) | 100MB | 12-18M | ~180 MB | 0.12-0.18s |
| WGS_GZ (compressed) | 100MB | 11-17M | ~180 MB | 0.20-0.30s |
| PacBio (long reads) | 100MB | 25-35M | ~120 MB | 0.06-0.10s |
| PacBio_GZ | 100MB | 20-30M | ~120 MB | 0.15-0.25s |
| RefGenome | 100MB | 15-25M | ~150 MB | 0.06-0.10s |

#### Large Datasets (Production Scale)
| Dataset Type | File Size | Expected OPS/s | Memory | Load Time |
|--------------|-----------|----------------|--------|-----------|
| WGS | 500MB | 10-15M | ~900 MB | 0.60-0.90s |
| PacBio | 500MB | 20-30M | ~600 MB | 0.30-0.50s |
| RefGenome | 1GB | 12-20M | ~1.8 GB | 0.60-1.00s |

**Note:** Performance scales with CPU cache hierarchy. Expect 3-4x degradation when dataset exceeds L3 cache (~16-32MB).

---

## 🐛 Troubleshooting

### Issue 1: Getting 6M OPS/s instead of 40M+

**Symptom:** Benchmark shows 5-10M OPS/s on 10MB dataset (expected: 40-55M).

**Diagnosis:**
```bash
# Check build mode
grep CMAKE_BUILD_TYPE build/CMakeCache.txt
# Expected: CMAKE_BUILD_TYPE:STRING=Release

# Check Robin Hood detection
grep "robin_hood" build/CMakeFiles/*.log
# Expected: Found robin_hood.h
```

**Fix:**
```bash
rm -rf build
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

---

### Issue 2: "Compile driver first!" Error

**Symptom:** `benchmark_runner.py` fails with "Compile driver first!"

**Cause:** `traceon_driver` binary not found at `../build/traceon_driver`

**Fix:**
```bash
# From benchmarks/ directory
cd ..
cmake --build build --target traceon_driver -j
cd benchmarks
python benchmark_runner.py
```

---

### Issue 3: PyFastX Shows Impossibly High Numbers

**Symptom:** PyFastX shows 100M+ OPS/s with `?` marker.

**Explanation:** PyFastX uses disk caching with indices. After first run, subsequent lookups may be cached. This is not directly comparable to TracEon's memory-resident approach.

**Interpretation:**
- **TracEon:** Pure memory-resident (predictable performance)
- **PyFastX:** Disk cache + memory (variable performance)
- **BioPython:** Pure memory (OOM on large files)

---

### Issue 4: BioPython Crashes with "(BioPy ERR)"

**Symptom:** BioPython shows `(BioPy ERR)` in memory comparison.

**Causes:**
1. **Out of Memory:** BioPython loaded entire dataset into RAM
2. **Exception:** Malformed records or encoding issues

**Expected behavior:** BioPython will OOM on datasets > 200MB on 16GB systems. This is normal and demonstrates TracEon's efficiency.

---

### Issue 5: High Variance Between Runs

**Symptom:** Performance varies 2-3x between runs (e.g., 20M → 8M → 25M OPS/s).

**Causes:**
1. Background processes consuming RAM/CPU
2. Thermal throttling (CPU overheating)
3. OS page cache effects

**Fix:**
```bash
# 1. Close background applications
# 2. Run multiple iterations, take median
for i in {1..5}; do
  python benchmark_runner.py --sizes 100 --scenario WGS
done | grep "WGS " | sort -n -k 7 | sed -n '3p'  # Median result
```

---

## 🔄 Regression Testing (CI/CD)

### Automated Performance Checks

**GitHub Actions example:**
```yaml
name: Performance Regression

on: [pull_request]

jobs:
  benchmark:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      
      - name: Install Dependencies
        run: pip install psutil tqdm
      
      - name: Build TracEon
        run: |
          cmake -B build -DCMAKE_BUILD_TYPE=Release
          cmake --build build --target traceon_driver -j
      
      - name: Run Regression Check
        run: |
          cd benchmarks
          python benchmark_runner.py --sizes 10,100 --quick
          python check_regression.py --baseline v1.0.0
```

**Exit codes:**
- `0`: Performance within 10% of baseline (pass)
- `1`: Performance degraded 10-20% (warning)
- `2`: Performance degraded > 20% (fail)

---

## 📁 File Descriptions

### Scripts

**`benchmark_runner.py`** (Main benchmark suite)
- Generates synthetic FASTA/FASTQ data
- Tests TracEon, BioPython, PyFastX, SeqKit
- Measures parsing, lookup throughput, memory usage
- Output: Formatted table with comparisons

**`validate_real_data.py`** (Real-world validation)
- Downloads actual genomic datasets from NCBI/SRA
- Tests against production data (E. coli, C. elegans, human chr21)
- Validates correctness, not just performance

**`check_regression.py`** (CI/CD integration)
- Compares current results against baseline
- Enforces performance thresholds
- Exit codes for automated testing

### Configuration

**`traceon_benchmark.yml`** (Conda environment)
- Reproducible Python environment
- Pins all dependency versions

---

## 🎓 Best Practices

### 1. Consistent Environment

```bash
# Use performance CPU governor (Linux)
sudo cpupower frequency-set -g performance

# Disable CPU boost (reduces variance)
echo 0 | sudo tee /sys/devices/system/cpu/cpufreq/boost

# Pin to specific cores (optional)
taskset -c 0-7 python benchmark_runner.py
```

### 2. Multiple Iterations

```bash
# Run 5 times, take median (reduces noise)
for i in {1..5}; do
  python benchmark_runner.py --sizes 100 >> results.txt
done
```

### 3. System Monitoring

```bash
# Monitor during benchmark
watch -n 1 'free -h; sensors | grep Core'
```

### 4. Result Archiving

```bash
# Save results with metadata
python benchmark_runner.py | tee results_$(date +%Y%m%d_%H%M%S).txt

# Add system info
uname -a >> results_*.txt
lscpu | grep "Model name" >> results_*.txt
```

---

## 📊 Understanding Performance Tiers

### The L3 Cache Effect

**Why performance drops with dataset size:**
```
10MB Dataset:
  CPU → L3 Cache (10ns) → Data
  Result: 50M OPS/s ✅

100MB Dataset:
  CPU → L3 Miss → RAM (100ns) → Data
  Result: 15M OPS/s (3x slower) ⚠️
```

**This is hardware physics, not a bug:**
- L3 cache: ~10ns latency
- Main RAM: ~100ns latency (10x slower)
- No software can eliminate this gap

### Short vs Long Reads

**Why PacBio is faster than WGS:**
- **WGS:** 666K sequences, 150bp each → High index overhead
- **PacBio:** 3K sequences, 15Kbp each → Low index overhead
- Hash map lookup time dominates for WGS
- Data retrieval dominates for PacBio

---

## 🔗 References

### Comparison Tools
- **BioPython:** https://biopython.org/
- **PyFastX:** https://github.com/lmdu/pyfastx
- **SeqKit:** https://bioinf.shenwei.me/seqkit/
- **SeqTK:** https://github.com/lh3/seqtk

### Datasets
- **NCBI SRA:** https://www.ncbi.nlm.nih.gov/sra
- **Ensembl:** http://www.ensembl.org/
- **UCSC Genome Browser:** https://genome.ucsc.edu/

### Performance Analysis
- [What Every Programmer Should Know About Memory](https://people.freebsd.org/~lstewart/articles/cpumemory.pdf)
- [Brendan Gregg's Performance Page](http://www.brendangregg.com/linuxperf.html)

---

## 🤝 Contributing Benchmarks

### Adding New Scenarios

Edit `benchmark_runner.py`:
```python
DATA_TYPES = [
    # Format: (label, extension, avg_read_length)
    ("MyScenario", "fastq", 250),  # Add here
]
```

### Adding New Comparisons

Add new script section:
```python
MY_TOOL_SCRIPT = """
import sys
# Your tool's benchmark code
"""
```

---

## 📞 Support

**Questions about benchmarks?**
- 📝 Check [Performance Profile](../docs/performance-profile.md)
- 💬 Open a [GitHub Discussion](https://github.com/woosflex/TracEon/discussions)
- 🐛 Report benchmark issues: [GitHub Issues](https://github.com/woosflex/TracEon/issues)

---

**Last Updated:** December 16, 2025  
**Version:** 1.0.0 "Avalon"  
**Platform:** Linux, macOS, Windows (with Python)

*"Trace On" - Benchmarking legendary performance.* 📊