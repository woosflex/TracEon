# TracEon Benchmarking Suite

## Quick Start

### Prerequisites
```bash
pip install psutil tqdm requests
```
### Running Matrix Benchmark
Tests TracEon across file sizes (10MB - 1GB) and biological scenarios (WGS, WES, PacBio, Genome).
```{bash}
# From project root
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release
cmake --build build --target traceon_driver
python benchmarks/benchmark_runner.py
```
### Running Real-Data Validation
Downloads and tests against actual genomic datasets from NCBI/SRA:
```{bash}
python benchmarks/validate_real_data.py
```
## Interpreting Results
### Throughput Metrics
- **30M+ OPS/s**: Excellent (L3 cache-resident data). Typical for datasets < 50MB.
- **13M-25M OPS/s**: Good (RAM latency wall). Typical for datasets > 100MB.
- **< 10M OPS/s**: Investigate. Check if you compiled in DEBUG mode.
### Comparison Flags
- **(?)**: PyFastX result may be cached/invalid (impossibly high numbers for RAM access).
- **CRASH**: Competitor (BioPython) ran out of memory.
## Performance Expectations (Release Build)
| Dataset Type | File Size | Expected OPS/s |
| --- | --- | --- |
|WGS (short reads)|10MB|40M - 55M|
|WGS (short reads)|100MB|12M - 18M|
|PacBio (long reads)|100MB|45M - 60M|
|Reference Genome|100MB|15M - 25M|

## Troubleshooting

### **Getting 6M OPS/s instead of 40M?**
- Check build mode: Must use `-DCMAKE_BUILD_TYPE=Release`.
- Check CMake output to ensure `robin_hood.h` was found.