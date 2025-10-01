# TracEon (v0.0.1)
TracEon is a high-performance, in-memory caching library for bioinformatics, engineered in C++ to dramatically accelerate genomic data access by eliminating repetitive file parsing. 

This project was developed as part of the MSc Bioinformatics program to bridge the gap between biological domain knowledge and low-level systems programming. The name is a nod to the "Trace On" command from Fate/stay night, re-contextualized for tracing biological data across eons.

## Key FeaturesHigh-Performance 
- **Multithreading**: Utilizes a robust producer-consumer model to parse large, uncompressed FASTA and FASTQ files in parallel, maximizing throughput on multi-core systems.
- **Robust Gzip Support**: Employs a single-threaded streaming parser for gzipped (.gz) files to ensure correctness, as compressed streams are not seekable.
- **Optimized In-Memory Storage**: Built around std::unordered_map for O(1) average-case data retrieval, ensuring the fastest possible random access.
- **Uncompressed Binary Cache**: The save() and restore() operations use a fast, uncompressed binary format, prioritizing speed to make cache persistence an I/O-bound operation rather than a CPU-bound one.
- **Cross-Platform**: Developed and tested on WSL2 (Ubuntu) and macOS, built with modern C++20 and CMake.
## Performance Benchmark
Benchmarks were run on a Google Cloud Platform VM (AMD EPYC 7B12, 8 cores) to test performance against large, real-world datasets: a 7.5 GB FASTQ file and the human genome assembly (FASTA), both compressed and uncompressed.

The results demonstrate the effectiveness of the multithreaded parser for plain text and the dramatic speedup achieved by using the uncompressed binary cache on subsequent runs.

| Metric              | FASTQ (7.5 GB) | FASTA (2.9 GB)  | FASTA.gz (841 MB)  |
|:--------------------|:--------------:|:---------------:|:------------------:|
| Number of Sequences |  23.6 million  |       194       |        194         |
| Initial Load Time   |    132.7 s     |     **30.6 s**      |       44.8 s       |
|Cache Restore Time|     63.7 s     |      2.8 s      |       **1.9 s**        |
|Restore Speedup|      2.1×      |      10.9×      |       **24.0×**        |
|Memory Footprint|    ~16.0 GB    |     ~2.9 GB     |      ~2.9 GB       |
|Lookup Throughput|    491k/sec    |    1.69M/sec    |     **1.76M/sec**      |

## Benchmark Analysis
- Multithreading is Highly Effective: The load time for the 2.9 GB uncompressed FASTA file (30.6s) is significantly faster than the smaller, but gzipped, 841 MB version (44.8s), showcasing the power of the parallel parser on uncompressed data.
- Uncompressed Caching is Key: The most dramatic result is the restore time. Restoring the human genome cache takes just 2-3 seconds, a speedup of 10-24x over the initial load, proving the effectiveness of the uncompressed binary format.
- FASTQ Complexity: The FASTQ file, with its millions of small records, presents a greater challenge, resulting in slower load times and lower lookup throughput due to the overhead of managing a much larger hash map.

## Getting Started
For end-users, pre-compiled binaries are available on the GitHub Releases page. It is recommended to download the appropriate binary for your system.
For developers who wish to build from source, follow the steps below.
### Building from Source

**Prerequisites**
- A C++20 compliant compiler (e.g., g++)
- CMake (version 3.20 or higher)
- Git

**Build Steps**
1. Clone the repository:
   ```shell
   git clone [https://github.com/woosflex/TracEon.git](https://github.com/woosflex/TracEon.git)
   cd TracEon
   ```
2. Configure with CMake:
    ```shell
    cmake -S . -B build
    ```
3. Compile the library and tests:
    ```shell
    cmake --build build
    ```
4. Run the unit tests (optional but recommended):
    ```shell
    cd build
    ctest
    ```
## Usage Example (C++)
```c++
#include <iostream>
#include "Cache.h" // Include the main TracEon header

int main() {
// 1. Create a cache instance
TracEon::Cache cache;

    // 2. Load a FASTA or FASTQ file (multithreaded for plain text)
    cache.loadFile("path/to/your/genome.fastq");

    // 3. Save the fast, uncompressed binary cache
    cache.save("genome_cache.bin");

    // 4. On subsequent runs, restore the cache instantly
    // TracEon::Cache new_cache;
    // new_cache.restore("genome_cache.bin");

    // 5. Access sequences with microsecond latency
    std::string my_sequence = cache.get("sequence_id_123");
    if (!my_sequence.empty()) {
        std::cout << "Found sequence: " << my_sequence.substr(0, 50) << "..." << std::endl;
    }

    return 0;
}
```
## Future Work
The project's focus is on maximizing speed for bioinformatics workflows. The next development phase will focus on expanding its accessibility.
- [ ] Develop a C API: Create a stable, C-compatible Application Binary Interface (ABI) to serve as the foundation for cross-language support.
- [ ] Create Python and R Bindings: Build user-friendly Python and R libraries on top of the C API, making TracEon accessible to a broader scientific community.

## License
This project is licensed under the MIT License. See the LICENSE file for details.