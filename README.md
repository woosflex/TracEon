# TracEon
**TracEon** is an exploratory C++ implementation of a high-performance, in-memory cache designed to accelerate genomic data access by eliminating repetitive file parsing.

This project was developed as part of the MSc Bioinformatics program to bridge the gap between biological domain knowledge and low-level systems programming. The name is a nod to the "Trace On" command from Fate/stay night, re-contextualized for tracing biological data across eons.

## Key Features
- **In-Memory Key-Value Store:** Built around `std::unordered_map` for fast, O(1) average-case data retrieval.

- **Pluggable Compression:** Uses a Strategy design pattern to allow for different encoding methods. Currently implemented:

  - **`PlainTextStrategy:`** A simple pass-through for baseline testing.

  - **`TwoBitDnaStrategy:`** Compresses DNA sequences by 4x using 2-bit encoding.

- **FASTA Parser:** Efficiently loads and parses standard FASTA files.

- **Persistent Cache:** Ability to `save()` the in-memory state to a compact binary file and `restore()` it on subsequent runs, avoiding the need to reparse large source files.

- **Cross-Platform:** Developed and tested on WSL2 (Ubuntu) and built with modern C++20 and CMake.

## Performance Benchmark
Benchmarks were run against a large, generated FASTA file (50,000 sequences of 1,000 bases each) to compare traditional disk parsing with TracEon's in-memory retrieval.

| Operation                       |    Time Taken     |
|:--------------------------------|:-----------------:|
| **Disk-Based Retrieval**        |   **83.83 ms**    |
| TracEon Initial Load (once)     |     174.47 ms     |
| **TracEon In-Memory Retrieval** |   **0.0013 ms**   |

**Conclusion:** Accessing a sequence from the cache is over **64,000 times faster** than parsing from disk. The initial load time is paid back after just 3 retrievals.

## Building from Source
- Prerequisites 
1. A C++20 compliant compiler (e.g., g++)
2. CMake (version 3.20 or higher)
3. Git

- Build Steps
1. Clone the repository:
```shell
git clone https://github.com/woosflex/TracEon.git
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
Here's a simple example of how to use the TracEon library in your own C++ project.

```c++
#include <iostream>
#include "Cache.h" // Include the main TracEon header

int main() {
// Create a cache instance
TracEon::Cache cache;

    // Use 2-bit compression for DNA
    cache.setStrategy(TracEon::StrategyType::TwoBit);

    // Load a FASTA file once
    cache.loadFasta("path/to/your/genome.fasta");

    // Save the processed cache for future runs
    cache.save("genome_cache.bin");

    // Access sequences instantly from memory
    std::string my_sequence = cache.get("gene_xyz_id");
    if (!my_sequence.empty()) {
        std::cout << "Found sequence: " << my_sequence.substr(0, 50) << "..." << std::endl;
    }

    return 0;
}
```

## Future Work
- [ ] FASTQ file support with quality score compression.
- [ ] Python and R bindings via pybind11 and Rcpp.
- [ ] Multi-threaded support for concurrent data access.
- [ ] Domain-specific query language (e.g., get_gc_content).

## License
This project is licensed under the MIT License. See the LICENSE file for details.