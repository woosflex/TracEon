/**
 * TracEon Simple Usage Example
 * * Demonstrates the recommended workflow for high-performance genomic data access.
 */

#include <iostream>
#include <chrono>
#include <TracEon.h>

int main() {
    std::cout << "╔════════════════════════════════════════╗\n";
    std::cout << "║   TracEon v" << TracEon::VERSION << " \"" << TracEon::CODENAME << "\"      ║\n";
    std::cout << "║   High-Performance Genomic Cache       ║\n";
    std::cout << "╚════════════════════════════════════════╝\n\n";
    
    // Step 1: Create cache instance
    TracEon::Cache cache;
    std::cout << "Created cache instance\n";
    
    // Step 2: Load FASTA/FASTQ file (multithreaded parsing for files > 10MB)
    // Note: Assuming running from build/ directory, adjust path as needed
    std::string p = "../test_data/simple.fasta";
    std::cout << "Loading genome file: " << p << "...\n";
    
    auto start = std::chrono::high_resolution_clock::now();
    try {
        cache.loadFile(p);
    } catch (...) {
        std::cerr << "File not found. Please run from build directory.\n";
        return 1;
    }
    
    auto load_time = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - start
    ).count();
    std::cout << "  Loaded " << cache.size() << " sequences in " << load_time << "ms\n";
    
    // Step 3: Save binary cache for instant loading next time
    std::cout << "\nSaving binary cache...\n";
    cache.save("example_cache.bin");
    std::cout << "  Cache saved\n";
    
    // Step 4: Restore from binary (memory-mapped, zero-copy)
    std::cout << "\nRestoring from binary cache...\n";
    TracEon::Cache fast_cache;
    start = std::chrono::high_resolution_clock::now();
    fast_cache.restore("example_cache.bin");
    auto restore_time = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - start
    ).count();
    std::cout << "  Restored " << fast_cache.size() << " sequences in " << restore_time << "ms\n";
    
    // Step 5: Zero-copy access (recommended for hot paths)
    std::cout << "\nAccessing sequences:\n";
    std::string_view seq = fast_cache.getView("seq1");
    if (!seq.empty()) {
        std::cout << "  Found seq1: " << seq.substr(0, 50) 
                  << (seq.size() > 50 ? "..." : "") << "\n";
        std::cout << "  Length: " << seq.size() << " bp\n";
        std::cout << "  Memory: Zero-copy view (no allocation)\n";
    }
    
    // Step 6: Legacy API (creates copy, slower but compatible)
    std::string seq_copy = fast_cache.get("seq2");
    if (!seq_copy.empty()) {
        std::cout << "\n  Found seq2 (copy): " << seq_copy.substr(0, 50) 
                  << (seq_copy.size() > 50 ? "..." : "") << "\n";
    }
    
    std::cout << "\n=== Example Complete ===\n";
    return 0;
}