/**
 * @file TracEon.h
 * @brief TracEon: High-Performance Genomic Data Cache
 * @version 1.0.0-rc1
 * @date December 2025
 * 
 * @section DESCRIPTION
 * TracEon is a C++20 library for indexing and caching genomic data (FASTA/FASTQ)
 * with zero-copy architecture and lock-free concurrency.
 * 
 * @section FEATURES
 * - Zero-copy string views (no allocation overhead)
 * - Lock-free concurrent reads (2x throughput)
 * - Memory-mapped binary cache (instant restore)
 * - Multithreaded parsing (8-core scaling)
 * - 40M+ OPS/s random access on modern hardware
 * 
 * @section USAGE
 * @code{.cpp}
 * #include <TracEon.h>
 * 
 * TracEon::Cache cache;
 * cache.loadFile("genome.fastq");          // Parse text file
 * cache.save("genome.bin");                // Save binary cache
 * 
 * // Later session:
 * cache.restore("genome.bin");             // Instant load via mmap
 * std::string_view seq = cache.getView("read_001"); // Zero-copy access
 * @endcode
 * 
 * @section PERFORMANCE
 * Benchmarked on AMD EPYC 7B12 (8 cores):
 * - Short reads (100MB FASTQ): 13M OPS/s
 * - Long reads (100MB FASTQ): 57M OPS/s
 * - Memory usage: 180MB for 100MB dataset
 * 
 * @section LICENSE
 * MIT License - See LICENSE file
 * 
 * @section ARCHITECTURE
 * - Cache: High-level API wrapper
 * - SmartStrategy: Core engine (arena allocation, lock-free reads)
 * - MapDefs: Conditional robin_hood vs std::unordered_map
 * 
 * For architectural details, see docs/architecture/ADR-001-lock-free-reads.md
 */

#ifndef TRACEON_H
#define TRACEON_H

// Public API
#include "Cache.h"
#include "DecodedRecordTypes.h"
#include "RecordTypes.h"

// Advanced users may need direct access
#include "SmartStrategy.h"

/**
 * @namespace TracEon
 * @brief Main namespace for all TracEon functionality
 */
namespace TracEon {
    /**
     * @brief Library version string
     * Format: MAJOR.MINOR.PATCH
     * Codename: Caliburn (v1.4.0)
     */
    constexpr const char* VERSION = "1.4.0";
    constexpr const char* CODENAME = "Caliburn";
    
    /**
     * @brief Check if lock-free optimization is active
     * @return true if using atomic acquire/release pattern
     */
    inline bool hasLockFreeReads() {
        return true; // Always enabled in RC1+
    }
    
    /**
     * @brief Check if Robin Hood hashing is available
     * @return true if third_party/robin_hood.h was found at compile time
     */
    inline bool hasRobinHoodOptimization() {
#ifdef TRACEON_USE_ROBIN_HOOD
        return true;
#else
        return false;
#endif
    }
}

#endif // TRACEON_H