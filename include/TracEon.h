/**
 * @file TracEon.h
 * @brief TracEon: High-Performance Genomic Data Cache
 * @version 2.0.0
 * @date 2026
 * 
 * @section DESCRIPTION
 * TracEon is a C++20 library for indexing and caching genomic data (FASTA/FASTQ)
 * with zero-copy architecture and lock-free concurrency.
 * 
 * @section FEATURES
 * - Zero-copy string views (no allocation overhead)
 * - Lock-free concurrent reads (2x throughput)
 * - `.traceon` v4 binary cache (LZ4 Frame + whole-payload CRC32C, streamed)
 * - Multithreaded parsing (8-core scaling)
 * - 40M+ OPS/s random access on modern hardware
 * 
 * @section USAGE
 * @code{.cpp}
 * #include <TracEon.h>
 * 
 * TracEon::Cache cache;
 * cache.loadFile("genome.fastq");          // Parse text file
 * cache.save("genome.bin");                // Save binary cache (v4)
 * 
 * // Later session:
 * cache.restore("genome.bin");             // Instant load via mmap
 * std::string_view seq = cache.getView("read_001"); // Zero-copy access
 * @endcode
 * 
 * @section LIFECYCLE
 * getView() returns a NON-OWNING std::string_view that is valid only while
 * the same loaded snapshot stays installed. It becomes dangling once
 * clearCache(), a reload (loadFile/loadGzipFile/restore), or destruction
 * begins — stop using views from a snapshot before clearing/reloading. See
 * docs/architecture/ADR-001-lock-free-reads.md.
 *
 * @section MIGRATION (1.x -> 2.0.0)
 * - `.traceon` binary caches are v4 only; v1/v2/v3 files are rejected and
 *   must be regenerated with 2.0.0 (cache.save()).
 * - Corrupt/truncated GZIP input throws std::runtime_error.
 * - set()/addEntry() after a load throws std::logic_error.
 * - No ABI stability is promised across 1.x -> 2.x; recompile consumers.
 * - The k-mer C API ships as an opt-in separate target (`traceon_kmer`) with a
 *   hardened C boundary (noexcept + status codes, caller-owned iterators,
 *   freeze semantics) — not part of the core umbrella.
 * 
 * @section LICENSE
 * MIT License - See LICENSE file
 * 
 * @section ARCHITECTURE
 * - Cache: High-level API wrapper
 * - SmartStrategy: Core engine (arena allocation, lock-free reads)
 * - MapDefs: ankerl::unordered_dense with transparent wyhash-style hashing
 * - Crc32c: CRC-32C with SSE4.2/AArch64/table dispatch (v4 checksum)
 * 
 * For architectural details, see docs/architecture/ADR-001-lock-free-reads.md
 */

#ifndef TRACEON_H
#define TRACEON_H

// Public API
#include "Cache.h"
#include "Crc32c.h"
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
     * Codename: Excalibur (v2.1.0) — the v2.0.0 "Durandal" codename is retired.
     */
    constexpr const char* VERSION = "2.1.0";
    constexpr const char* CODENAME = "Excalibur";
    
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