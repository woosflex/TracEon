/**
 * @file TracEon.h
 * @brief TracEon genomic snapshot cache
 * @version 3.0.0-dev
 * @date 2026
 * 
 * @section DESCRIPTION
 * TracEon is a C++20 library for indexing and caching genomic data (FASTA/FASTQ)
 * with zero-copy architecture and lock-free concurrency.
 * 
 * @section FEATURES
 * - Borrowed string views for local sequence-payload access
 * - Immutable-phase concurrent reads without a read-side mutex
 * - `.traceon` v5 portable binary cache (v4 migration reader; LZ4 Frame + CRC32C)
 * - Multithreaded parsing for supported large inputs
 * 
 * @section USAGE
 * @code{.cpp}
 * #include <TracEon.h>
 * 
 * TracEon::Cache cache;
 * cache.loadFile("genome.fastq");          // Parse text file
 * cache.save("genome.bin");                // Save portable binary cache (v5)
 * 
 * // Later session:
 * cache.restore("genome.bin");             // Decompress, validate, then publish
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
 * @section MIGRATION
 * - `save()` writes portable v5; `restore()` reads v4 and v5. v1/v2/v3 are
 *   rejected and must be regenerated.
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
 * - Crc32c: CRC-32C with SSE4.2/AArch64/table dispatch (v4/v5 checksum)
 * 
 * For architectural details, see docs/architecture/ADR-001-lock-free-reads.md
 */

#ifndef TRACEON_H
#define TRACEON_H

// Public API
#include "Cache.h"
#include "Crc32c.h"
#include "DecodedRecordTypes.h"

// Advanced users may need direct access
#include "SmartStrategy.h"

/**
 * @namespace TracEon
 * @brief Main namespace for all TracEon functionality
 */
namespace TracEon {
    /**
     * @brief Library version string
     * Development identifier. The latest published release remains v2.3.0.
     */
    constexpr const char* VERSION = "3.0.0-dev";
    constexpr const char* CODENAME = "Unreleased";
}

#endif // TRACEON_H
