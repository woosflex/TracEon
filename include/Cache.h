#ifndef TRACEON_CACHE_H
#define TRACEON_CACHE_H

#include <string>
#include <string_view>
#include <vector>
#include <memory>
#include <optional>

#include "IEncodingStrategy.h"
#include "DecodedRecordTypes.h"
#include "SmartStrategy.h"

namespace TracEon {

    /**
     * @class Cache
     * @brief High-level interface for the TracEon genomic data store.
     * * Provides a unified API for loading, saving, and accessing genomic sequences.
     * Internally routes requests to the zero-copy SmartStrategy (for file-backed data)
     * or a manual fallback store (for runtime-added data).
     */
    class Cache {
    public:
        explicit Cache(IndexMode mode = IndexMode::GENOME);
        ~Cache();

        // --- Data Access ---
        
        /**
         * @brief Retrieve a sequence by ID (Legacy API).
         * @return A copy of the sequence string, or empty string if not found.
         */
        std::string get(const std::string& key);

        /**
         * @brief Retrieve a zero-copy view of a sequence (High Performance).
         * @return A std::string_view pointing directly to the immutable arena.
         *
         * @warning LIFECYCLE CONTRACT (reader quiescence — see ADR-001 and the
         *          README "Lifecycle contract" section): the returned view is
         *          NON-OWNING and valid only while the same loaded snapshot
         *          stays installed. It becomes DANGLING the moment
         *          clearCache(), a reload (loadFile/restore), or destruction
         *          begins. Concurrent reads against a live snapshot are safe
         *          and lock-free, but callers MUST stop using every view from
         *          a snapshot before clearing/reloading/destroying the cache.
         *          The write-side mutex serializes writers only; it does NOT
         *          protect a view after the call returns.
         */
        std::string_view getView(const std::string& key);

        /**
         * @brief Retrieve both sequence and quality scores.
         */
        std::optional<DecodedFastqRecord> getFastqRecord(const std::string& key);

        /**
         * @brief Check if a key exists without retrieving data.
         */
        bool hasSequence(const std::string& key);

        // --- I/O Operations ---
        
        /**
         * @brief Parse a FASTA/FASTQ file into the cache.
         * @note Uses multithreaded parsing for files > 10MB.
         * @throws std::runtime_error if the file is corrupt or truncated —
         *         includes truncated GZIP streams, trailing garbage after the
         *         last GZIP member, and OOM-guard rejections.
         */
        void loadFile(const std::string& filepath);

        /**
         * @brief Explicitly load a GZIP-compressed FASTA/FASTQ file.
         *
         * Mirrors the README-documented example (`cache.loadGzipFile(...)`)
         * and the lifecycle contract above, which already name loadGzipFile()
         * as a reload path — this method completes the public Cache facade
         * (previously only SmartStrategy exposed it, so the documented example
         * did not compile).
         * @throws std::runtime_error if the file is not valid GZIP, or if the
         *         stream is truncated / has trailing data after the last member
         *         (data-integrity rejection — partial data is never served).
         * @note Failure atomicity: a pre-flight rejection (unreadable file or
         *       non-GZIP magic bytes) throws before teardown and PRESERVES the
         *       previously loaded snapshot; a mid-stream failure (truncated /
         *       corrupt GZIP stream, trailing garbage, OOM guard) throws after
         *       teardown has begun and leaves the cache EMPTY — never partial.
         *       See SmartStrategy::loadGzipFile() for the full contract.
         */
        void loadGzipFile(const std::string& filepath);
        
        /**
         * @brief Save the current index to a binary cache file (v4 format).
         *
         * v2.0.0 format: "TRO\x04" header (codec flags, index mode, logical
         * length, frame length, CRC32C) + streamed LZ4 Frame. The CRC32C
         * checksum is computed incrementally over the uncompressed payload as
         * it passes through the compressor — no second full-payload pass.
         * @note Entries added via set() are included in the serialized payload.
         * @throws std::runtime_error if the file cannot be written.
         */
        void save(const std::string& filepath);
        
        /**
         * @brief Instantaneously map a binary cache file into memory (v4 only).
         *
         * v1/v2/v3 caches are REJECTED ("unsupported cache version;
         * regenerate with v2.0.0") — v4 is the only readable format.
         * @throws std::runtime_error on malformed/implausible cache files:
         *         wrong magic/version, unsupported codec/mode, truncated data,
         *         wrong logical length, checksum mismatch, implausible size/
         *         count fields, or OOM-guard rejections. A failed load never
         *         publishes partial state (failure atomicity).
         */
        void restore(const std::string& filepath);

        /**
         * @brief Manually add a key-value pair (Heap storage).
         * @throws std::logic_error if called after a load (immutable-after-load
         *         contract, ADR-001). Legal only before any loadFile()/restore()
         *         or after clearCache().
         */
        void set(const std::string& key, const std::string& value);

        /**
         * @brief Release the current snapshot and reset the cache to empty.
         *
         * Exposes the documented lifecycle operation on the public facade
         * (see the README "Lifecycle contract" section and ADR-001). After
         * clearCache(), the cache is empty and set()/loadFile()/restore()
         * may be used again.
         * @warning Every std::string_view obtained from getView() before this
         *          call becomes DANGLING the moment clearCache() begins — stop
         *          using all views from the snapshot first (reader quiescence).
         */
        void clearCache();
        
        size_t size() const;

        /**
         * @brief Bytes held by the loaded snapshot's arena stores (text
         * arena + manual-entry bytes). Forwarded to SmartStrategy::getArenaBytes().
         * Used by the remote STATS message (v2.3.0 remote access).
         */
        size_t getArenaBytes() const;

        IndexMode getIndexMode() const;

        /**
         * @brief Retrieve the quality string for a key (FASTQ only).
         * @return A copy of the quality line, or empty string if not found /
         *         FASTA (no quality stored).
         */
        std::string getQuality(const std::string& key) const;

        /**
         * @brief Enumerate all keys currently in the cache (copies).
         */
        std::vector<std::string> getAllKeys() const;

        /**
         * @brief Detected format of the loaded data (UNKNOWN when empty).
         */
        FileFormat getDetectedFormat() const;

        // --- Internal / Diagnostics ---
        IEncodingStrategy* getStrategy() const;

    private:
        std::unique_ptr<SmartStrategy> m_strategy;
    };

} // namespace TracEon

#endif //TRACEON_CACHE_H