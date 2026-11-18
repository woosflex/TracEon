#ifndef TRACEON_CACHE_H
#define TRACEON_CACHE_H

#include <string>
#include <string_view>
#include <vector>
#include <memory>
#include <optional>

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
        std::string get(std::string_view key);

        /**
         * @brief Retrieve a borrowed view of the sequence payload.
         * @return A std::string_view pointing directly to the immutable arena.
         *         FASTA/FASTQ header descriptions and FASTQ plus-line metadata
         *         are not retained; this method never returns quality scores.
         *         Empty also represents a valid zero-length FASTQ sequence;
         *         use hasSequence() when absence must be distinguished.
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
        std::string_view getView(std::string_view key);

        /**
         * @brief Retrieve non-owning sequence and quality views for a record in a single lookup.
         * For FASTA records, quality is empty.
         */
        std::optional<FastqRecordView> getRecordView(std::string_view key);

        /**
         * @brief Retrieve owned sequence and quality payloads for a FASTQ ID.
         * Header descriptions and plus-line metadata are not retained.
         */
        std::optional<DecodedFastqRecord> getFastqRecord(std::string_view key);

        /**
         * @brief Check if a key exists without retrieving data.
         */
        bool hasSequence(std::string_view key);

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
         * @brief Save the current index to a portable v5 binary cache file.
         *
         * v5 uses explicit little-endian fields, stores full identifiers (not
         * implementation-specific hash values), preserves detected format, and
         * wraps the payload in an integrity-checked streamed LZ4 Frame.
         * @note Entries added via set() are included in the serialized payload.
         * @throws std::runtime_error if the file cannot be written.
         */
        void save(const std::string& filepath);

        // Explicit mapped persistence shares read-only payload pages across processes.
        // The backing file must remain immutable; replace it atomically, never in place.
        // restore() detects storage automatically. The one-argument save stays LZ4.
        void save(const std::string& filepath, CacheStorage storage);
        
        /**
         * @brief Restore a v4 or v5 binary cache with validated decoded or mapped storage.
         *
         * LZ4 files are mapped temporarily, then fully decoded and validated.
         * v5 mapped payloads retain the read-only mapping after validation.
         * Never modify/truncate an open mapped file; replace it atomically.
         * Both paths rebuild a private index before publication. v5 is canonical and
         * portable; v4 remains readable for migration. v1/v2/v3 are rejected.
         * Physical bytes after the declared frame are rejected.
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
         * @brief Backing bytes: text arena + manual entries + retained mapped
         * file (including its header). This is not private RSS or PSS. Forwarded to SmartStrategy::getArenaBytes().
         * Used by the remote STATS message (v2.3.0 remote access).
         */
        size_t getArenaBytes() const;

        IndexMode getIndexMode() const;

        /**
         * @brief Retrieve the quality string for a key (FASTQ only).
         * @return A copy of the quality line, or empty string if not found /
         *         FASTA (no quality stored).
         */
        std::string getQuality(std::string_view key) const;

        /**
         * @brief Enumerate all keys currently in the cache (copies).
         */
        std::vector<std::string> getAllKeys() const;

        /**
         * @brief Detected format of the loaded data (UNKNOWN when empty).
         */
        FileFormat getDetectedFormat() const;

    private:
        std::unique_ptr<SmartStrategy> m_strategy;
    };

} // namespace TracEon

#endif //TRACEON_CACHE_H
