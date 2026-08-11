#ifndef TRACEON_SMARTSTRATEGY_H
#define TRACEON_SMARTSTRATEGY_H

#include "IEncodingStrategy.h"
#include "MapDefs.h"

#include <string>
#include <string_view>
#include <vector>
#include <deque>
#include <memory>
#include <shared_mutex>
#include <variant>
#include <filesystem>
#include <atomic>
#include <functional>

namespace TracEon {

struct MMapHandle;

enum class FileFormat : uint8_t {
    DNA_FASTA = 0x00, RNA_FASTA = 0x01, PROTEIN_FASTA = 0x02,
    DNA_FASTQ = 0x03, RNA_FASTQ = 0x04, PROTEIN_FASTQ = 0x05,
    UNKNOWN = 0xFF
};

enum class CompressionMode {
    LZ4Default, ///< LZ4_compress_default() — fast, general case
    LZ4HC,      ///< LZ4_compress_HC() — higher ratio for large nucleotide data
};

struct SequenceView {
    // NOTE:
    // - The Map Key (std::string_view) is NON-OWNING: it points into either
    //   the text_arena_ buffer (parse and binary v4 restore) or manual_store_
    //   (addEntry()). Both backing stores are stable for the map's lifetime
    //   (arena is immutable-after-load; deque push_back never invalidates).
    //   No std::string is constructed on any insert path.
    // - SequenceView.id, sequence, and quality are NON-OWNING views pointing
    //   into the same stable buffers. For entries added via addEntry(), they
    //   point into manual_store_.
    std::string_view id;       
    std::string_view sequence; 
    std::string_view quality;  
};

// GenomeIndex: string_view keys (zero-copy, arena-backed — no std::string
// allocation on the insert path) with transparent lookup (find by
// string_view / const char* / std::string, no allocation on the read path).
using GenomeIndex = StringViewHashMap<SequenceView>;
using NGSIndex = HashMap<uint64_t, SequenceView>;

enum class IndexMode { GENOME, NGS };

class SmartStrategy : public IEncodingStrategy {
public:
    explicit SmartStrategy(IndexMode mode = IndexMode::GENOME);
    virtual ~SmartStrategy();

    std::vector<unsigned char> encode(const std::string& data, DataTypeHint hint = DataTypeHint::Generic) const override;
    std::string decode(const std::vector<unsigned char>& data) const override;

    /**
     * @brief Load file with automatic format detection.
     * Detects .gz extension and GZIP magic bytes.
     * @throws std::runtime_error on corrupt/truncated input — truncated GZIP
     *         streams, trailing garbage after the last GZIP member, truncated
     *         tail members in concatenated GZIP, and OOM-guard rejections.
     */
    void loadFile(const std::string& filepath);
    
    /**
     * @brief Explicitly load a GZIP-compressed file.
     * Throws if file is not valid GZIP.
     * @throws std::runtime_error if the file is not valid GZIP, or if the
     *         stream is truncated / has trailing data after the last member
     *         (data-integrity rejection — partial data is never served).
     */
    void loadGzipFile(const std::string& filepath);

    /**
     * @brief Load a binary cache file (.traceon v4 only).
     *
     * v2.0.0 removed the v1/v2/v3 readers — the v4 format ("TRO\x04", LZ4
     * Frame + whole-payload CRC32C) is the ONLY readable binary format. Old
     * caches are rejected with a clear "regenerate" error.
     * @throws std::runtime_error on malformed/implausible cache files: wrong
     *         magic or version, unsupported codec/mode, truncated payloads or
     *         frames, wrong logical length, checksum mismatch, record counts
     *         exceeding the payload bytes, or OOM-guard rejections. A failed
     *         load never publishes partial state (failure atomicity).
     */
    void loadBinary(const std::string& filepath);
    /**
     * @brief Write the current cache as a v4 binary cache file.
     *
     * v4 = "TRO\x04" header (codec flags, index mode, logical length, frame
     * length, CRC32C) + a streamed LZ4 Frame. The CRC32C checksum is computed
     * incrementally over the uncompressed payload as it passes through the
     * compressor — no second full-payload pass.
     * @throws std::runtime_error if the file cannot be written or LZ4 Frame
     *         compression fails.
     */
    void saveBinary(const std::string& filepath) const;

    // --- Manual Entry Support ---
    /**
     * @brief Add a record during the build phase.
     * @throws std::logic_error if called after a load (immutable-after-load
     *         contract, ADR-001): once data_loaded_ is set the index is frozen.
     *         Legal only before any loadFile()/loadBinary() or after
     *         clearCache().
     */
    void addEntry(const std::string& id, const std::string& seq, const std::string& qual);

    /**
     * @brief Zero-copy view of a sequence's data.
     *
     * @return A std::string_view pointing directly into the immutable arena
     *         (text_arena_/manual_store_/mmap), or an empty view if not found.
     *
     * @warning LIFECYCLE CONTRACT (reader quiescence — see ADR-001 and the
     *          README "Lifecycle contract" section): the returned view is
     *          NON-OWNING and is valid only while the same loaded snapshot
     *          stays installed. It becomes DANGLING the moment clearCache(),
     *          a reload (loadFile/loadGzipFile/loadBinary), or destruction
     *          begins. Reads are safe concurrently with other reads, but the
     *          application MUST stop using all views from a snapshot before
     *          clearing/reloading/destroying the cache. The write-side mutex
     *          serializes writers only; it does NOT reclaim memory from
     *          readers and is deliberately NOT taken on this read path.
     *
     * @note TRACEON_DEBUG_LIFECYCLE (opt-in) adds a reader counter that warns
     *       in clearInternal() if a getView() lookup overlaps a teardown — a
     *       diagnostic for coordinated misuse, not a synchronization mechanism.
     */
    std::string_view getView(const std::string_view& sequence_id) const;
    std::string getSequence(const std::string& sequence_id) const;
    std::string getQuality(const std::string& sequence_id) const;
    bool hasSequence(const std::string& sequence_id) const;
    
    size_t getFileCacheSize() const;
    std::vector<std::string> getAllKeys() const;
    void clearCache();
    
    FileFormat getDetectedFormat() const { return detected_format_; }
    IndexMode getIndexMode() const; 

private:
    std::vector<char> text_arena_;
    std::unique_ptr<MMapHandle> mmap_handle_;
    std::variant<GenomeIndex, NGSIndex> file_cache_;
    
    // Stable string storage for entries added via addEntry().
    // std::deque provides reference stability: push_back never invalidates
    // references/pointers to existing elements, so string_views into it are safe.
    std::deque<std::string> manual_store_;

    // Running total of bytes held in manual_store_. Checked against
    // getAvailableMemory() in addEntry() so manual entries — the only
    // remaining unbounded growth path in this class — fail loudly instead
    // of growing without limit.
    std::atomic<size_t> manual_store_bytes_{0};

    // Running total of the serialized (uncompressed) binary-cache payload
    // size, updated incrementally at every record insertion point instead
    // of being recomputed by walking the whole map at saveBinary() time.
    // Starts at sizeof(uint64_t): serializePayload() always writes an 8-byte
    // record-count header up front, even with zero records, so this base
    // must be present from construction (not just after a parse/restore
    // that calls refreshPayloadEstimate()) for addEntry()-only usage to
    // produce a correct estimate.
    mutable std::atomic<size_t> serialized_size_estimate_{sizeof(uint64_t)};

    FileFormat detected_format_;
    IndexMode index_mode_;
    mutable std::shared_mutex cache_mutex_;

    /**
     * Lock-Free Read Signal
     * True only when all data is loaded and immutable.
     */
    std::atomic<bool> data_ready_{false};

#ifdef TRACEON_DEBUG_LIFECYCLE
    /**
     * Reader-quiescence diagnostic (opt-in via -DTRACEON_DEBUG_LIFECYCLE).
     *
     * Incremented for the duration of each getView() call; clearInternal()
     * warns (std::cerr) if any readers are still inside getView() when the
     * snapshot is torn down. LIMITATION: catches a reader whose LOOKUP
     * overlaps a clear; it CANNOT detect a caller that retained the returned
     * std::string_view and dereferences it after clearCache()/reload/
     * destruction completes. Diagnostic only — NOT synchronization (see
     * ADR-001 reader-quiescence section).
     */
    mutable std::atomic<size_t> active_readers_{0};
#endif

    /**
     * Immutable-after-load signal (Bug 3 fix, ADR-001).
     *
     * Distinct from data_ready_: data_ready_ is also set by addEntry() so
     * lock-free reads can see manually-added entries during a build phase,
     * which would make it useless as an "is the cache frozen?" test. This
     * flag is set true ONLY by the load paths (loadFile/loadGzipFile via
     * parseArena(), and loadBinary()) and reset by clearInternal()
     * (clearCache()). Once true, the index is truly immutable: addEntry()
     * throws std::logic_error instead of mutating the map under concurrent
     * lock-free readers.
     */
    std::atomic<bool> data_loaded_{false};

    void determine_format_from_cache();
    bool isNucleotideSequence(std::string_view data) const;
    bool hasRNA(std::string_view data) const;
    uint64_t hash_key(std::string_view key) const;

    void clearInternal();

    /**
     * @brief Set by load paths (loadFile/loadGzipFile/loadBinary) so
     * addEntry() can enforce the immutable-after-load contract.
     */
    void markDataLoaded();

    /**
     * @brief Serialize payload (count + all records), pushing bytes through
     * `sink` as they're produced rather than materializing a full in-memory
     * copy. Used by saveBinary() (v3 streaming format) so peak memory during
     * save stays bounded regardless of dataset size.
     */
    void serializePayload(const std::function<void(const char*, size_t)>& sink) const;

    /**
     * @brief Choose a compression algorithm based on payload size and detected format.
     *
     * Rules (in priority order):
     *  1. payload_size > 10 MiB AND format is DNA/RNA → LZ4HC  (high repetition)
     *  2. everything else                              → LZ4Default
     */
    CompressionMode selectCompressionStrategy(size_t payload_size) const;

    /**
     * @brief Serialized (uncompressed) payload size, used to pick LZ4 vs
     * LZ4HC before streaming starts (selectCompressionStrategy() needs a
     * size up front). Backed by serialized_size_estimate_, which is updated
     * incrementally at every record insertion rather than recomputed by
     * walking the whole map here.
     */
    size_t estimatePayloadSize() const;

    /**
     * @brief Recompute serialized_size_estimate_ by walking file_cache_.
     * Called once after a full (re)parse or binary-cache restore where the
     * exact size isn't already known for free; addEntry() maintains the
     * estimate incrementally after this point without needing another walk.
     */
    void refreshPayloadEstimate();

    /**
     * @brief Scan a compressed GZIP file for stream boundaries.
     * Returns byte offsets of each GZIP stream header found.
     * Size == 1 means single-stream; size > 1 means concatenated streams.
     */
    std::vector<size_t> scanGzipStreams(const std::string& filepath) const;

    /**
     * @brief Decompress concatenated GZIP streams in parallel into text_arena_.
     * Each stream is decompressed by a separate thread using zlib-ng inflate().
     * @param filepath Path to the compressed file
     * @param stream_offsets Stream start byte offsets from scanGzipStreams()
     */
    void loadGzipParallel(const std::string& filepath,
                          const std::vector<size_t>& stream_offsets);

    /**
     * @brief Single-threaded GZIP decompression into text_arena_ (original path).
     */
    void loadGzipSingleStream(const std::string& filepath);

    /**
     * @brief Cross-platform available physical memory in bytes.
     * Uses /proc/meminfo (Linux), sysctl (macOS), GlobalMemoryStatusEx (Windows).
     * Returns 0 on failure (caller falls back to geometric growth).
     */
    static size_t getAvailableMemory();
    
    /**
     * @brief Normalize multi-line FASTA sequences in-place within text_arena_.
     * Removes embedded newlines from sequence regions so parsers always see
     * single-line sequences. Header lines are preserved verbatim.
     * Only call for FASTA (first char '>'). No-op for empty arenas.
     */
    void normalizeFastaArena();

    /**
     * @brief Internal GZIP loader (used by both public methods)
     * Reads compressed stream into text_arena_ without parsing.
     */
    void loadGzipInternal(const std::string& filepath);
    
    // Helper to centralize parsing logic
    void parseArena();

    template <typename MapType> void parseFastaInternal(std::string_view content, MapType& map);
    template <typename MapType> void parseFastqInternal(std::string_view content, MapType& map);
    template <typename MapType> void parseFastaMultithreadedTemplate(std::string_view content, MapType& dest_map);
    template <typename MapType> void parseFastqMultithreadedTemplate(std::string_view content, MapType& dest_map);

};

} // namespace TracEon

#endif // TRACEON_SMARTSTRATEGY_H