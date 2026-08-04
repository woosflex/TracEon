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
    // - The Map Key (std::string) owns the normalized identifier used for hashing/lookups.
    // - SequenceView.id, sequence, and quality are NON-OWNING views pointing
    //   into either the mmap buffer (binary mode) or the text_arena_ buffer (load mode).
    //   For entries added via addEntry(), they point into manual_store_.
    std::string_view id;       
    std::string_view sequence; 
    std::string_view quality;  
};

// GenomeIndex: string keys with transparent lookup (find by string_view, no allocation).
using GenomeIndex = StringHashMap<SequenceView>;
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
     */
    void loadFile(const std::string& filepath);
    
    /**
     * @brief Explicitly load a GZIP-compressed file.
     * Throws if file is not valid GZIP.
     */
    void loadGzipFile(const std::string& filepath);

    void loadBinary(const std::string& filepath);
    void saveBinary(const std::string& filepath) const;

    // --- Manual Entry Support ---
    void addEntry(const std::string& id, const std::string& seq, const std::string& qual);

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

    void determine_format_from_cache();
    bool isNucleotideSequence(std::string_view data) const;
    bool hasRNA(std::string_view data) const;
    uint64_t hash_key(std::string_view key) const;

    void clearInternal();

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