#ifndef TRACEON_SMARTSTRATEGY_H
#define TRACEON_SMARTSTRATEGY_H

#include "IEncodingStrategy.h"
#include "MapDefs.h"
#include "RecordTypes.h"

#include <string>
#include <string_view>
#include <vector>
#include <deque>
#include <memory>
#include <shared_mutex>
#include <variant>
#include <filesystem>
#include <atomic> 

namespace TracEon {

struct MMapHandle;

enum class FileFormat : uint8_t {
    DNA_FASTA = 0x00, RNA_FASTA = 0x01, PROTEIN_FASTA = 0x02,
    DNA_FASTQ = 0x03, RNA_FASTQ = 0x04, PROTEIN_FASTQ = 0x05,
    UNKNOWN = 0xFF
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
    SmartStrategy();
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

    FileFormat detected_format_;
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
     * @brief Serialize payload (count + all records) into a buffer.
     * Used by saveBinary() for both v1 and v2 formats.
     */
    void serializePayload(std::vector<char>& buf) const;

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

    // Wrappers
    void parseFastaOptimized(std::string_view content);
    void parseFastqOptimized(std::string_view content);
    void parseFastaMultithreaded(std::string_view content);
    void parseFastqMultithreaded(std::string_view content);
};

} // namespace TracEon

#endif // TRACEON_SMARTSTRATEGY_H