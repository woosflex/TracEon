#ifndef TRACEON_SMARTSTRATEGY_H
#define TRACEON_SMARTSTRATEGY_H

#include "IEncodingStrategy.h"
#include "MapDefs.h"
#include "RecordTypes.h"

#include <string>
#include <string_view>
#include <vector>
#include <memory>
#include <shared_mutex>
#include <variant>
#include <filesystem>

namespace TracEon {

// Forward declaration
struct MMapHandle;

enum class FileFormat : uint8_t {
    DNA_FASTA = 0x00, RNA_FASTA = 0x01, PROTEIN_FASTA = 0x02,
    DNA_FASTQ = 0x03, RNA_FASTQ = 0x04, PROTEIN_FASTQ = 0x05,
    UNKNOWN = 0xFF
};

// Zero-Copy Value Structure
struct SequenceView {
    // NOTE:
    // - The Map Key (std::string) owns the normalized identifier used for hashing/lookups.
    // - SequenceView.id, sequence, and quality are NON-OWNING views pointing
    //   into either the mmap buffer (binary mode) or the text_arena_ buffer (load mode).
    //   Their lifetime is guaranteed by SmartStrategy and valid until destruction.
    std::string_view id;       
    std::string_view sequence; 
    std::string_view quality;  
};

// --- ADAPTIVE INDEX TYPES ---
// Mode A: Genome/Long-Read (Safe, Owned Keys)
using GenomeIndex = HashMap<std::string, SequenceView>;

// Mode B: NGS/Short-Read (Fast, Hash Keys)
using NGSIndex = HashMap<uint64_t, SequenceView>;

enum class IndexMode {
    GENOME, // Hybrid (std::string key)
    NGS     // Hash-Only (uint64_t key)
};

class SmartStrategy : public IEncodingStrategy {
public:
    SmartStrategy();
    virtual ~SmartStrategy();

    // IEncodingStrategy interface
    std::vector<unsigned char> encode(const std::string& data, DataTypeHint hint = DataTypeHint::Generic) const override;
    std::string decode(const std::vector<unsigned char>& data) const override;

    // --- Core Loading Operations ---
    void loadFile(const std::string& filepath);
    void loadBinary(const std::string& filepath);
    void saveBinary(const std::string& filepath) const;

    // --- Access Methods ---
    std::string_view getView(const std::string_view& sequence_id) const;

    // Legacy Wrappers
    std::string getSequence(const std::string& sequence_id) const;
    std::string getQuality(const std::string& sequence_id) const;
    bool hasSequence(const std::string& sequence_id) const;
    size_t getFileCacheSize() const;
    std::vector<std::string> getAllKeys() const;

    // Utility
    void clearCache();
    FileFormat getDetectedFormat() const { return detected_format_; }
    IndexMode getIndexMode() const; 

private:
    // 1. Stable Memory
    std::vector<char> text_arena_;
    std::unique_ptr<MMapHandle> mmap_handle_;

    // 2. POLYMORPHIC CACHE (The Adaptive Core)
    std::variant<GenomeIndex, NGSIndex> file_cache_;
    
    FileFormat detected_format_;
    mutable std::shared_mutex cache_mutex_;

    // Helpers
    void determine_format_from_cache();
    bool isNucleotideSequence(std::string_view data) const;
    bool hasRNA(std::string_view data) const;
    uint64_t hash_key(std::string_view key) const; 
    
    // Templated Parsers (Internal)
    template <typename MapType>
    void parseFastaInternal(std::string_view content, MapType& map);
    
    template <typename MapType>
    void parseFastqInternal(std::string_view content, MapType& map);

    // Templated Multithreaded Parsers
    template <typename MapType>
    void parseFastaMultithreadedTemplate(std::string_view content, MapType& dest_map);

    template <typename MapType>
    void parseFastqMultithreadedTemplate(std::string_view content, MapType& dest_map);

    // --- Fix: Added missing declarations for wrappers defined in .cpp ---
    void parseFastaOptimized(std::string_view content);
    void parseFastqOptimized(std::string_view content);
    void parseFastaMultithreaded(std::string_view content);
    void parseFastqMultithreaded(std::string_view content);
};

} // namespace TracEon

#endif // TRACEON_SMARTSTRATEGY_H