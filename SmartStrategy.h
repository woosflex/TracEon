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

// Forward declaration for the PIMPL-style handling of OS-specific MMAP resources
struct MMapHandle;

// Enum for identifying the type of biological data
enum class FileFormat : uint8_t {
    DNA_FASTA = 0x00, RNA_FASTA = 0x01, PROTEIN_FASTA = 0x02,
    DNA_FASTQ = 0x03, RNA_FASTQ = 0x04, PROTEIN_FASTQ = 0x05,
    UNKNOWN = 0xFF
};

// Zero-Copy Data Structure
struct SequenceView {
    std::string_view id;
    std::string_view sequence;
    std::string_view quality; // Empty for FASTA
};

class SmartStrategy : public IEncodingStrategy {
public:
    SmartStrategy();
    virtual ~SmartStrategy();

    // IEncodingStrategy interface (Legacy wrappers)
    std::vector<unsigned char> encode(const std::string& data, DataTypeHint hint = DataTypeHint::Generic) const override;
    std::string decode(const std::vector<unsigned char>& data) const override;

    // --- Core Loading Operations ---

    // Parses text file into the Arena (In-Memory)
    void loadFile(const std::string& filepath);

    // Maps a binary file directly into memory (Zero-Copy)
    void loadBinary(const std::string& filepath);

    // Serializes current cache to the custom binary format
    void saveBinary(const std::string& filepath) const;

    // --- Access Methods ---

    // New Zero-Copy Accessor
    std::string_view getView(const std::string_view& sequence_id) const;

    // Legacy String Accessors (Constructs a string from the view)
    std::string getSequence(const std::string& sequence_id) const;
    std::string getQuality(const std::string& sequence_id) const;
    bool hasSequence(const std::string& sequence_id) const;
    size_t getFileCacheSize() const;
    std::vector<std::string> getAllKeys() const;

    // Utility
    void clearCache();
    FileFormat getDetectedFormat() const { return detected_format_; }

private:
    // --- Data Storage ---

    // 1. Arena for text-based loading. Stores the full file content contiguously.
    std::vector<char> text_arena_;

    // 2. Handle for Memory Mapped files (OS specific resource).
    std::unique_ptr<MMapHandle> mmap_handle_;

    // 3. CRITICAL FIX: Use std::string keys (owned) with SequenceView values (views)
    // This enables hash caching and fast lookups while maintaining zero-copy for values
    HashMap<std::string, SequenceView> file_cache_;

    FileFormat detected_format_;
    
    // Use shared_mutex for concurrent reads
    mutable std::shared_mutex cache_mutex_;

    // Helpers
    void determine_format_from_cache();
    bool isNucleotideSequence(std::string_view data) const;
    bool hasRNA(std::string_view data) const;
    
    // Optimized parsers
    void parseFastaOptimized(std::string_view content);
    void parseFastqOptimized(std::string_view content);
};

} // namespace TracEon

#endif // TRACEON_SMARTSTRATEGY_H