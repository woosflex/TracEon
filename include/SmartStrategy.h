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
#include <atomic> 

namespace TracEon {

struct MMapHandle;

enum class FileFormat : uint8_t {
    DNA_FASTA = 0x00, RNA_FASTA = 0x01, PROTEIN_FASTA = 0x02,
    DNA_FASTQ = 0x03, RNA_FASTQ = 0x04, PROTEIN_FASTQ = 0x05,
    UNKNOWN = 0xFF
};

/**
 * @struct SequenceView
 * @brief Lightweight pointer triplet for genomic records.
 * * Designed for Zero-Copy architecture.
 * - Key ownership is handled by the HashMap (std::string).
 * - View lifetime is tied to the SmartStrategy's Arena/MMAP.
 */
struct SequenceView {
    std::string_view id;       
    std::string_view sequence; 
    std::string_view quality;  
};

using GenomeIndex = HashMap<std::string, SequenceView>; 
using NGSIndex = HashMap<uint64_t, SequenceView>;       

enum class IndexMode {
    GENOME, 
    NGS     
};

/**
 * @class SmartStrategy
 * @brief Core engine implementing the Hybrid Architecture.
 * * Features:
 * - Arena Allocation for text data
 * - Memory Mapping for binary restoration
 * - Lock-Free atomic read path
 * - Templated Multithreaded Parsers
 */
class SmartStrategy : public IEncodingStrategy {
public:
    SmartStrategy();
    virtual ~SmartStrategy();

    std::vector<unsigned char> encode(const std::string& data, DataTypeHint hint = DataTypeHint::Generic) const override;
    std::string decode(const std::vector<unsigned char>& data) const override;

    void loadFile(const std::string& filepath);
    void loadBinary(const std::string& filepath);
    void saveBinary(const std::string& filepath) const;

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
    
    FileFormat detected_format_;
    mutable std::shared_mutex cache_mutex_;

    /**
     * Lock-Free Read Signal
     * True only when all data is loaded and immutable.
     * Establishes Happens-Before relationship for readers.
     */
    std::atomic<bool> data_ready_{false};

    void determine_format_from_cache();
    bool isNucleotideSequence(std::string_view data) const;
    bool hasRNA(std::string_view data) const;
    uint64_t hash_key(std::string_view key) const; 
    
    // Internal cleanup helper (No Lock)
    void clearInternal();

    // Templated Parsers
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