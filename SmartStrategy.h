#ifndef TRACEON_SMARTSTRATEGY_H
#define TRACEON_SMARTSTRATEGY_H

#include "IEncodingStrategy.h"
#include <string>
#include <vector>
#include <unordered_map>
#include <memory>
#include <mutex>
#include <variant>

#include "RecordTypes.h"

namespace TracEon {

// Forward-declare FileReader to avoid including its full header here
class FileReader;

// Enum for identifying the type of biological data in the file
enum class FileFormat : uint8_t {
    DNA_FASTA = 0x00, RNA_FASTA = 0x01, PROTEIN_FASTA = 0x02,
    DNA_FASTQ = 0x03, RNA_FASTQ = 0x04, PROTEIN_FASTQ = 0x05
};

class SmartStrategy : public IEncodingStrategy {
public:
    SmartStrategy();
    virtual ~SmartStrategy();

    // IEncodingStrategy interface
    std::vector<unsigned char> encode(const std::string& data, DataTypeHint hint = DataTypeHint::Generic) const override;
    std::string decode(const std::vector<unsigned char>& data) const override;

    // Main file-based operations (now multithreaded)
    void loadFile(const std::string& filepath);
    void saveBinary(const std::string& binary_filepath);
    void loadBinary(const std::string& binary_filepath);

    // File cache access methods (thread-safe)
    std::string getSequence(const std::string& sequence_id) const;
    std::string getQuality(const std::string& sequence_id) const;
    bool hasSequence(const std::string& sequence_id) const;
    size_t getFileCacheSize() const;

    void mergeFileCacheInto(
        std::unordered_map<std::string, std::variant<std::vector<unsigned char>, FastqRecord>>& store,
        const IEncodingStrategy& encoder
    );
    // Utility methods
    FileFormat getDetectedFormat() const { return detected_format_; }
    void clearFileCache();

private:
    // Internal struct to hold sequence data
    struct SequenceData {
        std::string id;
        std::string sequence;
        std::string quality;
    };

    // Type aliases for multithreading logic
    using SequenceCache = std::unordered_map<std::string, SequenceData>;
    using Chunk = std::vector<char>;
    using ParsedChunk = std::vector<SequenceData>;

    // Member variables
    SequenceCache file_cache_;
    FileFormat detected_format_;
    mutable std::mutex cache_mutex_; // Protects access to file_cache_

    // Private helper methods for multithreading
    ParsedChunk parse_chunk_fasta(const Chunk& buffer) const;
    ParsedChunk parse_chunk_fastq(const Chunk& buffer) const;
    void determine_format_from_cache();

    // Private helper methods for content analysis
    bool isNucleotideSequence(const std::string& data) const;
    bool hasRNA(const std::string& data) const;

    // Private helper methods for encoding/decoding
    std::vector<unsigned char> encode_nucleotide(const std::string& data) const;
    std::string decode_nucleotide(const std::vector<unsigned char>& data) const;
    std::vector<unsigned char> encode_rle(const std::string& data) const;
    std::string decode_rle(const std::vector<unsigned char>& data) const;
    std::vector<unsigned char> encode_plain(const std::string& data) const;
    std::string decode_plain(const std::vector<unsigned char>& data) const;
};

} // namespace TracEon

#endif // TRACEON_SMARTSTRATEGY_H

