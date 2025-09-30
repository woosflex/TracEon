#ifndef TRACEON_SMARTSTRATEGY_H
#define TRACEON_SMARTSTRATEGY_H

#include "IEncodingStrategy.h"
#include <string>
#include <vector>
#include <unordered_map>
#include <memory>
#include <fstream>
namespace TracEon {

class FileReader;

// File format identifiers for binary headers (your original spec)
enum class FileFormat : uint8_t {
    DNA_FASTA = 0x00,    // 0x000
    RNA_FASTA = 0x01,    // 0x001
    PROTEIN_FASTA = 0x02, // 0x010
    DNA_FASTQ = 0x03,    // 0x011
    RNA_FASTQ = 0x04,    // 0x100
    PROTEIN_FASTQ = 0x05  // 0x101
};

class SmartStrategy : public IEncodingStrategy {
private:
    // Internal cache for file-based operations
    struct SequenceData {
        std::string id;
        std::string sequence;
        std::string quality;
    };

    std::unordered_map<std::string, SequenceData> file_cache_;
    FileFormat detected_format_;

    // Detection methods
    bool isNucleotideSequence(const std::string& data) const;
    bool isQualityScore(const std::string& data) const;
    bool hasRNA(const std::string& data) const;

    // File parsing methods (missing declarations)
    void parseFASTA(FileReader& file, const std::string& first_line);
    void parseFASTQ(FileReader& file, const std::string& first_line);

    // Encoding helpers (your existing implementation)
    std::vector<unsigned char> encode_nucleotide(const std::string& data) const;
    std::string decode_nucleotide(const std::vector<unsigned char>& data) const;
    std::vector<unsigned char> encode_rle(const std::string& data) const;
    std::string decode_rle(const std::vector<unsigned char>& data) const;
    std::vector<unsigned char> encode_plain(const std::string& data) const;
    std::string decode_plain(const std::vector<unsigned char>& data) const;

public:
    SmartStrategy();
    virtual ~SmartStrategy();

    // IEncodingStrategy interface (compatibility with existing Cache)
    std::vector<unsigned char> encode(const std::string& data, DataTypeHint hint = DataTypeHint::Generic) const override;
    std::string decode(const std::vector<unsigned char>& data) const override;

    // File-based operations (your new requirements)
    void loadFile(const std::string& filepath);
    void saveBinary(const std::string& binary_filepath);
    void loadBinary(const std::string& binary_filepath);

    // File cache access
    std::string getSequence(const std::string& sequence_id) const;
    std::string getQuality(const std::string& sequence_id) const;
    bool hasSequence(const std::string& sequence_id) const;
    size_t getFileCacheSize() const { return file_cache_.size(); }

    // Utility
    FileFormat getDetectedFormat() const { return detected_format_; }
    void clearFileCache() { file_cache_.clear(); }
};
} // namespace TracEon
#endif //TRACEON_SMARTSTRATEGY_H