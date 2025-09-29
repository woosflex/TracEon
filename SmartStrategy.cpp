#include "SmartStrategy.h"
#include <cmath>
#include <cstdint>
#include <vector>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>

// Define our metadata header bytes for self-describing data (your existing system)
namespace {
    const uint8_t TYPE_DNA = 0x01;
    const uint8_t TYPE_RLE_QUALITY = 0x12;
    const uint8_t TYPE_PLAIN_PROTEIN = 0x21;

    // --- Helpers for 2-bit Nucleotide Encoding (your existing code) ---
    inline unsigned char baseToBits(char base) {
        switch (base) {
            case 'A': case 'a': return 0b00;
            case 'C': case 'c': return 0b01;
            case 'G': case 'g': return 0b10;
            case 'T': case 't': case 'U': case 'u': return 0b11;
            default:  return 0b00; // 'N' defaults to 'A'
        }
    }

    inline char bitsToBase(unsigned char bits) {
        switch (bits) {
            case 0b00: return 'A';
            case 0b01: return 'C';
            case 0b10: return 'G';
            case 0b11: return 'T';
            default:   return 'N';
        }
    }
}

SmartStrategy::SmartStrategy() : detected_format_(FileFormat::DNA_FASTA) {}

SmartStrategy::~SmartStrategy() = default;

// --- IEncodingStrategy interface (compatibility with existing Cache) ---

std::vector<unsigned char> SmartStrategy::encode(const std::string& data, DataTypeHint hint) const {
    if (data.empty()) return {};

    // If the caller gives us a hint (e.g., from a FASTQ parser), we trust it.
    if (hint == DataTypeHint::QualityScore) {
        auto encoded = encode_rle(data);
        encoded.insert(encoded.begin(), TYPE_RLE_QUALITY);
        return encoded;
    }

    // Otherwise, perform content sniffing.
    if (isNucleotideSequence(data)) {
        auto encoded = encode_nucleotide(data);
        encoded.insert(encoded.begin(), TYPE_DNA);
        return encoded;
    }

    // If it's not a nucleotide, it's probably protein.
    auto encoded = encode_plain(data);
    encoded.insert(encoded.begin(), TYPE_PLAIN_PROTEIN);
    return encoded;
}

std::string SmartStrategy::decode(const std::vector<unsigned char>& data) const {
    if (data.empty()) return "";

    uint8_t type_id = data[0];
    auto payload = std::vector<unsigned char>(data.begin() + 1, data.end());

    switch (type_id) {
        case TYPE_DNA:
            return decode_nucleotide(payload);
        case TYPE_RLE_QUALITY:
            return decode_rle(payload);
        case TYPE_PLAIN_PROTEIN:
            return decode_plain(payload);
        default:
            return "";
    }
}

// --- Detection Methods ---

bool SmartStrategy::isNucleotideSequence(const std::string& data) const {
    if (data.empty()) return false;

    int nucleotide_count = 0;
    int total_count = 0;

    for (char c : data) {
        if (std::isalpha(c)) {
            total_count++;
            char upper_c = std::toupper(c);
            if (upper_c == 'A' || upper_c == 'T' || upper_c == 'G' ||
                upper_c == 'C' || upper_c == 'U' || upper_c == 'N') {
                nucleotide_count++;
            }
        }
    }

    return total_count > 0 && (static_cast<double>(nucleotide_count) / total_count) > 0.8;
}

bool SmartStrategy::isQualityScore(const std::string& data) const {
    // Quality scores are typically ASCII 33-126
    for (char c : data) {
        if (c < 33 || c > 126) return false;
    }
    return !data.empty();
}

bool SmartStrategy::hasRNA(const std::string& data) const {
    return data.find('U') != std::string::npos || data.find('u') != std::string::npos;
}

// --- Encoding Helpers (your existing implementation) ---

std::vector<unsigned char> SmartStrategy::encode_nucleotide(const std::string& data) const {
    uint32_t original_length = data.length();
    std::vector<uint32_t> n_positions;
    n_positions.reserve(original_length / 10);

    for (uint32_t i = 0; i < original_length; ++i) {
        if (std::toupper(data[i]) == 'N') {
            n_positions.push_back(i);
        }
    }
    uint32_t n_count = n_positions.size();

    size_t packed_data_size = std::ceil(static_cast<double>(original_length) / 4.0);
    size_t n_pos_size = n_count * sizeof(uint32_t);

    std::vector<unsigned char> encoded_data(8 + packed_data_size + n_pos_size, 0);

    // Write header
    encoded_data[0] = (original_length >> 24) & 0xFF;
    encoded_data[1] = (original_length >> 16) & 0xFF;
    encoded_data[2] = (original_length >> 8) & 0xFF;
    encoded_data[3] = original_length & 0xFF;
    encoded_data[4] = (n_count >> 24) & 0xFF;
    encoded_data[5] = (n_count >> 16) & 0xFF;
    encoded_data[6] = (n_count >> 8) & 0xFF;
    encoded_data[7] = n_count & 0xFF;

    // Pack DNA data
    for (size_t i = 0; i < original_length; ++i) {
        unsigned char two_bit_value = baseToBits(data[i]);
        size_t byte_index = 8 + (i / 4);
        int bit_shift = (3 - (i % 4)) * 2;
        encoded_data[byte_index] |= (two_bit_value << bit_shift);
    }

    // Write N positions
    if (n_count > 0) {
        memcpy(&encoded_data[8 + packed_data_size], n_positions.data(), n_pos_size);
    }

    return encoded_data;
}

std::string SmartStrategy::decode_nucleotide(const std::vector<unsigned char>& data) const {
    if (data.size() < 8) return "";

    // Read header
    uint32_t original_length = (static_cast<uint32_t>(data[0]) << 24) |
                              (static_cast<uint32_t>(data[1]) << 16) |
                              (static_cast<uint32_t>(data[2]) << 8)  |
                              static_cast<uint32_t>(data[3]);
    uint32_t n_count = (static_cast<uint32_t>(data[4]) << 24) |
                      (static_cast<uint32_t>(data[5]) << 16) |
                      (static_cast<uint32_t>(data[6]) << 8)  |
                      static_cast<uint32_t>(data[7]);

    std::string decoded_dna;
    decoded_dna.reserve(original_length);
    size_t packed_data_size = std::ceil(static_cast<double>(original_length) / 4.0);

    // Unpack DNA data
    for (size_t i = 0; i < original_length; ++i) {
        size_t byte_index = 8 + (i / 4);
        int bit_shift = (3 - (i % 4)) * 2;
        unsigned char two_bit_value = (data[byte_index] >> bit_shift) & 0b11;
        decoded_dna += bitsToBase(two_bit_value);
    }

    // Reconstruct N's
    if (n_count > 0) {
        std::vector<uint32_t> n_positions(n_count);
        memcpy(n_positions.data(), &data[8 + packed_data_size], n_count * sizeof(uint32_t));
        for (uint32_t pos : n_positions) {
            if (pos < decoded_dna.length()) {
                decoded_dna[pos] = 'N';
            }
        }
    }

    return decoded_dna;
}

std::vector<unsigned char> SmartStrategy::encode_rle(const std::string& data) const {
    if (data.empty()) return {};
    std::vector<unsigned char> encoded_data;
    for (size_t i = 0; i < data.length(); ++i) {
        unsigned char count = 1;
        while (i + 1 < data.length() && data[i + 1] == data[i] && count < 255) {
            count++; i++;
        }
        encoded_data.push_back(count);
        encoded_data.push_back(data[i]);
    }
    return encoded_data;
}

std::string SmartStrategy::decode_rle(const std::vector<unsigned char>& data) const {
    std::string decoded_string;
    for (size_t i = 0; i < data.size(); i += 2) {
        if (i + 1 < data.size()) {
            decoded_string.append(data[i], data[i + 1]);
        }
    }
    return decoded_string;
}

std::vector<unsigned char> SmartStrategy::encode_plain(const std::string& data) const {
    return {data.begin(), data.end()};
}

std::string SmartStrategy::decode_plain(const std::vector<unsigned char>& data) const {
    return {data.begin(), data.end()};
}

// --- File-based Operations (your new requirements) ---

void SmartStrategy::loadFile(const std::string& filepath) {
    file_cache_.clear();

    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filepath);
    }

    std::string first_line;
    std::getline(file, first_line);
    file.seekg(0); // Reset to beginning

    if (first_line.empty()) {
        throw std::runtime_error("Empty file: " + filepath);
    }

    bool is_fastq = (first_line[0] == '@');
    bool is_fasta = (first_line[0] == '>');

    if (!is_fasta && !is_fastq) {
        throw std::runtime_error("Unknown file format: " + filepath);
    }

    // Parse based on format
    if (is_fasta) {
        parseFASTA(file);
    } else {
        parseFASTQ(file);
    }

    // Detect sequence type and set format
    if (!file_cache_.empty()) {
        auto first_seq = file_cache_.begin()->second.sequence;
        bool has_rna = hasRNA(first_seq);
        bool is_nucleotide = isNucleotideSequence(first_seq);

        if (is_fasta) {
            if (has_rna) {
                detected_format_ = FileFormat::RNA_FASTA;
            } else if (is_nucleotide) {
                detected_format_ = FileFormat::DNA_FASTA;
            } else {
                detected_format_ = FileFormat::PROTEIN_FASTA;
            }
        } else {
            if (has_rna) {
                detected_format_ = FileFormat::RNA_FASTQ;
            } else if (is_nucleotide) {
                detected_format_ = FileFormat::DNA_FASTQ;
            } else {
                detected_format_ = FileFormat::PROTEIN_FASTQ;
            }
        }
    }

    std::cout << "Loaded " << file_cache_.size() << " sequences (format: "
              << static_cast<int>(detected_format_) << ")" << std::endl;
}

void SmartStrategy::parseFASTA(std::ifstream& file) {
    std::string line, current_id, current_sequence;

    while (std::getline(file, line)) {
        if (line.empty()) continue;

        // Remove carriage return if present
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }

        if (line[0] == '>') {
            // Save previous sequence
            if (!current_id.empty()) {
                SequenceData seq_data;
                seq_data.id = current_id;
                seq_data.sequence = current_sequence;
                file_cache_[current_id] = seq_data;
            }

            // Start new sequence - extract ID (everything before first space)
            size_t first_space = line.find(' ');
            current_id = line.substr(1, first_space - 1);
            current_sequence.clear();
        } else {
            current_sequence += line;
        }
    }

    // Add last sequence
    if (!current_id.empty()) {
        SequenceData seq_data;
        seq_data.id = current_id;
        seq_data.sequence = current_sequence;
        file_cache_[current_id] = seq_data;
    }
}

void SmartStrategy::parseFASTQ(std::ifstream& file) {
    std::string header, sequence, plus, quality;

    while (std::getline(file, header) &&
           std::getline(file, sequence) &&
           std::getline(file, plus) &&
           std::getline(file, quality)) {

        // Remove carriage returns
        if (!header.empty() && header.back() == '\r') header.pop_back();
        if (!sequence.empty() && sequence.back() == '\r') sequence.pop_back();
        if (!quality.empty() && quality.back() == '\r') quality.pop_back();

        if (header.empty() || header[0] != '@') continue;

        // Extract ID (everything before first space)
        size_t first_space = header.find(' ');
        std::string key = header.substr(1, first_space - 1);

        SequenceData seq_data;
        seq_data.id = key;
        seq_data.sequence = sequence;
        seq_data.quality = quality;
        file_cache_[key] = seq_data;
    }
}

void SmartStrategy::saveBinary(const std::string& binary_filepath) {
    std::ofstream file(binary_filepath, std::ios::binary);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot create binary file: " + binary_filepath);
    }

    // Write format header
    uint8_t format_byte = static_cast<uint8_t>(detected_format_);
    file.write(reinterpret_cast<const char*>(&format_byte), sizeof(format_byte));

    // Write number of sequences
    uint32_t num_sequences = file_cache_.size();
    file.write(reinterpret_cast<const char*>(&num_sequences), sizeof(num_sequences));

    // Write each sequence
    for (const auto& [id, seq_data] : file_cache_) {
        // Write ID
        uint32_t id_len = id.length();
        file.write(reinterpret_cast<const char*>(&id_len), sizeof(id_len));
        file.write(id.c_str(), id_len);

        // Encode and write sequence
        auto encoded_seq = encode(seq_data.sequence);
        uint32_t seq_len = encoded_seq.size();
        file.write(reinterpret_cast<const char*>(&seq_len), sizeof(seq_len));
        file.write(reinterpret_cast<const char*>(encoded_seq.data()), seq_len);

        // Write quality if present
        if (!seq_data.quality.empty()) {
            auto encoded_qual = encode(seq_data.quality, DataTypeHint::QualityScore);
            uint32_t qual_len = encoded_qual.size();
            file.write(reinterpret_cast<const char*>(&qual_len), sizeof(qual_len));
            file.write(reinterpret_cast<const char*>(encoded_qual.data()), qual_len);
        } else {
            uint32_t qual_len = 0;
            file.write(reinterpret_cast<const char*>(&qual_len), sizeof(qual_len));
        }
    }

    std::cout << "Saved " << file_cache_.size() << " sequences to: " << binary_filepath << std::endl;
}

void SmartStrategy::loadBinary(const std::string& binary_filepath) {
    std::ifstream file(binary_filepath, std::ios::binary);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open binary file: " + binary_filepath);
    }

    // Check if this is a "TRAC" format file (wrong format for SmartStrategy)
    char magic[4];
    file.read(magic, 4);
    if (std::string(magic, 4) == "TRAC") {
        throw std::runtime_error("File is in TRAC format, not SmartStrategy format");
    }

    // Rewind to beginning
    file.seekg(0);

    file_cache_.clear();

    // Read format header
    uint8_t format_byte;
    file.read(reinterpret_cast<char*>(&format_byte), sizeof(format_byte));
    detected_format_ = static_cast<FileFormat>(format_byte);

    // Read number of sequences
    uint32_t num_sequences;
    file.read(reinterpret_cast<char*>(&num_sequences), sizeof(num_sequences));

    // Read each sequence
    for (uint32_t i = 0; i < num_sequences; ++i) {
        // Read ID
        uint32_t id_len;
        file.read(reinterpret_cast<char*>(&id_len), sizeof(id_len));
        std::string id(id_len, '\0');
        file.read(&id[0], id_len);

        // Read encoded sequence
        uint32_t seq_len;
        file.read(reinterpret_cast<char*>(&seq_len), sizeof(seq_len));
        std::vector<unsigned char> encoded_seq(seq_len);
        file.read(reinterpret_cast<char*>(encoded_seq.data()), seq_len);

        // Read quality if present
        uint32_t qual_len;
        file.read(reinterpret_cast<char*>(&qual_len), sizeof(qual_len));
        std::string quality;
        if (qual_len > 0) {
            std::vector<unsigned char> encoded_qual(qual_len);
            file.read(reinterpret_cast<char*>(encoded_qual.data()), qual_len);
            quality = decode(encoded_qual);
        }

        // Store decoded data
        SequenceData seq_data;
        seq_data.id = id;
        seq_data.sequence = decode(encoded_seq);
        seq_data.quality = quality;
        file_cache_[id] = seq_data;
    }

    std::cout << "Loaded " << file_cache_.size() << " sequences from: " << binary_filepath << std::endl;
}

// --- File cache access methods ---

std::string SmartStrategy::getSequence(const std::string& sequence_id) const {
    auto it = file_cache_.find(sequence_id);
    return (it != file_cache_.end()) ? it->second.sequence : "";
}

std::string SmartStrategy::getQuality(const std::string& sequence_id) const {
    auto it = file_cache_.find(sequence_id);
    return (it != file_cache_.end()) ? it->second.quality : "";
}

bool SmartStrategy::hasSequence(const std::string& sequence_id) const {
    return file_cache_.find(sequence_id) != file_cache_.end();
}