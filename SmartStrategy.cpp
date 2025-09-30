#include "SmartStrategy.h"
#include "FileReader.h"
#include <cmath>
#include <cstdint>
#include <vector>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <thread>
#include <future>
#include <filesystem>
#include <atomic>

namespace fs = std::filesystem;

namespace TracEon {

// --- Anonymous namespace for private helpers ---
namespace {
    const uint8_t TYPE_DNA = 0x01;
    const uint8_t TYPE_RLE_QUALITY = 0x12;
    const uint8_t TYPE_PLAIN_PROTEIN = 0x21;

    inline unsigned char baseToBits(char base) {
        switch (base) {
            case 'A': case 'a': return 0b00;
            case 'C': case 'c': return 0b01;
            case 'G': case 'g': return 0b10;
            case 'T': case 't': case 'U': case 'u': return 0b11;
            default: return 0b00;
        }
    }

    inline char bitsToBase(unsigned char bits) {
        switch (bits) {
            case 0b00: return 'A';
            case 0b01: return 'C';
            case 0b10: return 'G';
            case 0b11: return 'T';
            default: return 'N';
        }
    }
}


// --- Constructor / Destructor ---
SmartStrategy::SmartStrategy() : detected_format_(FileFormat::DNA_FASTA) {}
SmartStrategy::~SmartStrategy() = default;


// --- Main Public Methods (Multithreaded) ---

void SmartStrategy::loadFile(const std::string& filepath) {
    clearFileCache();

    // --- 1. Initial Setup & File Format Detection ---
    const size_t file_size = fs::file_size(filepath);
    if (file_size == 0) {
        throw std::runtime_error("Input file is empty: " + filepath);
    }
    const unsigned int num_threads = std::max(1u, std::thread::hardware_concurrency());
    const size_t chunk_size = file_size / num_threads;

    // Read the first line to determine the file format (@ for FASTQ, > for FASTA)
    FileReader first_line_reader(filepath);
    std::string first_line;
    if (!first_line_reader.getline(first_line) || first_line.empty()) {
        throw std::runtime_error("Cannot read from file: " + filepath);
    }
    const bool is_fastq = (first_line[0] == '@');
    const bool is_fasta = (first_line[0] == '>');
    if (!is_fasta && !is_fastq) {
        throw std::runtime_error("Unknown file format: " + filepath);
    }

    // --- 2. Pre-calculate Record-Aligned Chunk Boundaries ---

    // This helper lambda scans from an approximate position to find the true start of a record.
    auto find_next_record_start = [&](std::ifstream& file, size_t approx_pos) -> size_t {
        file.seekg(approx_pos);
        if (file.tellg() >= file_size) return file_size;

        // Read ahead in small blocks to find the start of the next record.
        const size_t buffer_size = 4096;
        std::vector<char> buffer(buffer_size);
        const char record_start_char = is_fastq ? '@' : '>';

        while (!file.eof()) {
            std::streampos block_start_pos = file.tellg();
            file.read(buffer.data(), buffer_size);
            std::streamsize bytes_read = file.gcount();

            for (std::streamsize i = 1; i < bytes_read; ++i) {
                // We look for a newline followed by '@' or '>' to signify a new record.
                if (buffer[i] == record_start_char && buffer[i - 1] == '\n') {
                    return block_start_pos + std::streampos(i);
                }
            }
            // If the record start is split across blocks, seek back slightly to catch it.
            if (!file.eof()) file.seekg(block_start_pos + std::streampos(bytes_read - 1));
        }
        return file_size; // If no new record is found, go to the end.
    };

    std::vector<size_t> boundaries;
    boundaries.push_back(0); // The first chunk always starts at the beginning.

    std::ifstream boundary_finder(filepath, std::ios::binary);
    for (unsigned int i = 1; i < num_threads; ++i) {
        boundaries.push_back(find_next_record_start(boundary_finder, i * chunk_size));
    }
    boundaries.push_back(file_size); // The last chunk ends at the file size.
    boundary_finder.close();


    // --- 3. Launch Threads with Precise Boundaries ---
    std::vector<std::future<ParsedChunk>> futures;
    for (unsigned int i = 0; i < num_threads; ++i) {
        futures.push_back(std::async(std::launch::async, [=, this] {
            const size_t start = boundaries[i];
            const size_t end = boundaries[i+1];

            if (start >= end) return ParsedChunk{}; // Skip empty chunks.

            std::ifstream file(filepath, std::ios::binary);
            Chunk buffer(end - start);
            file.seekg(start);
            file.read(buffer.data(), buffer.size());

            // The old, problematic logic of skipping a line is no longer needed.
            // Each chunk is now guaranteed to start with a full record.
            if (is_fasta) return parse_chunk_fasta(buffer);
            else return parse_chunk_fastq(buffer);
        }));
    }

    // --- 4. Merge Results ---
    for (auto& future : futures) {
        ParsedChunk parsed_chunk = future.get();
        std::lock_guard<std::mutex> lock(cache_mutex_);
        for (auto& seq_data : parsed_chunk) {
            file_cache_[seq_data.id] = std::move(seq_data);
        }
    }

    determine_format_from_cache();

    std::cout << "Loaded " << getFileCacheSize() << " sequences (format: "
              << (is_fastq ? "FASTQ" : "FASTA") << ") using " << num_threads << " threads." << std::endl;
}

void SmartStrategy::saveBinary(const std::string& binary_filepath) {
    std::ofstream file(binary_filepath, std::ios::binary);
    if (!file.is_open()) throw std::runtime_error("Cannot create binary file: " + binary_filepath);

    uint8_t format_byte = static_cast<uint8_t>(detected_format_);
    file.write(reinterpret_cast<const char*>(&format_byte), sizeof(format_byte));
    uint32_t num_sequences = getFileCacheSize();
    file.write(reinterpret_cast<const char*>(&num_sequences), sizeof(num_sequences));

    std::vector<SequenceData> records;
    {
        std::lock_guard<std::mutex> lock(cache_mutex_);
        records.reserve(file_cache_.size());
        for(auto const& [key, val] : file_cache_) {
            records.push_back(val);
        }
    }

    const unsigned int num_threads = std::max(1u, std::thread::hardware_concurrency());
    std::vector<std::future<std::vector<char>>> futures;

    auto process_slice = [&](size_t start, size_t end) {
        std::vector<char> buffer;
        for (size_t i = start; i < end; ++i) {
            const auto& seq_data = records[i];

            uint32_t id_len = seq_data.id.length();
            buffer.insert(buffer.end(), reinterpret_cast<const char*>(&id_len), reinterpret_cast<const char*>(&id_len) + sizeof(id_len));
            buffer.insert(buffer.end(), seq_data.id.begin(), seq_data.id.end());

            auto encoded_seq = encode(seq_data.sequence);
            uint32_t seq_len = encoded_seq.size();
            buffer.insert(buffer.end(), reinterpret_cast<const char*>(&seq_len), reinterpret_cast<const char*>(&seq_len) + sizeof(seq_len));
            buffer.insert(buffer.end(), encoded_seq.begin(), encoded_seq.end());

            auto encoded_qual = encode(seq_data.quality, DataTypeHint::QualityScore);
            uint32_t qual_len = encoded_qual.size();
            buffer.insert(buffer.end(), reinterpret_cast<const char*>(&qual_len), reinterpret_cast<const char*>(&qual_len) + sizeof(qual_len));
            if(qual_len > 0) buffer.insert(buffer.end(), encoded_qual.begin(), encoded_qual.end());
        }
        return buffer;
    };

    size_t slice_size = records.size() / num_threads;
    for (unsigned int i = 0; i < num_threads; ++i) {
        size_t start = i * slice_size;
        size_t end = (i == num_threads - 1) ? records.size() : (i + 1) * slice_size;
        if(start < end) {
            futures.push_back(std::async(std::launch::async, process_slice, start, end));
        }
    }

    for (auto& future : futures) {
        auto buffer = future.get();
        file.write(buffer.data(), buffer.size());
    }

    std::cout << "Saved " << getFileCacheSize() << " sequences to: " << binary_filepath << std::endl;
}

void SmartStrategy::loadBinary(const std::string& binary_filepath) {
    std::ifstream file(binary_filepath, std::ios::binary | std::ios::ate);
    if (!file.is_open()) throw std::runtime_error("Cannot open binary file: " + binary_filepath);

    size_t file_size = file.tellg();
    file.seekg(0, std::ios::beg);

    std::vector<char> buffer(file_size);
    file.read(buffer.data(), file_size);
    file.close();

    clearFileCache();

    size_t cursor = 0;
    detected_format_ = static_cast<FileFormat>(buffer[cursor++]);
    uint32_t num_sequences;
    memcpy(&num_sequences, &buffer[cursor], sizeof(uint32_t));
    cursor += sizeof(uint32_t);

    const unsigned int num_threads = std::max(1u, std::thread::hardware_concurrency());
    std::vector<std::future<ParsedChunk>> futures;

    std::atomic<size_t> atomic_cursor(cursor);

    for (unsigned int i = 0; i < num_threads; ++i) {
        futures.push_back(std::async(std::launch::async, [&] {
            ParsedChunk local_chunk;
            while(true) {
                size_t current_pos = atomic_cursor.load();
                if (current_pos >= buffer.size()) break;

                uint32_t id_len;
                memcpy(&id_len, &buffer[current_pos], sizeof(uint32_t));
                std::string id(&buffer[current_pos + sizeof(uint32_t)], id_len);

                uint32_t seq_len;
                memcpy(&seq_len, &buffer[current_pos + sizeof(uint32_t) + id_len], sizeof(uint32_t));
                auto seq_start = reinterpret_cast<const unsigned char*>(&buffer[current_pos + sizeof(uint32_t) * 2 + id_len]);
                std::vector<unsigned char> encoded_seq(seq_start, seq_start + seq_len);

                uint32_t qual_len;
                memcpy(&qual_len, &buffer[current_pos + sizeof(uint32_t) * 2 + id_len + seq_len], sizeof(uint32_t));

                size_t next_pos = current_pos + sizeof(uint32_t) * 3 + id_len + seq_len + qual_len;

                if (!atomic_cursor.compare_exchange_strong(current_pos, next_pos)) continue;

                std::string quality;
                if(qual_len > 0) {
                    auto qual_start = reinterpret_cast<const unsigned char*>(&buffer[current_pos + sizeof(uint32_t) * 3 + id_len + seq_len]);
                    std::vector<unsigned char> encoded_qual(qual_start, qual_start + qual_len);
                    quality = decode(encoded_qual);
                }

                local_chunk.push_back({id, decode(encoded_seq), quality});
            }
            return local_chunk;
        }));
    }

    for (auto& future : futures) {
        ParsedChunk parsed_chunk = future.get();
        std::lock_guard<std::mutex> lock(cache_mutex_);
        for (auto& seq_data : parsed_chunk) {
            file_cache_[seq_data.id] = std::move(seq_data);
        }
    }

    std::cout << "Loaded " << getFileCacheSize() << " sequences from: " << binary_filepath << std::endl;
}

// --- Private Chunk Parsers ---

SmartStrategy::ParsedChunk SmartStrategy::parse_chunk_fasta(const Chunk& buffer) const {
    ParsedChunk records;
    std::string_view content(buffer.data(), buffer.size());
    std::string current_id, current_sequence;

    size_t start_pos = 0;
    while(start_pos < content.size()) {
        size_t end_pos = content.find('\n', start_pos);
        if (end_pos == std::string_view::npos) end_pos = content.size();

        std::string_view line = content.substr(start_pos, end_pos - start_pos);

        if(!line.empty() && line.back() == '\r') line.remove_suffix(1);

        if (!line.empty()) {
            if (line[0] == '>') {
                if (!current_id.empty()) {
                    records.push_back({current_id, current_sequence, ""});
                }
                size_t first_space = line.find(' ');
                current_id = std::string(line.substr(1, first_space != std::string_view::npos ? first_space - 1 : std::string_view::npos));
                current_sequence.clear();
            } else {
                current_sequence.append(line);
            }
        }
        start_pos = end_pos + 1;
    }

    if (!current_id.empty()) {
        records.push_back({current_id, current_sequence, ""});
    }
    return records;
}

SmartStrategy::ParsedChunk SmartStrategy::parse_chunk_fastq(const Chunk& buffer) const {
    ParsedChunk records;
    std::string_view content(buffer.data(), buffer.size());

    size_t start_pos = 0;
    while(start_pos < content.size()) {
        size_t p1 = content.find('\n', start_pos);
        if(p1 == std::string_view::npos) break;
        size_t p2 = content.find('\n', p1 + 1);
        if(p2 == std::string_view::npos) break;
        size_t p3 = content.find('\n', p2 + 1);
        if(p3 == std::string_view::npos) break;
        size_t p4 = content.find('\n', p3 + 1);
        if(p4 == std::string_view::npos) p4 = content.size();

        std::string_view header = content.substr(start_pos, p1-start_pos);
        std::string_view sequence = content.substr(p1+1, p2-(p1+1));
        std::string_view quality = content.substr(p3+1, p4-(p3+1));

        // Trim trailing newline and carriage return characters from the views
        auto trim_view = [](std::string_view sv) {
            while (!sv.empty() && (sv.back() == '\n' || sv.back() == '\r')) {
                sv.remove_suffix(1);
            }
            return sv;
        };

        header = trim_view(header);
        sequence = trim_view(sequence);
        quality = trim_view(quality);

        if (!header.empty() && header[0] == '@') {
            size_t first_space = header.find(' ');
            std::string key(header.substr(1, first_space != std::string_view::npos ? first_space - 1 : std::string_view::npos));
            records.push_back({key, std::string(sequence), std::string(quality)});
        }
        start_pos = p4 + 1;
    }
    return records;
}

void SmartStrategy::determine_format_from_cache() {
    std::lock_guard<std::mutex> lock(cache_mutex_);
    if (file_cache_.empty()) return;

    const auto& first_seq_data = file_cache_.begin()->second;
    const auto& first_seq = first_seq_data.sequence;
    bool is_rna = hasRNA(first_seq);
    bool is_nuc = isNucleotideSequence(first_seq);
    bool is_fastq = !first_seq_data.quality.empty();

    if (!is_fastq) {
        detected_format_ = is_rna ? FileFormat::RNA_FASTA :
                          is_nuc ? FileFormat::DNA_FASTA :
                          FileFormat::PROTEIN_FASTA;
    } else {
        detected_format_ = is_rna ? FileFormat::RNA_FASTQ :
                          is_nuc ? FileFormat::DNA_FASTQ :
                          FileFormat::PROTEIN_FASTQ;
    }
}

// --- Other Public and Private Methods ---
    void SmartStrategy::mergeFileCacheInto(
     std::unordered_map<std::string, std::variant<std::vector<unsigned char>, FastqRecord>>& store,
     const IEncodingStrategy& encoder)
{
    // This function can now access member variables
    std::lock_guard<std::mutex> lock(cache_mutex_);

    for (const auto& [key, data] : file_cache_) {
        if (!data.quality.empty()) { // It's a FASTQ record
            FastqRecord record;
            record.compressed_sequence = encoder.encode(data.sequence, DataTypeHint::Generic);
            record.compressed_quality = encoder.encode(data.quality, DataTypeHint::QualityScore);
            store[key] = record;
        } else { // It's a FASTA record
            store[key] = encoder.encode(data.sequence);
        }
    }
}

size_t SmartStrategy::getFileCacheSize() const {
    std::lock_guard<std::mutex> lock(cache_mutex_);
    return file_cache_.size();
}

void SmartStrategy::clearFileCache() {
    std::lock_guard<std::mutex> lock(cache_mutex_);
    file_cache_.clear();
}

std::string SmartStrategy::getSequence(const std::string& sequence_id) const {
    std::lock_guard<std::mutex> lock(cache_mutex_);
    auto it = file_cache_.find(sequence_id);
    return (it != file_cache_.end()) ? it->second.sequence : "";
}

std::string SmartStrategy::getQuality(const std::string& sequence_id) const {
    std::lock_guard<std::mutex> lock(cache_mutex_);
    auto it = file_cache_.find(sequence_id);
    return (it != file_cache_.end()) ? it->second.quality : "";
}

bool SmartStrategy::hasSequence(const std::string& sequence_id) const {
    std::lock_guard<std::mutex> lock(cache_mutex_);
    return file_cache_.find(sequence_id) != file_cache_.end();
}

std::vector<unsigned char> SmartStrategy::encode(const std::string& data, DataTypeHint hint) const {
    if (data.empty()) return {};
    if (hint == DataTypeHint::QualityScore) {
        auto encoded = encode_rle(data);
        encoded.insert(encoded.begin(), TYPE_RLE_QUALITY);
        return encoded;
    }
    if (isNucleotideSequence(data)) {
        auto encoded = encode_nucleotide(data);
        encoded.insert(encoded.begin(), TYPE_DNA);
        return encoded;
    }
    auto encoded = encode_plain(data);
    encoded.insert(encoded.begin(), TYPE_PLAIN_PROTEIN);
    return encoded;
}

std::string SmartStrategy::decode(const std::vector<unsigned char>& data) const {
    if (data.empty()) return "";
    uint8_t type_id = data[0];
    auto payload = std::vector<unsigned char>(data.begin() + 1, data.end());
    switch (type_id) {
        case TYPE_DNA: return decode_nucleotide(payload);
        case TYPE_RLE_QUALITY: return decode_rle(payload);
        case TYPE_PLAIN_PROTEIN: return decode_plain(payload);
        default: return "";
    }
}

bool SmartStrategy::isNucleotideSequence(const std::string& data) const {
    if (data.empty()) return false;
    int nucleotide_count = 0;
    int total_count = 0;
    for (char c : data) {
        if (std::isalpha(c)) {
            total_count++;
            char upper_c = std::toupper(c);
            if (upper_c == 'A' || upper_c == 'T' || upper_c == 'G' || upper_c == 'C' || upper_c == 'U' || upper_c == 'N') {
                nucleotide_count++;
            }
        }
    }
    return total_count > 0 && (static_cast<double>(nucleotide_count) / total_count) > 0.8;
}

bool SmartStrategy::hasRNA(const std::string& data) const {
    return data.find('U') != std::string::npos || data.find('u') != std::string::npos;
}

std::vector<unsigned char> SmartStrategy::encode_nucleotide(const std::string& data) const {
    uint32_t original_length = data.length();
    std::vector<uint32_t> n_positions;
    if(original_length > 0) n_positions.reserve(original_length / 10);
    for (uint32_t i = 0; i < original_length; ++i) {
        if (std::toupper(data[i]) == 'N') n_positions.push_back(i);
    }
    uint32_t n_count = n_positions.size();
    size_t packed_data_size = (original_length + 3) / 4;
    size_t n_pos_size = n_count * sizeof(uint32_t);
    std::vector<unsigned char> encoded_data(8 + packed_data_size + n_pos_size, 0);
    encoded_data[0] = (original_length >> 24) & 0xFF;
    encoded_data[1] = (original_length >> 16) & 0xFF;
    encoded_data[2] = (original_length >> 8) & 0xFF;
    encoded_data[3] = original_length & 0xFF;
    encoded_data[4] = (n_count >> 24) & 0xFF;
    encoded_data[5] = (n_count >> 16) & 0xFF;
    encoded_data[6] = (n_count >> 8) & 0xFF;
    encoded_data[7] = n_count & 0xFF;
    for (size_t i = 0; i < original_length; ++i) {
        encoded_data[8 + (i / 4)] |= (baseToBits(data[i]) << ((3 - (i % 4)) * 2));
    }
    if (n_count > 0) {
        memcpy(&encoded_data[8 + packed_data_size], n_positions.data(), n_pos_size);
    }
    return encoded_data;
}

std::string SmartStrategy::decode_nucleotide(const std::vector<unsigned char>& data) const {
    if (data.size() < 8) return "";
    uint32_t original_length = (uint32_t(data[0]) << 24) | (uint32_t(data[1]) << 16) | (uint32_t(data[2]) << 8) | uint32_t(data[3]);
    uint32_t n_count = (uint32_t(data[4]) << 24) | (uint32_t(data[5]) << 16) | (uint32_t(data[6]) << 8) | uint32_t(data[7]);
    std::string decoded_dna;
    decoded_dna.reserve(original_length);
    size_t packed_data_size = (original_length + 3) / 4;
    for (size_t i = 0; i < original_length; ++i) {
        decoded_dna += bitsToBase((data[8 + (i / 4)] >> ((3 - (i % 4)) * 2)) & 0b11);
    }
    if (n_count > 0) {
        std::vector<uint32_t> n_positions(n_count);
        memcpy(n_positions.data(), &data[8 + packed_data_size], n_count * sizeof(uint32_t));
        for (uint32_t pos : n_positions) {
            if (pos < decoded_dna.length()) decoded_dna[pos] = 'N';
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
    decoded_string.reserve(data.size()); // Heuristic
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

} // namespace TracEon

