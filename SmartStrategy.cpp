#include "SmartStrategy.h"
#include "FileReader.h"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <thread>
#include <future>
#include <filesystem>

namespace fs = std::filesystem;
namespace TracEon {

SmartStrategy::SmartStrategy() : detected_format_(FileFormat::DNA_FASTA) {}
SmartStrategy::~SmartStrategy() = default;

// Simple pass-through - no compression
std::vector<unsigned char> SmartStrategy::encode(const std::string& data, DataTypeHint hint) const {
    return {data.begin(), data.end()};
}

std::string SmartStrategy::decode(const std::vector<unsigned char>& data) const {
    return {data.begin(), data.end()};
}

void SmartStrategy::loadFile(const std::string& filepath) {
    clearFileCache();

    const size_t file_size = fs::file_size(filepath);
    if (file_size == 0) return;

    const bool is_gzipped = (filepath.size() > 3 && filepath.substr(filepath.size() - 3) == ".gz");
    char record_start_char;

    {
        FileReader first_line_reader(filepath);
        std::string first_line;
        if (!first_line_reader.getline(first_line) || first_line.empty()) {
            throw std::runtime_error("Cannot read from file: " + filepath);
        }
        record_start_char = first_line[0];
    }
    const bool is_fastq = (record_start_char == '@');

    // Use single-threaded for small files (< 1MB) or gzipped files
    const bool use_single_threaded = (file_size < 1024 * 1024) || is_gzipped;
    const unsigned int num_threads = use_single_threaded ? 1 : std::max(1u, std::thread::hardware_concurrency());

    if (use_single_threaded) {
        // Single-threaded parsing
        FileReader reader(filepath);
        std::string line;

        if (is_fastq) {
            while (reader.getline(line)) {
                if (line.empty() || line[0] != '@') continue;
                std::string s, p, q;
                if (!reader.getline(s) || !reader.getline(p) || !reader.getline(q)) break;

                size_t first_space = line.find(' ');
                std::string id = line.substr(1, first_space != std::string::npos ? first_space - 1 : std::string::npos);
                file_cache_[id] = {id, s, q};
            }
        } else {
            std::string id, seq_buffer;
            while(reader.getline(line)) {
                if (line.empty()) continue;
                if (line[0] == '>') {
                    if (!id.empty()) file_cache_[id] = {id, seq_buffer, ""};
                    seq_buffer.clear();
                    size_t first_space = line.find(' ');
                    id = line.substr(1, first_space != std::string::npos ? first_space - 1 : std::string::npos);
                } else {
                    seq_buffer += line;
                }
            }
            if (!id.empty()) file_cache_[id] = {id, seq_buffer, ""};
        }
    } else {
        // Multithreaded parsing with pre-scan
        std::vector<size_t> record_starts;

        // Phase 1: Pre-scan to build index of verified record start positions
        {
            std::ifstream scan_file(filepath, std::ios::binary);
            std::string line;

            if (is_fastq) {
                // For FASTQ: validate 4-line structure before adding to index
                while (std::getline(scan_file, line)) {
                    size_t header_pos = scan_file.tellg();
                    header_pos = header_pos - line.length() - 1;

                    if (!line.empty() && line[0] == '@') {
                        // Verify this is a true header by checking next 3 lines exist
                        std::string seq, plus, qual;
                        size_t saved_pos = scan_file.tellg();

                        if (std::getline(scan_file, seq) &&
                            std::getline(scan_file, plus) &&
                            std::getline(scan_file, qual) &&
                            !plus.empty() && plus[0] == '+') {
                            // Valid 4-line FASTQ record
                            record_starts.push_back(header_pos);
                        } else {
                            // Not a valid record, rewind and continue
                            scan_file.seekg(saved_pos);
                        }
                    }
                }
            } else {
                // For FASTA: simpler - just look for '>'
                while (std::getline(scan_file, line)) {
                    if (!line.empty() && line[0] == '>') {
                        size_t pos = scan_file.tellg();
                        record_starts.push_back(pos - line.length() - 1);
                    }
                }
            }
        }

        if (record_starts.empty()) {
            determine_format_from_cache();
            return;
        }

        // Phase 2: Divide records among threads
        size_t total_records = record_starts.size();
        size_t records_per_thread = (total_records + num_threads - 1) / num_threads;

        std::vector<std::future<std::vector<SequenceData>>> futures;

        for (unsigned int i = 0; i < num_threads; ++i) {
            size_t start_idx = i * records_per_thread;
            if (start_idx >= total_records) break;

            size_t end_idx = std::min(start_idx + records_per_thread, total_records);

            futures.push_back(std::async(std::launch::async, [=, this]() {
                std::vector<SequenceData> local_results;
                std::ifstream thread_file(filepath, std::ios::binary);

                for (size_t idx = start_idx; idx < end_idx; ++idx) {
                    thread_file.seekg(record_starts[idx]);

                    if (is_fastq) {
                        std::string header, sequence, plus, quality;

                        if (!std::getline(thread_file, header)) continue;
                        if (!std::getline(thread_file, sequence)) continue;
                        if (!std::getline(thread_file, plus)) continue;
                        if (!std::getline(thread_file, quality)) continue;

                        if (header.empty() || header[0] != '@') continue;

                        size_t first_space = header.find(' ');
                        std::string id = header.substr(1, first_space != std::string::npos ? first_space - 1 : std::string::npos);

                        local_results.push_back({id, sequence, quality});
                    } else {
                        // FASTA
                        std::string header, seq_buffer, line;

                        if (!std::getline(thread_file, header)) continue;

                        size_t first_space = header.find(' ');
                        std::string id = header.substr(1, first_space != std::string::npos ? first_space - 1 : std::string::npos);

                        // Read sequence lines until next record or end
                        while (std::getline(thread_file, line)) {
                            if (!line.empty() && line[0] == '>') break;
                            seq_buffer += line;
                        }

                        local_results.push_back({id, seq_buffer, ""});
                    }
                }

                return local_results;
            }));
        }

        // Phase 3: Collect results
        for (auto& future : futures) {
            auto chunk_results = future.get();
            for (const auto& data : chunk_results) {
                file_cache_[data.id] = data;
            }
        }
    }

    determine_format_from_cache();
    std::cout << "Loaded " << getFileCacheSize() << " sequences (format: "
              << (is_fastq ? "FASTQ" : "FASTA") << ") using "
              << num_threads << " parsing thread(s)." << std::endl;
}

void SmartStrategy::saveBinary(const std::string& binary_filepath) {
    std::ofstream file(binary_filepath, std::ios::binary);
    if (!file.is_open()) throw std::runtime_error("Cannot create binary file: " + binary_filepath);

    std::lock_guard<std::mutex> lock(cache_mutex_);

    // Write magic number for SmartStrategy format
    const char magic[] = "SMRT";
    file.write(magic, 4);

    uint8_t format_byte = static_cast<uint8_t>(detected_format_);
    file.write(reinterpret_cast<const char*>(&format_byte), sizeof(format_byte));

    uint64_t num_sequences = file_cache_.size();
    file.write(reinterpret_cast<const char*>(&num_sequences), sizeof(num_sequences));

    for (const auto& [key, data] : file_cache_) {
        uint32_t id_len = data.id.length();
        file.write(reinterpret_cast<const char*>(&id_len), sizeof(id_len));
        file.write(data.id.c_str(), id_len);

        uint32_t seq_len = data.sequence.length();
        file.write(reinterpret_cast<const char*>(&seq_len), sizeof(seq_len));
        file.write(data.sequence.c_str(), seq_len);

        uint32_t qual_len = data.quality.length();
        file.write(reinterpret_cast<const char*>(&qual_len), sizeof(qual_len));
        if (qual_len > 0) {
            file.write(data.quality.c_str(), qual_len);
        }
    }

    std::cout << "Saved " << file_cache_.size() << " sequences to: " << binary_filepath << std::endl;
}

void SmartStrategy::loadBinary(const std::string& binary_filepath) {
    std::ifstream file(binary_filepath, std::ios::binary);
    if (!file.is_open()) throw std::runtime_error("Cannot open binary file: " + binary_filepath);

    // Check magic number
    char magic[4];
    file.read(magic, 4);
    if (std::string(magic, 4) != "SMRT") {
        throw std::runtime_error("Invalid SmartStrategy binary file format");
    }

    clearFileCache();

    uint8_t format_byte;
    file.read(reinterpret_cast<char*>(&format_byte), sizeof(format_byte));
    detected_format_ = static_cast<FileFormat>(format_byte);

    uint64_t num_sequences;
    file.read(reinterpret_cast<char*>(&num_sequences), sizeof(num_sequences));

    for (uint64_t i = 0; i < num_sequences; ++i) {
        SequenceData data;
        uint32_t len;

        file.read(reinterpret_cast<char*>(&len), sizeof(len));
        data.id.resize(len);
        file.read(&data.id[0], len);

        file.read(reinterpret_cast<char*>(&len), sizeof(len));
        data.sequence.resize(len);
        file.read(&data.sequence[0], len);

        file.read(reinterpret_cast<char*>(&len), sizeof(len));
        if (len > 0) {
            data.quality.resize(len);
            file.read(&data.quality[0], len);
        }

        file_cache_[data.id] = std::move(data);
    }

    std::cout << "Loaded " << file_cache_.size() << " sequences from: " << binary_filepath << std::endl;
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

bool SmartStrategy::isNucleotideSequence(const std::string& data) const {
    if (data.empty()) return false;
    int nucleotide_count = 0, total_count = 0;
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

} // namespace TracEon