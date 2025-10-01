#include "Cache.h"
#include "SmartStrategy.h"
#include "FileReader.h"
#include <fstream>
#include <sstream>
#include <ios>
#include <variant>

namespace TracEon {

Cache::Cache() {
    // We only have one strategy now!
    m_strategy = std::make_unique<SmartStrategy>();
}

Cache::~Cache() = default;

// --- Data Interaction ---

std::string Cache::get(const std::string& key) {
    // Check if it's in our regular store first
    if (m_store.count(key)) {
        const auto& record_variant = m_store.at(key);
        if (const auto* data = std::get_if<FastaRecordData>(&record_variant)) {
            return m_strategy->decode(*data);
        }
    }

    // Check the SmartStrategy file cache
    if (auto* smart_strategy = dynamic_cast<SmartStrategy*>(m_strategy.get())) {
        return smart_strategy->getSequence(key);
    }

    return "";
}

std::optional<DecodedFastqRecord> Cache::getFastqRecord(const std::string& key) {
    // Check regular store first
    if (m_store.count(key)) {
        const auto& record_variant = m_store.at(key);
        if (const auto* data = std::get_if<FastqRecord>(&record_variant)) {
            return DecodedFastqRecord{
                m_strategy->decode(data->compressed_sequence),
                m_strategy->decode(data->compressed_quality)
            };
        }
    }

    // Check SmartStrategy file cache
    if (auto* smart_strategy = dynamic_cast<SmartStrategy*>(m_strategy.get())) {
        if (smart_strategy->hasSequence(key)) {
            return DecodedFastqRecord{
                smart_strategy->getSequence(key),
                smart_strategy->getQuality(key)
            };
        }
    }

    return std::nullopt;
}

// --- Status & Inspection ---

size_t Cache::size() const {
    size_t total = m_store.size();
    if (auto* smart_strategy = dynamic_cast<SmartStrategy*>(m_strategy.get())) {
        total += smart_strategy->getFileCacheSize();
    }
    return total;
}

size_t Cache::getStoredSize(const std::string& key) const {
    if (m_store.count(key)) {
        const auto& record_variant = m_store.at(key);
        return std::visit([](const auto& value) {
            using T = std::decay_t<decltype(value)>;
            if constexpr (std::is_same_v<T, FastaRecordData>) {
                return value.size();
            } else if constexpr (std::is_same_v<T, FastqRecord>) {
                return value.compressed_sequence.size() + value.compressed_quality.size();
            }
            return size_t(0);
        }, record_variant);
    }
    return 0;
}

void Cache::set(const std::string& key, const std::string& value) {
    m_store[key] = m_strategy->encode(value);
}
// --- File I/O ---

void Cache::loadFile(const std::string& filepath) {
    // Always delegate loading to the SmartStrategy
    if (auto* smart_strategy = dynamic_cast<SmartStrategy*>(m_strategy.get())) {
        smart_strategy->loadFile(filepath);
    }
}

void Cache::loadFastaFromReader(FileReader& file, const std::string& first_line) {
    std::string line = first_line;
    std::string current_key;
    std::stringstream current_value_stream;

    if (!line.empty() && line[0] == '>') {
        size_t first_space = line.find(' ');
        current_key = line.substr(1, first_space != std::string::npos ? first_space - 1 : std::string::npos);
    }

    while (file.getline(line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!current_key.empty()) m_store[current_key] = m_strategy->encode(current_value_stream.str());
            size_t first_space = line.find(' ');
            current_key = line.substr(1, first_space != std::string::npos ? first_space - 1 : std::string::npos);
            current_value_stream.str("");
            current_value_stream.clear();
        } else {
            current_value_stream << line;
        }
    }
    if (!current_key.empty()) m_store[current_key] = m_strategy->encode(current_value_stream.str());
}

void Cache::loadFastqFromReader(FileReader& file, const std::string& first_line) {
    std::string header = first_line, sequence, plus, quality;
    while (true) {
        if (!file.getline(sequence) || !file.getline(plus) || !file.getline(quality)) break;
        if (header.empty() || header[0] != '@') continue;

        size_t first_space = header.find(' ');
        std::string key = header.substr(1, first_space != std::string::npos ? first_space - 1 : std::string::npos);

        FastqRecord record;
        record.compressed_sequence = m_strategy->encode(sequence, DataTypeHint::Generic);
        record.compressed_quality = m_strategy->encode(quality, DataTypeHint::QualityScore);
        m_store[key] = record;

        if(!file.getline(header)) break;
    }
}

void Cache::save(const std::string& filepath) {
    // If m_store has data (from `set` calls), use the old TRAC format.
    if (!m_store.empty()) {
        std::ofstream file(filepath, std::ios::binary);
        if (!file.is_open()) return;
        const char magic[] = "TRAC";
        uint8_t version = 2;
        uint64_t record_count = m_store.size();
        file.write(magic, 4);
        file.write(reinterpret_cast<const char*>(&version), sizeof(version));
        file.write(reinterpret_cast<const char*>(&record_count), sizeof(record_count));

        for (const auto& [key, record_variant] : m_store) {
            uint32_t key_len = key.length();
            file.write(reinterpret_cast<const char*>(&key_len), sizeof(key_len));
            file.write(key.c_str(), key_len);

            if (const auto* fasta_data = std::get_if<FastaRecordData>(&record_variant)) {
                uint8_t record_type = 0; // FASTA
                file.write(reinterpret_cast<const char*>(&record_type), sizeof(record_type));
                uint32_t data_len = fasta_data->size();
                file.write(reinterpret_cast<const char*>(&data_len), sizeof(data_len));
                file.write(reinterpret_cast<const char*>(fasta_data->data()), data_len);
            } else if (const auto* fastq_data = std::get_if<FastqRecord>(&record_variant)) {
                uint8_t record_type = 1; // FASTQ
                file.write(reinterpret_cast<const char*>(&record_type), sizeof(record_type));
                uint32_t seq_len = fastq_data->compressed_sequence.size();
                file.write(reinterpret_cast<const char*>(&seq_len), sizeof(seq_len));
                file.write(reinterpret_cast<const char*>(fastq_data->compressed_sequence.data()), seq_len);
                uint32_t qual_len = fastq_data->compressed_quality.size();
                file.write(reinterpret_cast<const char*>(&qual_len), sizeof(qual_len));
                file.write(reinterpret_cast<const char*>(fastq_data->compressed_quality.data()), qual_len);
            }
        }
    }
    // Otherwise, if the data was loaded from a file, use the SmartStrategy's saver.
    else if (auto* smart_strategy = dynamic_cast<SmartStrategy*>(m_strategy.get())) {
        if (smart_strategy->getFileCacheSize() > 0) {
            smart_strategy->saveBinary(filepath);
        }
    }
}

void Cache::restore(const std::string& filepath) {
    std::ifstream file(filepath, std::ios::binary);
    if (!file.is_open()) throw std::runtime_error("Cannot open cache file: " + filepath);

    char magic_check[4] = {0};
    file.read(magic_check, 4);

    if (file.gcount() < 4) { // File is too small, assume it's the new format
        file.close();
        if (auto* smart_strategy = dynamic_cast<SmartStrategy*>(m_strategy.get())) {
            smart_strategy->loadBinary(filepath);
        }
        return;
    }

    std::string magic_str(magic_check, 4);

    if (magic_str == "TRAC") {
        // --- It's the OLD format, use the manual parser ---
        uint8_t version;
        file.read(reinterpret_cast<char*>(&version), sizeof(version));
        if (version != 2) { file.close(); return; }

        uint64_t record_count;
        file.read(reinterpret_cast<char*>(&record_count), sizeof(record_count));
        m_store.clear();

        for (uint64_t i = 0; i < record_count; ++i) {
            uint32_t key_len;
            file.read(reinterpret_cast<char*>(&key_len), sizeof(key_len));
            std::string key(key_len, '\0');
            file.read(&key[0], key_len);

            uint8_t record_type;
            file.read(reinterpret_cast<char*>(&record_type), sizeof(record_type));

            if (record_type == 0) { // FASTA
                uint32_t data_len;
                file.read(reinterpret_cast<char*>(&data_len), sizeof(data_len));
                FastaRecordData data(data_len);
                file.read(reinterpret_cast<char*>(data.data()), data_len);
                m_store[key] = data;
            } else if (record_type == 1) { // FASTQ
                FastqRecord record;
                uint32_t seq_len;
                file.read(reinterpret_cast<char*>(&seq_len), sizeof(seq_len));
                record.compressed_sequence.resize(seq_len);
                file.read(reinterpret_cast<char*>(record.compressed_sequence.data()), seq_len);
                uint32_t qual_len;
                file.read(reinterpret_cast<char*>(&qual_len), sizeof(qual_len));
                record.compressed_quality.resize(qual_len);
                file.read(reinterpret_cast<char*>(record.compressed_quality.data()), qual_len);
                m_store[key] = record;
            }
        }
    } else {
        // --- Otherwise, assume it's the NEW SmartStrategy format ---
        file.close();
        if (auto* smart_strategy = dynamic_cast<SmartStrategy*>(m_strategy.get())) {
            smart_strategy->loadBinary(filepath);
        }
    }
}

// Add method for direct access to SmartStrategy file operations
void Cache::loadSmartFile(const std::string& filepath) {
    if (auto* smart_strategy = dynamic_cast<SmartStrategy*>(m_strategy.get())) {
        smart_strategy->loadFile(filepath);
    } else {
        throw std::runtime_error("SmartStrategy not available");
    }
}

void Cache::saveSmartBinary(const std::string& filepath) {
    if (auto* smart_strategy = dynamic_cast<SmartStrategy*>(m_strategy.get())) {
        smart_strategy->saveBinary(filepath);
    } else {
        throw std::runtime_error("SmartStrategy not available");
    }
}

void Cache::loadSmartBinary(const std::string& filepath) {
    if (auto* smart_strategy = dynamic_cast<SmartStrategy*>(m_strategy.get())) {
        smart_strategy->loadBinary(filepath);
    } else {
        throw std::runtime_error("SmartStrategy not available");
    }
}

} // namespace TracEon