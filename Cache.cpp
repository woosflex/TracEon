#include "Cache.h"
#include "SmartStrategy.h"
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

// --- File I/O ---

void Cache::loadFile(const std::string& filepath) {
    // Use SmartStrategy's file loading capability
    if (auto* smart_strategy = dynamic_cast<SmartStrategy*>(m_strategy.get())) {
        smart_strategy->loadFile(filepath);
        return;
    }

    // Fallback to old method if not SmartStrategy
    std::ifstream file(filepath);
    if (!file.is_open()) return;

    // Read the first line to sniff the file type.
    std::string first_line;
    std::getline(file, first_line);

    // Rewind the file to the beginning so the parsers can read the whole thing.
    file.seekg(0);

    // Sniff the file type based on the first character.
    if (!first_line.empty()) {
        if (first_line[0] == '>') {
            loadFasta(file);
        } else if (first_line[0] == '@') {
            loadFastq(file);
        }
    }
}

void Cache::loadFasta(std::ifstream& file) {
    std::string line, current_key;
    std::stringstream current_value_stream;

    while (std::getline(file, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;

        if (line[0] == '>') {
            if (!current_key.empty()) {
                m_store[current_key] = m_strategy->encode(current_value_stream.str());
            }
            size_t first_space = line.find(' ');
            current_key = line.substr(1, first_space - 1);
            current_value_stream.str("");
            current_value_stream.clear();
        } else {
            current_value_stream << line;
        }
    }
    if (!current_key.empty()) {
        m_store[current_key] = m_strategy->encode(current_value_stream.str());
    }
}

void Cache::loadFastq(std::ifstream& file) {
    std::string header, sequence, plus, quality;

    while (std::getline(file, header) &&
           std::getline(file, sequence) &&
           std::getline(file, plus) &&
           std::getline(file, quality))
    {
        if (!header.empty() && header.back() == '\r') header.pop_back();
        if (!sequence.empty() && sequence.back() == '\r') sequence.pop_back();
        if (!quality.empty() && quality.back() == '\r') quality.pop_back();

        if (header.empty() || header[0] != '@') continue;

        size_t first_space = header.find(' ');
        std::string key = header.substr(1, first_space - 1);

        FastqRecord record;
        record.compressed_sequence = m_strategy->encode(sequence, DataTypeHint::Generic);
        record.compressed_quality = m_strategy->encode(quality, DataTypeHint::QualityScore);
        m_store[key] = record;
    }
}

void Cache::save(const std::string& filepath) {
    // Check if we can use SmartStrategy's binary save
    if (auto* smart_strategy = dynamic_cast<SmartStrategy*>(m_strategy.get())) {
        if (smart_strategy->getFileCacheSize() > 0) {
            smart_strategy->saveBinary(filepath);
            return;
        }
    }

    // Fallback to original save format
    std::ofstream file(filepath, std::ios::binary);
    if (!file.is_open()) return;

    // Write Header
    const char magic[] = "TRAC";
    uint8_t version = 2; // New version for the SmartStrategy format
    uint64_t record_count = m_store.size();
    file.write(magic, 4);
    file.write(reinterpret_cast<const char*>(&version), sizeof(version));
    file.write(reinterpret_cast<const char*>(&record_count), sizeof(record_count));

    // Write Records (simplified for now)
    for (const auto& [key, record_variant] : m_store) {
        // Write key
        uint32_t key_len = key.length();
        file.write(reinterpret_cast<const char*>(&key_len), sizeof(key_len));
        file.write(key.c_str(), key_len);

        // Write record type and data
        if (const auto* fasta_data = std::get_if<FastaRecordData>(&record_variant)) {
            uint8_t record_type = 0; // FASTA
            file.write(reinterpret_cast<const char*>(&record_type), sizeof(record_type));

            uint32_t data_len = fasta_data->size();
            file.write(reinterpret_cast<const char*>(&data_len), sizeof(data_len));
            file.write(reinterpret_cast<const char*>(fasta_data->data()), data_len);
        } else if (const auto* fastq_data = std::get_if<FastqRecord>(&record_variant)) {
            uint8_t record_type = 1; // FASTQ
            file.write(reinterpret_cast<const char*>(&record_type), sizeof(record_type));

            // Write sequence
            uint32_t seq_len = fastq_data->compressed_sequence.size();
            file.write(reinterpret_cast<const char*>(&seq_len), sizeof(seq_len));
            file.write(reinterpret_cast<const char*>(fastq_data->compressed_sequence.data()), seq_len);

            // Write quality
            uint32_t qual_len = fastq_data->compressed_quality.size();
            file.write(reinterpret_cast<const char*>(&qual_len), sizeof(qual_len));
            file.write(reinterpret_cast<const char*>(fastq_data->compressed_quality.data()), qual_len);
        }
    }
}

void Cache::restore(const std::string& filepath) {
    // Try SmartStrategy binary load first
    if (auto* smart_strategy = dynamic_cast<SmartStrategy*>(m_strategy.get())) {
        try {
            smart_strategy->loadBinary(filepath);
            return;
        } catch (const std::exception&) {
            // If SmartStrategy load fails, fall back to original format
        }
    }

    // Fallback to original restore format
    std::ifstream file(filepath, std::ios::binary);
    if (!file.is_open()) return;

    // Read and verify header
    char magic[4];
    file.read(magic, 4);
    if (std::string(magic, 4) != "TRAC") return;

    uint8_t version;
    file.read(reinterpret_cast<char*>(&version), sizeof(version));
    if (version != 2) return;

    uint64_t record_count;
    file.read(reinterpret_cast<char*>(&record_count), sizeof(record_count));

    m_store.clear();

    // Read records
    for (uint64_t i = 0; i < record_count; ++i) {
        // Read key
        uint32_t key_len;
        file.read(reinterpret_cast<char*>(&key_len), sizeof(key_len));
        std::string key(key_len, '\0');
        file.read(&key[0], key_len);

        // Read record type
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

            // Read sequence
            uint32_t seq_len;
            file.read(reinterpret_cast<char*>(&seq_len), sizeof(seq_len));
            record.compressed_sequence.resize(seq_len);
            file.read(reinterpret_cast<char*>(record.compressed_sequence.data()), seq_len);

            // Read quality
            uint32_t qual_len;
            file.read(reinterpret_cast<char*>(&qual_len), sizeof(qual_len));
            record.compressed_quality.resize(qual_len);
            file.read(reinterpret_cast<char*>(record.compressed_quality.data()), qual_len);

            m_store[key] = record;
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