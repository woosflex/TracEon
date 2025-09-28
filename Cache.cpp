#include "Cache.h"
#include "PlainTextStrategy.h"
#include "TwoBitDnaStrategy.h"
#include <fstream>
#include <sstream>
#include <ios>
#include <variant>

namespace TracEon {

Cache::Cache() {
    // Default to plain text for both strategies on startup.
    setStrategy(StrategyType::PlainText);
    m_quality_strategy = std::make_unique<PlainTextStrategy>();
}

Cache::~Cache() = default;

// --- Configuration ---

void Cache::setStrategy(StrategyType type) {
    m_current_strategy_type = type;
    switch (type) {
        case StrategyType::TwoBit:
            m_dna_strategy = std::make_unique<TwoBitDnaStrategy>();
            break;
        case StrategyType::PlainText:
        default:
            m_dna_strategy = std::make_unique<PlainTextStrategy>();
            break;
    }
}

// --- Data Interaction ---

void Cache::set(const std::string& key, const std::string& value) {
    FastaRecordData data = m_dna_strategy->encode(value);
    m_store[key] = data;
}

std::string Cache::get(const std::string& key) {
    if (m_store.count(key)) {
        const auto& record_variant = m_store.at(key);
        if (const auto* data = std::get_if<FastaRecordData>(&record_variant)) {
            return m_dna_strategy->decode(*data);
        }
    }
    return "";
}

std::optional<DecodedFastqRecord> Cache::getFastqRecord(const std::string& key) {
    if (m_store.count(key)) {
        const auto& record_variant = m_store.at(key);
        if (const auto* data = std::get_if<FastqRecord>(&record_variant)) {
            // It's a FASTQ record, so decode both parts.
            return DecodedFastqRecord{
                m_dna_strategy->decode(data->compressed_sequence),
                m_quality_strategy->decode(data->compressed_quality)
            };
        }
    }
    // Key not found or it wasn't a FASTQ record.
    return std::nullopt;
}

// --- Status & Inspection ---

size_t Cache::size() const {
    return m_store.size();
}

size_t Cache::getStoredSize(const std::string& key) const {
    if (m_store.count(key)) {
        const auto& record_variant = m_store.at(key);
        // Use std::visit to handle the different types in the variant.
        return std::visit([](const auto& value) {
            using T = std::decay_t<decltype(value)>;
            if constexpr (std::is_same_v<T, FastaRecordData>) {
                return value.size(); // It's a vector, return its size.
            } else if constexpr (std::is_same_v<T, FastqRecord>) {
                // It's a struct, return the combined size of its vectors.
                return value.compressed_sequence.size() + value.compressed_quality.size();
            }
            return size_t(0);
        }, record_variant);
    }
    return 0;
}

// --- File I/O ---

void Cache::loadFasta(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file.is_open()) return;

    std::string line, current_key;
    std::stringstream current_value_stream;

    while (std::getline(file, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;

        if (line[0] == '>') {
            if (!current_key.empty()) {
                set(current_key, current_value_stream.str());
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
        set(current_key, current_value_stream.str());
    }
}

void Cache::loadFastq(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file.is_open()) return;

    std::string header, sequence, plus, quality;

    // Read the file in chunks of 4 lines.
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

        // Create a FastqRecord, encode its parts, and store it.
        FastqRecord record;
        record.compressed_sequence = m_dna_strategy->encode(sequence);
        record.compressed_quality = m_quality_strategy->encode(quality);
        m_store[key] = record;
    }
}

void Cache::save(const std::string& filepath) {
    std::ofstream file(filepath, std::ios::binary);
    if (!file.is_open()) return;

    // Write Header
    const char magic[] = "TRAC";
    uint8_t version = 1;
    uint8_t strategy_id = static_cast<uint8_t>(m_current_strategy_type);
    uint64_t record_count = m_store.size();
    file.write(magic, 4);
    file.write(reinterpret_cast<const char*>(&version), sizeof(version));
    file.write(reinterpret_cast<const char*>(&strategy_id), sizeof(strategy_id));
    file.write(reinterpret_cast<const char*>(&record_count), sizeof(record_count));

    // Write Records
    for (const auto& pair : m_store) {
        // Write the key first
        const std::string& key = pair.first;
        uint64_t key_len = key.length();
        file.write(reinterpret_cast<const char*>(&key_len), sizeof(key_len));
        file.write(key.c_str(), key_len);

        // Now visit the variant to write the correct data type
        std::visit([&file](const auto& value){
            using T = std::decay_t<decltype(value)>;
            if constexpr (std::is_same_v<T, FastaRecordData>) {
                uint8_t type_id = 0; // 0 for FASTA
                file.write(reinterpret_cast<const char*>(&type_id), sizeof(type_id));
                uint64_t val_len = value.size();
                file.write(reinterpret_cast<const char*>(&val_len), sizeof(val_len));
                file.write(reinterpret_cast<const char*>(value.data()), val_len);
            } else if constexpr (std::is_same_v<T, FastqRecord>) {
                uint8_t type_id = 1; // 1 for FASTQ
                file.write(reinterpret_cast<const char*>(&type_id), sizeof(type_id));
                // Write sequence
                uint64_t seq_len = value.compressed_sequence.size();
                file.write(reinterpret_cast<const char*>(&seq_len), sizeof(seq_len));
                file.write(reinterpret_cast<const char*>(value.compressed_sequence.data()), seq_len);
                // Write quality
                uint64_t qual_len = value.compressed_quality.size();
                file.write(reinterpret_cast<const char*>(&qual_len), sizeof(qual_len));
                file.write(reinterpret_cast<const char*>(value.compressed_quality.data()), qual_len);
            }
        }, pair.second);
    }
}

void Cache::restore(const std::string& filepath) {
    std::ifstream file(filepath, std::ios::binary);
    if (!file.is_open()) return;

    m_store.clear();

    // Read and Verify Header
    char magic[4];
    uint8_t version, strategy_id;
    uint64_t record_count;
    file.read(magic, 4);
    if (std::string(magic, 4) != "TRAC") return;
    file.read(reinterpret_cast<char*>(&version), sizeof(version));
    file.read(reinterpret_cast<char*>(&strategy_id), sizeof(strategy_id));
    file.read(reinterpret_cast<char*>(&record_count), sizeof(record_count));
    setStrategy(static_cast<StrategyType>(strategy_id));

    // Read Records
    for (uint64_t i = 0; i < record_count; ++i) {
        uint64_t key_len;
        file.read(reinterpret_cast<char*>(&key_len), sizeof(key_len));
        std::string key(key_len, '\0');
        file.read(&key[0], key_len);

        uint8_t type_id;
        file.read(reinterpret_cast<char*>(&type_id), sizeof(type_id));

        if (type_id == 0) { // FASTA
            uint64_t val_len;
            file.read(reinterpret_cast<char*>(&val_len), sizeof(val_len));
            FastaRecordData value(val_len);
            file.read(reinterpret_cast<char*>(value.data()), val_len);
            m_store[key] = value;
        } else if (type_id == 1) { // FASTQ
            FastqRecord record;
            uint64_t seq_len, qual_len;
            // Read sequence
            file.read(reinterpret_cast<char*>(&seq_len), sizeof(seq_len));
            record.compressed_sequence.resize(seq_len);
            file.read(reinterpret_cast<char*>(record.compressed_sequence.data()), seq_len);
            // Read quality
            file.read(reinterpret_cast<char*>(&qual_len), sizeof(qual_len));
            record.compressed_quality.resize(qual_len);
            file.read(reinterpret_cast<char*>(record.compressed_quality.data()), qual_len);
            m_store[key] = record;
        }
    }
}

} // namespace TracEon