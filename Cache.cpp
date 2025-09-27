#include "Cache.h"
#include "PlainTextStrategy.h"
#include "TwoBitDnaStrategy.h"
#include <fstream>
#include <sstream>
#include <ios>

namespace TracEon {

Cache::Cache() {
    // Default to the plain text strategy on startup.
    setStrategy(StrategyType::PlainText);
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
    // Let the current strategy handle the encoding.
    FastaRecordData data = m_dna_strategy->encode(value);
    m_store[key] = data;
}

std::string Cache::get(const std::string& key) {
    if (m_store.count(key)) {
        const auto& record_variant = m_store.at(key);
        // Make sure we're getting the data type we expect.
        if (const auto* data = std::get_if<FastaRecordData>(&record_variant)) {
            // Let the current strategy handle decoding.
            return m_dna_strategy->decode(*data);
        }
    }
    // Didn't find the key, so return an empty string.
    return "";
}

// --- Status & Inspection ---

size_t Cache::size() const {
    return m_store.size();
}

size_t Cache::getStoredSize(const std::string& key) const {
    if (m_store.count(key)) {
        const auto& record_variant = m_store.at(key);
        if (const auto* data = std::get_if<FastaRecordData>(&record_variant)) {
            // Return the size of the raw stored data in bytes.
            return data->size();
        }
    }
    return 0;
}

// --- File I/O ---

void Cache::loadFasta(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        // A real-world version might throw an error here.
        return;
    }

    std::string line;
    std::string current_key;
    std::stringstream current_value_stream;

    while (std::getline(file, line)) {
        // Handle Windows-style line endings (\r\n)
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (line.empty()) {
            continue;
        }

        if (line[0] == '>') { // New sequence header
            // If we have a sequence buffered, save it.
            if (!current_key.empty()) {
                set(current_key, current_value_stream.str());
            }

            // Parse the new key (the first word after '>')
            size_t first_space = line.find(' ');
            current_key = line.substr(1, first_space - 1);

            // Reset for the new sequence.
            current_value_stream.str("");
            current_value_stream.clear();
        } else { // Sequence data
            current_value_stream << line;
        }
    }

    // Save the very last sequence in the file.
    if (!current_key.empty()) {
        set(current_key, current_value_stream.str());
    }
}

void Cache::save(const std::string& filepath) {
    std::ofstream file(filepath, std::ios::binary);
    if (!file.is_open()) return;

    // --- Write Header ---
    const char magic[] = "TRAC";
    uint8_t version = 1;
    uint8_t strategy_id = static_cast<uint8_t>(m_current_strategy_type);
    uint64_t record_count = m_store.size();

    file.write(magic, 4);
    file.write(reinterpret_cast<const char*>(&version), sizeof(version));
    file.write(reinterpret_cast<const char*>(&strategy_id), sizeof(strategy_id));
    file.write(reinterpret_cast<const char*>(&record_count), sizeof(record_count));

    // --- Write Records ---
    for (const auto& pair : m_store) {
        // Write key (length-prefixed)
        const std::string& key = pair.first;
        uint64_t key_len = key.length();
        file.write(reinterpret_cast<const char*>(&key_len), sizeof(key_len));
        file.write(key.c_str(), key_len);

        // Write value (length-prefixed)
        if (const auto* data = std::get_if<FastaRecordData>(&pair.second)) {
            uint64_t val_len = data->size();
            file.write(reinterpret_cast<const char*>(&val_len), sizeof(val_len));
            file.write(reinterpret_cast<const char*>(data->data()), val_len);
        }
    }
}

void Cache::restore(const std::string& filepath) {
    std::ifstream file(filepath, std::ios::binary);
    if (!file.is_open()) return;

    m_store.clear(); // Start with a fresh cache

    // --- Read and Verify Header ---
    char magic[4];
    uint8_t version;
    uint8_t strategy_id;
    uint64_t record_count;

    file.read(magic, 4);
    if (std::string(magic, 4) != "TRAC") return; // Not a TracEon file

    file.read(reinterpret_cast<char*>(&version), sizeof(version));
    file.read(reinterpret_cast<char*>(&strategy_id), sizeof(strategy_id));
    file.read(reinterpret_cast<char*>(&record_count), sizeof(record_count));

    // Set our strategy to match the one from the file.
    setStrategy(static_cast<StrategyType>(strategy_id));

    // --- Read Records ---
    for (uint64_t i = 0; i < record_count; ++i) {
        // Read key
        uint64_t key_len;
        file.read(reinterpret_cast<char*>(&key_len), sizeof(key_len));
        std::string key(key_len, '\0');
        file.read(&key[0], key_len);

        // Read value
        uint64_t val_len;
        file.read(reinterpret_cast<char*>(&val_len), sizeof(val_len));
        FastaRecordData value(val_len);
        file.read(reinterpret_cast<char*>(value.data()), val_len);

        // Insert directly to avoid re-encoding.
        m_store[key] = value;
    }
}

} // namespace TracEon