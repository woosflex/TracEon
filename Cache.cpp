#include "Cache.h"
#include "PlainTextStrategy.h"
#include "TwoBitDnaStrategy.h"
#include <fstream>
#include <sstream>

namespace TracEon {

Cache::Cache() {
    // We'll default to the plain text strategy on startup.
    setStrategy(StrategyType::PlainText);
}

Cache::~Cache() = default;

void Cache::setStrategy(StrategyType type) {
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

size_t Cache::getStoredSize(const std::string& key) const {
    if (m_store.count(key)) {
        const auto& record_variant = m_store.at(key);
        if (const auto* data = std::get_if<FastaRecordData>(&record_variant)) {
            // Just return the size of the raw stored data.
            return data->size();
        }
    }
    return 0;
}

size_t Cache::size() const {
    return m_store.size();
}

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

void Cache::loadFasta(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        // Silently fail if we can't open the file.
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

        if (line[0] == '>') { // New sequence header found
            // If we have a sequence buffered, save it first.
            if (!current_key.empty()) {
                set(current_key, current_value_stream.str());
            }

            // Parse the new key (the first word after '>')
            size_t first_space = line.find(' ');
            current_key = line.substr(1, first_space - 1);

            // Reset the stream for the new sequence data.
            current_value_stream.str("");
            current_value_stream.clear();
        } else { // Sequence data line
            current_value_stream << line;
        }
    }

    // After the loop, make sure to save the very last sequence.
    if (!current_key.empty()) {
        set(current_key, current_value_stream.str());
    }
}

} // namespace TracEon