#include "Cache.h"
#include "PlainTextStrategy.h"

namespace TracEon {

    Cache::Cache() {
        // For now, our cache defaults to using the plain text strategy.
        m_dna_strategy = std::make_unique<PlainTextStrategy>();
    }

    Cache::~Cache() = default;

    // Just returns the current number of key-value pairs in the store.
    size_t Cache::size() const {
        return m_store.size();
    }

    void Cache::set(const std::string& key, const std::string& value) {
        // Delegate the encoding work to whatever strategy we're using.
        FastaRecordData data = m_dna_strategy->encode(value);
        m_store[key] = data;
    }

    std::string Cache::get(const std::string& key) {
        // See if the key even exists.
        if (m_store.count(key)) {
            const auto& record_variant = m_store.at(key);
            // Make sure the value is the type we expect (a simple byte vector for now).
            if (const FastaRecordData* data = std::get_if<FastaRecordData>(&record_variant)) {
                // Delegate the decoding work to the strategy and send it back.
                return m_dna_strategy->decode(*data);
            }
        }
        // If we didn't find the key, just return an empty string.
        return "";
    }

} // namespace TracEon