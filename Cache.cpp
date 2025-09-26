#include "Cache.h"

// We need to define the namespace we declared in the header file.
namespace TracEon {

    // The constructor. For now, it doesn't need to do anything.
    Cache::Cache() {
        // We'll initialize our strategies here later.
    }

    // The destructor. Also doesn't need to do anything yet.
    Cache::~Cache() = default;

    // The implementation for the size() method.
    // It just returns the current number of elements in our map.
    size_t Cache::size() const {
        return m_store.size();
    }

    void Cache::set(const std::string& key, const std::string& value) {
        // For now, we'll just treat the string's raw data as a vector of bytes
        // and store it as a FastaRecordData variant.
        FastaRecordData data(value.begin(), value.end());
        m_store[key] = data;
    }

    std::string Cache::get(const std::string& key) {
        // Check if the key exists in our map
        if (m_store.count(key)) {
            // Get the variant associated with the key
            const auto& record_variant = m_store.at(key);
            // We expect it to be a FastaRecordData (a vector of bytes)
            if (const FastaRecordData* data = std::get_if<FastaRecordData>(&record_variant)) {
                // Convert the vector of bytes back into a string and return it
                return {data->begin(), data->end()};
            }
        }
        // If the key is not found or it's the wrong type, return an empty string.
        return "";
    }

} // namespace TracEon
//
// Created by Adnan Raza on 9/26/2025.
//