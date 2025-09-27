#ifndef TRACEON_CACHE_H
#define TRACEON_CACHE_H

#include <string>
#include <unordered_map>
#include <memory>
#include <variant>
#include <cstddef>

#include "IEncodingStrategy.h"
#include "RecordTypes.h"
#include "DecodedRecordTypes.h"

namespace TracEon {

    // The main class for the TracEon in-memory cache.
    class Cache {
    public:
        Cache();
        ~Cache();

        // --- Configuration ---
        // Sets the encoding/compression strategy for subsequent 'set' operations.
        void setStrategy(StrategyType type);

        // --- Data Interaction ---
        // Adds or updates a key-value pair in the cache, encoding the value.
        void set(const std::string& key, const std::string& value);
        // Retrieves and decodes a value by its key. Returns empty string if not found.
        std::string get(const std::string& key);
        std::optional<DecodedFastqRecord> getFastqRecord(const std::string& key);

        // --- Status & Inspection ---
        // Returns the total number of items in the cache.
        size_t size() const;
        // Returns the size in bytes of the stored (and possibly compressed) data for a key.
        size_t getStoredSize(const std::string& key) const;

        // --- File I/O ---
        // Parses a FASTA file and loads its contents into the cache.
        void loadFasta(const std::string& filepath);
        void loadFastq(const std::string& filepath);
        // Writes the entire cache's state to a binary file.
        void save(const std::string& filepath);
        // Loads the cache's state from a binary file.
        void restore(const std::string& filepath);
        // void loadFastq(const std::string& filepath); // Future work

    private:
        StrategyType m_current_strategy_type;
        std::unique_ptr<IEncodingStrategy> m_dna_strategy;
        // We'll use this later for FASTQ quality scores.
        std::unique_ptr<IEncodingStrategy> m_quality_strategy;

        using FastaRecordData = std::vector<unsigned char>;
        using RecordData = std::variant<FastaRecordData, FastqRecord>;
        std::unordered_map<std::string, RecordData> m_store;
    };

} // namespace TracEon

#endif //TRACEON_CACHE_H