#ifndef TRACEON_CACHE_H
#define TRACEON_CACHE_H

#include <string>
#include <unordered_map>
#include <memory>
#include <variant>
#include <cstddef>
#include <optional>

#include "IEncodingStrategy.h"
#include "RecordTypes.h"
#include "DecodedRecordTypes.h"

namespace TracEon {

    // The main class for the TracEon in-memory cache.
    // It now uses a single, smart strategy for all encoding.
    class Cache {
    public:
        Cache();
        ~Cache();

        // Data Interaction
        std::string get(const std::string& key);
        std::optional<DecodedFastqRecord> getFastqRecord(const std::string& key);

        // Status & Inspection
        size_t size() const;
        size_t getStoredSize(const std::string& key) const;

        // File I/O (existing interface)
        void loadFile(const std::string& filepath); // <-- The unified loader
        void save(const std::string& filepath);
        void restore(const std::string& filepath);

        // SmartStrategy specific file operations (NEW - these were missing!)
        void loadSmartFile(const std::string& filepath);
        void saveSmartBinary(const std::string& filepath);
        void loadSmartBinary(const std::string& filepath);

    private:
        // Private helpers for loading specific file types
        void loadFasta(std::ifstream& file);
        void loadFastq(std::ifstream& file);

        // A single, smart strategy for all data types.
        std::unique_ptr<IEncodingStrategy> m_strategy;

        using FastaRecordData = std::vector<unsigned char>;
        using RecordData = std::variant<FastaRecordData, FastqRecord>;
        std::unordered_map<std::string, RecordData> m_store;
    };

} // namespace TracEon

#endif //TRACEON_CACHE_H