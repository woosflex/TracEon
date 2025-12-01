#ifndef TRACEON_CACHE_H
#define TRACEON_CACHE_H

#include <string>
#include <string_view>
#include <vector>
#include <memory>
#include <optional>
#include <unordered_map>

#include "IEncodingStrategy.h"
#include "RecordTypes.h"
#include "DecodedRecordTypes.h"
#include "MapDefs.h" 

namespace TracEon {
    class SmartStrategy; // Forward declaration

    class Cache {
    public:
        Cache();
        ~Cache();

        // --- Data Access ---
        std::string get(const std::string& key);
        std::string_view getView(const std::string& key);
        std::optional<DecodedFastqRecord> getFastqRecord(const std::string& key);

        // --- NEW: Existence Check ---
        bool hasSequence(const std::string& key);

        // --- I/O Operations ---
        void loadFile(const std::string& filepath);
        void save(const std::string& filepath);
        void restore(const std::string& filepath);
        void set(const std::string& key, const std::string& value);
        
        size_t size() const;

        // --- Benchmarking Accessor ---
        IEncodingStrategy* getStrategy() const;

    private:
        std::unique_ptr<SmartStrategy> m_strategy;
        std::unordered_map<std::string, std::string> m_manual_store;
    };

} // namespace TracEon

#endif //TRACEON_CACHE_H