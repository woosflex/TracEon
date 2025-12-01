#ifndef TRACEON_CACHE_H
#define TRACEON_CACHE_H

#include <string>
#include <string_view>
#include <vector>
#include <memory>
#include <optional>

#include "IEncodingStrategy.h"
#include "RecordTypes.h"
#include "DecodedRecordTypes.h"
#include "MapDefs.h" // For HashMap definition if needed elsewhere

namespace TracEon {
    class SmartStrategy; // Forward declare to avoid circular dependency issues if any

    class Cache {
    public:
        Cache();
        ~Cache();

        // --- Data Access ---
        
        // Legacy: Returns a copy of the string (compatibility)
        std::string get(const std::string& key);

        // V2.0: Returns a Zero-Copy view into the arena/mmap
        std::string_view getView(const std::string& key);

        std::optional<DecodedFastqRecord> getFastqRecord(const std::string& key);

        // --- I/O Operations ---
        
        void loadFile(const std::string& filepath);
        
        // Optimized Binary I/O
        void save(const std::string& filepath);
        void restore(const std::string& filepath);

        // Manual Insertion (Legacy support)
        void set(const std::string& key, const std::string& value);
        
        size_t size() const;

    private:
        // Use unique_ptr to abstract the strategy
        std::unique_ptr<SmartStrategy> m_strategy;

        // Legacy fallback store for manually set items (not used in high-perf path)
        std::unordered_map<std::string, std::string> m_manual_store;
    };

} // namespace TracEon

#endif //TRACEON_CACHE_H