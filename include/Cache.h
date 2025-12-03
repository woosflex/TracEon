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
    class SmartStrategy;

    /**
     * @class Cache
     * @brief High-level interface for the TracEon genomic data store.
     * * Provides a unified API for loading, saving, and accessing genomic sequences.
     * Internally routes requests to the zero-copy SmartStrategy (for file-backed data)
     * or a manual fallback store (for runtime-added data).
     */
    class Cache {
    public:
        Cache();
        ~Cache();

        // --- Data Access ---
        
        /**
         * @brief Retrieve a sequence by ID (Legacy API).
         * @return A copy of the sequence string, or empty string if not found.
         */
        std::string get(const std::string& key);

        /**
         * @brief Retrieve a zero-copy view of a sequence (High Performance).
         * @return A std::string_view pointing directly to the memory-mapped region.
         * @note Safe to use concurrently if data is loaded.
         */
        std::string_view getView(const std::string& key);

        /**
         * @brief Retrieve both sequence and quality scores.
         */
        std::optional<DecodedFastqRecord> getFastqRecord(const std::string& key);

        /**
         * @brief Check if a key exists without retrieving data.
         */
        bool hasSequence(const std::string& key);

        // --- I/O Operations ---
        
        /**
         * @brief Parse a FASTA/FASTQ file into the cache.
         * @note Uses multithreaded parsing for files > 10MB.
         */
        void loadFile(const std::string& filepath);
        
        /**
         * @brief Save the current FILE-BACKED index to a binary cache file.
         * @note Does NOT persist items added via set() in v1.0.
         */
        void save(const std::string& filepath);
        
        /**
         * @brief Instantaneously map a binary cache file into memory.
         */
        void restore(const std::string& filepath);

        /**
         * @brief Manually add a key-value pair (Heap storage).
         */
        void set(const std::string& key, const std::string& value);
        
        size_t size() const;

        // --- Internal / Diagnostics ---
        IEncodingStrategy* getStrategy() const;

    private:
        std::unique_ptr<SmartStrategy> m_strategy;
        std::unordered_map<std::string, std::string> m_manual_store;
    };

} // namespace TracEon

#endif //TRACEON_CACHE_H