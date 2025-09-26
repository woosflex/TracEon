#ifndef TRACEON_CACHE_H
#define TRACEON_CACHE_H

#include <string>
#include <unordered_map>
#include <memory>
#include <variant>
#include <cstddef>

#include "IEncodingStrategy.h"
#include "RecordTypes.h"

namespace TracEon {

    class Cache {
    public:
        Cache();
        ~Cache();

        // --- Configuration ---
        void setStrategy(StrategyType type);

        // --- Data Interaction ---
        void set(const std::string& key, const std::string& value);
        std::string get(const std::string& key);

        // --- Status & Inspection ---
        size_t size() const;
        size_t getStoredSize(const std::string& key) const;

        // --- Unimplemented file I/O ---
        void loadFasta(const std::string& filepath);
        void loadFastq(const std::string& filepath);
        void save(const std::string& filepath);
        void restore(const std::string& filepath);

    private:
        StrategyType m_current_strategy_type;
        std::unique_ptr<IEncodingStrategy> m_dna_strategy;
        // We'll use this later for FASTQ files.
        std::unique_ptr<IEncodingStrategy> m_quality_strategy;

        using FastaRecordData = std::vector<unsigned char>;
        using RecordData = std::variant<FastaRecordData, FastqRecord>;
        std::unordered_map<std::string, RecordData> m_store;
    };

} // namespace TracEon

#endif //TRACEON_CACHE_H