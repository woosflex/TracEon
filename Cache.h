//
// Created by Adnan Raza on 9/26/2025.
//

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
        void loadFasta(const std::string& filepath);
        void loadFastq(const std::string& filepath);
        void save(const std::string& filepath);
        void restore(const std::string& filepath);
        std::string getFastaSequence(const std::string& key);
        size_t size() const;
        void set(const std::string& key, const std::string& value);
        std::string get(const std::string& key);

    private:
        std::unique_ptr<IEncodingStrategy> m_dna_strategy;
        std::unique_ptr<IEncodingStrategy> m_quality_strategy;
        using FastaRecordData = std::vector<unsigned char>;
        using RecordData = std::variant<FastaRecordData, FastqRecord>;
        std::unordered_map<std::string, RecordData> m_store;
    };

} // namespace TracEon

#endif //TRACEON_CACHE_H