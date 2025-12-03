#include "../include/Cache.h"
#include "../include/SmartStrategy.h"
#include <stdexcept>

namespace TracEon {

    Cache::Cache() : m_strategy(std::make_unique<SmartStrategy>()) {}
    Cache::~Cache() = default;

    // --- Data Access ---

    std::string Cache::get(const std::string& key) {
        // Check Heap Store first (Write-Through)
        if (m_manual_store.count(key)) return m_manual_store.at(key);
        // Check Zero-Copy Store
        return m_strategy->getSequence(key);
    }

    std::string_view Cache::getView(const std::string& key) {
        // Optimization: Prioritize Zero-Copy strategy for views
        std::string_view v = m_strategy->getView(key);
        if (!v.empty()) return v;
        
        // Fallback: Heap store (pointer stability not guaranteed if map rehashes, be careful)
        if (m_manual_store.count(key)) return m_manual_store.at(key);
        return {};
    }

    std::optional<DecodedFastqRecord> Cache::getFastqRecord(const std::string& key) {
        if (m_strategy->hasSequence(key)) {
            return DecodedFastqRecord{
                m_strategy->getSequence(key),
                m_strategy->getQuality(key)
            };
        }
        return std::nullopt;
    }

    // --- Fix: Implemented missing method ---
    bool Cache::hasSequence(const std::string& key) {
        if (m_strategy->hasSequence(key)) return true;
        return m_manual_store.count(key) > 0;
    }

    void Cache::loadFile(const std::string& filepath) {
        m_strategy->loadFile(filepath);
    }

    void Cache::save(const std::string& filepath) {
        // Note: In v1.0, save() only persists the SmartStrategy (file-backed) index.
        // Manual entries in m_manual_store are NOT merged into the binary format yet.
        m_strategy->saveBinary(filepath);
    }

    void Cache::restore(const std::string& filepath) {
        m_strategy->loadBinary(filepath);
    }

    void Cache::set(const std::string& key, const std::string& value) {
        m_manual_store[key] = value;
    }

    size_t Cache::size() const {
        return m_strategy->getFileCacheSize() + m_manual_store.size();
    }

    IEncodingStrategy* Cache::getStrategy() const {
        return m_strategy.get();
    }

} // namespace TracEon