#include "Cache.h"
#include "SmartStrategy.h"
#include <stdexcept>

namespace TracEon {

    Cache::Cache() : m_strategy(std::make_unique<SmartStrategy>()) {}
    Cache::~Cache() = default;

    std::string Cache::get(const std::string& key) {
        // 1. Try manual store
        if (m_manual_store.count(key)) return m_manual_store.at(key);

        // 2. Delegate to SmartStrategy
        return m_strategy->getSequence(key);
    }

    std::string_view Cache::getView(const std::string& key) {
        // 1. Check SmartStrategy (Primary High-Perf Path)
        // Note: Manual store items are std::string, so we can return view to them,
        // but the pointer is only valid as long as m_manual_store doesn't rehash.
        // For safety/simplicity in V2, we prioritize the Strategy.
        
        std::string_view v = m_strategy->getView(key);
        if (!v.empty()) return v;

        // 2. Fallback to manual store
        if (m_manual_store.count(key)) {
            return m_manual_store.at(key);
        }

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

    void Cache::loadFile(const std::string& filepath) {
        m_strategy->loadFile(filepath);
    }

    void Cache::save(const std::string& filepath) {
        // V2.0 always uses the optimized binary format
        m_strategy->saveBinary(filepath);
    }

    void Cache::restore(const std::string& filepath) {
        // V2.0 uses MMAP loading
        m_strategy->loadBinary(filepath);
    }

    void Cache::set(const std::string& key, const std::string& value) {
        // Items manually added go into a separate fallback map
        // because they aren't backed by the Arena/MMAP
        m_manual_store[key] = value;
    }

    size_t Cache::size() const {
        return m_strategy->getFileCacheSize() + m_manual_store.size();
    }

} // namespace TracEon