#include "Cache.h"
#include "SmartStrategy.h"

namespace TracEon {

    Cache::Cache(IndexMode mode) : m_strategy(std::make_unique<SmartStrategy>(mode)) {}
    Cache::~Cache() = default;

    // --- Data Access ---

    std::string Cache::get(const std::string& key) {
        return m_strategy->getSequence(key);
    }

    std::string_view Cache::getView(const std::string& key) {
        return m_strategy->getView(key);
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

    bool Cache::hasSequence(const std::string& key) {
        return m_strategy->hasSequence(key);
    }

    void Cache::loadFile(const std::string& filepath) {
        m_strategy->loadFile(filepath);
    }

    void Cache::loadGzipFile(const std::string& filepath) {
        m_strategy->loadGzipFile(filepath);
    }

    void Cache::save(const std::string& filepath) {
        m_strategy->saveBinary(filepath);
    }

    void Cache::restore(const std::string& filepath) {
        m_strategy->loadBinary(filepath);
    }

    void Cache::set(const std::string& key, const std::string& value) {
        m_strategy->addEntry(key, value, "");
    }

    void Cache::clearCache() {
        m_strategy->clearCache();
    }

    size_t Cache::size() const {
        return m_strategy->getFileCacheSize();
    }

    size_t Cache::getArenaBytes() const {
        return m_strategy->getArenaBytes();
    }

    IndexMode Cache::getIndexMode() const {
        return m_strategy->getIndexMode();
    }

    std::string Cache::getQuality(const std::string& key) const {
        return m_strategy->getQuality(key);
    }

    std::vector<std::string> Cache::getAllKeys() const {
        return m_strategy->getAllKeys();
    }

    FileFormat Cache::getDetectedFormat() const {
        return m_strategy->getDetectedFormat();
    }

    IEncodingStrategy* Cache::getStrategy() const {
        return m_strategy.get();
    }

} // namespace TracEon
