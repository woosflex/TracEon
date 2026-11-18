#include "Cache.h"
#include "SmartStrategy.h"

namespace TracEon {

    Cache::Cache(IndexMode mode) : m_strategy(std::make_unique<SmartStrategy>(mode)) {}
    Cache::~Cache() = default;

    // --- Data Access ---

    std::string Cache::get(std::string_view key) {
        return m_strategy->getSequence(key);
    }

    std::string_view Cache::getView(std::string_view key) {
        return m_strategy->getView(key);
    }

    std::optional<FastqRecordView> Cache::getRecordView(std::string_view key) {
        auto rec = m_strategy->findRecord(key);
        if (!rec) return std::nullopt;
        return FastqRecordView{rec->sequence, rec->quality};
    }

    std::optional<DecodedFastqRecord> Cache::getFastqRecord(std::string_view key) {
        auto rec = m_strategy->findRecord(key);
        if (!rec) return std::nullopt;
        return DecodedFastqRecord{
            std::string(rec->sequence),
            std::string(rec->quality)
        };
    }

    bool Cache::hasSequence(std::string_view key) {
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

    void Cache::save(const std::string& filepath, CacheStorage storage) {
        m_strategy->saveBinary(filepath, storage);
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

    std::string Cache::getQuality(std::string_view key) const {
        return m_strategy->getQuality(key);
    }

    std::vector<std::string> Cache::getAllKeys() const {
        return m_strategy->getAllKeys();
    }

    FileFormat Cache::getDetectedFormat() const {
        return m_strategy->getDetectedFormat();
    }

} // namespace TracEon
