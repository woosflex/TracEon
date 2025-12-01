#include "SmartStrategy.h"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <cstring>
#include <stdexcept>
#include <mutex>
#include <shared_mutex>
#include <thread>
#include <vector>
#include <cassert>

#ifdef _WIN32
    #define WIN32_LEAN_AND_MEAN
    #include <windows.h>
#else
    #include <sys/mman.h>
    #include <sys/stat.h>
    #include <fcntl.h>
    #include <unistd.h>
#endif

namespace TracEon {

struct MMapHandle {
    void* data = nullptr;
    size_t size = 0;
#ifdef _WIN32
    HANDLE hFile = INVALID_HANDLE_VALUE;
    HANDLE hMap = NULL;
#else
    int fd = -1;
#endif

    MMapHandle() = default;
    ~MMapHandle() { cleanup(); }

    void cleanup() {
#ifdef _WIN32
        if (data) UnmapViewOfFile(data);
        if (hMap) CloseHandle(hMap);
        if (hFile != INVALID_HANDLE_VALUE) CloseHandle(hFile);
        data = nullptr; hMap = NULL; hFile = INVALID_HANDLE_VALUE;
#else
        if (data && data != MAP_FAILED) munmap(data, size);
        if (fd != -1) close(fd);
        data = nullptr; fd = -1;
#endif
    }
};

static const char MAGIC_BYTES[] = "MMAP";

SmartStrategy::SmartStrategy() : detected_format_(FileFormat::UNKNOWN) {}
SmartStrategy::~SmartStrategy() { clearCache(); }

std::vector<unsigned char> SmartStrategy::encode(const std::string& data, DataTypeHint hint) const {
    return {data.begin(), data.end()};
}

std::string SmartStrategy::decode(const std::vector<unsigned char>& data) const {
    return {data.begin(), data.end()};
}

// --- Multithreaded Parsers ---

void SmartStrategy::parseFastaMultithreaded(std::string_view content) {
    const size_t num_threads = std::thread::hardware_concurrency();
    const size_t chunk_size = content.size() / num_threads;
    
    if (chunk_size < 1024 * 1024) { 
        parseFastaOptimized(content);
        return;
    }
    
    std::vector<std::thread> threads;
    std::vector<HashMap<std::string, SequenceView>> thread_caches(num_threads);
    
    auto worker = [&](size_t thread_id, size_t start, size_t end) {
        const char* ptr = content.data() + start;
        const char* chunk_end = content.data() + end;
        const char* global_end = content.data() + content.size();
        
        if (thread_id > 0) {
            while (ptr > content.data() && *(ptr - 1) != '\n') --ptr;
            while (ptr < chunk_end && *ptr != '>') ++ptr;
        }
        
        while (ptr < chunk_end) {
            while (ptr < chunk_end && (*ptr == '\n' || *ptr == '\r' || *ptr == ' ')) ++ptr;
            if (ptr >= chunk_end || *ptr != '>') break;
            
            ++ptr; 
            const char* id_start = ptr;
            while (ptr < global_end && *ptr != ' ' && *ptr != '\t' && *ptr != '\n' && *ptr != '\r') ++ptr;
            const char* id_end = ptr;
            
            while (ptr < global_end && *ptr != '\n' && *ptr != '\r') ++ptr;
            while (ptr < global_end && (*ptr == '\n' || *ptr == '\r')) ++ptr;
            
            const char* seq_start = ptr;
            const char* seq_end = seq_start;
            while (ptr < global_end) {
                if (*ptr == '>') {
                    seq_end = ptr;
                    while (seq_end > seq_start && (*(seq_end-1) == '\n' || *(seq_end-1) == '\r' || *(seq_end-1) == ' ')) --seq_end;
                    break;
                }
                ++ptr;
            }
            if (ptr >= global_end) {
                seq_end = global_end;
                while (seq_end > seq_start && (*(seq_end-1) == '\n' || *(seq_end-1) == '\r' || *(seq_end-1) == ' ')) --seq_end;
            }
            
            std::string_view id(id_start, id_end - id_start);
            std::string_view seq(seq_start, seq_end - seq_start);
            
            if (!id.empty() && !seq.empty()) {
                #ifndef NDEBUG
                assert(id.find('\n') == std::string_view::npos);
                assert(id.find('\r') == std::string_view::npos);
                #endif
                thread_caches[thread_id].emplace(std::string(id), SequenceView{id, seq, {}});
            }
        }
    };
    
    for (size_t i = 0; i < num_threads; ++i) {
        size_t start = i * chunk_size;
        size_t end = (i == num_threads - 1) ? content.size() : (i + 1) * chunk_size;
        threads.emplace_back(worker, i, start, end);
    }
    for (auto& t : threads) t.join();
    
    size_t total_size = 0;
    for (const auto& cache : thread_caches) total_size += cache.size();
    file_cache_.reserve(total_size);
    
    for (auto& cache : thread_caches) {
        for (auto& entry : cache) {
            file_cache_.insert(std::move(entry)); 
        }
    }
}

void SmartStrategy::parseFastqMultithreaded(std::string_view content) {
    const size_t num_threads = std::thread::hardware_concurrency();
    const size_t chunk_size = content.size() / num_threads;
    
    if (chunk_size < 1024 * 1024) {
        parseFastqOptimized(content);
        return;
    }
    
    std::vector<std::thread> threads;
    std::vector<HashMap<std::string, SequenceView>> thread_caches(num_threads);
    
    auto worker = [&](size_t thread_id, size_t start, size_t end) {
        const char* ptr = content.data() + start;
        const char* chunk_end = content.data() + end;
        const char* global_end = content.data() + content.size();
        
        if (thread_id > 0) {
            while (ptr > content.data() && *(ptr - 1) != '\n') --ptr;
            while (ptr < chunk_end) {
                if (*ptr == '@') {
                    const char* check = ptr - 1;
                    while (check > content.data() && (*check == '\n' || *check == '\r')) --check;
                    while (check > content.data() && *check != '\n' && *check != '\r') --check;
                    if (check > content.data()) ++check;
                    if (*check != '+') break;
                }
                ++ptr;
            }
        }
        
        while (ptr < chunk_end) {
            while (ptr < chunk_end && (*ptr == '\n' || *ptr == '\r' || *ptr == ' ')) ++ptr;
            if (ptr >= chunk_end || *ptr != '@') break;
            
            ++ptr;
            const char* id_start = ptr;
            while (ptr < global_end && *ptr != ' ' && *ptr != '\t' && *ptr != '\n' && *ptr != '\r') ++ptr;
            const char* id_end = ptr;
            
            while (ptr < global_end && *ptr != '\n' && *ptr != '\r') ++ptr;
            while (ptr < global_end && (*ptr == '\n' || *ptr == '\r')) ++ptr;
            
            const char* seq_start = ptr;
            while (ptr < global_end && *ptr != '\n' && *ptr != '\r') ++ptr;
            const char* seq_end = ptr;
            while (ptr < global_end && (*ptr == '\n' || *ptr == '\r')) ++ptr;
            while (ptr < global_end && *ptr != '\n' && *ptr != '\r') ++ptr;
            while (ptr < global_end && (*ptr == '\n' || *ptr == '\r')) ++ptr;
            
            const char* qual_start = ptr;
            while (ptr < global_end && *ptr != '\n' && *ptr != '\r') ++ptr;
            const char* qual_end = ptr;
            while (ptr < global_end && (*ptr == '\n' || *ptr == '\r')) ++ptr;
            
            std::string_view id(id_start, id_end - id_start);
            std::string_view seq(seq_start, seq_end - seq_start);
            std::string_view qual(qual_start, qual_end - qual_start);
            
            if (!id.empty() && !seq.empty()) {
                #ifndef NDEBUG
                assert(id.find('\n') == std::string_view::npos);
                assert(id.find('\r') == std::string_view::npos);
                #endif
                thread_caches[thread_id].emplace(std::string(id), SequenceView{id, seq, qual});
            }
        }
    };
    
    for (size_t i = 0; i < num_threads; ++i) {
        size_t start = i * chunk_size;
        size_t end = (i == num_threads - 1) ? content.size() : (i + 1) * chunk_size;
        threads.emplace_back(worker, i, start, end);
    }
    for (auto& t : threads) t.join();
    
    size_t total_size = 0;
    for (const auto& cache : thread_caches) total_size += cache.size();
    file_cache_.reserve(total_size);
    for (auto& cache : thread_caches) {
        for (auto& entry : cache) file_cache_.insert(std::move(entry));
    }
}

// --- Optimized Parsers ---

void SmartStrategy::parseFastaOptimized(std::string_view content) {
    const char* ptr = content.data();
    const char* end = ptr + content.size();
    
    while (ptr < end) {
        while (ptr < end && (*ptr == '\n' || *ptr == '\r' || *ptr == ' ')) ++ptr;
        if (ptr >= end || *ptr != '>') break;
        
        ++ptr;
        const char* id_start = ptr;
        while (ptr < end && *ptr != ' ' && *ptr != '\t' && *ptr != '\n' && *ptr != '\r') ++ptr;
        const char* id_end = ptr;
        
        while (ptr < end && *ptr != '\n' && *ptr != '\r') ++ptr;
        while (ptr < end && (*ptr == '\n' || *ptr == '\r')) ++ptr;
        
        const char* seq_start = ptr;
        const char* seq_end = seq_start;
        while (ptr < end) {
            if (*ptr == '>') {
                seq_end = ptr;
                while (seq_end > seq_start && (*(seq_end-1) == '\n' || *(seq_end-1) == '\r' || *(seq_end-1) == ' ')) --seq_end;
                break;
            }
            ++ptr;
        }
        if (ptr >= end) {
            seq_end = end;
            while (seq_end > seq_start && (*(seq_end-1) == '\n' || *(seq_end-1) == '\r' || *(seq_end-1) == ' ')) --seq_end;
        }
        
        std::string_view id(id_start, id_end - id_start);
        std::string_view seq(seq_start, seq_end - seq_start);
        if (!id.empty() && !seq.empty()) {
            #ifndef NDEBUG
            assert(id.find('\n') == std::string_view::npos);
            assert(id.find('\r') == std::string_view::npos);
            #endif
            file_cache_.emplace(std::string(id), SequenceView{id, seq, {}});
        }
    }
}

void SmartStrategy::parseFastqOptimized(std::string_view content) {
    const char* ptr = content.data();
    const char* end = ptr + content.size();
    
    while (ptr < end) {
        while (ptr < end && (*ptr == '\n' || *ptr == '\r' || *ptr == ' ')) ++ptr;
        if (ptr >= end || *ptr != '@') break;
        
        ++ptr;
        const char* id_start = ptr;
        while (ptr < end && *ptr != ' ' && *ptr != '\t' && *ptr != '\n' && *ptr != '\r') ++ptr;
        const char* id_end = ptr;
        
        while (ptr < end && *ptr != '\n' && *ptr != '\r') ++ptr;
        while (ptr < end && (*ptr == '\n' || *ptr == '\r')) ++ptr;
        
        const char* seq_start = ptr;
        while (ptr < end && *ptr != '\n' && *ptr != '\r') ++ptr;
        const char* seq_end = ptr;
        while (ptr < end && (*ptr == '\n' || *ptr == '\r')) ++ptr;
        
        while (ptr < end && *ptr != '\n' && *ptr != '\r') ++ptr;
        while (ptr < end && (*ptr == '\n' || *ptr == '\r')) ++ptr;
        
        const char* qual_start = ptr;
        while (ptr < end && *ptr != '\n' && *ptr != '\r') ++ptr;
        const char* qual_end = ptr;
        while (ptr < end && (*ptr == '\n' || *ptr == '\r')) ++ptr;
        
        std::string_view id(id_start, id_end - id_start);
        std::string_view seq(seq_start, seq_end - seq_start);
        std::string_view qual(qual_start, qual_end - qual_start);
        
        if (!id.empty() && !seq.empty()) {
            #ifndef NDEBUG
            assert(id.find('\n') == std::string_view::npos);
            assert(id.find('\r') == std::string_view::npos);
            #endif
            file_cache_.emplace(std::string(id), SequenceView{id, seq, qual});
        }
    }
}

// --- Loading / Saving ---

void SmartStrategy::loadFile(const std::string& filepath) {
    std::unique_lock<std::shared_mutex> lock(cache_mutex_);
    clearCache();

    std::ifstream file(filepath, std::ios::binary | std::ios::ate);
    if (!file) throw std::runtime_error("Cannot open file: " + filepath);
    size_t fileSize = file.tellg();
    file.seekg(0, std::ios::beg);

    text_arena_.resize(fileSize);
    if (!file.read(text_arena_.data(), fileSize)) throw std::runtime_error("Read failed");

    std::string_view content(text_arena_.data(), fileSize);
    if (content.empty()) return;
    
    bool isFastq = (content[0] == '@');
    file_cache_.reserve(fileSize / 100);
    
    const size_t MULTITHREAD_THRESHOLD = 10 * 1024 * 1024;
    if (fileSize > MULTITHREAD_THRESHOLD) {
        if (isFastq) parseFastqMultithreaded(content);
        else parseFastaMultithreaded(content);
    } else {
        if (isFastq) parseFastqOptimized(content);
        else parseFastaOptimized(content);
    }
    determine_format_from_cache();
}

void SmartStrategy::saveBinary(const std::string& filepath) const {
    std::shared_lock<std::shared_mutex> lock(cache_mutex_);
    std::ofstream out(filepath, std::ios::binary);
    if (!out) throw std::runtime_error("Cannot write: " + filepath);

    out.write(MAGIC_BYTES, 4);
    uint64_t count = file_cache_.size();
    out.write(reinterpret_cast<const char*>(&count), sizeof(count));

    for (const auto& [key, view] : file_cache_) {
        // Key is std::string now
        uint32_t len = static_cast<uint32_t>(key.size());
        out.write(reinterpret_cast<const char*>(&len), sizeof(len));
        out.write(key.data(), len);

        // Sequence/Qual are views
        len = static_cast<uint32_t>(view.sequence.size());
        out.write(reinterpret_cast<const char*>(&len), sizeof(len));
        out.write(view.sequence.data(), len);

        len = static_cast<uint32_t>(view.quality.size());
        out.write(reinterpret_cast<const char*>(&len), sizeof(len));
        if (len > 0) out.write(view.quality.data(), len);
    }
}

void SmartStrategy::loadBinary(const std::string& filepath) {
    std::unique_lock<std::shared_mutex> lock(cache_mutex_);
    clearCache();
    
    mmap_handle_ = std::make_unique<MMapHandle>();
    // ... MMAP Setup ...
#ifdef _WIN32
    mmap_handle_->hFile = CreateFileA(filepath.c_str(), GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
    GetFileSizeEx(mmap_handle_->hFile, (PLARGE_INTEGER)&mmap_handle_->size);
    mmap_handle_->hMap = CreateFileMappingA(mmap_handle_->hFile, NULL, PAGE_READONLY, 0, 0, NULL);
    mmap_handle_->data = MapViewOfFile(mmap_handle_->hMap, FILE_MAP_READ, 0, 0, 0);
#else
    mmap_handle_->fd = open(filepath.c_str(), O_RDONLY);
    struct stat sb; fstat(mmap_handle_->fd, &sb);
    mmap_handle_->size = sb.st_size;
    mmap_handle_->data = mmap(NULL, mmap_handle_->size, PROT_READ, MAP_PRIVATE, mmap_handle_->fd, 0);
#endif

    const char* ptr = static_cast<const char*>(mmap_handle_->data);
    const char* end = ptr + mmap_handle_->size;

    if (mmap_handle_->size < 4 || std::strncmp(ptr, MAGIC_BYTES, 4) != 0) throw std::runtime_error("Invalid binary");
    ptr += 4;

    uint64_t count = *reinterpret_cast<const uint64_t*>(ptr);
    ptr += sizeof(uint64_t);
    file_cache_.reserve(count);

    for (uint64_t i = 0; i < count; ++i) {
        if (ptr >= end) break;
        uint32_t len = *reinterpret_cast<const uint32_t*>(ptr); ptr += 4;
        std::string_view id_view(ptr, len); ptr += len;

        len = *reinterpret_cast<const uint32_t*>(ptr); ptr += 4;
        std::string_view seq(ptr, len); ptr += len;

        len = *reinterpret_cast<const uint32_t*>(ptr); ptr += 4;
        std::string_view qual;
        if (len > 0) { qual = std::string_view(ptr, len); ptr += len; }

        // HYBRID ARCHITECTURE: 
        // 1. Copy ID into std::string (KEY)
        // 2. Keep Sequence as string_view (VALUE) pointing to MMAP
        // NOTE: We deliberately allow SequenceView.id to point to MMAP (id_view) 
        // to ensure it points to stable memory, while the map key is an owned copy.
        file_cache_.emplace(std::string(id_view), SequenceView{id_view, seq, qual});
    }
    determine_format_from_cache();
}

void SmartStrategy::clearCache() {
    file_cache_.clear();
    text_arena_.clear();
    text_arena_.shrink_to_fit();
    mmap_handle_.reset();
}

// --- ACCESS ---

std::string_view SmartStrategy::getView(const std::string_view& sequence_id) const {
    std::shared_lock<std::shared_mutex> lock(cache_mutex_);
    auto it = file_cache_.find(std::string(sequence_id)); 
    if (it != file_cache_.end()) return it->second.sequence;
    return {};
}

std::string SmartStrategy::getSequence(const std::string& sequence_id) const {
    return std::string(getView(sequence_id));
}
std::string SmartStrategy::getQuality(const std::string& sequence_id) const {
    std::shared_lock<std::shared_mutex> lock(cache_mutex_);
    auto it = file_cache_.find(sequence_id);
    return (it != file_cache_.end()) ? std::string(it->second.quality) : "";
}
bool SmartStrategy::hasSequence(const std::string& sequence_id) const {
    std::shared_lock<std::shared_mutex> lock(cache_mutex_);
    return file_cache_.find(sequence_id) != file_cache_.end();
}
size_t SmartStrategy::getFileCacheSize() const {
    std::shared_lock<std::shared_mutex> lock(cache_mutex_);
    return file_cache_.size();
}
std::vector<std::string> SmartStrategy::getAllKeys() const {
    std::shared_lock<std::shared_mutex> lock(cache_mutex_);
    std::vector<std::string> keys;
    keys.reserve(file_cache_.size());
    for (const auto& pair : file_cache_) keys.emplace_back(pair.first);
    return keys;
}
void SmartStrategy::determine_format_from_cache() { /* ... */ }
bool SmartStrategy::isNucleotideSequence(std::string_view data) const { return !data.empty(); }
bool SmartStrategy::hasRNA(std::string_view data) const { return data.find('U') != std::string_view::npos; }

} // namespace