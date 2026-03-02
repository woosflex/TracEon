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
#include <cmath>
#include <functional> 
#include <cctype>
#include <zlib.h>

#ifdef _WIN32
    #define WIN32_LEAN_AND_MEAN
    #define NOMINMAX
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

static const char MAGIC_BYTES[] = "TRO\x01"; // v1.1+ format (was "MMAP" in v1.0)

SmartStrategy::SmartStrategy() : detected_format_(FileFormat::UNKNOWN), file_cache_(GenomeIndex{}) {}
SmartStrategy::~SmartStrategy() { clearCache(); }

std::vector<unsigned char> SmartStrategy::encode(const std::string& data, DataTypeHint hint) const {
    return {data.begin(), data.end()};
}
std::string SmartStrategy::decode(const std::vector<unsigned char>& data) const {
    return {data.begin(), data.end()};
}
inline uint64_t SmartStrategy::hash_key(std::string_view key) const {
    return std::hash<std::string_view>{}(key);
}

void SmartStrategy::clearInternal() {
    data_ready_.store(false, std::memory_order_release);
    if (std::holds_alternative<GenomeIndex>(file_cache_)) std::get<GenomeIndex>(file_cache_).clear();
    else std::get<NGSIndex>(file_cache_).clear();
    text_arena_.clear();
    text_arena_.shrink_to_fit();
    manual_store_.clear();
    mmap_handle_.reset();
}

void SmartStrategy::clearCache() {
    std::unique_lock<std::shared_mutex> lock(cache_mutex_);
    clearInternal();
}

void SmartStrategy::addEntry(const std::string& id, const std::string& seq, const std::string& qual) {
    std::unique_lock<std::shared_mutex> lock(cache_mutex_);
    data_ready_.store(false, std::memory_order_release);

    // Store copies in manual_store_ (std::deque: push_back never invalidates
    // references to existing elements, so these string_views stay valid).
    manual_store_.push_back(id);
    std::string_view id_view(manual_store_.back());
    manual_store_.push_back(seq);
    std::string_view seq_view(manual_store_.back());
    manual_store_.push_back(qual);
    std::string_view qual_view(manual_store_.back());

    if (std::holds_alternative<GenomeIndex>(file_cache_)) {
        std::get<GenomeIndex>(file_cache_).emplace(id, SequenceView{id_view, seq_view, qual_view});
    } else {
        std::get<NGSIndex>(file_cache_).emplace(hash_key(id), SequenceView{id_view, seq_view, qual_view});
    }
    data_ready_.store(true, std::memory_order_release);
}

// --- GZIP Loader ---
void SmartStrategy::loadGzipInternal(const std::string& filepath) {
    gzFile file = gzopen(filepath.c_str(), "rb");
    if (!file) throw std::runtime_error("Cannot open GZIP file: " + filepath);

    const size_t CHUNK_SIZE = 1024 * 1024; // 1MB
    std::vector<char> temp_buffer;
    temp_buffer.reserve(CHUNK_SIZE * 4); // modest initial reservation; grows geometrically below

    auto chunk = std::make_unique<char[]>(CHUNK_SIZE);
    while (true) {
        int bytes_read = gzread(file, chunk.get(), CHUNK_SIZE);
        if (bytes_read < 0) {
            int err;
            const char* error_msg = gzerror(file, &err);
            gzclose(file);
            throw std::runtime_error("GZIP read error: " + std::string(error_msg));
        }
        if (bytes_read == 0) break;
        // Geometric growth: double capacity when needed to minimise reallocations.
        if (temp_buffer.size() + bytes_read > temp_buffer.capacity()) {
            temp_buffer.reserve(std::max(temp_buffer.capacity() * 2, temp_buffer.size() + bytes_read));
        }
        size_t old_size = temp_buffer.size();
        temp_buffer.resize(old_size + bytes_read);
        std::memcpy(temp_buffer.data() + old_size, chunk.get(), bytes_read);
    }
    gzclose(file);
    text_arena_ = std::move(temp_buffer);
}

void SmartStrategy::normalizeFastaArena() {
    if (text_arena_.empty()) return;
    char* write = text_arena_.data();
    const char* read = text_arena_.data();
    const char* end = read + text_arena_.size();

    // Skip any leading blank lines
    while (read < end && (*read == '\n' || *read == '\r')) ++read;

    while (read < end) {
        if (*read != '>') break; // Not FASTA, stop
        // Copy header line verbatim, including the terminating '\n'
        while (read < end && *read != '\n' && *read != '\r') *write++ = *read++;
        if (read < end && *read == '\r') ++read; // strip \r
        if (read < end && *read == '\n') *write++ = *read++; // keep \n

        // Copy sequence chars, stripping all newlines, until the next '>' or end
        while (read < end && *read != '>') {
            char c = *read++;
            if (c != '\n' && c != '\r') *write++ = c;
        }
        *write++ = '\n'; // terminate the (now single-line) sequence
    }
    text_arena_.resize(static_cast<size_t>(write - text_arena_.data()));
}

void SmartStrategy::parseArena() {
    std::string_view content(text_arena_.data(), text_arena_.size());
    if (content.empty()) {
        data_ready_.store(true, std::memory_order_release);
        return;
    }
    bool isFastq = (content[0] == '@');
    bool use_ngs_mode = false; // Force Hybrid
    const size_t MULTITHREAD_THRESHOLD = 10 * 1024 * 1024;

    // Normalize multi-line FASTA in-place before any parser runs.
    // Ensures sequence string_views never contain embedded newlines.
    if (!isFastq) {
        normalizeFastaArena();
        content = std::string_view(text_arena_.data(), text_arena_.size());
    }

    if (use_ngs_mode) {
        file_cache_ = NGSIndex{};
        auto& map = std::get<NGSIndex>(file_cache_);
        if (content.size() > MULTITHREAD_THRESHOLD) {
            if (isFastq) parseFastqMultithreadedTemplate(content, map);
            else parseFastaMultithreadedTemplate(content, map);
        } else {
            if (isFastq) parseFastqInternal(content, map);
            else parseFastaInternal(content, map);
        }
    } else {
        file_cache_ = GenomeIndex{};
        auto& map = std::get<GenomeIndex>(file_cache_);
        if (content.size() > MULTITHREAD_THRESHOLD) {
            if (isFastq) parseFastqMultithreadedTemplate(content, map);
            else parseFastaMultithreadedTemplate(content, map);
        } else {
            if (isFastq) parseFastqInternal(content, map);
            else parseFastaInternal(content, map);
        }
    }
    determine_format_from_cache();
    data_ready_.store(true, std::memory_order_release);
}

void SmartStrategy::loadGzipFile(const std::string& filepath) {
    std::unique_lock<std::shared_mutex> lock(cache_mutex_);
    clearInternal();
    loadGzipInternal(filepath);
    parseArena();
}

void SmartStrategy::loadFile(const std::string& filepath) {
    bool is_gzip = false;
    if (filepath.size() > 3 && filepath.substr(filepath.size() - 3) == ".gz") is_gzip = true;
    if (!is_gzip) {
        std::ifstream check(filepath, std::ios::binary);
        if (check) {
            unsigned char magic[2];
            check.read(reinterpret_cast<char*>(magic), 2);
            if (check.gcount() == 2 && magic[0] == 0x1f && magic[1] == 0x8b) is_gzip = true;
        }
    }
    if (is_gzip) { loadGzipFile(filepath); return; }

    std::unique_lock<std::shared_mutex> lock(cache_mutex_);
    clearInternal();
    std::ifstream file(filepath, std::ios::binary | std::ios::ate);
    if (!file) throw std::runtime_error("Cannot open file: " + filepath);
    size_t fileSize = file.tellg();
    file.seekg(0, std::ios::beg);
    text_arena_.resize(fileSize);
    if (!file.read(text_arena_.data(), fileSize)) throw std::runtime_error("Read failed");
    parseArena();
}

// --- Templated Parsers ---
template <typename MapType> void SmartStrategy::parseFastaMultithreadedTemplate(std::string_view content, MapType& dest_map) {
    const size_t num_threads = std::max(1u, std::thread::hardware_concurrency());
    const size_t chunk_size = content.size() / num_threads;
    std::vector<std::thread> threads;
    std::vector<MapType> thread_caches(num_threads);
    auto worker = [&](size_t thread_id, size_t start, size_t end) {
        const char* ptr = content.data() + start;
        const char* chunk_end = content.data() + end;
        const char* global_end = content.data() + content.size();
        if (thread_id > 0) { while (ptr > content.data() && *(ptr - 1) != '\n') --ptr; while (ptr < chunk_end && *ptr != '>') ++ptr; }
        while (ptr < chunk_end) {
            while (ptr < chunk_end && (*ptr == '\n' || *ptr == '\r' || *ptr == ' ')) ++ptr;
            if (ptr >= chunk_end || *ptr != '>') break;
            ++ptr; const char* id_start = ptr;
            while (ptr < global_end && *ptr != ' ' && *ptr != '\t' && *ptr != '\n' && *ptr != '\r') ++ptr;
            const char* id_end = ptr;
            while (ptr < global_end && *ptr != '\n' && *ptr != '\r') ++ptr;
            while (ptr < global_end && (*ptr == '\n' || *ptr == '\r')) ++ptr;
            const char* seq_start = ptr; const char* seq_end = seq_start;
            while (ptr < global_end) {
                if (*ptr == '>') { seq_end = ptr; while (seq_end > seq_start && (*(seq_end-1) == '\n' || *(seq_end-1) == '\r' || *(seq_end-1) == ' ')) --seq_end; break; }
                ++ptr;
            }
            if (ptr >= global_end) { seq_end = global_end; while (seq_end > seq_start && (*(seq_end-1) == '\n' || *(seq_end-1) == '\r' || *(seq_end-1) == ' ')) --seq_end; }
            std::string_view id(id_start, id_end - id_start);
            std::string_view seq(seq_start, seq_end - seq_start);
            if (!id.empty() && !seq.empty()) {
                if constexpr (std::is_same_v<MapType, GenomeIndex>) { thread_caches[thread_id].emplace(std::string(id), SequenceView{id, seq, {}}); }
                else { thread_caches[thread_id].emplace(std::hash<std::string_view>{}(id), SequenceView{id, seq, {}}); }
            }
        }
    };
    for (size_t i = 0; i < num_threads; ++i) { size_t start = i * chunk_size; size_t end = (i == num_threads - 1) ? content.size() : (i + 1) * chunk_size; threads.emplace_back(worker, i, start, end); }
    for (auto& t : threads) t.join();
    size_t total_size = 0;
    for (const auto& cache : thread_caches) total_size += cache.size();
    dest_map.reserve(total_size);
    for (auto& cache : thread_caches) { for (auto& entry : cache) dest_map.insert(std::move(entry)); }
}

template <typename MapType> void SmartStrategy::parseFastqMultithreadedTemplate(std::string_view content, MapType& dest_map) {
    const size_t num_threads = std::max(1u, std::thread::hardware_concurrency());
    const size_t chunk_size = content.size() / num_threads;
    std::vector<std::thread> threads;
    std::vector<MapType> thread_caches(num_threads);
    auto worker = [&](size_t thread_id, size_t start, size_t end) {
        const char* ptr = content.data() + start;
        const char* chunk_end = content.data() + end;
        const char* global_end = content.data() + content.size();
        if (thread_id > 0) { while (ptr > content.data() && *(ptr - 1) != '\n') --ptr; while (ptr < chunk_end) { if (*ptr == '@') { const char* check = ptr - 1; while (check > content.data() && (*check == '\n' || *check == '\r')) --check; while (check > content.data() && *check != '\n' && *check != '\r') --check; if (check > content.data()) ++check; if (*check != '+') break; } ++ptr; } }
        while (ptr < chunk_end) {
            while (ptr < chunk_end && (*ptr == '\n' || *ptr == '\r' || *ptr == ' ')) ++ptr;
            if (ptr >= chunk_end || *ptr != '@') break;
            ++ptr; const char* id_start = ptr;
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
                if constexpr (std::is_same_v<MapType, GenomeIndex>) { thread_caches[thread_id].emplace(std::string(id), SequenceView{id, seq, qual}); }
                else { thread_caches[thread_id].emplace(std::hash<std::string_view>{}(id), SequenceView{id, seq, qual}); }
            }
        }
    };
    for (size_t i = 0; i < num_threads; ++i) { size_t start = i * chunk_size; size_t end = (i == num_threads - 1) ? content.size() : (i + 1) * chunk_size; threads.emplace_back(worker, i, start, end); }
    for (auto& t : threads) t.join();
    size_t total_size = 0;
    for (const auto& cache : thread_caches) total_size += cache.size();
    dest_map.reserve(total_size);
    for (auto& cache : thread_caches) { for (auto& entry : cache) dest_map.insert(std::move(entry)); }
}

template <typename MapType> void SmartStrategy::parseFastaInternal(std::string_view content, MapType& map) {
    const char* ptr = content.data(); const char* end = ptr + content.size();
    size_t estimated_records = (content.size() < 50 * 1024 * 1024) ? content.size() / 100 : content.size() / 200;
    map.reserve(static_cast<size_t>(estimated_records * 1.25));
    while (ptr < end) {
        while (ptr < end && (*ptr == '\n' || *ptr == '\r' || *ptr == ' ')) ++ptr;
        if (ptr >= end || *ptr != '>') break;
        ++ptr; const char* id_start = ptr; while (ptr < end && *ptr != ' ' && *ptr != '\t' && *ptr != '\n' && *ptr != '\r') ++ptr; const char* id_end = ptr;
        while (ptr < end && *ptr != '\n' && *ptr != '\r') ++ptr; while (ptr < end && (*ptr == '\n' || *ptr == '\r')) ++ptr;
        const char* seq_start = ptr; const char* seq_end = seq_start;
        while (ptr < end) { if (*ptr == '>') { seq_end = ptr; while (seq_end > seq_start && (*(seq_end-1) == '\n' || *(seq_end-1) == '\r' || *(seq_end-1) == ' ')) --seq_end; break; } ++ptr; }
        if (ptr >= end) { seq_end = end; while (seq_end > seq_start && (*(seq_end-1) == '\n' || *(seq_end-1) == '\r' || *(seq_end-1) == ' ')) --seq_end; }
        std::string_view id(id_start, id_end - id_start); std::string_view seq(seq_start, seq_end - seq_start);
        if (!id.empty() && !seq.empty()) {
            if constexpr (std::is_same_v<MapType, GenomeIndex>) { map.emplace(std::string(id), SequenceView{id, seq, {}}); }
            else { map.emplace(hash_key(id), SequenceView{id, seq, {}}); }
        }
    }
}

template <typename MapType> void SmartStrategy::parseFastqInternal(std::string_view content, MapType& map) {
    const char* ptr = content.data(); const char* end = ptr + content.size();
    size_t estimated_records = (content.size() < 50 * 1024 * 1024) ? content.size() / 150 : content.size() / 200;
    map.reserve(static_cast<size_t>(estimated_records * 1.25));
    while (ptr < end) {
        while (ptr < end && (*ptr == '\n' || *ptr == '\r' || *ptr == ' ')) ++ptr;
        if (ptr >= end || *ptr != '@') break;
        ++ptr; const char* id_start = ptr; while (ptr < end && *ptr != ' ' && *ptr != '\t' && *ptr != '\n' && *ptr != '\r') ++ptr; const char* id_end = ptr;
        while (ptr < end && *ptr != '\n' && *ptr != '\r') ++ptr; while (ptr < end && (*ptr == '\n' || *ptr == '\r')) ++ptr;
        const char* seq_start = ptr; while (ptr < end && *ptr != '\n' && *ptr != '\r') ++ptr; const char* seq_end = ptr;
        while (ptr < end && (*ptr == '\n' || *ptr == '\r')) ++ptr; while (ptr < end && *ptr != '\n' && *ptr != '\r') ++ptr; while (ptr < end && (*ptr == '\n' || *ptr == '\r')) ++ptr;
        const char* qual_start = ptr; while (ptr < end && *ptr != '\n' && *ptr != '\r') ++ptr; const char* qual_end = ptr; while (ptr < end && (*ptr == '\n' || *ptr == '\r')) ++ptr;
        std::string_view id(id_start, id_end - id_start); std::string_view seq(seq_start, seq_end - seq_start); std::string_view qual(qual_start, qual_end - qual_start);
        if (!id.empty() && !seq.empty()) {
            if constexpr (std::is_same_v<MapType, GenomeIndex>) { map.emplace(std::string(id), SequenceView{id, seq, qual}); }
            else { map.emplace(hash_key(id), SequenceView{id, seq, qual}); }
        }
    }
}

void SmartStrategy::saveBinary(const std::string& filepath) const {
    std::shared_lock<std::shared_mutex> lock(cache_mutex_);
    std::ofstream out(filepath, std::ios::binary);
    if (!out) throw std::runtime_error("Cannot write: " + filepath);
    out.write(MAGIC_BYTES, 4);
    uint8_t mode = std::holds_alternative<NGSIndex>(file_cache_) ? 1 : 0;
    out.write(reinterpret_cast<const char*>(&mode), 1);
    if (mode == 0) {
        const auto& map = std::get<GenomeIndex>(file_cache_);
        uint64_t count = map.size();
        out.write(reinterpret_cast<const char*>(&count), sizeof(count));
        for (const auto& [key, view] : map) {
            uint32_t len = static_cast<uint32_t>(key.size());
            out.write(reinterpret_cast<const char*>(&len), sizeof(len));
            out.write(key.data(), len);
            len = static_cast<uint32_t>(view.sequence.size());
            out.write(reinterpret_cast<const char*>(&len), sizeof(len));
            out.write(view.sequence.data(), len);
            len = static_cast<uint32_t>(view.quality.size());
            out.write(reinterpret_cast<const char*>(&len), sizeof(len));
            if (len > 0) out.write(view.quality.data(), len);
        }
    } else {
        const auto& map = std::get<NGSIndex>(file_cache_);
        uint64_t count = map.size();
        out.write(reinterpret_cast<const char*>(&count), sizeof(count));
        for (const auto& [key, view] : map) {
            out.write(reinterpret_cast<const char*>(&key), sizeof(key));
            uint32_t len = static_cast<uint32_t>(view.id.size());
            out.write(reinterpret_cast<const char*>(&len), sizeof(len));
            out.write(view.id.data(), len);
            len = static_cast<uint32_t>(view.sequence.size());
            out.write(reinterpret_cast<const char*>(&len), sizeof(len));
            out.write(view.sequence.data(), len);
            len = static_cast<uint32_t>(view.quality.size());
            out.write(reinterpret_cast<const char*>(&len), sizeof(len));
            if (len > 0) out.write(view.quality.data(), len);
        }
    }
}

void SmartStrategy::loadBinary(const std::string& filepath) {
    std::unique_lock<std::shared_mutex> lock(cache_mutex_);
    clearInternal();
    mmap_handle_ = std::make_unique<MMapHandle>();
#ifdef _WIN32
    mmap_handle_->hFile = CreateFileA(filepath.c_str(), GENERIC_READ, FILE_SHARE_READ,
                                       NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
    if (mmap_handle_->hFile == INVALID_HANDLE_VALUE)
        throw std::runtime_error("Cannot open binary cache: " + filepath);
    LARGE_INTEGER file_size{};
    if (!GetFileSizeEx(mmap_handle_->hFile, &file_size))
        throw std::runtime_error("GetFileSizeEx failed: " + filepath);
    mmap_handle_->size = static_cast<size_t>(file_size.QuadPart);
    mmap_handle_->hMap = CreateFileMappingA(mmap_handle_->hFile, NULL, PAGE_READONLY, 0, 0, NULL);
    if (!mmap_handle_->hMap)
        throw std::runtime_error("CreateFileMapping failed: " + filepath);
    mmap_handle_->data = MapViewOfFile(mmap_handle_->hMap, FILE_MAP_READ, 0, 0, 0);
    if (!mmap_handle_->data)
        throw std::runtime_error("MapViewOfFile failed: " + filepath);
#else
    mmap_handle_->fd = open(filepath.c_str(), O_RDONLY);
    if (mmap_handle_->fd == -1)
        throw std::runtime_error("Cannot open binary cache: " + filepath);
    struct stat sb{};
    if (fstat(mmap_handle_->fd, &sb) == -1)
        throw std::runtime_error("fstat failed: " + filepath);
    mmap_handle_->size = static_cast<size_t>(sb.st_size);
    mmap_handle_->data = mmap(NULL, mmap_handle_->size, PROT_READ, MAP_PRIVATE, mmap_handle_->fd, 0);
    if (mmap_handle_->data == MAP_FAILED)
        throw std::runtime_error("mmap failed: " + filepath);
#endif

    const char* ptr = static_cast<const char*>(mmap_handle_->data);
    const char* end = ptr + mmap_handle_->size;

    // Bounds-checked read helper. Uses memcpy to avoid unaligned-access UB on ARM.
    auto safe_advance = [&](size_t n) -> const char* {
        if (ptr + n > end)
            throw std::runtime_error("Binary cache is truncated or corrupt: " + filepath);
        const char* result = ptr;
        ptr += n;
        return result;
    };

    // Validate magic / format version
    if (mmap_handle_->size < 5 || std::strncmp(safe_advance(4), MAGIC_BYTES, 4) != 0) {
        throw std::runtime_error(
            "Invalid or outdated binary cache (re-generate with save()): " + filepath);
    }

    uint8_t mode;
    std::memcpy(&mode, safe_advance(1), 1);
    uint64_t count;
    std::memcpy(&count, safe_advance(sizeof(uint64_t)), sizeof(uint64_t));

    // Sanity-cap count to avoid looping on a corrupt large value.
    constexpr uint64_t MAX_RECORDS = 1'000'000'000ULL;
    if (count > MAX_RECORDS)
        throw std::runtime_error("Binary cache record count implausible, file may be corrupt: " + filepath);

    if (mode == 0) {
        file_cache_ = GenomeIndex{};
        auto& map = std::get<GenomeIndex>(file_cache_);
        map.reserve(count);
        for (uint64_t i = 0; i < count; ++i) {
            uint32_t len; std::memcpy(&len, safe_advance(4), 4);
            std::string_view id_view(safe_advance(len), len);
            uint32_t seq_len; std::memcpy(&seq_len, safe_advance(4), 4);
            std::string_view seq(safe_advance(seq_len), seq_len);
            uint32_t qual_len; std::memcpy(&qual_len, safe_advance(4), 4);
            std::string_view qual;
            if (qual_len > 0) qual = std::string_view(safe_advance(qual_len), qual_len);
            map.emplace(std::string(id_view), SequenceView{id_view, seq, qual});
        }
    } else {
        file_cache_ = NGSIndex{};
        auto& map = std::get<NGSIndex>(file_cache_);
        map.reserve(count);
        for (uint64_t i = 0; i < count; ++i) {
            uint64_t hash; std::memcpy(&hash, safe_advance(8), 8);
            uint32_t len; std::memcpy(&len, safe_advance(4), 4);
            std::string_view id_view(safe_advance(len), len);
            uint32_t seq_len; std::memcpy(&seq_len, safe_advance(4), 4);
            std::string_view seq(safe_advance(seq_len), seq_len);
            uint32_t qual_len; std::memcpy(&qual_len, safe_advance(4), 4);
            std::string_view qual;
            if (qual_len > 0) qual = std::string_view(safe_advance(qual_len), qual_len);
            map.emplace(hash, SequenceView{id_view, seq, qual});
        }
    }
    determine_format_from_cache();
    data_ready_.store(true, std::memory_order_release);
}

// --- ACCESS ---

std::string_view SmartStrategy::getView(const std::string_view& sequence_id) const {
    if (data_ready_.load(std::memory_order_acquire)) {
        if (std::holds_alternative<GenomeIndex>(file_cache_)) {
            const auto& map = std::get<GenomeIndex>(file_cache_);
            auto it = map.find(sequence_id); // heterogeneous lookup: no std::string allocation
            if (it != map.end()) return it->second.sequence;
        } else {
            const auto& map = std::get<NGSIndex>(file_cache_);
            uint64_t h = hash_key(sequence_id); 
            auto it = map.find(h);
            if (it != map.end() && it->second.id == sequence_id) return it->second.sequence;
        }
        return {};
    }
    return {};
}

std::string SmartStrategy::getSequence(const std::string& sequence_id) const { return std::string(getView(sequence_id)); }

std::string SmartStrategy::getQuality(const std::string& sequence_id) const { 
    if (data_ready_.load(std::memory_order_acquire)) {
        if (std::holds_alternative<GenomeIndex>(file_cache_)) {
            auto& map = std::get<GenomeIndex>(file_cache_);
            auto it = map.find(sequence_id);
            return (it != map.end()) ? std::string(it->second.quality) : "";
        } else {
            auto& map = std::get<NGSIndex>(file_cache_);
            auto it = map.find(hash_key(sequence_id));
            return (it != map.end() && it->second.id == sequence_id) ? std::string(it->second.quality) : "";
        }
    }
    return "";
}

bool SmartStrategy::hasSequence(const std::string& sequence_id) const {
    if (data_ready_.load(std::memory_order_acquire)) {
        if (std::holds_alternative<GenomeIndex>(file_cache_)) {
            return std::get<GenomeIndex>(file_cache_).count(sequence_id) > 0;
        } else {
            auto& map = std::get<NGSIndex>(file_cache_);
            auto it = map.find(hash_key(sequence_id));
            return (it != map.end() && it->second.id == sequence_id);
        }
    }
    return false;
}

size_t SmartStrategy::getFileCacheSize() const {
    std::shared_lock<std::shared_mutex> lock(cache_mutex_);
    if (std::holds_alternative<GenomeIndex>(file_cache_)) return std::get<GenomeIndex>(file_cache_).size();
    return std::get<NGSIndex>(file_cache_).size();
}

std::vector<std::string> SmartStrategy::getAllKeys() const {
    std::shared_lock<std::shared_mutex> lock(cache_mutex_);
    std::vector<std::string> keys;
    if (std::holds_alternative<GenomeIndex>(file_cache_)) {
        auto& map = std::get<GenomeIndex>(file_cache_);
        keys.reserve(map.size());
        for (const auto& pair : map) keys.emplace_back(pair.first);
    } else {
        auto& map = std::get<NGSIndex>(file_cache_);
        keys.reserve(map.size());
        for (const auto& pair : map) keys.emplace_back(std::string(pair.second.id));
    }
    return keys;
}

IndexMode SmartStrategy::getIndexMode() const {
    return std::holds_alternative<GenomeIndex>(file_cache_) ? IndexMode::GENOME : IndexMode::NGS;
}

void SmartStrategy::determine_format_from_cache() {
    bool empty = false;
    SequenceView first_view;
    if (std::holds_alternative<GenomeIndex>(file_cache_)) {
        const auto& map = std::get<GenomeIndex>(file_cache_);
        if (map.empty()) empty = true;
        else first_view = map.begin()->second;
    } else {
        const auto& map = std::get<NGSIndex>(file_cache_);
        if (map.empty()) empty = true;
        else first_view = map.begin()->second;
    }

    if (empty) { detected_format_ = FileFormat::UNKNOWN; return; }
    bool is_fastq = !first_view.quality.empty();
    bool is_rna = hasRNA(first_view.sequence);
    bool is_nuc = isNucleotideSequence(first_view.sequence);
    if (!is_fastq) {
        detected_format_ = is_rna ? FileFormat::RNA_FASTA : is_nuc ? FileFormat::DNA_FASTA : FileFormat::PROTEIN_FASTA;
    } else {
        detected_format_ = is_rna ? FileFormat::RNA_FASTQ : is_nuc ? FileFormat::DNA_FASTQ : FileFormat::PROTEIN_FASTQ;
    }
}

bool SmartStrategy::isNucleotideSequence(std::string_view data) const {
    if (data.empty()) return false;
    size_t nucleotide_count = 0; size_t total_count = 0;
    size_t scan_limit = std::min(data.size(), (size_t)1000);
    for (size_t i = 0; i < scan_limit; ++i) {
        char c = data[i];
        if (std::isalpha(c)) {
            total_count++;
            char upper_c = std::toupper(c);
            if (upper_c == 'A' || upper_c == 'T' || upper_c == 'G' || upper_c == 'C' || upper_c == 'U' || upper_c == 'N') nucleotide_count++;
        }
    }
    return total_count > 0 && (static_cast<double>(nucleotide_count) / total_count) > 0.8;
}

bool SmartStrategy::hasRNA(std::string_view data) const {
    size_t scan_limit = std::min(data.size(), (size_t)1000);
    for (size_t i = 0; i < scan_limit; ++i) {
        if (data[i] == 'U' || data[i] == 'u') return true;
    }
    return false;
}

// Dummy parsers
void SmartStrategy::parseFastaOptimized(std::string_view content) { 
    if (std::holds_alternative<GenomeIndex>(file_cache_)) parseFastaInternal(content, std::get<GenomeIndex>(file_cache_));
    else parseFastaInternal(content, std::get<NGSIndex>(file_cache_));
}
void SmartStrategy::parseFastqOptimized(std::string_view content) {
    if (std::holds_alternative<GenomeIndex>(file_cache_)) parseFastqInternal(content, std::get<GenomeIndex>(file_cache_));
    else parseFastqInternal(content, std::get<NGSIndex>(file_cache_));
}
void SmartStrategy::parseFastaMultithreaded(std::string_view content) {
    if (std::holds_alternative<GenomeIndex>(file_cache_)) parseFastaMultithreadedTemplate(content, std::get<GenomeIndex>(file_cache_));
    else parseFastaMultithreadedTemplate(content, std::get<NGSIndex>(file_cache_));
}
void SmartStrategy::parseFastqMultithreaded(std::string_view content) {
    if (std::holds_alternative<GenomeIndex>(file_cache_)) parseFastqMultithreadedTemplate(content, std::get<GenomeIndex>(file_cache_));
    else parseFastqMultithreadedTemplate(content, std::get<NGSIndex>(file_cache_));
}

} // namespace