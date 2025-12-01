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

// --- Platform Specific Includes for MMAP ---
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

// --- MMapHandle Implementation ---

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

    ~MMapHandle() {
        cleanup();
    }

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

// --- Constants ---
static const char MAGIC_BYTES[] = "MMAP";

// --- Constructor / Destructor ---

SmartStrategy::SmartStrategy() : detected_format_(FileFormat::UNKNOWN) {}
SmartStrategy::~SmartStrategy() {
    clearCache(); // Ensure views are cleared before destroying arena/mmap
}

// --- IEncodingStrategy (Legacy/Pass-through) ---

std::vector<unsigned char> SmartStrategy::encode(const std::string& data, DataTypeHint hint) const {
    return {data.begin(), data.end()};
}

std::string SmartStrategy::decode(const std::vector<unsigned char>& data) const {
    return {data.begin(), data.end()};
}

// --- Multithreaded FASTA Parser ---

void SmartStrategy::parseFastaMultithreaded(std::string_view content) {
    const size_t num_threads = std::thread::hardware_concurrency();
    const size_t chunk_size = content.size() / num_threads;
    
    if (chunk_size < 1024 * 1024) { // If chunks too small (< 1MB), use single-threaded
        parseFastaOptimized(content);
        return;
    }
    
    std::vector<std::thread> threads;
    std::vector<HashMap<std::string, SequenceView>> thread_caches(num_threads);
    
    auto worker = [&](size_t thread_id, size_t start, size_t end) {
        const char* ptr = content.data() + start;
        const char* chunk_end = content.data() + end;
        const char* global_end = content.data() + content.size();
        
        // Align to record boundary: scan backward to find last '>' before our chunk
        if (thread_id > 0) {
            while (ptr > content.data() && *(ptr - 1) != '\n') {
                --ptr;
            }
            // Now find the next '>'
            while (ptr < chunk_end && *ptr != '>') {
                ++ptr;
            }
        }
        
        // Parse records in this chunk
        while (ptr < chunk_end) {
            // Skip whitespace
            while (ptr < chunk_end && (*ptr == '\n' || *ptr == '\r' || *ptr == ' ')) {
                ++ptr;
            }
            
            if (ptr >= chunk_end || *ptr != '>') break;
            
            // Parse header
            ++ptr; // Skip '>'
            const char* id_start = ptr;
            
            while (ptr < global_end && *ptr != ' ' && *ptr != '\t' && *ptr != '\n' && *ptr != '\r') {
                ++ptr;
            }
            const char* id_end = ptr;
            
            // Skip rest of header
            while (ptr < global_end && *ptr != '\n' && *ptr != '\r') {
                ++ptr;
            }
            while (ptr < global_end && (*ptr == '\n' || *ptr == '\r')) {
                ++ptr;
            }
            
            // Collect sequence
            const char* seq_start = ptr;
            const char* seq_end = seq_start;
            
            while (ptr < global_end) {
                if (*ptr == '>') {
                    seq_end = ptr;
                    while (seq_end > seq_start && (*(seq_end-1) == '\n' || *(seq_end-1) == '\r' || *(seq_end-1) == ' ')) {
                        --seq_end;
                    }
                    break;
                }
                ++ptr;
            }
            
            if (ptr >= global_end) {
                seq_end = global_end;
                while (seq_end > seq_start && (*(seq_end-1) == '\n' || *(seq_end-1) == '\r' || *(seq_end-1) == ' ')) {
                    --seq_end;
                }
            }
            
            // Store in thread-local cache
            std::string_view id(id_start, id_end - id_start);
            std::string_view seq(seq_start, seq_end - seq_start);
            
            if (!id.empty() && !seq.empty()) {
                thread_caches[thread_id].emplace(std::string(id), SequenceView{id, seq, {}});
            }
        }
    };
    
    // Launch threads
    for (size_t i = 0; i < num_threads; ++i) {
        size_t start = i * chunk_size;
        size_t end = (i == num_threads - 1) ? content.size() : (i + 1) * chunk_size;
        threads.emplace_back(worker, i, start, end);
    }
    
    // Wait for completion
    for (auto& t : threads) {
        t.join();
    }
    
    // Merge thread-local caches
    size_t total_size = 0;
    for (const auto& cache : thread_caches) {
        total_size += cache.size();
    }
    file_cache_.reserve(total_size);
    
    for (auto& cache : thread_caches) {
        for (auto& entry : cache) {
            file_cache_.insert(std::move(entry));
        }
    }
}

// --- Multithreaded FASTQ Parser ---

void SmartStrategy::parseFastqMultithreaded(std::string_view content) {
    const size_t num_threads = std::thread::hardware_concurrency();
    const size_t chunk_size = content.size() / num_threads;
    
    if (chunk_size < 1024 * 1024) { // If chunks too small (< 1MB), use single-threaded
        parseFastqOptimized(content);
        return;
    }
    
    std::vector<std::thread> threads;
    std::vector<HashMap<std::string, SequenceView>> thread_caches(num_threads);
    
    auto worker = [&](size_t thread_id, size_t start, size_t end) {
        const char* ptr = content.data() + start;
        const char* chunk_end = content.data() + end;
        const char* global_end = content.data() + content.size();
        
        // Align to record boundary: FASTQ records are 4 lines, so find next '@' at line start
        if (thread_id > 0) {
            while (ptr > content.data() && *(ptr - 1) != '\n') {
                --ptr;
            }
            // Scan forward to find '@' at beginning of line
            while (ptr < chunk_end) {
                if (*ptr == '@') {
                    // Verify it's not quality line by checking if previous line was '+'
                    const char* check = ptr - 1;
                    while (check > content.data() && (*check == '\n' || *check == '\r')) --check;
                    while (check > content.data() && *check != '\n' && *check != '\r') --check;
                    if (check > content.data()) ++check;
                    
                    // If previous line doesn't start with '+', this is a valid header
                    if (*check != '+') break;
                }
                ++ptr;
            }
        }
        
        // Parse FASTQ records
        while (ptr < chunk_end) {
            // Skip whitespace
            while (ptr < chunk_end && (*ptr == '\n' || *ptr == '\r' || *ptr == ' ')) {
                ++ptr;
            }
            
            if (ptr >= chunk_end || *ptr != '@') break;
            
            // Line 1: @Header
            ++ptr;
            const char* id_start = ptr;
            
            while (ptr < global_end && *ptr != ' ' && *ptr != '\t' && *ptr != '\n' && *ptr != '\r') {
                ++ptr;
            }
            const char* id_end = ptr;
            
            while (ptr < global_end && *ptr != '\n' && *ptr != '\r') {
                ++ptr;
            }
            while (ptr < global_end && (*ptr == '\n' || *ptr == '\r')) {
                ++ptr;
            }
            
            // Line 2: Sequence
            const char* seq_start = ptr;
            while (ptr < global_end && *ptr != '\n' && *ptr != '\r') {
                ++ptr;
            }
            const char* seq_end = ptr;
            
            while (ptr < global_end && (*ptr == '\n' || *ptr == '\r')) {
                ++ptr;
            }
            
            // Line 3: +
            while (ptr < global_end && *ptr != '\n' && *ptr != '\r') {
                ++ptr;
            }
            while (ptr < global_end && (*ptr == '\n' || *ptr == '\r')) {
                ++ptr;
            }
            
            // Line 4: Quality
            const char* qual_start = ptr;
            while (ptr < global_end && *ptr != '\n' && *ptr != '\r') {
                ++ptr;
            }
            const char* qual_end = ptr;
            
            while (ptr < global_end && (*ptr == '\n' || *ptr == '\r')) {
                ++ptr;
            }
            
            // Store in thread-local cache
            std::string_view id(id_start, id_end - id_start);
            std::string_view seq(seq_start, seq_end - seq_start);
            std::string_view qual(qual_start, qual_end - qual_start);
            
            if (!id.empty() && !seq.empty()) {
                thread_caches[thread_id].emplace(std::string(id), SequenceView{id, seq, qual});
            }
        }
    };
    
    // Launch threads
    for (size_t i = 0; i < num_threads; ++i) {
        size_t start = i * chunk_size;
        size_t end = (i == num_threads - 1) ? content.size() : (i + 1) * chunk_size;
        threads.emplace_back(worker, i, start, end);
    }
    
    // Wait for completion
    for (auto& t : threads) {
        t.join();
    }
    
    // Merge thread-local caches
    size_t total_size = 0;
    for (const auto& cache : thread_caches) {
        total_size += cache.size();
    }
    file_cache_.reserve(total_size);
    
    for (auto& cache : thread_caches) {
        for (auto& entry : cache) {
            file_cache_.insert(std::move(entry));
        }
    }
}

// --- Optimized FASTA Parser (O(N) Single Pass) ---

void SmartStrategy::parseFastaOptimized(std::string_view content) {
    const char* ptr = content.data();
    const char* end = ptr + content.size();
    
    while (ptr < end) {
        // Skip whitespace/empty lines
        while (ptr < end && (*ptr == '\n' || *ptr == '\r' || *ptr == ' ')) {
            ++ptr;
        }
        
        if (ptr >= end || *ptr != '>') break;
        
        // Parse header line
        ++ptr; // Skip '>'
        const char* id_start = ptr;
        
        // ID is from here until first space or newline
        while (ptr < end && *ptr != ' ' && *ptr != '\t' && *ptr != '\n' && *ptr != '\r') {
            ++ptr;
        }
        const char* id_end = ptr;
        
        // Skip rest of header line
        while (ptr < end && *ptr != '\n' && *ptr != '\r') {
            ++ptr;
        }
        while (ptr < end && (*ptr == '\n' || *ptr == '\r')) {
            ++ptr;
        }
        
        // Collect sequence until next '>' or EOF
        const char* seq_start = ptr;
        const char* seq_end = seq_start;
        
        // Scan for next record marker or EOF
        while (ptr < end) {
            if (*ptr == '>') {
                // Found next record, backtrack to remove trailing whitespace
                seq_end = ptr;
                while (seq_end > seq_start && (*(seq_end-1) == '\n' || *(seq_end-1) == '\r' || *(seq_end-1) == ' ')) {
                    --seq_end;
                }
                break;
            }
            ++ptr;
        }
        
        if (ptr >= end) {
            // Reached EOF, trim trailing whitespace
            seq_end = end;
            while (seq_end > seq_start && (*(seq_end-1) == '\n' || *(seq_end-1) == '\r' || *(seq_end-1) == ' ')) {
                --seq_end;
            }
        }
        
        // Create views and insert
        std::string_view id(id_start, id_end - id_start);
        std::string_view seq(seq_start, seq_end - seq_start);
        
        if (!id.empty() && !seq.empty()) {
            // CRITICAL: Use std::string for key (owned), string_view for values (zero-copy)
            file_cache_.emplace(std::string(id), SequenceView{id, seq, {}});
        }
    }
}

// --- Optimized FASTQ Parser (O(N) Single Pass) ---

void SmartStrategy::parseFastqOptimized(std::string_view content) {
    const char* ptr = content.data();
    const char* end = ptr + content.size();
    
    while (ptr < end) {
        // Skip whitespace
        while (ptr < end && (*ptr == '\n' || *ptr == '\r' || *ptr == ' ')) {
            ++ptr;
        }
        
        if (ptr >= end || *ptr != '@') break;
        
        // Line 1: @Header
        ++ptr; // Skip '@'
        const char* id_start = ptr;
        
        while (ptr < end && *ptr != ' ' && *ptr != '\t' && *ptr != '\n' && *ptr != '\r') {
            ++ptr;
        }
        const char* id_end = ptr;
        
        // Skip rest of header
        while (ptr < end && *ptr != '\n' && *ptr != '\r') {
            ++ptr;
        }
        while (ptr < end && (*ptr == '\n' || *ptr == '\r')) {
            ++ptr;
        }
        
        // Line 2: Sequence
        const char* seq_start = ptr;
        while (ptr < end && *ptr != '\n' && *ptr != '\r') {
            ++ptr;
        }
        const char* seq_end = ptr;
        
        while (ptr < end && (*ptr == '\n' || *ptr == '\r')) {
            ++ptr;
        }
        
        // Line 3: + (skip entire line)
        while (ptr < end && *ptr != '\n' && *ptr != '\r') {
            ++ptr;
        }
        while (ptr < end && (*ptr == '\n' || *ptr == '\r')) {
            ++ptr;
        }
        
        // Line 4: Quality
        const char* qual_start = ptr;
        while (ptr < end && *ptr != '\n' && *ptr != '\r') {
            ++ptr;
        }
        const char* qual_end = ptr;
        
        while (ptr < end && (*ptr == '\n' || *ptr == '\r')) {
            ++ptr;
        }
        
        // Create views and insert
        std::string_view id(id_start, id_end - id_start);
        std::string_view seq(seq_start, seq_end - seq_start);
        std::string_view qual(qual_start, qual_end - qual_start);
        
        if (!id.empty() && !seq.empty()) {
            file_cache_.emplace(std::string(id), SequenceView{id, seq, qual});
        }
    }
}

// --- Core Logic: Text Loading (Arena) ---

void SmartStrategy::loadFile(const std::string& filepath) {
    std::unique_lock<std::shared_mutex> lock(cache_mutex_);
    clearCache();

    std::ifstream file(filepath, std::ios::binary | std::ios::ate);
    if (!file) throw std::runtime_error("Cannot open file: " + filepath);

    size_t fileSize = file.tellg();
    file.seekg(0, std::ios::beg);

    // 1. Allocate Arena
    text_arena_.resize(fileSize);
    
    // 2. Read entire file into Arena
    if (!file.read(text_arena_.data(), fileSize)) {
        throw std::runtime_error("Failed to read file into arena");
    }

    // 3. Parse In-Place (Zero Copy view creation)
    std::string_view content(text_arena_.data(), fileSize);
    
    // Determine format from first character
    if (content.empty()) return;
    
    bool isFastq = (content[0] == '@');
    
    // Reserve space for expected entries (heuristic: assume avg 100 bytes per record)
    file_cache_.reserve(fileSize / 100);
    
    // Use multithreaded parsers for large files, single-threaded for small
    const size_t MULTITHREAD_THRESHOLD = 10 * 1024 * 1024; // 10MB
    
    if (fileSize > MULTITHREAD_THRESHOLD) {
        if (isFastq) {
            parseFastqMultithreaded(content);
        } else {
            parseFastaMultithreaded(content);
        }
    } else {
        // For smaller files, single-threaded is faster (no overhead)
        if (isFastq) {
            parseFastqOptimized(content);
        } else {
            parseFastaOptimized(content);
        }
    }
    
    determine_format_from_cache();
}

// --- Core Logic: Binary Save (Serialization) ---

void SmartStrategy::saveBinary(const std::string& filepath) const {
    std::shared_lock<std::shared_mutex> lock(cache_mutex_);
    std::ofstream out(filepath, std::ios::binary);
    if (!out) throw std::runtime_error("Cannot create file: " + filepath);

    // 1. Header
    out.write(MAGIC_BYTES, 4);
    
    // 2. Size
    uint64_t count = file_cache_.size();
    out.write(reinterpret_cast<const char*>(&count), sizeof(count));

    // 3. Data Loop
    for (const auto& [key, view] : file_cache_) {
        uint32_t len;
        
        // ID
        len = static_cast<uint32_t>(view.id.size());
        out.write(reinterpret_cast<const char*>(&len), sizeof(len));
        out.write(view.id.data(), len);

        // Sequence
        len = static_cast<uint32_t>(view.sequence.size());
        out.write(reinterpret_cast<const char*>(&len), sizeof(len));
        out.write(view.sequence.data(), len);

        // Quality
        len = static_cast<uint32_t>(view.quality.size());
        out.write(reinterpret_cast<const char*>(&len), sizeof(len));
        if (len > 0) out.write(view.quality.data(), len);
    }
}

// --- Core Logic: Binary Load (MMAP / Zero-Copy) ---

void SmartStrategy::loadBinary(const std::string& filepath) {
    std::unique_lock<std::shared_mutex> lock(cache_mutex_);
    clearCache();
    
    mmap_handle_ = std::make_unique<MMapHandle>();

#ifdef _WIN32
    // Windows MMAP
    mmap_handle_->hFile = CreateFileA(filepath.c_str(), GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
    if (mmap_handle_->hFile == INVALID_HANDLE_VALUE) throw std::runtime_error("Win32: Open file failed");

    LARGE_INTEGER size;
    GetFileSizeEx(mmap_handle_->hFile, &size);
    mmap_handle_->size = static_cast<size_t>(size.QuadPart);

    mmap_handle_->hMap = CreateFileMappingA(mmap_handle_->hFile, NULL, PAGE_READONLY, 0, 0, NULL);
    if (!mmap_handle_->hMap) throw std::runtime_error("Win32: CreateFileMapping failed");

    mmap_handle_->data = MapViewOfFile(mmap_handle_->hMap, FILE_MAP_READ, 0, 0, 0);
    if (!mmap_handle_->data) throw std::runtime_error("Win32: MapViewOfFile failed");
#else
    // POSIX MMAP
    mmap_handle_->fd = open(filepath.c_str(), O_RDONLY);
    if (mmap_handle_->fd == -1) throw std::runtime_error("Posix: Open failed");

    struct stat sb;
    if (fstat(mmap_handle_->fd, &sb) == -1) throw std::runtime_error("Posix: fstat failed");
    mmap_handle_->size = sb.st_size;

    mmap_handle_->data = mmap(NULL, mmap_handle_->size, PROT_READ, MAP_PRIVATE, mmap_handle_->fd, 0);
    if (mmap_handle_->data == MAP_FAILED) throw std::runtime_error("Posix: mmap failed");
#endif

    // --- Pointer Arithmetic for Zero-Copy Reconstruction ---
    
    const char* ptr = static_cast<const char*>(mmap_handle_->data);
    const char* end = ptr + mmap_handle_->size;

    // Check Magic
    if (mmap_handle_->size < 4 || std::strncmp(ptr, MAGIC_BYTES, 4) != 0) {
        throw std::runtime_error("Invalid binary format: Missing Magic Bytes");
    }
    ptr += 4;

    // Read Count
    uint64_t count = *reinterpret_cast<const uint64_t*>(ptr);
    ptr += sizeof(uint64_t);

    file_cache_.reserve(count);

    // O(N) Loop
    for (uint64_t i = 0; i < count; ++i) {
        if (ptr >= end) break;

        // ID
        uint32_t id_len = *reinterpret_cast<const uint32_t*>(ptr);
        ptr += sizeof(uint32_t);
        std::string_view id_view(ptr, id_len);
        ptr += id_len;

        // Sequence
        uint32_t seq_len = *reinterpret_cast<const uint32_t*>(ptr);
        ptr += sizeof(uint32_t);
        std::string_view seq_view(ptr, seq_len);
        ptr += seq_len;

        // Quality
        uint32_t qual_len = *reinterpret_cast<const uint32_t*>(ptr);
        ptr += sizeof(uint32_t);
        std::string_view qual_view;
        if (qual_len > 0) {
            qual_view = std::string_view(ptr, qual_len);
            ptr += qual_len;
        }

        // Insert into map - KEY IS OWNED, VALUES ARE VIEWS
        file_cache_.emplace(std::string(id_view), SequenceView{id_view, seq_view, qual_view});
    }

    determine_format_from_cache();
}

// --- Utilities ---

void SmartStrategy::clearCache() {
    file_cache_.clear();
    text_arena_.clear();
    text_arena_.shrink_to_fit();
    mmap_handle_.reset(); // Unmaps memory
}

std::string_view SmartStrategy::getView(const std::string_view& sequence_id) const {
    std::shared_lock<std::shared_mutex> lock(cache_mutex_); // Read lock (allows concurrent reads)
    auto it = file_cache_.find(std::string(sequence_id)); // Temporary string for lookup
    if (it != file_cache_.end()) {
        return it->second.sequence;
    }
    return {};
}

std::string SmartStrategy::getSequence(const std::string& sequence_id) const {
    std::string_view v = getView(sequence_id);
    return std::string(v);
}

std::string SmartStrategy::getQuality(const std::string& sequence_id) const {
    std::shared_lock<std::shared_mutex> lock(cache_mutex_);
    auto it = file_cache_.find(sequence_id);
    if (it != file_cache_.end()) {
        return std::string(it->second.quality);
    }
    return "";
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
    for (const auto& pair : file_cache_) {
        keys.push_back(pair.first);
    }
    return keys;
}

void SmartStrategy::determine_format_from_cache() {
    if (file_cache_.empty()) {
        detected_format_ = FileFormat::UNKNOWN;
        return;
    }
    const auto& first = file_cache_.begin()->second;
    bool hasQ = !first.quality.empty();
    bool isRNA = hasRNA(first.sequence);
    
    if (hasQ) detected_format_ = isRNA ? FileFormat::RNA_FASTQ : FileFormat::DNA_FASTQ;
    else detected_format_ = isRNA ? FileFormat::RNA_FASTA : FileFormat::DNA_FASTA;
}

bool SmartStrategy::isNucleotideSequence(std::string_view data) const {
    if (data.empty()) return false;
    return true; 
}

bool SmartStrategy::hasRNA(std::string_view data) const {
    return data.find('U') != std::string_view::npos || data.find('u') != std::string_view::npos;
}

} // namespace TracEon