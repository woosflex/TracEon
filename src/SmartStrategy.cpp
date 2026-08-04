#include "SmartStrategy.h"
#include "SimdUtils.h"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <cstring>
#include <cstdio>
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
#include <lz4.h>
#include <lz4hc.h>
#include <lz4frame.h>
#include <new>
#include <atomic>

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

#ifdef __APPLE__
    #include <sys/sysctl.h>
    #include <sys/types.h>
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

static const char MAGIC_BYTES_V1[] = "TRO\x01"; // v1.1 format (uncompressed)
static const char MAGIC_BYTES_V2[] = "TRO\x02"; // v1.3 format (LZ4-compressed, single block)
static const char MAGIC_BYTES_V3[] = "TRO\x03"; // v1.5 format (LZ4 Frame, streamed)

// Streaming (de)compression chunk size for the v3 binary cache format.
// Bounds peak memory during saveBinary()/loadBinary() to a small constant
// regardless of dataset size (see ADR-004 follow-up: streaming binary cache).
static constexpr size_t STREAM_CHUNK_SIZE = 1024 * 1024; // 1MB

SmartStrategy::SmartStrategy(IndexMode mode)
    : detected_format_(FileFormat::UNKNOWN),
      index_mode_(mode),
      file_cache_(mode == IndexMode::NGS ? std::variant<GenomeIndex, NGSIndex>(NGSIndex{}) : std::variant<GenomeIndex, NGSIndex>(GenomeIndex{})) {}
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
    manual_store_bytes_.store(0, std::memory_order_relaxed);
    // sizeof(uint64_t): serializePayload()'s record-count header is always
    // present, even for an empty cache (see serialized_size_estimate_ decl).
    serialized_size_estimate_.store(sizeof(uint64_t), std::memory_order_relaxed);
    mmap_handle_.reset();
}

void SmartStrategy::clearCache() {
    std::unique_lock<std::shared_mutex> lock(cache_mutex_);
    clearInternal();
}

void SmartStrategy::addEntry(const std::string& id, const std::string& seq, const std::string& qual) {
    std::unique_lock<std::shared_mutex> lock(cache_mutex_);

    // OOM guard: manual_store_ is the only remaining unbounded growth path
    // in this class (every load path already guards its own allocations).
    // Fail loudly instead of growing without limit.
    const size_t entry_bytes = id.size() + seq.size() + qual.size();
    const size_t projected = manual_store_bytes_.load(std::memory_order_relaxed) + entry_bytes;
    const size_t avail_mem = getAvailableMemory();
    if (avail_mem > 0 && projected > avail_mem / 2) {
        throw std::runtime_error(
            "SmartStrategy OOM guard: addEntry() would grow manual_store_ to " +
            std::to_string(projected) + " bytes (available memory ~" +
            std::to_string(avail_mem) + " bytes)");
    }

    // Do NOT flip data_ready_ to false here. getView() checks data_ready_ without
    // holding a lock (lock-free fast path), so a false→true toggle would cause
    // concurrent readers to see empty results during the window. The unique_lock
    // already prevents map corruption; just emit true at the end.

    // Store copies in manual_store_ (std::deque: push_back never invalidates
    // references to existing elements, so these string_views stay valid).
    // The proactive check above uses system-wide available memory, which a
    // ulimit/cgroup-constrained process can still exceed sooner — catch
    // bad_alloc here too so that case throws the same descriptive error
    // instead of an unguarded std::bad_alloc.
    std::string_view id_view, seq_view, qual_view;
    try {
        manual_store_.push_back(id);
        id_view = manual_store_.back();
        manual_store_.push_back(seq);
        seq_view = manual_store_.back();
        manual_store_.push_back(qual);
        qual_view = manual_store_.back();
    } catch (const std::bad_alloc&) {
        std::cerr << "SmartStrategy OOM: failed to allocate " << entry_bytes
                  << " bytes for addEntry()" << std::endl;
        throw;
    }
    manual_store_bytes_.fetch_add(entry_bytes, std::memory_order_relaxed);

    if (std::holds_alternative<GenomeIndex>(file_cache_)) {
        std::get<GenomeIndex>(file_cache_).emplace(id, SequenceView{id_view, seq_view, qual_view});
        serialized_size_estimate_.fetch_add(
            3 * sizeof(uint32_t) + id.size() + seq.size() + qual.size(), std::memory_order_relaxed);
    } else {
        std::get<NGSIndex>(file_cache_).emplace(hash_key(id), SequenceView{id_view, seq_view, qual_view});
        serialized_size_estimate_.fetch_add(
            sizeof(uint64_t) + 3 * sizeof(uint32_t) + id.size() + seq.size() + qual.size(),
            std::memory_order_relaxed);
    }
    data_ready_.store(true, std::memory_order_release);
}

// ── Cross-platform available memory ───────────────────────────────────────────
size_t SmartStrategy::getAvailableMemory() {
#if defined(_WIN32)
    MEMORYSTATUSEX status;
    status.dwLength = sizeof(status);
    if (GlobalMemoryStatusEx(&status)) {
        return static_cast<size_t>(status.ullAvailPhys);
    }
    return 0;
#elif defined(__APPLE__)
    uint64_t free_pages = 0;
    uint64_t page_size = 0;
    size_t len = sizeof(free_pages);
    if (sysctlbyname("vm.page_free_count", &free_pages, &len, nullptr, 0) == 0) {
        len = sizeof(page_size);
        if (sysctlbyname("hw.pagesize", &page_size, &len, nullptr, 0) != 0 || page_size == 0) {
            page_size = 4096;
        }
        return static_cast<size_t>(free_pages) * static_cast<size_t>(page_size);
    }
    // Fallback when vm.page_free_count is unavailable: use a conservative
    // quarter of physical memory to preserve reserve-cap behavior.
    uint64_t phys_mem = 0;
    len = sizeof(phys_mem);
    if (sysctlbyname("hw.memsize", &phys_mem, &len, nullptr, 0) == 0) {
        return static_cast<size_t>(phys_mem / 4);
    }
    return 0;
#else
    // Linux: parse MemAvailable from /proc/meminfo
    std::ifstream meminfo("/proc/meminfo");
    if (meminfo) {
        std::string line;
        while (std::getline(meminfo, line)) {
            if (line.find("MemAvailable:") == 0) {
                unsigned long long kb = 0;
                if (std::sscanf(line.c_str(), "MemAvailable: %llu kB", &kb) == 1) {
                    return static_cast<size_t>(kb * 1024);
                }
            }
        }
    }
    return 0;
#endif
}

// --- GZIP Loaders ---

// Minimum compressed size to try parallel decompression (avoid overhead for tiny files)
static constexpr size_t PARALLEL_GZIP_THRESHOLD = 1024 * 1024; // 1MB

std::vector<size_t> SmartStrategy::scanGzipStreams(const std::string& filepath) const {
    std::vector<size_t> offsets;

    // mmap instead of reading the whole compressed file into a heap buffer:
    // loadGzipParallel() mmaps this same file again immediately afterwards,
    // so a full heap read here was a redundant full-size copy of the
    // compressed input held in memory purely to scan for magic bytes.
    MMapHandle mmap_scan;
#ifdef _WIN32
    mmap_scan.hFile = CreateFileA(filepath.c_str(), GENERIC_READ, FILE_SHARE_READ,
                                   NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
    if (mmap_scan.hFile == INVALID_HANDLE_VALUE) return offsets;
    LARGE_INTEGER fs{};
    if (!GetFileSizeEx(mmap_scan.hFile, &fs)) return offsets;
    mmap_scan.size = static_cast<size_t>(fs.QuadPart);
    if (mmap_scan.size < 18) return offsets; // Too small to contain even one GZIP stream
    mmap_scan.hMap = CreateFileMappingA(mmap_scan.hFile, NULL, PAGE_READONLY, 0, 0, NULL);
    if (!mmap_scan.hMap) return offsets;
    mmap_scan.data = MapViewOfFile(mmap_scan.hMap, FILE_MAP_READ, 0, 0, 0);
    if (!mmap_scan.data) return offsets;
#else
    mmap_scan.fd = open(filepath.c_str(), O_RDONLY);
    if (mmap_scan.fd == -1) return offsets;
    struct stat sb{};
    if (fstat(mmap_scan.fd, &sb) == -1) return offsets;
    mmap_scan.size = static_cast<size_t>(sb.st_size);
    if (mmap_scan.size < 18) return offsets; // Too small to contain even one GZIP stream
    mmap_scan.data = ::mmap(NULL, mmap_scan.size, PROT_READ, MAP_PRIVATE, mmap_scan.fd, 0);
    if (mmap_scan.data == MAP_FAILED) return offsets;
#endif

    const unsigned char* compressed = static_cast<const unsigned char*>(mmap_scan.data);
    size_t file_size = mmap_scan.size;

    // Find real stream boundaries by actually decoding each member with zlib
    // rather than scanning for the 0x1f 0x8b 0x08 magic bytes in isolation.
    // Highly-compressible genomic payloads make that byte sequence common
    // enough to appear *inside* a single member's compressed data purely by
    // chance (observed in practice: NCBI single-member FASTA.gz files with
    // 2-3 "matches" that were never real stream headers). Splitting at a
    // bogus offset feeds inflate() a truncated/misaligned deflate stream,
    // which can decode into gigabytes of garbage via LZ77 back-references —
    // not a real OOM, but indistinguishable from one without this fix.
    //
    // Output is discarded into a small reusable scratch buffer: we only need
    // to know where each member's compressed data ends (via inflate's
    // avail_in bookkeeping and Z_STREAM_END), not the decompressed content
    // itself, so this dry run stays cheap and constant-memory regardless of
    // file size.
    std::vector<char> scratch(64 * 1024);
    size_t offset = 0;
    while (offset + 3 <= file_size &&
           compressed[offset] == 0x1f && compressed[offset + 1] == 0x8b &&
           compressed[offset + 2] == 0x08) {
        offsets.push_back(offset);

        z_stream strm = {};
        if (inflateInit2(&strm, 15 + 16) != Z_OK) break;
        strm.next_in  = const_cast<Bytef*>(compressed + offset);
        strm.avail_in = static_cast<uInt>(file_size - offset);

        int ret = Z_OK;
        while (ret != Z_STREAM_END) {
            strm.next_out  = reinterpret_cast<Bytef*>(scratch.data());
            strm.avail_out = static_cast<uInt>(scratch.size());
            ret = inflate(&strm, Z_NO_FLUSH);
            if (ret == Z_STREAM_ERROR || ret == Z_DATA_ERROR || ret == Z_MEM_ERROR || ret == Z_BUF_ERROR) {
                // Not a genuine member boundary (or truncated stream) —
                // stop scanning; loadGzipParallel/loadGzipSingleStream will
                // surface the real error when they attempt full decoding.
                break;
            }
        }

        size_t consumed = (file_size - offset) - strm.avail_in;
        inflateEnd(&strm);
        if (ret != Z_STREAM_END || consumed == 0) break;
        offset += consumed;
    }

    // Fallback: if no valid header found (shouldn't happen for valid GZIP), return offset 0
    if (offsets.empty()) offsets.push_back(0);

    return offsets;
}

void SmartStrategy::loadGzipParallel(const std::string& filepath,
                                      const std::vector<size_t>& stream_offsets) {
    // mmap the compressed file for zero-copy read access across threads
    auto mmap = std::make_unique<MMapHandle>();
#ifdef _WIN32
    mmap->hFile = CreateFileA(filepath.c_str(), GENERIC_READ, FILE_SHARE_READ,
                              NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
    if (mmap->hFile == INVALID_HANDLE_VALUE)
        throw std::runtime_error("Cannot open GZIP file: " + filepath);
    LARGE_INTEGER fs{};
    GetFileSizeEx(mmap->hFile, &fs);
    mmap->size = static_cast<size_t>(fs.QuadPart);
    mmap->hMap = CreateFileMappingA(mmap->hFile, NULL, PAGE_READONLY, 0, 0, NULL);
    mmap->data = MapViewOfFile(mmap->hMap, FILE_MAP_READ, 0, 0, 0);
    if (!mmap->data) throw std::runtime_error("mmap failed: " + filepath);
#else
    mmap->fd = open(filepath.c_str(), O_RDONLY);
    if (mmap->fd == -1) throw std::runtime_error("Cannot open GZIP file: " + filepath);
    struct stat sb{};
    fstat(mmap->fd, &sb);
    mmap->size = static_cast<size_t>(sb.st_size);
    mmap->data = ::mmap(NULL, mmap->size, PROT_READ, MAP_PRIVATE, mmap->fd, 0);
    if (mmap->data == MAP_FAILED) throw std::runtime_error("mmap failed: " + filepath);
#endif

    const unsigned char* compressed_data = static_cast<const unsigned char*>(mmap->data);
    size_t file_size = mmap->size;

    size_t num_streams = stream_offsets.size();

    // Build per-stream compressed sizes
    std::vector<size_t> stream_compressed_sizes(num_streams);
    for (size_t i = 0; i < num_streams; ++i) {
        size_t stream_end = (i + 1 < num_streams) ? stream_offsets[i + 1] : file_size;
        stream_compressed_sizes[i] = stream_end - stream_offsets[i];
    }

    // Decompress each stream in parallel into its own buffer
    size_t num_threads = std::min(num_streams, static_cast<size_t>(std::thread::hardware_concurrency()));
    if (num_threads == 0) num_threads = 1;

    std::vector<std::vector<char>> decompressed_chunks(num_streams);
    std::vector<std::string> thread_errors(num_streams);

    // Shared budget so concurrently-growing per-stream buffers within the
    // *current batch* can't collectively exceed available memory. Streams
    // from earlier batches are merged into text_arena_ and freed as soon as
    // their batch joins (see batch loop below), so this only ever needs to
    // account for the up-to-`num_threads` streams live right now — not every
    // stream in the whole file.
    std::atomic<size_t> total_reserved_bytes{0};
    // Refreshed at the start of every batch so the guard reacts to memory
    // the OS has reclaimed since the previous batch's buffers were freed,
    // rather than working off one stale snapshot for the whole file.
    size_t avail_mem = getAvailableMemory();

    auto decompress_stream = [&](size_t stream_idx) {
        const unsigned char* src = compressed_data + stream_offsets[stream_idx];
        size_t src_size = stream_compressed_sizes[stream_idx];

        size_t estimated = src_size * 3;
        if (avail_mem > 0)
            estimated = std::min(estimated, avail_mem / 4 / num_threads);
        size_t reserve_size = estimated > 0 ? estimated : src_size * 3;

        auto& out = decompressed_chunks[stream_idx];
        try {
            out.resize(reserve_size);
        } catch (const std::bad_alloc&) {
            thread_errors[stream_idx] = "OOM: failed to allocate " +
                std::to_string(reserve_size) + " bytes for stream decompression";
            return;
        }
        total_reserved_bytes.fetch_add(reserve_size, std::memory_order_relaxed);

        z_stream strm = {};
        if (inflateInit2(&strm, 15 + 16) != Z_OK) {
            thread_errors[stream_idx] = "inflateInit2 failed";
            return;
        }

        strm.next_in  = const_cast<Bytef*>(src);
        strm.avail_in = static_cast<uInt>(src_size);

        size_t write_pos = 0;
        int ret = Z_OK;

        // Grows `out` geometrically, guarding against both a single
        // allocation failure and the combined budget across all threads.
        auto grow = [&]() -> bool {
            size_t old_size = out.size();
            size_t new_size = old_size * 2;
            if (avail_mem > 0 &&
                total_reserved_bytes.load(std::memory_order_relaxed) + (new_size - old_size) > avail_mem / 2) {
                thread_errors[stream_idx] = "OOM guard: aggregate GZIP stream buffers "
                    "would exceed available memory budget";
                return false;
            }
            try {
                out.resize(new_size);
            } catch (const std::bad_alloc&) {
                thread_errors[stream_idx] = "OOM: failed to grow stream buffer to " +
                    std::to_string(new_size) + " bytes";
                return false;
            }
            total_reserved_bytes.fetch_add(new_size - old_size, std::memory_order_relaxed);
            return true;
        };

        while (ret != Z_STREAM_END) {
            // Geometric growth if buffer is full
            if (write_pos == out.size()) {
                if (!grow()) { inflateEnd(&strm); return; }
            }

            strm.next_out  = reinterpret_cast<Bytef*>(out.data() + write_pos);
            strm.avail_out = static_cast<uInt>(out.size() - write_pos);

            ret = inflate(&strm, Z_NO_FLUSH);

            if (ret == Z_STREAM_ERROR || ret == Z_DATA_ERROR || ret == Z_MEM_ERROR) {
                inflateEnd(&strm);
                thread_errors[stream_idx] = "inflate error: " + std::to_string(ret);
                return;
            }

            write_pos = out.size() - strm.avail_out;

            if (ret == Z_BUF_ERROR) {
                // Output full but not done — grow and continue
                if (!grow()) { inflateEnd(&strm); return; }
                continue;
            }
        }

        inflateEnd(&strm);
        out.resize(write_pos);
        out.shrink_to_fit();
    };

    // Spawn threads in batches of num_threads. Streams are processed in
    // index order (batch_start increases monotonically, indices within a
    // batch are in order), so merging into text_arena_ immediately after
    // each batch joins — rather than once at the very end — preserves
    // stream order while keeping at most num_threads decompressed buffers
    // alive at any given time. This is what lets the OOM guard's
    // `total_reserved_bytes` reflect real concurrent memory use instead of
    // accumulating across every stream in the file.
    text_arena_.clear();
    for (size_t batch_start = 0; batch_start < num_streams; batch_start += num_threads) {
        size_t batch_end = std::min(batch_start + num_threads, num_streams);

        // Refresh the available-memory reading for this batch so the guard
        // reacts to memory freed by previously-merged batches instead of
        // working off a snapshot taken before any decompression happened.
        avail_mem = getAvailableMemory();
        total_reserved_bytes.store(0, std::memory_order_relaxed);

        std::vector<std::thread> threads;
        for (size_t i = batch_start; i < batch_end; ++i) {
            threads.emplace_back(decompress_stream, i);
        }
        for (auto& t : threads) t.join();

        // Check for per-thread errors before merging this batch.
        for (size_t i = batch_start; i < batch_end; ++i) {
            if (!thread_errors[i].empty()) {
                mmap.reset();
                throw std::runtime_error("Parallel GZIP stream " + std::to_string(i) +
                                         " failed: " + thread_errors[i]);
            }
        }

        // Merge this batch's chunks into text_arena_ and free them
        // immediately so the next batch's memory budget starts fresh.
        for (size_t i = batch_start; i < batch_end; ++i) {
            auto& dc = decompressed_chunks[i];
            size_t pos = text_arena_.size();
            text_arena_.resize(pos + dc.size());
            std::memcpy(text_arena_.data() + pos, dc.data(), dc.size());
            dc.clear();
            dc.shrink_to_fit();
        }
    }

    // Release mmap — compressed data no longer needed
    mmap.reset();
}

void SmartStrategy::loadGzipSingleStream(const std::string& filepath) {
    const size_t CHUNK_SIZE = 1024 * 1024; // 1MB

    text_arena_.clear();

    // Pre-size text_arena_ with 3x heuristic and OOM guard
    {
        std::error_code ec;
        const uintmax_t raw_size = std::filesystem::file_size(filepath, ec);
        if (!ec && raw_size > 0) {
            const size_t estimated_size = static_cast<size_t>(raw_size * 3);
            const size_t avail_mem = getAvailableMemory();
            const size_t reserve_size = (avail_mem > 0)
                ? std::min(estimated_size, avail_mem / 4)
                : estimated_size;
            try {
                text_arena_.reserve(reserve_size);
            } catch (const std::bad_alloc&) {
                std::cerr << "SmartStrategy OOM: failed to reserve "
                          << reserve_size << " bytes for " << filepath
                          << " (compressed=" << raw_size << ", estimated="
                          << estimated_size << ")" << std::endl;
                throw;
            }
        } else {
            text_arena_.reserve(CHUNK_SIZE * 4);
        }
    }

    gzFile file = gzopen(filepath.c_str(), "rb");
    if (!file) throw std::runtime_error("Cannot open GZIP file: " + filepath);

    auto chunk = std::make_unique<char[]>(CHUNK_SIZE);
    size_t write_pos = 0;

    try {
        while (true) {
            size_t required = write_pos + CHUNK_SIZE;
            if (required > text_arena_.capacity()) {
                text_arena_.reserve(std::max(text_arena_.capacity() * 2, required));
            }
            text_arena_.resize(required);

            int bytes_read = gzread(file, chunk.get(), CHUNK_SIZE);
            if (bytes_read < 0) {
                int err;
                const char* error_msg = gzerror(file, &err);
                text_arena_.clear();
                throw std::runtime_error("GZIP read error: " + std::string(error_msg));
            }
            if (bytes_read == 0) break;
            std::memcpy(text_arena_.data() + write_pos, chunk.get(), bytes_read);
            write_pos += static_cast<size_t>(bytes_read);
        }
    } catch (...) {
        gzclose(file);
        text_arena_.clear();
        throw;
    }
    gzclose(file);

    text_arena_.resize(write_pos);
    text_arena_.shrink_to_fit();
}

void SmartStrategy::loadGzipInternal(const std::string& filepath) {
    std::error_code ec;
    const uintmax_t file_size = std::filesystem::file_size(filepath, ec);
    const bool large_enough = (!ec && file_size >= PARALLEL_GZIP_THRESHOLD);

    if (large_enough) {
        auto streams = scanGzipStreams(filepath);
        if (streams.size() > 1) {
            loadGzipParallel(filepath, streams);
            return;
        }
    }

    loadGzipSingleStream(filepath);
}

void SmartStrategy::normalizeFastaArena() {
    if (text_arena_.empty()) return;
    const size_t original_size = text_arena_.size();
    // Grow vector by 1 to safely write trailing '\n' past original end.
    // resize() value-initializes the new element, but we overwrite it immediately.
    text_arena_.resize(original_size + 1);
    char*       write = text_arena_.data();
    const char* read  = text_arena_.data();
    const char* end   = read + original_size;  // Read boundary is original content

    // Skip any leading blank lines
    while (read < end && (*read == '\n' || *read == '\r')) ++read;

    while (read < end) {
        if (*read != '>') break;

        // ── Copy header line up to (but not including) \r\n ─────────────────
        // simd_find_char scans 32 bytes/iter (AVX2) or 16 bytes/iter (NEON).
        const char* nl      = simd_find_char(read, end, '\n');
        const char* hdr_end = (nl > read && *(nl - 1) == '\r') ? nl - 1 : nl;
        size_t hdr_len      = static_cast<size_t>(hdr_end - read);
        // write <= read always (compaction), so memmove handles overlap safely
        std::memmove(write, read, hdr_len);
        write += hdr_len;
        if (nl < end) *write++ = '\n'; // keep the terminating newline
        read = (nl < end) ? nl + 1 : end; // advance past \n

        // ── Locate the end of this record: next '>' or end of arena ──────────
        // Biggest win on large genomes: 32 bytes/iter instead of 1 byte/iter
        // for the O(sequence_length) scan that formerly drove the inner loop.
        const char* rec_end = simd_find_char(read, end, '>');

        // ── Copy sequence bytes, stripping all \r and \n ─────────────────────
        // Use simd_find_char to jump to the next \n, then bulk-copy the
        // non-newline run with memmove.  Tail bytes (<32 or <16) handled by
        // simd_find_char's scalar fall-through — no special case needed here.
        while (read < rec_end) {
            const char* line_nl  = simd_find_char(read, rec_end, '\n');
            // Strip a trailing \r if it immediately precedes the \n (or rec_end)
            const char* copy_end = (line_nl > read && *(line_nl - 1) == '\r')
                                   ? line_nl - 1 : line_nl;
            size_t run = static_cast<size_t>(copy_end - read);
            if (run) std::memmove(write, read, run);
            write += run;
            // Advance past the \n; step to rec_end if no \n was found
            read = (line_nl < rec_end) ? line_nl + 1 : rec_end;
        }
        *write++ = '\n'; // terminate the now-single-line sequence
        // read == rec_end (pointing at '>' of next record, or end)
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
    bool use_ngs_mode = (index_mode_ == IndexMode::NGS);
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
    // One-time O(n) walk to seed serialized_size_estimate_ so saveBinary()
    // doesn't need to re-walk the map on every call; addEntry() maintains
    // this incrementally after this point.
    refreshPayloadEstimate();
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

    // OOM guard: mirrors loadGzipSingleStream's reserve-with-catch pattern.
    // Uncompressed loads have no size-reduction headroom (unlike GZIP, where
    // the compressed size hints at a smaller footprint), so we only guard
    // against an implausible allocation relative to available memory.
    const size_t avail_mem = getAvailableMemory();
    if (avail_mem > 0 && fileSize > avail_mem / 2) {
        throw std::runtime_error(
            "SmartStrategy OOM guard: refusing to load " + std::to_string(fileSize) +
            " bytes from " + filepath + " (available memory ~" +
            std::to_string(avail_mem) + " bytes)");
    }
    try {
        text_arena_.resize(fileSize);
    } catch (const std::bad_alloc&) {
        std::cerr << "SmartStrategy OOM: failed to allocate " << fileSize
                  << " bytes for " << filepath << std::endl;
        throw;
    }
    if (!file.read(text_arena_.data(), fileSize)) throw std::runtime_error("Read failed");
    parseArena();
}

// --- Templated Parsers ---
template <typename MapType> void SmartStrategy::parseFastaMultithreadedTemplate(std::string_view content, MapType& dest_map) {
    const size_t num_threads = std::max(1u, std::thread::hardware_concurrency());
    const size_t chunk_size = content.size() / num_threads;
    std::vector<std::thread> threads;
    std::vector<MapType> thread_caches(num_threads);
    // Pre-reserve thread-local maps to avoid mid-parse rehashing
    // Use conservative heuristic: assume avg record ~500 bytes (header + sequence + newlines)
    // This avoids over-reservation with large records while still reducing rehashing
    {
        const size_t est_per_thread = chunk_size / 500;
        for (auto& cache : thread_caches) {
            cache.reserve(static_cast<size_t>(est_per_thread * 1.25));
        }
    }
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
            ptr = simd_find_char(ptr, global_end, '\n');
            while (ptr < global_end && (*ptr == '\n' || *ptr == '\r')) ++ptr;
            const char* seq_start = ptr; const char* seq_end = seq_start;
            const char* next_record = simd_find_char(ptr, global_end, '>');
            if (next_record < global_end) { seq_end = next_record; while (seq_end > seq_start && (*(seq_end-1) == '\n' || *(seq_end-1) == '\r' || *(seq_end-1) == ' ')) --seq_end; ptr = next_record; }
            else { seq_end = global_end; while (seq_end > seq_start && (*(seq_end-1) == '\n' || *(seq_end-1) == '\r' || *(seq_end-1) == ' ')) --seq_end; ptr = global_end; }
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
    // Pre-reserve thread-local maps to avoid mid-parse rehashing
    // Use conservative heuristic: assume avg record ~600 bytes (header + seq + qual + newlines)
    // This avoids over-reservation with large records while still reducing rehashing
    {
        const size_t est_per_thread = chunk_size / 600;
        for (auto& cache : thread_caches) {
            cache.reserve(static_cast<size_t>(est_per_thread * 1.25));
        }
    }
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
            ptr = simd_find_char(ptr, global_end, '\n');
            while (ptr < global_end && (*ptr == '\n' || *ptr == '\r')) ++ptr;
            const char* seq_start = ptr;
            ptr = simd_find_char(ptr, global_end, '\n');
            const char* seq_end = ptr;
            while (seq_end > seq_start && *(seq_end - 1) == '\r') --seq_end;
            while (ptr < global_end && (*ptr == '\n' || *ptr == '\r')) ++ptr;
            ptr = simd_find_char(ptr, global_end, '\n');
            while (ptr < global_end && (*ptr == '\n' || *ptr == '\r')) ++ptr;
            const char* qual_start = ptr;
            ptr = simd_find_char(ptr, global_end, '\n');
            const char* qual_end = ptr;
            while (qual_end > qual_start && *(qual_end - 1) == '\r') --qual_end;
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
        ptr = simd_find_char(ptr, end, '\n'); if (ptr < end) ++ptr; while (ptr < end && (*ptr == '\n' || *ptr == '\r')) ++ptr;
        const char* seq_start = ptr; const char* seq_end = seq_start;
        const char* next_record = simd_find_char(ptr, end, '>'); if (next_record < end) { seq_end = next_record; while (seq_end > seq_start && (*(seq_end-1) == '\n' || *(seq_end-1) == '\r' || *(seq_end-1) == ' ')) --seq_end; ptr = next_record; } else { seq_end = end; while (seq_end > seq_start && (*(seq_end-1) == '\n' || *(seq_end-1) == '\r' || *(seq_end-1) == ' ')) --seq_end; ptr = end; }
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
        ptr = simd_find_char(ptr, end, '\n'); while (ptr < end && (*ptr == '\n' || *ptr == '\r')) ++ptr;
        const char* seq_start = ptr; ptr = simd_find_char(ptr, end, '\n'); const char* seq_end = ptr;
        while (seq_end > seq_start && *(seq_end - 1) == '\r') --seq_end;
        while (ptr < end && (*ptr == '\n' || *ptr == '\r')) ++ptr; ptr = simd_find_char(ptr, end, '\n'); while (ptr < end && (*ptr == '\n' || *ptr == '\r')) ++ptr;
        const char* qual_start = ptr; ptr = simd_find_char(ptr, end, '\n'); const char* qual_end = ptr;
        while (qual_end > qual_start && *(qual_end - 1) == '\r') --qual_end;
        while (ptr < end && (*ptr == '\n' || *ptr == '\r')) ++ptr;
        std::string_view id(id_start, id_end - id_start); std::string_view seq(seq_start, seq_end - seq_start); std::string_view qual(qual_start, qual_end - qual_start);
        if (!id.empty() && !seq.empty()) {
            if constexpr (std::is_same_v<MapType, GenomeIndex>) { map.emplace(std::string(id), SequenceView{id, seq, qual}); }
            else { map.emplace(hash_key(id), SequenceView{id, seq, qual}); }
        }
    }
}

void SmartStrategy::serializePayload(const std::function<void(const char*, size_t)>& sink) const {
    uint8_t mode = std::holds_alternative<NGSIndex>(file_cache_) ? 1 : 0;

    auto emit = [&sink](const void* p, size_t n) {
        sink(reinterpret_cast<const char*>(p), n);
    };

    if (mode == 0) {
        const auto& map = std::get<GenomeIndex>(file_cache_);
        uint64_t count = map.size();
        emit(&count, sizeof(count));

        for (const auto& [key, view] : map) {
            uint32_t len = static_cast<uint32_t>(key.size());
            emit(&len, sizeof(len));
            emit(key.data(), len);

            len = static_cast<uint32_t>(view.sequence.size());
            emit(&len, sizeof(len));
            emit(view.sequence.data(), len);

            len = static_cast<uint32_t>(view.quality.size());
            emit(&len, sizeof(len));
            if (len > 0) emit(view.quality.data(), len);
        }
    } else {
        const auto& map = std::get<NGSIndex>(file_cache_);
        uint64_t count = map.size();
        emit(&count, sizeof(count));

        for (const auto& [key, view] : map) {
            emit(&key, sizeof(key));

            uint32_t len = static_cast<uint32_t>(view.id.size());
            emit(&len, sizeof(len));
            emit(view.id.data(), len);

            len = static_cast<uint32_t>(view.sequence.size());
            emit(&len, sizeof(len));
            emit(view.sequence.data(), len);

            len = static_cast<uint32_t>(view.quality.size());
            emit(&len, sizeof(len));
            if (len > 0) emit(view.quality.data(), len);
        }
    }
}

size_t SmartStrategy::estimatePayloadSize() const {
    // O(1): backed by serialized_size_estimate_, kept exact and up to date
    // by refreshPayloadEstimate() (called once after a fresh parse/restore)
    // and addEntry() (incremental updates thereafter). This must stay byte-
    // exact with serializePayload()'s actual output — loadBinary() relies on
    // it to pre-size text_arena_ before streaming decompression.
    return serialized_size_estimate_.load(std::memory_order_relaxed);
}

void SmartStrategy::refreshPayloadEstimate() {
    // Mirrors serializePayload()'s byte layout without writing anything.
    // Only called after a full (re)parse/restore — NOT on saveBinary()'s hot
    // path — so repeated saveBinary() calls don't re-walk the whole map.
    size_t total = sizeof(uint64_t); // record count
    if (std::holds_alternative<GenomeIndex>(file_cache_)) {
        const auto& map = std::get<GenomeIndex>(file_cache_);
        for (const auto& [key, view] : map) {
            total += 3 * sizeof(uint32_t) + key.size() + view.sequence.size() + view.quality.size();
        }
    } else {
        const auto& map = std::get<NGSIndex>(file_cache_);
        for (const auto& [key, view] : map) {
            total += sizeof(key) + 3 * sizeof(uint32_t) + view.id.size() + view.sequence.size() + view.quality.size();
        }
    }
    serialized_size_estimate_.store(total, std::memory_order_relaxed);
}

CompressionMode
SmartStrategy::selectCompressionStrategy(size_t payload_size) const {
    constexpr size_t LARGE_THRESHOLD = 10u * 1024u * 1024u; // 10 MiB

    if (payload_size > LARGE_THRESHOLD) {
        switch (detected_format_) {
            case FileFormat::DNA_FASTA:
            case FileFormat::RNA_FASTA:
            case FileFormat::DNA_FASTQ:
            case FileFormat::RNA_FASTQ:
                return CompressionMode::LZ4HC;
            default:
                break;
        }
    }
    return CompressionMode::LZ4Default;
}

void SmartStrategy::saveBinary(const std::string& filepath) const {
    std::shared_lock<std::shared_mutex> lock(cache_mutex_);
    std::ofstream out(filepath, std::ios::binary);
    if (!out) throw std::runtime_error("Cannot write: " + filepath);

    // v3 format: stream-serialize + stream-compress via LZ4 Frame (LZ4F)
    // instead of materializing a full serialized "payload" copy plus a full
    // "compressed" copy of the whole dataset (the old v2 path's ~3x peak
    // memory: text_arena_ + payload + compressed, all resident at once).
    // Peak extra memory here is bounded by STREAM_CHUNK_SIZE regardless of
    // dataset size.
    const size_t estimated_size = estimatePayloadSize();
    const CompressionMode comp_mode = selectCompressionStrategy(estimated_size);

    out.write(MAGIC_BYTES_V3, 4);
    uint8_t index_mode = std::holds_alternative<NGSIndex>(file_cache_) ? 1 : 0;
    out.write(reinterpret_cast<const char*>(&index_mode), 1);
    uint64_t original_size = static_cast<uint64_t>(estimated_size);
    out.write(reinterpret_cast<const char*>(&original_size), sizeof(original_size));

    LZ4F_preferences_t prefs;
    std::memset(&prefs, 0, sizeof(prefs));
    prefs.compressionLevel = (comp_mode == CompressionMode::LZ4HC) ? LZ4HC_CLEVEL_DEFAULT : 0;

    LZ4F_cctx* cctx = nullptr;
    if (LZ4F_isError(LZ4F_createCompressionContext(&cctx, LZ4F_VERSION)))
        throw std::runtime_error("LZ4F_createCompressionContext failed: " + filepath);

    struct CtxGuard {
        LZ4F_cctx* c;
        ~CtxGuard() { if (c) LZ4F_freeCompressionContext(c); }
    } guard{cctx};

    std::vector<char> in_buf;
    in_buf.reserve(STREAM_CHUNK_SIZE);
    std::vector<char> out_buf(std::max(LZ4F_compressBound(STREAM_CHUNK_SIZE, &prefs),
                                        static_cast<size_t>(64 * 1024)));

    size_t header_size = LZ4F_compressBegin(cctx, out_buf.data(), out_buf.size(), &prefs);
    if (LZ4F_isError(header_size))
        throw std::runtime_error("LZ4F_compressBegin failed: " + filepath);
    out.write(out_buf.data(), static_cast<std::streamsize>(header_size));

    auto flush_chunk = [&](const char* data, size_t len) {
        size_t bound = LZ4F_compressBound(len, &prefs);
        if (out_buf.size() < bound) out_buf.resize(bound);
        size_t written = LZ4F_compressUpdate(cctx, out_buf.data(), out_buf.size(),
                                              data, len, nullptr);
        if (LZ4F_isError(written))
            throw std::runtime_error("LZ4F_compressUpdate failed: " + filepath);
        if (written > 0) out.write(out_buf.data(), static_cast<std::streamsize>(written));
    };

    serializePayload([&](const char* data, size_t len) {
        in_buf.insert(in_buf.end(), data, data + len);
        if (in_buf.size() >= STREAM_CHUNK_SIZE) {
            flush_chunk(in_buf.data(), in_buf.size());
            in_buf.clear();
        }
    });
    if (!in_buf.empty()) {
        flush_chunk(in_buf.data(), in_buf.size());
        in_buf.clear();
    }

    size_t end_size = LZ4F_compressEnd(cctx, out_buf.data(), out_buf.size(), nullptr);
    if (LZ4F_isError(end_size))
        throw std::runtime_error("LZ4F_compressEnd failed: " + filepath);
    out.write(out_buf.data(), static_cast<std::streamsize>(end_size));

    if (!out) throw std::runtime_error("Binary cache write failed: " + filepath);
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

    // Check magic bytes and determine format version
    if (mmap_handle_->size < 5)
        throw std::runtime_error("Binary cache too small to contain header: " + filepath);

    const char* magic = safe_advance(4);
    uint8_t format_version = magic[3];

    if (std::strncmp(magic, "TRO", 3) != 0)
        throw std::runtime_error("Invalid binary cache magic bytes: " + filepath);

    uint8_t mode;
    std::memcpy(&mode, safe_advance(1), 1);
    // Sync index_mode_ with what's stored in the binary so that subsequent
    // loadFile() / parseArena() calls use the same index type as restored data.
    index_mode_ = (mode == 1) ? IndexMode::NGS : IndexMode::GENOME;

    // --- V2 Format (LZ4-compressed) ---
    if (format_version == 0x02) {
        uint64_t original_size, compressed_size;
        std::memcpy(&original_size, safe_advance(sizeof(uint64_t)), sizeof(uint64_t));
        std::memcpy(&compressed_size, safe_advance(sizeof(uint64_t)), sizeof(uint64_t));

        // Bounds check for compressed data
        if (ptr + compressed_size > end)
            throw std::runtime_error("Binary cache v2 is truncated: " + filepath);

        const char* compressed_data = safe_advance(compressed_size);

        // Decompress into text_arena_
        text_arena_.resize(original_size);
        int decompressed_size = LZ4_decompress_safe(compressed_data, text_arena_.data(),
                                                     static_cast<int>(compressed_size),
                                                     static_cast<int>(original_size));
        if (decompressed_size < 0)
            throw std::runtime_error("LZ4 decompression failed for binary cache: " + filepath);
        if (decompressed_size != static_cast<int>(original_size))
            throw std::runtime_error("LZ4 decompressed size mismatch: " + filepath);

        // Release mmap handle since we now own the data in text_arena_
        mmap_handle_.reset();

        // Parse from text_arena_
        const char* arena_ptr = text_arena_.data();
        const char* arena_end = arena_ptr + text_arena_.size();

        auto arena_safe_advance = [&](size_t n) -> const char* {
            if (arena_ptr + n > arena_end)
                throw std::runtime_error("Binary cache v2 payload is corrupt: " + filepath);
            const char* result = arena_ptr;
            arena_ptr += n;
            return result;
        };

        uint64_t count;
        std::memcpy(&count, arena_safe_advance(sizeof(uint64_t)), sizeof(uint64_t));

        constexpr uint64_t MAX_RECORDS = 1'000'000'000ULL;
        if (count > MAX_RECORDS)
            throw std::runtime_error("Binary cache record count implausible: " + filepath);

        if (mode == 0) {
            file_cache_ = GenomeIndex{};
            auto& map = std::get<GenomeIndex>(file_cache_);
            map.reserve(count);
            for (uint64_t i = 0; i < count; ++i) {
                uint32_t len; std::memcpy(&len, arena_safe_advance(4), 4);
                std::string_view id_view(arena_safe_advance(len), len);
                uint32_t seq_len; std::memcpy(&seq_len, arena_safe_advance(4), 4);
                std::string_view seq(arena_safe_advance(seq_len), seq_len);
                uint32_t qual_len; std::memcpy(&qual_len, arena_safe_advance(4), 4);
                std::string_view qual;
                if (qual_len > 0) qual = std::string_view(arena_safe_advance(qual_len), qual_len);
                map.emplace(std::string(id_view), SequenceView{id_view, seq, qual});
            }
        } else {
            file_cache_ = NGSIndex{};
            auto& map = std::get<NGSIndex>(file_cache_);
            map.reserve(count);
            for (uint64_t i = 0; i < count; ++i) {
                uint64_t hash; std::memcpy(&hash, arena_safe_advance(8), 8);
                uint32_t len; std::memcpy(&len, arena_safe_advance(4), 4);
                std::string_view id_view(arena_safe_advance(len), len);
                uint32_t seq_len; std::memcpy(&seq_len, arena_safe_advance(4), 4);
                std::string_view seq(arena_safe_advance(seq_len), seq_len);
                uint32_t qual_len; std::memcpy(&qual_len, arena_safe_advance(4), 4);
                std::string_view qual;
                if (qual_len > 0) qual = std::string_view(arena_safe_advance(qual_len), qual_len);
                map.emplace(hash, SequenceView{id_view, seq, qual});
            }
        }
    }
    // --- V3 Format (LZ4 Frame, streamed) ---
    else if (format_version == 0x03) {
        uint64_t original_size;
        std::memcpy(&original_size, safe_advance(sizeof(uint64_t)), sizeof(uint64_t));

        const char* frame_data = ptr; // remaining mmap'd bytes are the LZ4F frame
        const size_t frame_size = static_cast<size_t>(end - ptr);

        // OOM guard: mirrors the pattern used by loadFile()/loadGzipSingleStream().
        const size_t avail_mem = getAvailableMemory();
        if (avail_mem > 0 && original_size > avail_mem / 2) {
            throw std::runtime_error(
                "SmartStrategy OOM guard: refusing to decompress " +
                std::to_string(original_size) + " bytes from " + filepath +
                " (available memory ~" + std::to_string(avail_mem) + " bytes)");
        }
        try {
            text_arena_.resize(original_size);
        } catch (const std::bad_alloc&) {
            std::cerr << "SmartStrategy OOM: failed to allocate " << original_size
                      << " bytes for " << filepath << std::endl;
            throw;
        }

        LZ4F_dctx* dctx = nullptr;
        if (LZ4F_isError(LZ4F_createDecompressionContext(&dctx, LZ4F_VERSION)))
            throw std::runtime_error("LZ4F_createDecompressionContext failed: " + filepath);
        struct DctxGuard {
            LZ4F_dctx* d;
            ~DctxGuard() { if (d) LZ4F_freeDecompressionContext(d); }
        } dguard{dctx};

        const char* src = frame_data;
        size_t src_remaining = frame_size;
        char* dst = text_arena_.data();
        size_t dst_remaining = text_arena_.size();
        size_t hint = 1;

        while (hint != 0 && src_remaining > 0) {
            size_t dst_size = dst_remaining;
            size_t src_size = src_remaining;
            hint = LZ4F_decompress(dctx, dst, &dst_size, src, &src_size, nullptr);
            if (LZ4F_isError(hint))
                throw std::runtime_error("LZ4F_decompress failed for binary cache: " + filepath);
            dst += dst_size;
            dst_remaining -= dst_size;
            src += src_size;
            src_remaining -= src_size;
            if (dst_size == 0 && src_size == 0 && hint != 0) break; // no progress; malformed frame
        }
        if (dst_remaining != 0)
            throw std::runtime_error("Binary cache v3 decompressed size mismatch: " + filepath);

        // Release mmap handle since we now own the data in text_arena_
        mmap_handle_.reset();

        // Parse from text_arena_
        const char* arena_ptr = text_arena_.data();
        const char* arena_end = arena_ptr + text_arena_.size();

        auto arena_safe_advance = [&](size_t n) -> const char* {
            if (arena_ptr + n > arena_end)
                throw std::runtime_error("Binary cache v3 payload is corrupt: " + filepath);
            const char* result = arena_ptr;
            arena_ptr += n;
            return result;
        };

        uint64_t count;
        std::memcpy(&count, arena_safe_advance(sizeof(uint64_t)), sizeof(uint64_t));

        constexpr uint64_t MAX_RECORDS = 1'000'000'000ULL;
        if (count > MAX_RECORDS)
            throw std::runtime_error("Binary cache record count implausible: " + filepath);

        if (mode == 0) {
            file_cache_ = GenomeIndex{};
            auto& map = std::get<GenomeIndex>(file_cache_);
            map.reserve(count);
            for (uint64_t i = 0; i < count; ++i) {
                uint32_t len; std::memcpy(&len, arena_safe_advance(4), 4);
                std::string_view id_view(arena_safe_advance(len), len);
                uint32_t seq_len; std::memcpy(&seq_len, arena_safe_advance(4), 4);
                std::string_view seq(arena_safe_advance(seq_len), seq_len);
                uint32_t qual_len; std::memcpy(&qual_len, arena_safe_advance(4), 4);
                std::string_view qual;
                if (qual_len > 0) qual = std::string_view(arena_safe_advance(qual_len), qual_len);
                map.emplace(std::string(id_view), SequenceView{id_view, seq, qual});
            }
        } else {
            file_cache_ = NGSIndex{};
            auto& map = std::get<NGSIndex>(file_cache_);
            map.reserve(count);
            for (uint64_t i = 0; i < count; ++i) {
                uint64_t hash; std::memcpy(&hash, arena_safe_advance(8), 8);
                uint32_t len; std::memcpy(&len, arena_safe_advance(4), 4);
                std::string_view id_view(arena_safe_advance(len), len);
                uint32_t seq_len; std::memcpy(&seq_len, arena_safe_advance(4), 4);
                std::string_view seq(arena_safe_advance(seq_len), seq_len);
                uint32_t qual_len; std::memcpy(&qual_len, arena_safe_advance(4), 4);
                std::string_view qual;
                if (qual_len > 0) qual = std::string_view(arena_safe_advance(qual_len), qual_len);
                map.emplace(hash, SequenceView{id_view, seq, qual});
            }
        }
    }
    // --- V1 Format (uncompressed, mmap'd) ---
    else if (format_version == 0x01) {
        uint64_t count;
        std::memcpy(&count, safe_advance(sizeof(uint64_t)), sizeof(uint64_t));

        constexpr uint64_t MAX_RECORDS = 1'000'000'000ULL;
        if (count > MAX_RECORDS)
            throw std::runtime_error("Binary cache record count implausible: " + filepath);

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
    } else {
        throw std::runtime_error("Unknown binary cache format version: " + filepath);
    }

    determine_format_from_cache();
    // One-time O(n) walk to seed serialized_size_estimate_ for a subsequent
    // saveBinary() (e.g. re-saving after restoring + addEntry()'ing more
    // data); addEntry() maintains it incrementally after this point.
    refreshPayloadEstimate();
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