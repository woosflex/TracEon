#include "SmartStrategy.h"
#include "SimdUtils.h"
#include "Crc32c.h"
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
#include <functional> 
#include <cctype>
#include <zlib.h>
#include <lz4.h>
#include <lz4hc.h>
#include <lz4frame.h>
#include <new>
#include <atomic>
#include <climits>
#include <limits>

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

// ── `.traceon` v4 binary format ─────────────────────────────────────────────
// Header layout (all multi-byte fields little-endian):
//
//   offset 0  magic          "TRO\x04"           (4 bytes)
//   offset 4  codec flags    0x01 = LZ4 Frame;   (1 byte)
//                            bits 1..7 reserved, must be 0
//   offset 5  index mode     0 = GENOME, 1 = NGS (1 byte)
//   offset 6  logical length  uncompressed payload length (u64 LE)
//   offset 14 frame length    LZ4 Frame length    (u64 LE)
//   offset 22 CRC32C          Castagnoli, u32 LE, over the ENTIRE
//                            uncompressed logical payload followed by the
//                            canonical header bytes [0..22) (the checksum
//                            field itself excluded)
//   offset 26 LZ4 Frame       (frame length bytes)
//
// Total header size: 26 bytes. The checksum is a whole-payload integrity
// check (accidental-corruption detection, NOT authentication — a deliberate
// attacker can construct collisions). Truncation is detected by the exact
// logical-length / exact-frame-termination requirements on load, not by the
// checksum alone. The payload-first ordering is deliberate: it lets the
// checksum be computed incrementally on both sides — as serialized chunks
// pass through the LZ4F compressor on save, and as decompressed chunks are
// written to text_arena_ on load — with the 22 header bytes fed afterwards
// (frame_len is only known once the frame is complete). See
// docs/architecture/ADR-005-traceon-v4-binary-format.md and the design
// review at outputs/traceon-v2-design-review.md.
static constexpr char V4_MAGIC[4] = {'T', 'R', 'O', '\x04'};
static constexpr uint8_t V4_CODEC_LZ4F = 0x01;   // codec flags: LZ4 Frame
static constexpr uint8_t V4_CODEC_MASK = 0x01;   // only bit 0 is defined
static constexpr size_t V4_HEADER_LEN = 26;      // 4+1+1+8+8+4
static constexpr size_t V4_CRC_OFFSET = 22;      // header bytes [0..22) are checksummed

// Streaming (de)compression chunk size for the v4 binary cache format.
// Bounds peak memory during saveBinary()/loadBinary() to a small constant
// regardless of dataset size (see ADR-004 follow-up: streaming binary cache).
static constexpr size_t STREAM_CHUNK_SIZE = 1024 * 1024; // 1MB

// ── Explicit little-endian header field codecs ──────────────────────────────
// The v4 header fields are little-endian on the wire regardless of host byte
// order. The payload itself keeps the existing v1-v3 native-endian record
// layout (serializePayload()) — unchanged by the v4 format.
static void write_le64(char* dst, uint64_t v) noexcept {
    for (int i = 0; i < 8; ++i) dst[i] = static_cast<char>((v >> (8 * i)) & 0xFFu);
}
static void write_le32(char* dst, uint32_t v) noexcept {
    for (int i = 0; i < 4; ++i) dst[i] = static_cast<char>((v >> (8 * i)) & 0xFFu);
}
static uint64_t read_le64(const char* p) noexcept {
    uint64_t v = 0;
    for (int i = 7; i >= 0; --i) v = (v << 8) | static_cast<uint8_t>(p[i]);
    return v;
}
static uint32_t read_le32(const char* p) noexcept {
    uint32_t v = 0;
    for (int i = 3; i >= 0; --i) v = (v << 8) | static_cast<uint8_t>(p[i]);
    return v;
}

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
    // ── Reader-quiescence diagnostic (TRACEON_DEBUG_LIFECYCLE, opt-in) ──────
    // Detects a getView() call whose LOOKUP overlaps this teardown. It cannot
    // detect a caller that retained the returned std::string_view and
    // dereferences it AFTER clearCache()/reload/destruction completes — a
    // debug assertion cannot prove anything about post-return use. It is a
    // diagnostic for coordinated misuse, NOT a synchronization or reclamation
    // mechanism (see ADR-001 / README lifecycle contract).
#ifdef TRACEON_DEBUG_LIFECYCLE
    const size_t active = active_readers_.load(std::memory_order_relaxed);
    if (active != 0) {
        std::cerr << "TracEon lifecycle warning: clearCache()/reload while "
                  << active << " getView() call(s) are still in flight — "
                     "retained string_views into the old snapshot are now "
                     "dangling (reader quiescence contract, see ADR-001)"
                  << std::endl;
    }
#endif
    // ORDERING CONTRACT (string_view-keyed index): the map is cleared BEFORE
    // any backing store is released. GenomeIndex keys/values are non-owning
    // views into text_arena_, the mmap region, or manual_store_; all three
    // are reset/released below. Clearing the map first guarantees no
    // string_view ever outlives its backing buffer — a dangling view would
    // only matter if the map outlived the buffer, which this order forbids.
    data_ready_.store(false, std::memory_order_release);
    data_loaded_.store(false, std::memory_order_release);
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

void SmartStrategy::markDataLoaded() {
    data_loaded_.store(true, std::memory_order_release);
}

void SmartStrategy::addEntry(const std::string& id, const std::string& seq, const std::string& qual) {
    std::unique_lock<std::shared_mutex> lock(cache_mutex_);

    // Immutable-after-load contract (Bug 3 fix — see ADR-001). The lock-free
    // read path (getView()/getQuality()/hasSequence()) checks data_ready_
    // and then reads the map WITHOUT holding a lock, which is only safe if
    // the map is truly immutable once published. addEntry() is therefore
    // only legal while building a cache (before any loadFile()/loadBinary())
    // or after clearCache(); mutating a loaded cache would race with
    // concurrent lock-free readers (torn reads, and a segfault on rehash).
    // Note: data_ready_ alone can't gate this — addEntry() itself publishes
    // manual entries by setting data_ready_, so a separate data_loaded_ flag
    // records whether a LOAD has happened.
    if (data_loaded_.load(std::memory_order_acquire)) {
        throw std::logic_error(
            "TracEon: cannot add entries after data has been loaded "
            "(immutable-after-load contract; call clearCache() first)");
    }

    // Duplicate-key guard (Bug 2 fix). emplace() keeps the FIRST value for a
    // given key, so a duplicate set() is a silent no-op. Bail out before
    // touching manual_store_ or serialized_size_estimate_: an unconditional
    // estimate increment for a record that was never inserted would inflate
    // the binary-cache header's declared size, making restore() throw
    // "Binary cache v3 decompressed size mismatch" (the payload actually
    // holds fewer bytes than the header claims).
    if (std::holds_alternative<GenomeIndex>(file_cache_)) {
        if (std::get<GenomeIndex>(file_cache_).count(id) > 0) return;
    } else {
        if (std::get<NGSIndex>(file_cache_).count(hash_key(id)) > 0) return;
    }

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

    // Do NOT flip data_ready_ to false here. getView() checks data_ready_
    // without holding a lock (lock-free fast path), so a false→true toggle
    // would cause concurrent readers to see empty results during the window.
    // The unique_lock already prevents map corruption between concurrent
    // writers; data_loaded_ (checked above) enforces the stricter
    // immutable-after-load contract for readers that run after a load.
    // Build-phase set() calls are single-threaded by contract (see
    // ADR-001), so mutating the map here is safe; just emit true at the end.

    // Store copies in manual_store_ (std::deque: push_back never invalidates
    // references to existing elements, so these string_views stay valid).
    // The GenomeIndex key is emplace()'d as id_view — a string_view into this
    // deque element — NEVER as the caller's `id` std::string: the caller's
    // buffer dies when addEntry() returns, and GenomeIndex keys are now
    // non-owning string_views (zero-copy refactor). manual_store_ is the
    // stable backing for both the key view and the SequenceView fields; it is
    // only destroyed by clearInternal(), which clears the map FIRST.
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
        std::get<GenomeIndex>(file_cache_).emplace(id_view, SequenceView{id_view, seq_view, qual_view});
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
    bool last_member_complete = true; // false when a member decode failed
    while (offset + 3 <= file_size &&
           compressed[offset] == 0x1f && compressed[offset + 1] == 0x8b &&
           compressed[offset + 2] == 0x08) {
        offsets.push_back(offset);

        z_stream strm = {};
        if (inflateInit2(&strm, 15 + 16) != Z_OK) { last_member_complete = false; break; }
        strm.next_in  = const_cast<Bytef*>(compressed + offset);
        strm.avail_in = static_cast<uInt>(file_size - offset);

        int ret = Z_OK;
        while (ret != Z_STREAM_END) {
            strm.next_out  = reinterpret_cast<Bytef*>(scratch.data());
            strm.avail_out = static_cast<uInt>(scratch.size());
            ret = inflate(&strm, Z_NO_FLUSH);
            if (ret == Z_STREAM_ERROR || ret == Z_DATA_ERROR || ret == Z_MEM_ERROR || ret == Z_BUF_ERROR) {
                // Not a genuine member boundary (or truncated stream) —
                // stop scanning; the trailing-data/truncation check below
                // rejects the file with an accurate error.
                break;
            }
        }

        size_t consumed = (file_size - offset) - strm.avail_in;
        inflateEnd(&strm);
        if (ret != Z_STREAM_END || consumed == 0) {
            last_member_complete = false;
            break;
        }
        offset += consumed;
    }

    // Trailing data after the last member. zlib-ng's gz_look() silently
    // discards bytes that do not begin a new member, so a file like
    // 'valid.gz || garbage || valid2.gz' would otherwise decompress only the
    // first part and serve it as complete (P0 data-integrity bug). Raw
    // inflate knows the exact member end, so any leftover bytes that do not
    // begin a valid member are rejected here. A member decode that failed to
    // reach Z_STREAM_END means the stream is truncated/corrupt — rejected
    // with an accurate message instead of letting the parallel decoder
    // balloon memory on repeated Z_BUF_ERROR before failing with a
    // misleading OOM error.
    if (!offsets.empty() && offset < file_size) {
        if (last_member_complete) {
            throw std::runtime_error(
                "GZIP file has trailing data after the last gzip stream: " + filepath);
        }
        throw std::runtime_error(
            "GZIP stream truncated or corrupt at offset " +
            std::to_string(offsets.back()) + ": " + filepath);
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
                // Z_BUF_ERROR fires for two distinct conditions: the output
                // buffer is full, OR the input is exhausted mid-stream
                // (strm.avail_in == 0). The latter means the deflate stream
                // can never reach Z_STREAM_END — a truncated/corrupt member.
                // The old code treated it as output-full and grew the buffer
                // geometrically forever, so a truncated tail member of a
                // concatenated gzip ballooned memory (~3.4 GB from a 23-byte
                // tail) until the aggregate budget guard tripped with a
                // misleading OOM error naming an innocent earlier stream.
                if (strm.avail_in == 0) {
                    inflateEnd(&strm);
                    thread_errors[stream_idx] = "GZIP stream truncated: unexpected end of input";
                    return;
                }
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

    // Pre-size text_arena_ with 3x heuristic and OOM guard. `raw_size` is
    // kept in scope for the trailing-garbage reconciliation below.
    std::error_code size_ec;
    const uintmax_t raw_size = std::filesystem::file_size(filepath, size_ec);
    if (!size_ec && raw_size > 0) {
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

    gzFile file = gzopen(filepath.c_str(), "rb");
    if (!file) throw std::runtime_error("Cannot open GZIP file: " + filepath);

    auto chunk = std::make_unique<char[]>(CHUNK_SIZE);
    size_t write_pos = 0;

    try {
        while (true) {
            size_t required = write_pos + CHUNK_SIZE;
            if (required > text_arena_.capacity()) {
                // Zip-bomb guard: pre-check the new size against available
                // memory BEFORE reserving. A post-hoc try/catch around
                // reserve() is not enough — under Linux overcommit the
                // allocation can succeed and the process gets OOM-killed
                // during the later memcpy. The initial reserve above caps at
                // avail_mem/4, but the geometric growth here had no budget
                // check of its own, so a tiny highly-compressible file could
                // balloon text_arena_ to the full machine memory.
                const size_t avail_mem_gz = getAvailableMemory();
                if (avail_mem_gz > 0 && required > avail_mem_gz / 2) {
                    throw std::runtime_error(
                        "SmartStrategy OOM guard: refusing to grow GZIP buffer to " +
                        std::to_string(required) + " bytes from " + filepath +
                        " (available memory ~" + std::to_string(avail_mem_gz) + " bytes)");
                }
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
            if (bytes_read == 0) {
                // gzread() returns 0 both at a clean EOF and when the stream
                // is truncated mid-deflate. zlib-ng's zlib-compat gzread does
                // NOT return -1 for a truncated stream (verified: gzread
                // returns partial data, then 0, with gzerror == Z_BUF_ERROR
                // 'unexpected end of file'), so the old code treated a
                // truncated stream as clean EOF and silently served the
                // partial data as complete (P0 data-integrity bug). Consult
                // gzerror and reject anything that did not end cleanly.
                int err = Z_OK;
                const char* error_msg = gzerror(file, &err);
                if (err != Z_OK) {
                    text_arena_.clear();
                    throw std::runtime_error(
                        "GZIP stream truncated or corrupt: " +
                        std::string(error_msg ? error_msg : ""));
                }
                break;
            }
            std::memcpy(text_arena_.data() + write_pos, chunk.get(), bytes_read);
            write_pos += static_cast<size_t>(bytes_read);
        }

        // ── Trailing-garbage check ─────────────────────────────────────────
        // zlib-ng's gz_look() silently discards bytes after the last valid
        // member that do not begin the gzip magic, reporting a clean Z_OK
        // EOF — so gzerror alone cannot catch 'valid.gz || garbage ||
        // valid2.gz', which would otherwise decompress only up to the garbage
        // and serve that as the complete file. gzoffset() reports the
        // compressed-file position of the next byte to be consumed; for a
        // well-formed file this equals the file size exactly (verified for
        // single-member and cleanly-concatenated files). Trailing data larger
        // than zlib-ng's internal 128 KiB read buffer leaves gzoffset short
        // of the file size and is rejected here. A smaller tail is absorbed
        // by that buffer during the final read and is indistinguishable from
        // a clean EOF through the gzread() API; files ≥ 1 MiB additionally
        // get exact coverage from scanGzipStreams(), which decodes every
        // member with raw inflate and knows the precise member end.
        if (!size_ec && raw_size > 0) {
            const z_off_t consumed = gzoffset(file);
            if (consumed >= 0 && static_cast<uintmax_t>(consumed) < raw_size) {
                text_arena_.clear();
                throw std::runtime_error(
                    "GZIP file has trailing data after the last gzip stream: " +
                    filepath);
            }
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
        const char* const seq_out_start = write; // write position before the sequence
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
        // Terminate the now-single-line sequence ONLY if the sequence actually
        // produced output bytes. Writing unconditionally clobbers the next
        // record's '>' marker when the sequence ends exactly at rec_end —
        // either because it abuts the next header with no newline
        // ('>a\nACGT>b\n...') or because the record is header-only
        // ('>a\n>b\nACGT\n', emitted by seqkit seq -m / bioawk filters) — and
        // the outer loop then breaks at the destroyed marker, silently dropping
        // every later record (P0 data-integrity bug). `write == read` here
        // means the sequence consumed no separator bytes: it ran up to the
        // next '>' (abutting header — no terminator may be written) or up to
        // EOF without a trailing newline (write into the +1 resize slot).
        if (write != seq_out_start && (write != read || read == end)) {
            *write++ = '\n';
        }
        // read == rec_end (pointing at '>' of next record, or end)
    }
    text_arena_.resize(static_cast<size_t>(write - text_arena_.data()));
}

void SmartStrategy::parseArena() {
    std::string_view content(text_arena_.data(), text_arena_.size());
    if (content.empty()) {
        markDataLoaded();
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
    markDataLoaded();
    data_ready_.store(true, std::memory_order_release);
}

void SmartStrategy::loadGzipFile(const std::string& filepath) {
    // Enforce the documented contract ("Throws if file is not valid GZIP").
    // zlib's gzread() transparently reads non-gzip files as raw data
    // (verified: gzread on a plain FASTA returns the bytes with Z_OK), so
    // without an explicit magic-byte check an invalid file silently loads
    // as plain text — contradicting both the header contract and the
    // data-integrity model (partial/garbage data is never served).
    // Matches loadFile()'s own magic-byte detection (0x1f 0x8b).
    //
    // Deliberately BEFORE clearInternal(): a pre-flight validation failure
    // (unreadable file / not GZIP) must PRESERVE the currently loaded
    // snapshot — the failed load is a no-op on cache state. Mid-stream
    // failures (truncation, trailing garbage, OOM) happen after teardown
    // and leave the cache empty; see the failure-atomicity contract in
    // SmartStrategy.h and the pre-loaded-cache tests in CacheTests.cpp.
    {
        std::ifstream check(filepath, std::ios::binary);
        if (!check) throw std::runtime_error("Cannot open file: " + filepath);
        unsigned char magic[2];
        check.read(reinterpret_cast<char*>(magic), 2);
        if (check.gcount() != 2 || magic[0] != 0x1f || magic[1] != 0x8b)
            throw std::runtime_error("File is not valid GZIP: " + filepath);
    }
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

// ── Per-thread reserve estimation ───────────────────────────────────────────
// ankerl::unordered_dense rehashes (re-inserts every element, O(n)) whenever
// a map outgrows its reserved capacity. The old fixed heuristics
// (chunk/500 for FASTA, chunk/600 for FASTQ) undershot badly on short reads:
// 100 bp WGS records are ~208 B, so a ~166K-record chunk (1 GB decompressed /
// 18 cores) was reserved at only ~99K slots → every thread mid-parse rehashed.
//
// Estimate the average record size from a small sample of the ACTUAL file so
// the per-thread reserve matches the chunk's real record count:
//   est_records_in_chunk = chunk_byte_size / avg_bytes_per_record
//   reserve(est_records_in_chunk * 1.25)
//
//   - FASTQ: strict 4-line framing ⇒ records = newlines / 4. LF and CRLF both
//     carry exactly one '\n' per line; empty seq/qual lines still count, and
//     no header/quality classification is needed ('@' inside quality lines is
//     irrelevant to a newline count).
//   - FASTA: normalizeFastaArena() has already collapsed sequences to a single
//     line when this runs, and every header starts with '>' at a line start,
//     so records = count of '>' line-starts in the sample.
static size_t estimateAvgRecordBytes(std::string_view content, bool isFastq) noexcept {
    constexpr size_t kSampleBytes = 1u << 20; // 1 MiB sample (~0.1 ms scan)
    const size_t span = std::min(content.size(), kSampleBytes);
    if (span == 0) return 1;
    size_t records = 0;
    if (isFastq) {
        size_t newlines = 0;
        const char* p = content.data();
        const char* end = p + span;
        while (true) {
            const char* nl = simd_find_char(p, end, '\n');
            if (nl >= end) break;
            ++newlines;
            p = nl + 1;
        }
        records = newlines / 4; // strict 4-line cycle
    } else {
        const char* p = content.data();
        const char* end = p + span;
        if (*p == '>') ++records;
        while (true) {
            const char* nl = simd_find_char(p, end, '\n');
            if (nl >= end) break;
            p = nl + 1;
            if (p < end && *p == '>') ++records;
        }
    }
    // Degenerate sample (no complete record / no header in the first MiB):
    // fall back to the internal parsers' long-standing assumptions.
    if (records == 0) return isFastq ? 200 : 150;
    return (span + records / 2) / records; // rounded avg bytes/record
}

// --- Templated Parsers ---
template <typename MapType> void SmartStrategy::parseFastaMultithreadedTemplate(std::string_view content, MapType& dest_map) {
    const size_t num_threads = std::max(1u, std::thread::hardware_concurrency());
    const size_t chunk_size = content.size() / num_threads;
    std::vector<std::thread> threads;
    std::vector<MapType> thread_caches(num_threads);
    // Pre-reserve thread-local maps to avoid mid-parse rehashing. The reserve
    // is derived from the chunk's actual byte span and a measured average
    // record size (1.25x headroom) instead of a fixed bytes/record constant
    // that silently undershot on short records (see estimateAvgRecordBytes).
    {
        const size_t avg_record = estimateAvgRecordBytes(content, /*isFastq=*/false);
        const size_t est_per_thread = chunk_size / avg_record;
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
                if constexpr (std::is_same_v<MapType, GenomeIndex>) { thread_caches[thread_id].emplace(id, SequenceView{id, seq, {}}); }
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
    // Pre-reserve thread-local maps to avoid mid-parse rehashing. The reserve
    // is derived from the chunk's actual byte span and a measured average
    // record size (1.25x headroom) instead of a fixed bytes/record constant
    // that silently undershot on short records (see estimateAvgRecordBytes).
    {
        const size_t avg_record = estimateAvgRecordBytes(content, /*isFastq=*/true);
        const size_t est_per_thread = chunk_size / avg_record;
        for (auto& cache : thread_caches) {
            cache.reserve(static_cast<size_t>(est_per_thread * 1.25));
        }
    }
    auto worker = [&](size_t thread_id, size_t start, size_t end) {
        const char* ptr = content.data() + start;
        const char* chunk_end = content.data() + end;
        const char* global_end = content.data() + content.size();
        if (thread_id > 0) {
            // ── Chunk-boundary scan (Bug 1 fix) ─────────────────────────────
            // Find the first GENUINE FASTQ record header at/after the chunk
            // start using a FORWARD classification. A genuine header always
            // has its sequence on the next line and the '+' separator on the
            // line after that (4-line record: header / seq / '+' / qual), so
            // a candidate line-start '@' is a genuine header iff the line
            // TWO lines later begins with '+'. A quality line that happens
            // to start with '@' (Phred Q31+ scores) is the 4th line of its
            // record, so its line+1 is the NEXT record's header and its
            // line+2 is that record's sequence — which never starts with
            // '+' — so it is correctly skipped.
            //
            // The earlier classifier walked BACKWARD and tested the previous
            // line's first byte; a quality line that itself starts with '+'
            // (Phred+33 Q10, ASCII 43) made a genuine following header look
            // like a quality line, so every worker after chunk 0 silently
            // dropped its whole chunk (repro: 80k reads with '+'-leading
            // quality → 4,468 loaded, 94.4% loss). The forward test has no
            // such false negative: for a genuine header, seq / '+' / qual
            // always follow in a valid FASTQ file.
            //
            // Align to the start of the line containing the chunk boundary.
            while (ptr > content.data() && *(ptr - 1) != '\n') --ptr;
            while (ptr < chunk_end) {
                if (*ptr == '@') {
                    // Look two lines ahead (reads stay within the buffer:
                    // simd_find_char is bounded by global_end; the lookahead
                    // is never dereferenced past global_end). Exactly ONE line
                    // terminator is consumed per boundary ('\n', with the
                    // optional preceding '\r' already part of the previous
                    // line): collapsing runs of newlines here would skip an
                    // empty sequence/quality line and misclassify a genuine
                    // header (or a quality '@') whenever a zero-length read
                    // is involved.
                    const char* next = simd_find_char(ptr, global_end, '\n'); // end of '@' line
                    if (next < global_end) ++next; // consume exactly one terminator → start of line+1
                    const char* plus = simd_find_char(next, global_end, '\n'); // end of line+1
                    if (plus < global_end) ++plus; // consume exactly one terminator → start of line+2
                    if (plus < global_end && *plus == '+') break; // genuine header
                    // Line+2 is not '+' → this '@' is a quality line, not a
                    // header — keep scanning for the next candidate.
                }
                // Advance to the start of the next line.
                while (ptr < chunk_end && *ptr != '\n') ++ptr;
                if (ptr < chunk_end) ++ptr; // consume '\n'
                if (ptr < chunk_end && *ptr == '\n') ++ptr; // consume exactly one more terminator
            }
        }
        while (ptr < chunk_end) {
            // Strict 4-line FASTQ cycle: header, sequence, '+', quality.
            // Quality lines are consumed as the 4th line of the cycle and
            // never re-scanned, so '@' inside quality (Phred Q31+) cannot be
            // mistaken for a record header. The '+' separator is verified so
            // a malformed file can't shift the cycle either.
            while (ptr < chunk_end && (*ptr == '\n' || *ptr == '\r' || *ptr == ' ')) ++ptr;
            if (ptr >= chunk_end || *ptr != '@') break;
            ++ptr; const char* id_start = ptr;
            while (ptr < global_end && *ptr != ' ' && *ptr != '\t' && *ptr != '\n' && *ptr != '\r') ++ptr;
            const char* id_end = ptr;
            ptr = simd_find_char(ptr, global_end, '\n');
            if (ptr < global_end && *ptr == '\n') ++ptr; // exactly one terminator per boundary
            const char* seq_start = ptr;
            ptr = simd_find_char(ptr, global_end, '\n');
            const char* seq_end = ptr;
            while (seq_end > seq_start && *(seq_end - 1) == '\r') --seq_end;
            if (ptr < global_end && *ptr == '\n') ++ptr; // exactly one terminator per boundary
            if (ptr >= global_end || *ptr != '+') break; // missing '+' separator
            ptr = simd_find_char(ptr, global_end, '\n');
            if (ptr < global_end && *ptr == '\n') ++ptr; // exactly one terminator per boundary
            const char* qual_start = ptr;
            ptr = simd_find_char(ptr, global_end, '\n');
            const char* qual_end = ptr;
            while (qual_end > qual_start && *(qual_end - 1) == '\r') --qual_end;
            if (ptr < global_end && *ptr == '\n') ++ptr; // exactly one terminator per boundary
            std::string_view id(id_start, id_end - id_start);
            std::string_view seq(seq_start, seq_end - seq_start);
            std::string_view qual(qual_start, qual_end - qual_start);
            if (!id.empty() && !seq.empty()) {
                if constexpr (std::is_same_v<MapType, GenomeIndex>) { thread_caches[thread_id].emplace(id, SequenceView{id, seq, qual}); }
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
            if constexpr (std::is_same_v<MapType, GenomeIndex>) { map.emplace(id, SequenceView{id, seq, {}}); }
            else { map.emplace(hash_key(id), SequenceView{id, seq, {}}); }
        }
    }
}

template <typename MapType> void SmartStrategy::parseFastqInternal(std::string_view content, MapType& map) {
    const char* ptr = content.data(); const char* end = ptr + content.size();
    size_t estimated_records = (content.size() < 50 * 1024 * 1024) ? content.size() / 150 : content.size() / 200;
    map.reserve(static_cast<size_t>(estimated_records * 1.25));
    while (ptr < end) {
        // Strict 4-line FASTQ cycle (header, sequence, '+', quality). The
        // quality line is consumed as the 4th line of the cycle and never
        // re-scanned, so '@' inside quality (Phred Q31+ scores ≥ 64) cannot
        // be misparsed as a new record header — the root cause of the record
        // loss fixed for the multithreaded path. The '+' separator is
        // verified so a malformed file can't shift the cycle either.
        while (ptr < end && (*ptr == '\n' || *ptr == '\r' || *ptr == ' ')) ++ptr;
        if (ptr >= end || *ptr != '@') break;
        ++ptr; const char* id_start = ptr; while (ptr < end && *ptr != ' ' && *ptr != '\t' && *ptr != '\n' && *ptr != '\r') ++ptr; const char* id_end = ptr;
        ptr = simd_find_char(ptr, end, '\n'); if (ptr < end && *ptr == '\n') ++ptr; // exactly one terminator per boundary
        const char* seq_start = ptr; ptr = simd_find_char(ptr, end, '\n'); const char* seq_end = ptr;
        while (seq_end > seq_start && *(seq_end - 1) == '\r') --seq_end;
        if (ptr < end && *ptr == '\n') ++ptr; // exactly one terminator per boundary
        if (ptr >= end || *ptr != '+') break; // missing '+' separator
        ptr = simd_find_char(ptr, end, '\n'); if (ptr < end && *ptr == '\n') ++ptr; // exactly one terminator per boundary
        const char* qual_start = ptr; ptr = simd_find_char(ptr, end, '\n'); const char* qual_end = ptr;
        while (qual_end > qual_start && *(qual_end - 1) == '\r') --qual_end;
        if (ptr < end && *ptr == '\n') ++ptr; // exactly one terminator per boundary
        std::string_view id(id_start, id_end - id_start); std::string_view seq(seq_start, seq_end - seq_start); std::string_view qual(qual_start, qual_end - qual_start);
        if (!id.empty() && !seq.empty()) {
            if constexpr (std::is_same_v<MapType, GenomeIndex>) { map.emplace(id, SequenceView{id, seq, qual}); }
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

    // v4 format: stream-serialize + stream-compress via LZ4 Frame (LZ4F)
    // instead of materializing a full serialized "payload" copy plus a full
    // "compressed" copy of the whole dataset (the old v2 path's ~3x peak
    // memory: text_arena_ + payload + compressed, all resident at once).
    // Peak extra memory here is bounded by STREAM_CHUNK_SIZE regardless of
    // dataset size. The v4 header adds a whole-payload CRC32C: the checksum
    // accumulator is updated as payload chunks pass through the compressor
    // (below), so there is no second full-payload pass and no extra
    // allocation. The frame length and checksum are unknown until streaming
    // finishes, so the header is written with placeholders and patched in
    // place afterwards (seekp) — the frame itself is never buffered.
    const size_t estimated_size = estimatePayloadSize();
    const CompressionMode comp_mode = selectCompressionStrategy(estimated_size);

    char header[V4_HEADER_LEN];
    std::memcpy(header, V4_MAGIC, 4);
    header[4] = V4_CODEC_LZ4F;
    header[5] = std::holds_alternative<NGSIndex>(file_cache_) ? 1 : 0;
    write_le64(header + 6, static_cast<uint64_t>(estimated_size));
    write_le64(header + 14, 0); // placeholder: frame length (patched below)
    write_le32(header + 22, 0); // placeholder: CRC32C (patched below)
    out.write(header, V4_HEADER_LEN);

    // Streaming CRC-32C: ENTIRE uncompressed payload first (fed here, chunk
    // by chunk, as it passes through the compressor), then the canonical
    // header bytes [0..22) once the frame length is known (patched below).
    // The payload-first ordering is what makes this single-pass possible:
    // frame_len is only known after the frame is complete, so the header is
    // fed to the checksum last on BOTH save and load.
    Crc32c crc;

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

    // The CRC is updated here, over the UNCOMPRESSED payload bytes, as they
    // pass through on their way into the LZ4F compressor — the required
    // "streaming CRC" with no second pass over the payload.
    serializePayload([&](const char* data, size_t len) {
        crc.update(data, len);
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

    // Patch the placeholder header fields now that the frame is complete.
    // The checksum is finished with the canonical header bytes [0..22)
    // exactly as they appear on disk: frame_len is written into the in-memory
    // header first, then both header regions are fed to the accumulator.
    const std::streampos frame_end = out.tellp();
    if (frame_end < static_cast<std::streampos>(V4_HEADER_LEN))
        throw std::runtime_error("Binary cache write failed: " + filepath);
    const uint64_t frame_len = static_cast<uint64_t>(frame_end) - V4_HEADER_LEN;
    write_le64(header + 14, frame_len);
    crc.update(header, 14);     // magic + codec flags + index mode + logical_len
    crc.update(header + 14, 8); // frame_len, as stored on disk
    write_le32(header + 22, crc.finalize());
    out.seekp(static_cast<std::streampos>(14));
    out.write(header + 14, 12); // frame_len(8) + CRC32C(4)
    out.seekp(frame_end); // restore EOF position

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

    // ── Failure atomicity ───────────────────────────────────────────────────
    // Every malformed-cache path below throws. The whole reader body is
    // wrapped so an invalid load NEVER publishes data_ready_ and NEVER leaves
    // partial map/arena state behind: on any exception the cache is reset to
    // the pristine empty state (map cleared, arena released, mmap dropped,
    // flags cleared). A failed restore therefore looks exactly like a fresh
    // SmartStrategy, never like a half-loaded cache.
    try {
        // Bounds-checked read helper. Uses memcpy to avoid unaligned-access UB
        // on ARM. `n > end - ptr` instead of `ptr + n > end`: n is
        // attacker-controlled (up to uint64_t from the count/compressed-size
        // fields), and the addition could wrap near 2^64 — letting the check
        // pass and turning the subsequent `ptr += n` into pointer-arithmetic
        // UB. `end - ptr` is exact and cannot overflow.
        auto safe_advance = [&](size_t n) -> const char* {
            if (n > static_cast<size_t>(end - ptr))
                throw std::runtime_error("Binary cache is truncated or corrupt: " + filepath);
            const char* result = ptr;
            ptr += n;
            return result;
        };

        // ── Header (26 bytes, little-endian fields) ─────────────────────────
        if (mmap_handle_->size < V4_HEADER_LEN)
            throw std::runtime_error("Binary cache too small to contain header: " + filepath);

        const char* magic = safe_advance(4);
        if (std::strncmp(magic, "TRO", 3) != 0)
            throw std::runtime_error("Invalid binary cache magic bytes: " + filepath);
        // v2.0.0 ships v4 as the ONLY readable format (design review Q1/Q5:
        // no legacy readers, clean break). v1/v2/v3 caches must be regenerated.
        if (magic[3] != 0x04)
            throw std::runtime_error(
                "unsupported cache version; regenerate with v2.0.0: " + filepath);

        const uint8_t flags = static_cast<uint8_t>(*safe_advance(1));
        // Codec flags: only the LZ4 Frame bit is defined; flags == 0 (no
        // codec) and unknown future bits are both rejected.
        if ((flags & V4_CODEC_MASK) != V4_CODEC_LZ4F || (flags & ~V4_CODEC_MASK) != 0)
            throw std::runtime_error("Binary cache declares unsupported codec flags: " + filepath);

        const uint8_t mode = static_cast<uint8_t>(*safe_advance(1));
        if (mode > 1)
            throw std::runtime_error("Binary cache declares invalid index mode: " + filepath);
        // Sync index_mode_ with what's stored in the binary so that subsequent
        // loadFile() / parseArena() calls use the same index type as restored data.
        index_mode_ = (mode == 1) ? IndexMode::NGS : IndexMode::GENOME;

        const uint64_t logical_len = read_le64(safe_advance(8));
        const uint64_t frame_len = read_le64(safe_advance(8));
        const uint32_t stored_crc = read_le32(safe_advance(4));
        const char* const header_begin = static_cast<const char*>(mmap_handle_->data);

        // ── Hardening: implausible sizes rejected BEFORE any allocation ─────
        // Every serialized payload begins with the 8-byte record count, so a
        // smaller logical_len is always malformed. v4 uses the 64-bit-clean
        // LZ4F streaming API, so the v2 class of INT_MAX int-truncation bugs
        // (2^32+N -> N silently passing the size check) cannot occur here; the
        // equivalent protections are this floor check, the 32-bit size_t
        // guard, and the OOM guard below — all applied before resize().
        if (logical_len < sizeof(uint64_t))
            throw std::runtime_error("Binary cache declares implausible decompressed size: " + filepath);
        if (frame_len == 0)
            throw std::runtime_error("Binary cache declares empty compressed frame: " + filepath);
        if (logical_len > static_cast<uint64_t>(std::numeric_limits<size_t>::max()))
            throw std::runtime_error("Binary cache declares implausible decompressed size: " + filepath);

        // Subtraction-form bounds check for the compressed frame: no pointer
        // arithmetic overflow possible even for attacker-controlled frame_len
        // (up to 2^64-1).
        if (frame_len > static_cast<uint64_t>(end - ptr))
            throw std::runtime_error("Binary cache is truncated: " + filepath);

        // OOM guard: refuse to allocate a multi-GB arena from a small file
        // (mirrors loadGzipSingleStream()/loadFile()).
        const size_t avail_mem = getAvailableMemory();
        if (avail_mem > 0 && logical_len > avail_mem / 2) {
            throw std::runtime_error(
                "SmartStrategy OOM guard: refusing to decompress " +
                std::to_string(logical_len) + " bytes from " + filepath +
                " (available memory ~" + std::to_string(avail_mem) + " bytes)");
        }
        try {
            text_arena_.resize(logical_len);
        } catch (const std::bad_alloc&) {
            std::cerr << "SmartStrategy OOM: failed to allocate " << logical_len
                      << " bytes for " << filepath << std::endl;
            throw;
        }

        // ── Streaming LZ4 Frame decompression (bounded memory) ──────────────
        LZ4F_dctx* dctx = nullptr;
        if (LZ4F_isError(LZ4F_createDecompressionContext(&dctx, LZ4F_VERSION)))
            throw std::runtime_error("LZ4F_createDecompressionContext failed: " + filepath);
        struct DctxGuard {
            LZ4F_dctx* d;
            ~DctxGuard() { if (d) LZ4F_freeDecompressionContext(d); }
        } dguard{dctx};

        const char* src = ptr;
        size_t src_remaining = static_cast<size_t>(frame_len);
        char* dst = text_arena_.data();
        size_t dst_remaining = text_arena_.size();
        size_t hint = 1;

        // Streaming CRC32C: the ENTIRE decompressed payload first (each chunk
        // as it is written to text_arena_ — no second full-payload pass, no
        // extra allocation), then the canonical header bytes [0..22) once the
        // frame has been fully consumed. Payload-first ordering matches
        // saveBinary() (see the v4 format comment at the top of this file).
        Crc32c crc;

        while (hint != 0 && src_remaining > 0) {
            size_t dst_size = dst_remaining;
            size_t src_size = src_remaining;
            hint = LZ4F_decompress(dctx, dst, &dst_size, src, &src_size, nullptr);
            if (LZ4F_isError(hint))
                throw std::runtime_error("LZ4F_decompress failed for binary cache: " + filepath);
            crc.update(dst, dst_size);
            dst += dst_size;
            dst_remaining -= dst_size;
            src += src_size;
            src_remaining -= src_size;
            if (dst_size == 0 && src_size == 0 && hint != 0) break; // no progress; malformed frame
        }

        // Exact-length requirement: the declared logical length must equal the
        // actual decompressed size. A truncated file that loses its tail (or a
        // wrong logical length) fails here, not on the checksum — a digest
        // cannot know that a missing tail was supposed to exist.
        if (dst_remaining != 0)
            throw std::runtime_error("Binary cache decompressed size mismatch: " + filepath);
        // Complete-frame requirement: the LZ4 Frame must end EXACTLY at the
        // declared frame boundary. hint == 0 means the frame end was reached;
        // leftover src bytes mean the frame ended before the declared
        // boundary (trailing garbage inside the declared frame length).
        if (src_remaining != 0 || hint != 0)
            throw std::runtime_error(
                "Binary cache frame not terminated at declared boundary: " + filepath);

        // ── Whole-payload integrity check ───────────────────────────────────
        // CRC32C over the ENTIRE uncompressed logical payload + the canonical
        // header fields (checksum field excluded). A mismatch means the file
        // was modified/corrupted after writing (accidental-corruption
        // detection, NOT authentication).
        crc.update(header_begin, V4_CRC_OFFSET);
        if (crc.finalize() != stored_crc)
            throw std::runtime_error(
                "Binary cache checksum mismatch (corrupt or modified file): " + filepath);

        // Release mmap handle since we now own the data in text_arena_
        mmap_handle_.reset();

        // ── Parse from text_arena_ (payload layout unchanged from v1–v3) ───
        const char* arena_ptr = text_arena_.data();
        const char* arena_end = arena_ptr + text_arena_.size();

        auto arena_safe_advance = [&](size_t n) -> const char* {
            if (n > static_cast<size_t>(arena_end - arena_ptr))
                throw std::runtime_error("Binary cache v4 payload is corrupt: " + filepath);
            const char* result = arena_ptr;
            arena_ptr += n;
            return result;
        };

        uint64_t count;
        std::memcpy(&count, arena_safe_advance(sizeof(uint64_t)), sizeof(uint64_t));

        // Count-bounded reserve: every serialized record needs at least
        // MIN_RECORD_BYTES (GENOME: 3×u32 length fields; NGS: u64 hash +
        // 3×u32 — see serializePayload()), so a count exceeding remaining_bytes
        // / MIN_RECORD_BYTES is impossible in a genuine cache. Subsumes the
        // old 1e9 cap, which still allowed a ~80 GB reserve().
        const uint64_t remaining_bytes = static_cast<uint64_t>(arena_end - arena_ptr);
        const uint64_t min_record_bytes = (mode == 0) ? 12 : 20;
        if (count > remaining_bytes / min_record_bytes)
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
                map.emplace(id_view, SequenceView{id_view, seq, qual});
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
    } catch (...) {
        // Failure atomicity: an invalid load must NOT publish data_ready_ or
        // leave partial map state behind (design review P0 sprint item 4).
        clearInternal();
        throw;
    }

    determine_format_from_cache();
    // One-time O(n) walk to seed serialized_size_estimate_ for a subsequent
    // saveBinary() (e.g. re-saving after restoring + addEntry()'ing more
    // data); addEntry() maintains it incrementally after this point.
    refreshPayloadEstimate();
    markDataLoaded();
    data_ready_.store(true, std::memory_order_release);
}


// --- ACCESS ---

std::string_view SmartStrategy::getView(const std::string_view& sequence_id) const {
#ifdef TRACEON_DEBUG_LIFECYCLE
    // Reader-quiescence diagnostic (opt-in, debug builds). Tracks whether any
    // getView() lookup is in flight when a clear/reload tears the snapshot
    // down. LIMITATION: this catches a reader whose LOOKUP overlaps the clear;
    // it cannot detect a caller that retained the returned std::string_view
    // and dereferences it after clearCache()/reload/destruction completes.
    // Diagnostic only — NOT synchronization (see ADR-001 / README).
    struct ReaderScope {
        std::atomic<size_t>& c;
        explicit ReaderScope(std::atomic<size_t>& c_) : c(c_) {
            c.fetch_add(1, std::memory_order_relaxed);
        }
        ~ReaderScope() { c.fetch_sub(1, std::memory_order_relaxed); }
    } scope(active_readers_);
#endif
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

} // namespace