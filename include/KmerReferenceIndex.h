#ifndef TRACEON_KMER_REFERENCE_INDEX_H
#define TRACEON_KMER_REFERENCE_INDEX_H

// Multi-value, mmap-persistent k-mer reference index.
//
// The problem this solves: alignment/analysis tools query the SAME
// reference genome index over and over across runs, and rebuilding (or
// even just re-parsing/re-inserting into a fresh hash table) that index
// every run wastes enormous time -- Phase 2 measured ~47-52s just to
// build minimap2's minimizer index over the human genome. This class
// makes that a ONE-TIME cost: build once, save() once, and every later
// run mmap_open()s the file -- the on-disk byte layout IS the in-memory
// query layout, so loading is a single mmap() + pointer cast, with no
// parsing, no decompression, no rehashing, no per-record insertion.
//
// This is deliberately NOT built on ankerl::unordered_dense (used
// elsewhere in TracEon, e.g. KmerIndex.h) -- that map's internal
// representation isn't a stable, position-independent byte layout you
// can mmap back verbatim. Instead this is a hand-rolled flat
// open-addressing hash table (linear probing, power-of-two size) whose
// slot array is plain old data, safe to write to disk and map back as-is.
//
// ── TRKI on-disk format (v1) ────────────────────────────────────────────────
// All integers are LITTLE-ENDIAN on disk (explicit byte serialization --
// the format is endian-clean; a big-endian host decodes the tables into
// owned buffers instead of zero-copy pointer casts).
//
//   magic      "TRKI"                        (4 bytes)
//   version    u32 LE  == 1
//   k          u32 LE  k-mer length, [1, 32]
//   reserved   u32 LE  (0)
//   n_slots    u64 LE  slots table length (power of two, or 0 = empty)
//   n_keys     u64 LE  distinct stored k-mers (<= n_slots)
//   n_positions u64 LE positions array length
//   slots      n_slots x (key u64 LE, value u64 LE)   -- KmerSlot
//   positions  n_positions x u64 LE
//
// Total file size MUST equal 40 + n_slots*16 + n_positions*8 exactly;
// any mismatch (truncation OR trailing garbage) is rejected at load.
//
// ── Security contract (vuln-0006) ──────────────────────────────────────────
// mmap_open() validates every file-controlled field with overflow-safe
// arithmetic BEFORE computing any product or pointer: n_slots and
// n_positions are checked against the actual mapped length via the
// division form (n_slots > bytes_after_header / sizeof(KmerSlot) => reject),
// n_slots must be 0 or a power of two, n_keys <= n_slots, k in [1,32],
// the declared tables must match the file size EXACTLY, and every
// non-empty slot's (offset, count) must lie within n_positions (walked at
// load time with explicit LE reads). query() additionally re-validates the
// slot's (offset, count) before building the return span (defense in
// depth). Malformed/crafted files fail loudly with std::runtime_error --
// never an out-of-bounds read, never a wild pointer.

#include "KmerEncoding.h"

#include <cstdint>
#include <cstddef>
#include <cstring>
#include <string>
#include <string_view>
#include <vector>
#include <span>
#include <algorithm>
#include <memory>
#include <stdexcept>

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

namespace detail {

// Minimal mmap-a-whole-file-read-only wrapper, independent of
// SmartStrategy's private MMapHandle (src/SmartStrategy.cpp) so this
// header has no dependency on the rest of the library.
struct ReadOnlyMapping {
    void* data = nullptr;
    size_t size = 0;
#ifdef _WIN32
    HANDLE hFile = INVALID_HANDLE_VALUE;
    HANDLE hMap = nullptr;
#else
    int fd = -1;
#endif

    ReadOnlyMapping() = default;
    ReadOnlyMapping(const ReadOnlyMapping&) = delete;
    ReadOnlyMapping& operator=(const ReadOnlyMapping&) = delete;
    ReadOnlyMapping(ReadOnlyMapping&& other) noexcept { *this = std::move(other); }
    ReadOnlyMapping& operator=(ReadOnlyMapping&& other) noexcept {
        if (this != &other) {
            cleanup();
            data = other.data; size = other.size;
#ifdef _WIN32
            hFile = other.hFile; hMap = other.hMap;
            other.hFile = INVALID_HANDLE_VALUE; other.hMap = nullptr;
#else
            fd = other.fd;
            other.fd = -1;
#endif
            other.data = nullptr; other.size = 0;
        }
        return *this;
    }
    ~ReadOnlyMapping() { cleanup(); }

    void cleanup() {
#ifdef _WIN32
        if (data) UnmapViewOfFile(data);
        if (hMap) CloseHandle(hMap);
        if (hFile != INVALID_HANDLE_VALUE) CloseHandle(hFile);
        data = nullptr; hMap = nullptr; hFile = INVALID_HANDLE_VALUE;
#else
        if (data && data != MAP_FAILED) munmap(data, size);
        if (fd != -1) close(fd);
        data = nullptr; fd = -1;
#endif
    }

    static ReadOnlyMapping open(const std::string& path) {
        ReadOnlyMapping m;
#ifdef _WIN32
        m.hFile = CreateFileA(path.c_str(), GENERIC_READ, FILE_SHARE_READ, nullptr,
                              OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, nullptr);
        if (m.hFile == INVALID_HANDLE_VALUE)
            throw std::runtime_error("KmerReferenceIndex: cannot open " + path);
        LARGE_INTEGER sz;
        GetFileSizeEx(m.hFile, &sz);
        m.size = (size_t)sz.QuadPart;
        m.hMap = CreateFileMappingA(m.hFile, nullptr, PAGE_READONLY, 0, 0, nullptr);
        if (!m.hMap) throw std::runtime_error("KmerReferenceIndex: CreateFileMapping failed");
        m.data = MapViewOfFile(m.hMap, FILE_MAP_READ, 0, 0, m.size);
        if (!m.data) throw std::runtime_error("KmerReferenceIndex: MapViewOfFile failed");
#else
        m.fd = ::open(path.c_str(), O_RDONLY);
        if (m.fd == -1) throw std::runtime_error("KmerReferenceIndex: cannot open " + path);
        struct stat st{};
        if (fstat(m.fd, &st) != 0) throw std::runtime_error("KmerReferenceIndex: fstat failed");
        m.size = (size_t)st.st_size;
        m.data = mmap(nullptr, m.size, PROT_READ, MAP_PRIVATE, m.fd, 0);
        if (m.data == MAP_FAILED) throw std::runtime_error("KmerReferenceIndex: mmap failed");
#endif
        return m;
    }
};

// ── Explicit little-endian byte I/O ─────────────────────────────────────────
// Used for ALL TRKI header/table disk I/O so the on-disk format is
// endian-clean and validation never depends on host byte order.

constexpr size_t kTrkiHeaderSize = 40; // 4 magic + 3*u32 + 3*u64

inline uint32_t load_u32_le(const unsigned char* p) noexcept {
    return (uint32_t)p[0] | ((uint32_t)p[1] << 8) | ((uint32_t)p[2] << 16) | ((uint32_t)p[3] << 24);
}

inline uint64_t load_u64_le(const unsigned char* p) noexcept {
    return (uint64_t)p[0] | ((uint64_t)p[1] << 8) | ((uint64_t)p[2] << 16) |
           ((uint64_t)p[3] << 24) | ((uint64_t)p[4] << 32) | ((uint64_t)p[5] << 40) |
           ((uint64_t)p[6] << 48) | ((uint64_t)p[7] << 56);
}

inline void store_u32_le(unsigned char* p, uint32_t v) noexcept {
    p[0] = (unsigned char)(v & 0xFFu);
    p[1] = (unsigned char)((v >> 8) & 0xFFu);
    p[2] = (unsigned char)((v >> 16) & 0xFFu);
    p[3] = (unsigned char)((v >> 24) & 0xFFu);
}

inline void store_u64_le(unsigned char* p, uint64_t v) noexcept {
    for (int i = 0; i < 8; ++i) p[i] = (unsigned char)((v >> (8 * i)) & 0xFFu);
}

// Compile-time host endianness: little-endian hosts get zero-copy pointer
// casts into the mapping; big-endian hosts decode the LE tables into owned
// buffers (correct, at the cost of one copy at load).
constexpr bool is_little_endian_host() noexcept {
#if defined(__BYTE_ORDER__) && defined(__ORDER_BIG_ENDIAN__) && \
    (__BYTE_ORDER__ == __ORDER_BIG_ENDIAN__)
    return false;
#else
    return true;
#endif
}

} // namespace detail

// On-disk / in-memory slot layout for the flat hash table. Plain old
// data: no pointers, no padding-sensitive members beyond natural uint64
// alignment -- safe to mmap back verbatim on the same/compatible platform.
struct KmerSlot {
    uint64_t key;   // packed 2-bit-per-base k-mer; UINT64_MAX reserved as "empty"
    uint64_t value; // (offset << 32) | count -- index + length into the positions[] array
};

static constexpr uint64_t KMER_SLOT_EMPTY = ~0ULL;

// splitmix64 finalizer: mixes the packed k-mer's bits before using them as
// a bucket index. Required because raw 2-bit-packed k-mers are NOT
// uniformly distributed in their low bits when inserted in sorted order
// (adjacent k-mers along a sequence differ by a 2-bit shift, so sorted
// unique keys cluster on low bits) -- indexing directly by `key & mask`
// on sorted input causes pathological linear-probe pileups (observed:
// effectively O(n^2) insertion on a real chromosome). Must be applied
// identically at insert time and query time.
inline uint64_t kmer_slot_hash(uint64_t x) {
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    x = x ^ (x >> 31);
    return x;
}

// Documented in-memory header struct (native layout, informational).
// The on-disk representation is the explicit-LE layout in detail:: above.
struct KmerReferenceIndexHeader {
    char magic[4];      // "TRKI"
    uint32_t version;   // format version
    uint32_t k;         // k-mer length this index was built with
    uint32_t reserved;  // padding/future use
    uint64_t n_slots;   // size of the slots[] table (power of two)
    uint64_t n_keys;    // number of distinct k-mers actually stored
    uint64_t n_positions; // total length of the positions[] array
};
static_assert(sizeof(KmerReferenceIndexHeader) == detail::kTrkiHeaderSize,
              "TRKI header must be exactly 40 bytes");
static_assert(sizeof(size_t) >= 8,
              "KmerReferenceIndex: 64-bit size_t required (TRKI uses 64-bit offsets)");

// Multi-value packed-k-mer -> {positions...} index, backed by a flat
// open-addressing hash table (linear probing) whose byte layout is
// identical whether it's freshly built in RAM or mmap'd from disk.
class KmerReferenceIndex {
public:
    KmerReferenceIndex() = default;

    // --- Build phase (owns its storage) ---

    // Builds the index over every valid k-mer window across `seq`,
    // recording EVERY occurrence (multi-value), using the rolling 2-bit
    // encoder so build is O(1) amortized per position. k must be in
    // [1, 32]; invalid k (or a sequence shorter than k) yields an empty
    // index instead of undefined behavior.
    void build_from_sequence(std::string_view seq, int k) {
        k_ = k;
        if (k <= 0 || k > 32 || (int)seq.size() < k) { finalize_empty(); return; }

        const uint64_t mask = (k >= 32) ? ~0ULL : ((1ULL << (2 * k)) - 1);
        uint64_t code = 0;
        int valid_run = 0;

        std::vector<std::pair<uint64_t, uint64_t>> kmer_pos; // (kmer, position)
        kmer_pos.reserve(seq.size() >= (size_t)k ? seq.size() - k + 1 : 0);

        for (size_t i = 0; i < seq.size(); ++i) {
            int b = kmer_base_code(seq[i]);
            if (b < 0) { valid_run = 0; code = 0; continue; }
            code = ((code << 2) | (uint64_t)b) & mask;
            valid_run++;
            if (valid_run >= k) kmer_pos.emplace_back(code, i - (size_t)k + 1);
        }

        build_from_pairs(std::move(kmer_pos));
    }

    // --- Query (works identically whether owned or mmap'd) ---

    // Returns all recorded positions for `kmer` (empty span on miss).
    // Every file-supplied slot's (offset, count) is re-validated against
    // n_positions_ before the span is built (defense in depth on top of
    // the load-time validation).
    std::span<const uint64_t> query(uint64_t kmer) const {
        if (n_slots_ == 0) return {};
        size_t mask = n_slots_ - 1;
        size_t idx = (size_t)(kmer_slot_hash(kmer) & mask);
        for (size_t probes = 0; probes < n_slots_; ++probes) {
            const KmerSlot& s = slots_[idx];
            if (s.key == KMER_SLOT_EMPTY) return {};
            if (s.key == kmer) {
                uint32_t count = (uint32_t)(s.value & 0xFFFFFFFFu);
                uint64_t offset = s.value >> 32;
                if (count == 0 || offset > n_positions_ || count > n_positions_ - offset)
                    return {};
                return std::span<const uint64_t>(positions_ + offset, count);
            }
            idx = (idx + 1) & mask;
        }
        return {};
    }

    size_t n_keys() const { return n_keys_; }
    int k() const { return k_; }

    // --- Persistence: save() writes the exact in-memory layout; the
    // saved file can be mmap_open()'d back with zero parsing/copying.
    // The on-disk format is explicit little-endian (see format block at
    // the top of this file); on a little-endian host the bytes written
    // are identical to the native in-memory layout. ---

    void save(const std::string& path) const {
        FILE* f = std::fopen(path.c_str(), "wb");
        if (!f) throw std::runtime_error("KmerReferenceIndex: cannot write " + path);

        unsigned char hdr[detail::kTrkiHeaderSize];
        std::memcpy(hdr, "TRKI", 4);
        detail::store_u32_le(hdr + 4, 1);              // version
        detail::store_u32_le(hdr + 8, (uint32_t)k_);   // k
        detail::store_u32_le(hdr + 12, 0);             // reserved
        detail::store_u64_le(hdr + 16, (uint64_t)n_slots_);
        detail::store_u64_le(hdr + 24, (uint64_t)n_keys_);
        detail::store_u64_le(hdr + 32, (uint64_t)n_positions_);
        std::fwrite(hdr, 1, sizeof(hdr), f);

        for (size_t i = 0; i < n_slots_; ++i) {
            unsigned char slot[2 * sizeof(uint64_t)];
            detail::store_u64_le(slot, slots_[i].key);
            detail::store_u64_le(slot + 8, slots_[i].value);
            std::fwrite(slot, 1, sizeof(slot), f);
        }
        for (size_t i = 0; i < n_positions_; ++i) {
            unsigned char buf[sizeof(uint64_t)];
            detail::store_u64_le(buf, positions_[i]);
            std::fwrite(buf, 1, sizeof(buf), f);
        }
        std::fclose(f);
    }

    // Loads a TRKI file with full validation of every file-controlled
    // field (vuln-0006): overflow-safe table-size checks, power-of-two
    // n_slots, k/n_keys sanity, exact file-size match, and a full walk of
    // the slots table verifying each (offset, count) lies within
    // n_positions. Any violation throws std::runtime_error -- loads never
    // publish a partially-validated index and never read out of bounds.
    static KmerReferenceIndex mmap_open(const std::string& path) {
        KmerReferenceIndex idx;
        idx.mapping_ = std::make_shared<detail::ReadOnlyMapping>(detail::ReadOnlyMapping::open(path));

        const auto* base = static_cast<const unsigned char*>(idx.mapping_->data);
        const size_t file_size = idx.mapping_->size;

        if (file_size < detail::kTrkiHeaderSize)
            throw std::runtime_error("KmerReferenceIndex: file too small to contain a header");

        // Field-wise explicit-LE reads (no struct memcpy -> endian-clean).
        if (std::memcmp(base, "TRKI", 4) != 0)
            throw std::runtime_error("KmerReferenceIndex: bad magic bytes");
        const uint32_t version = detail::load_u32_le(base + 4);
        if (version != 1)
            throw std::runtime_error("KmerReferenceIndex: unsupported version");
        const uint32_t k = detail::load_u32_le(base + 8);
        if (k == 0 || k > 32)
            throw std::runtime_error("KmerReferenceIndex: invalid k-mer length");
        const uint64_t n_slots = detail::load_u64_le(base + 16);
        const uint64_t n_keys = detail::load_u64_le(base + 24);
        const uint64_t n_positions = detail::load_u64_le(base + 32);

        // n_slots must be 0 (empty index) or a power of two -- the hash
        // table masks bucket indices with (n_slots - 1).
        if (n_slots != 0 && (n_slots & (n_slots - 1)) != 0)
            throw std::runtime_error("KmerReferenceIndex: n_slots must be a power of two");
        if (n_keys > n_slots)
            throw std::runtime_error("KmerReferenceIndex: n_keys exceeds n_slots");

        // Overflow-safe size validation (division form -- a product like
        // n_slots * sizeof(KmerSlot) can wrap mod 2^64; the division
        // comparison cannot). Every check below is against the ACTUAL
        // mapped length, so a crafted n_slots=2^60 is rejected here
        // instead of wrapping to 0 and passing a truncation check.
        const size_t bytes_after_header = file_size - detail::kTrkiHeaderSize;
        if (n_slots > bytes_after_header / sizeof(KmerSlot))
            throw std::runtime_error("KmerReferenceIndex: slots table exceeds file size");
        const size_t slots_bytes = (size_t)n_slots * sizeof(KmerSlot); // safe: <= file_size
        if (n_positions > (bytes_after_header - slots_bytes) / sizeof(uint64_t))
            throw std::runtime_error("KmerReferenceIndex: positions table exceeds file size");
        const size_t positions_bytes = (size_t)n_positions * sizeof(uint64_t); // safe: <= file_size

        // Exact match: the declared tables must account for the WHOLE file.
        // A mismatch means either truncation or trailing garbage -- both
        // rejected (vuln-0006 hardening).
        if (bytes_after_header != slots_bytes + positions_bytes)
            throw std::runtime_error("KmerReferenceIndex: file size does not match declared tables");

        const unsigned char* slot_bytes = base + detail::kTrkiHeaderSize;
        const unsigned char* pos_bytes = slot_bytes + slots_bytes;

        // Load-time validation walk: every non-empty slot's (offset, count)
        // must lie inside the positions table. Uses explicit LE reads so
        // the validation itself is endian-correct on any host. Malformed
        // entries fail loudly at LOAD time, not at query time.
        for (uint64_t i = 0; i < n_slots; ++i) {
            const uint64_t key = detail::load_u64_le(slot_bytes + i * sizeof(KmerSlot));
            if (key == KMER_SLOT_EMPTY) continue;
            const uint64_t value = detail::load_u64_le(slot_bytes + i * sizeof(KmerSlot) + sizeof(uint64_t));
            const uint32_t count = (uint32_t)(value & 0xFFFFFFFFu);
            const uint64_t offset = value >> 32;
            if (count == 0 || offset > n_positions || count > n_positions - offset)
                throw std::runtime_error("KmerReferenceIndex: slot references out-of-range positions");
        }

        idx.k_ = (int)k;
        idx.n_slots_ = (size_t)n_slots;
        idx.n_keys_ = (size_t)n_keys;
        idx.n_positions_ = (size_t)n_positions;

        if (detail::is_little_endian_host()) {
            // Zero-copy fast path: the mapped bytes ARE the query layout.
            idx.slots_ = reinterpret_cast<const KmerSlot*>(slot_bytes);
            idx.positions_ = reinterpret_cast<const uint64_t*>(pos_bytes);
        } else {
            // Big-endian host: decode the LE tables into owned buffers so
            // native struct access is correct.
            idx.owned_slots_.resize(n_slots);
            for (uint64_t i = 0; i < n_slots; ++i) {
                idx.owned_slots_[i].key = detail::load_u64_le(slot_bytes + i * sizeof(KmerSlot));
                idx.owned_slots_[i].value = detail::load_u64_le(slot_bytes + i * sizeof(KmerSlot) + sizeof(uint64_t));
            }
            idx.owned_positions_.resize(n_positions);
            for (uint64_t i = 0; i < n_positions; ++i)
                idx.owned_positions_[i] = detail::load_u64_le(pos_bytes + i * sizeof(uint64_t));
            idx.slots_ = idx.owned_slots_.data();
            idx.positions_ = idx.owned_positions_.data();
        }
        return idx;
    }

private:
    void finalize_empty() {
        n_slots_ = 0; n_keys_ = 0; n_positions_ = 0;
        slots_ = nullptr; positions_ = nullptr;
    }

    void build_from_pairs(std::vector<std::pair<uint64_t, uint64_t>>&& kmer_pos) {
        std::sort(kmer_pos.begin(), kmer_pos.end());

        owned_positions_.clear();
        owned_positions_.reserve(kmer_pos.size());

        std::vector<std::pair<uint64_t, uint64_t>> key_offset_count; // (key, offset<<32|count)
        size_t i = 0;
        while (i < kmer_pos.size()) {
            uint64_t key = kmer_pos[i].first;
            uint64_t offset = owned_positions_.size();
            size_t j = i;
            while (j < kmer_pos.size() && kmer_pos[j].first == key) {
                owned_positions_.push_back(kmer_pos[j].second);
                ++j;
            }
            uint32_t count = (uint32_t)(j - i);
            key_offset_count.emplace_back(key, (offset << 32) | count);
            i = j;
        }

        // The (offset << 32) | count encoding caps the positions table at
        // 2^32 entries; reject anything beyond instead of silently
        // corrupting offset/count fields.
        if (owned_positions_.size() > 0xFFFFFFFFull)
            throw std::runtime_error("KmerReferenceIndex: too many positions for TRKI format (max 2^32)");

        n_keys_ = key_offset_count.size();
        n_positions_ = owned_positions_.size();

        size_t n_slots = 1;
        while (n_slots < n_keys_ * 2 + 1) n_slots <<= 1; // ~50% max load factor
        n_slots_ = n_slots;

        owned_slots_.assign(n_slots_, KmerSlot{KMER_SLOT_EMPTY, 0});
        size_t mask = n_slots_ - 1;
        for (auto& [key, val] : key_offset_count) {
            size_t idx = (size_t)(kmer_slot_hash(key) & mask);
            while (owned_slots_[idx].key != KMER_SLOT_EMPTY) idx = (idx + 1) & mask;
            owned_slots_[idx] = KmerSlot{key, val};
        }

        slots_ = owned_slots_.data();
        positions_ = owned_positions_.data();
    }

    int k_ = 0;
    size_t n_slots_ = 0;
    size_t n_keys_ = 0;
    size_t n_positions_ = 0;

    const KmerSlot* slots_ = nullptr;
    const uint64_t* positions_ = nullptr;

    // Owning storage when built in-process (or decoded on a big-endian
    // host at load); empty when mmap_open()'d on a little-endian host
    // (in which case slots_/positions_ point into mapping_ instead).
    std::vector<KmerSlot> owned_slots_;
    std::vector<uint64_t> owned_positions_;
    std::shared_ptr<detail::ReadOnlyMapping> mapping_;
};

} // namespace TracEon

#endif // TRACEON_KMER_REFERENCE_INDEX_H
