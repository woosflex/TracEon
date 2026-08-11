#ifndef TRACEON_MAPDEFS_H
#define TRACEON_MAPDEFS_H

#include <ankerl/unordered_dense.h>
#include <string>
#include <string_view>
#include <functional>
#include <cstring>

namespace TracEon {

// OPTIMIZATION: Use ankerl::unordered_dense for better cache locality and insert performance
template<typename Key, typename Value>
using HashMap = ankerl::unordered_dense::map<Key, Value>;

// Transparent hasher: enables heterogeneous lookup (find by string_view,
// const char*, or std::string — no std::string allocation on the hot read
// path). All three overloads hash the CONTENT, so a key stored as
// std::string_view and a lookup passed as std::string or const char* hash
// identically. The explicit const char* overload is required: with only the
// string_view and std::string overloads, a const char* argument is ambiguous
// (both are user-defined conversions); with it, string literals and C strings
// dispatch directly without constructing a temporary std::string.
// OPTIMIZATION (OPT-1): fast non-cryptographic string hashing for the map
// insert hotspot. The previous implementation used libstdc++
// std::hash<std::string_view> (_Hash_bytes, a Murmur-style hash) — slow for
// the short genomic identifiers this library stores ('read12345' etc.) — and
// because it did NOT declare is_avalanching, ankerl::unordered_dense wrapped
// its output in ANOTHER wyhash re-mix on every probe (double hashing).
//
// This is a self-contained wyhash-style hash (public domain, Wang Yi — the
// same MUM-mix construction ankerl itself vendors, but header-local with no
// dependency on ankerl internals): a bijective multiply-xor mix over 128-bit
// multiplies gives full 64-bit avalanche, so the same value is safe to feed
// ankerl directly. Declaring is_avalanching makes ankerl skip its extra
// re-mix entirely. The len<=16 fast path reads at most four dwords with no
// loop, which covers essentially every genomic ID.
struct TransparentStringHash {
    using is_transparent = void;
    // Tells ankerl::unordered_dense the hash output is already fully mixed:
    // mixed_hash() then uses it verbatim instead of applying wyhash::hash()
    // on top (one 128-bit multiply + xor saved per probe).
    using is_avalanching = void;

    // wyhash secret constants (public domain).
    static constexpr uint64_t kSeed    = 0xa0761d6478bd642fULL;
    static constexpr uint64_t kSecret1 = 0xe7037ed1a0b428dbULL;
    static constexpr uint64_t kSecret2 = 0x8ebc6af09c88c6e3ULL;
    static constexpr uint64_t kSecret3 = 0x589965cc75374cc3ULL;

    // 64x64 -> 128 multiply, xor the halves. Bijective, fully avalanching.
    static uint64_t mum(uint64_t a, uint64_t b) noexcept {
#if defined(__SIZEOF_INT128__)
        __uint128_t r = static_cast<__uint128_t>(a) * b;
        return static_cast<uint64_t>(r) ^ static_cast<uint64_t>(r >> 64U);
#else
        // Schoolbook 32-bit fallback (portable; identical result).
        uint64_t ha = a >> 32U, hb = b >> 32U;
        uint64_t la = static_cast<uint32_t>(a), lb = static_cast<uint32_t>(b);
        uint64_t rh = ha * hb, rm0 = ha * lb, rm1 = hb * la, rl = la * lb;
        uint64_t t = rl + (rm0 << 32U);
        auto c = static_cast<uint64_t>(t < rl);
        uint64_t lo = t + (rm1 << 32U);
        c += static_cast<uint64_t>(lo < t);
        uint64_t hi = rh + (rm0 >> 32U) + (rm1 >> 32U) + c;
        return lo ^ hi;
#endif
    }

    // multiply-and-xor mix (MUM).
    static uint64_t mix(uint64_t a, uint64_t b) noexcept { return mum(a, b); }

    static uint64_t load8(const char* p) noexcept {
        uint64_t v;
        std::memcpy(&v, p, sizeof(v)); // unaligned-safe; one load on x86/ARM
        return v;
    }
    static uint64_t load4(const char* p) noexcept {
        uint32_t v;
        std::memcpy(&v, p, sizeof(v));
        return v;
    }
    // reads the 1..3 remaining bytes, wyhash r3 style
    static uint64_t load3(const char* p, size_t k) noexcept {
        return (static_cast<uint64_t>(static_cast<unsigned char>(p[0])) << 16U) |
               (static_cast<uint64_t>(static_cast<unsigned char>(p[k >> 1U])) << 8U) |
               static_cast<uint64_t>(static_cast<unsigned char>(p[k - 1]));
    }

    // Content hash — identical for std::string_view, const char*, and
    // std::string arguments of the same bytes (all three overloads forward
    // here), so transparent lookup never depends on the argument type.
    size_t hash_bytes(const char* data, size_t len) const noexcept {
        const char* p = data;
        uint64_t seed = kSeed;
        uint64_t a = 0, b = 0;
        if (len <= 16) {
            // Fast path for short IDs: at most four dword loads, no loop.
            if (len >= 4) {
                a = (load4(p) << 32U) | load4(p + ((len >> 3U) << 2U));
                b = (load4(p + len - 4) << 32U) | load4(p + len - 4 - ((len >> 3U) << 2U));
            } else if (len > 0) {
                a = load3(p, len);
            }
        } else {
            size_t i = len;
            if (i > 48) {
                uint64_t see1 = seed, see2 = seed;
                do {
                    seed = mix(load8(p) ^ kSecret2, load8(p + 8) ^ seed);
                    see1 = mix(load8(p + 16) ^ kSecret3, load8(p + 24) ^ see1);
                    see2 = mix(load8(p + 32) ^ kSecret1, load8(p + 40) ^ see2);
                    p += 48;
                    i -= 48;
                } while (i > 48);
                seed ^= see1 ^ see2;
            }
            while (i > 16) {
                seed = mix(load8(p) ^ kSecret2, load8(p + 8) ^ seed);
                i -= 16;
                p += 16;
            }
            a = load8(p + i - 16);
            b = load8(p + i - 8);
        }
        return mix(kSecret1 ^ len, mix(a ^ kSecret1, b ^ seed));
    }

    size_t operator()(std::string_view sv) const noexcept {
        return hash_bytes(sv.data(), sv.size());
    }
    size_t operator()(const char* s) const noexcept {
        return hash_bytes(s, std::char_traits<char>::length(s));
    }
    size_t operator()(const std::string& s) const noexcept {
        return hash_bytes(s.data(), s.size());
    }
};

// StringHashMap: std::string keys, but find/count/contains accept string_view.
template<typename Value>
using StringHashMap = ankerl::unordered_dense::map<
    std::string, Value, TransparentStringHash, std::equal_to<>>;

// StringViewHashMap: std::string_view keys, transparent lookup. The KEY IS
// NON-OWNING — the caller MUST guarantee the backing buffer outlives the map
// and is never moved/reallocated while the map is alive. TracEon satisfies
// this for every insert path:
//   - parse paths: views point into text_arena_, which is immutable after
//     load (data_loaded_) and only rebuilt by clearInternal(), which clears
//     the map FIRST;
//   - loadBinary v1: views point into the mmap region (mmap_handle_), kept
//     alive until clearInternal() (map cleared first);
//   - addEntry(): views point into manual_store_ (std::deque<std::string> —
//     push_back never invalidates references to existing elements).
// Zero allocations on the insert path: the map stores ptr+len, never copies
// the identifier bytes.
template<typename Value>
using StringViewHashMap = ankerl::unordered_dense::map<
    std::string_view, Value, TransparentStringHash, std::equal_to<>>;

} // namespace TracEon

#endif // TRACEON_MAPDEFS_H