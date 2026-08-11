#pragma once
/**
 * @file Crc32c.h
 * @brief CRC-32C (Castagnoli) with runtime hardware dispatch.
 *
 * Implements the CRC-32C variant used by SCTP/iSCSI/ext4 (reflected in/out,
 * polynomial 0x1EDC6F41, init 0xFFFFFFFF, final XOR 0xFFFFFFFF) as the
 * `.traceon` v4 whole-payload integrity checksum. Three implementations:
 *
 *   - crc32c_x86_run   — x86 SSE4.2 CRC32 instruction (`_mm_crc32_u*`),
 *                        compiled only when TRACEON_HAS_AVX2 is defined
 *                        (CMake sets it on x86-64 when the compiler accepts
 *                        -mavx2; the function carries a target("sse4.2")
 *                        attribute so it compiles without global -msse4.2).
 *                        Selected at runtime via __builtin_cpu_supports("sse4.2").
 *   - crc32c_arm_run   — AArch64 crc32cx/crc32cw/crc32cb/crc32ch intrinsics,
 *                        compiled only when the compiler advertises
 *                        __ARM_FEATURE_CRC32 (AArch64 builds always define
 *                        TRACEON_HAS_NEON via CMake; the CRC32 instruction is
 *                        mandatory on ARMv8.1+ and universally present on the
 *                        AArch64 hosts TracEon targets).
 *   - crc32c_table_run — portable byte-wise table fallback (always available).
 *
 * All three share the same *raw running* semantics so they can be mixed
 * across chunk boundaries: the caller initialises the accumulator to
 * 0xFFFFFFFF, feeds chunks, and applies the final XOR 0xFFFFFFFF once at the
 * end (see Crc32c below). The dispatcher crc32c_dispatch_run() picks the
 * fastest compiled path for the current CPU.
 *
 * Equivalence between the fallback and hardware paths is tested in
 * tests/SmartStrategyTests.cpp (random inputs, many lengths).
 */

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <array>

#if defined(__x86_64__) || defined(_M_X64)
#  include <immintrin.h>
#endif

#if defined(__aarch64__) && defined(__ARM_FEATURE_CRC32)
#  include <arm_acle.h>
#endif

namespace TracEon {

namespace detail {

// Reflected CRC-32C table. Entry i is the CRC of byte i with the reflected
// polynomial 0x82F63B78 (bit-reversal of 0x1EDC6F41). The table is generated
// at compile time (constexpr) — no static-init cost, no mutable global.
inline constexpr uint32_t crc32c_make_table_entry(int i) noexcept {
    uint32_t c = static_cast<uint32_t>(i);
    for (int k = 0; k < 8; ++k)
        c = (c & 1u) ? (0x82F63B78u ^ (c >> 1u)) : (c >> 1u);
    return c;
}

inline constexpr std::array<uint32_t, 256> kCrc32cTable = [] {
    std::array<uint32_t, 256> t{};
    for (int i = 0; i < 256; ++i) t[static_cast<size_t>(i)] = crc32c_make_table_entry(i);
    return t;
}();

// Portable byte-wise fallback. `crc` is the raw running accumulator
// (init 0xFFFFFFFF); no final XOR is applied here.
[[nodiscard]] inline uint32_t crc32c_table_run(const uint8_t* data, size_t len,
                                               uint32_t crc) noexcept {
    for (size_t i = 0; i < len; ++i)
        crc = (crc >> 8u) ^ kCrc32cTable[(crc ^ data[i]) & 0xFFu];
    return crc;
}

// x86 SSE4.2 CRC32 instruction. Compiled with target("sse4.2") so the
// intrinsic is available even though the TU as a whole is not built with
// -msse4.2. `crc` is the raw running accumulator; the instruction processes
// 8/4/2/1 bytes per step and returns the updated raw accumulator.
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__)) && defined(TRACEON_HAS_AVX2)
[[nodiscard]] __attribute__((target("sse4.2"))) inline uint32_t
crc32c_x86_run(const uint8_t* data, size_t len, uint32_t crc) noexcept {
    while (len >= 8) {
        uint64_t v;
        std::memcpy(&v, data, sizeof(v)); // unaligned-safe
        crc = static_cast<uint32_t>(_mm_crc32_u64(crc, v));
        data += 8;
        len -= 8;
    }
    if (len >= 4) {
        uint32_t v;
        std::memcpy(&v, data, sizeof(v));
        crc = _mm_crc32_u32(crc, v);
        data += 4;
        len -= 4;
    }
    if (len >= 2) {
        uint16_t v;
        std::memcpy(&v, data, sizeof(v));
        crc = _mm_crc32_u16(crc, v);
        data += 2;
        len -= 2;
    }
    if (len >= 1) crc = _mm_crc32_u8(crc, *data);
    return crc;
}
#else
[[nodiscard]] inline uint32_t crc32c_x86_run(const uint8_t* data, size_t len,
                                             uint32_t crc) noexcept {
    return crc32c_table_run(data, len, crc);
}
#endif

// AArch64 crc32cx/crc32cw/crc32ch/crc32cb instructions. Same raw-running
// semantics as the x86 path and the table fallback.
// Apple Clang does not declare the ACLE __crc32c* intrinsics in its
// (minimal) <arm_acle.h>, so on Apple platforms we use the portable table
// fallback (correct and covered by the fallback-equivalence tests). Linux
// GCC/Clang aarch64 builds get the hardware path.
#if defined(__aarch64__) && defined(__ARM_FEATURE_CRC32) && !defined(__APPLE__)
#  include <arm_acle.h>
[[nodiscard]] inline uint32_t crc32c_arm_run(const uint8_t* data, size_t len,
                                             uint32_t crc) noexcept {
    while (len >= 8) {
        uint64_t v;
        std::memcpy(&v, data, sizeof(v));
        crc = __crc32cx(crc, v);
        data += 8;
        len -= 8;
    }
    if (len >= 4) {
        uint32_t v;
        std::memcpy(&v, data, sizeof(v));
        crc = __crc32cw(crc, v);
        data += 4;
        len -= 4;
    }
    if (len >= 2) {
        uint16_t v;
        std::memcpy(&v, data, sizeof(v));
        crc = __crc32ch(crc, v);
        data += 2;
        len -= 2;
    }
    if (len >= 1) crc = __crc32cb(crc, *data);
    return crc;
}
#else
[[nodiscard]] inline uint32_t crc32c_arm_run(const uint8_t* data, size_t len,
                                             uint32_t crc) noexcept {
    return crc32c_table_run(data, len, crc);
}
#endif

// Runtime dispatcher. On x86-64 (GCC/Clang + TRACEON_HAS_AVX2) the SSE4.2
// path is used when the CPU advertises sse4.2; on AArch64 with CRC32 the
// intrinsic path is used; otherwise the table fallback.
[[nodiscard]] inline uint32_t crc32c_dispatch_run(const uint8_t* data, size_t len,
                                                  uint32_t crc) noexcept {
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__)) && defined(TRACEON_HAS_AVX2)
    if (__builtin_cpu_supports("sse4.2")) return crc32c_x86_run(data, len, crc);
#elif defined(__aarch64__) && defined(__ARM_FEATURE_CRC32)
    return crc32c_arm_run(data, len, crc);
#endif
    (void)data; (void)len; (void)crc;
    return crc32c_table_run(data, len, crc);
}

[[nodiscard]] inline const char* crc32c_impl_name() noexcept {
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__)) && defined(TRACEON_HAS_AVX2)
    if (__builtin_cpu_supports("sse4.2")) return "sse4.2";
#elif defined(__aarch64__) && defined(__ARM_FEATURE_CRC32)
    return "aarch64-crc32";
#endif
    return "table";
}

} // namespace detail

/**
 * @brief Streaming CRC-32C (Castagnoli) accumulator.
 *
 * Feed arbitrary chunks with update(); call finalize() once for the standard
 * CRC-32C value (init 0xFFFFFFFF, final XOR 0xFFFFFFFF). One-shot callers
 * can use the free function crc32c().
 *
 * This is the checksum used by the `.traceon` v4 binary format: it covers
 * the canonical header fields (everything up to but excluding the checksum
 * field itself) plus the ENTIRE uncompressed logical payload. On the
 * save/load paths the accumulator is updated as chunks pass through the
 * LZ4F compressor/decompressor — there is deliberately no second full-payload
 * pass and no extra allocation.
 */
class Crc32c {
public:
    Crc32c() noexcept = default;

    void update(const void* data, size_t len) noexcept {
        crc_ = detail::crc32c_dispatch_run(static_cast<const uint8_t*>(data), len, crc_);
    }

    [[nodiscard]] uint32_t finalize() const noexcept { return crc_ ^ 0xFFFFFFFFu; }

    void reset() noexcept { crc_ = 0xFFFFFFFFu; }

private:
    uint32_t crc_ = 0xFFFFFFFFu;
};

/** @brief One-shot CRC-32C over [data, data+len). */
[[nodiscard]] inline uint32_t crc32c(const void* data, size_t len) noexcept {
    return detail::crc32c_dispatch_run(static_cast<const uint8_t*>(data), len,
                                       0xFFFFFFFFu) ^ 0xFFFFFFFFu;
}

/** @brief Force the portable table implementation (tests/benchmarks). */
[[nodiscard]] inline uint32_t crc32c_table_only(const void* data, size_t len) noexcept {
    return detail::crc32c_table_run(static_cast<const uint8_t*>(data), len,
                                    0xFFFFFFFFu) ^ 0xFFFFFFFFu;
}

/** @brief Force the hardware implementation when compiled; table otherwise. */
[[nodiscard]] inline uint32_t crc32c_hw_only(const void* data, size_t len) noexcept {
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__)) && defined(TRACEON_HAS_AVX2)
    return detail::crc32c_x86_run(static_cast<const uint8_t*>(data), len,
                                  0xFFFFFFFFu) ^ 0xFFFFFFFFu;
#elif defined(__aarch64__) && defined(__ARM_FEATURE_CRC32)
    return detail::crc32c_arm_run(static_cast<const uint8_t*>(data), len,
                                  0xFFFFFFFFu) ^ 0xFFFFFFFFu;
#else
    return detail::crc32c_table_run(static_cast<const uint8_t*>(data), len,
                                    0xFFFFFFFFu) ^ 0xFFFFFFFFu;
#endif
}

/** @brief Name of the implementation the dispatcher selected on this CPU. */
[[nodiscard]] inline const char* crc32c_active_impl() noexcept {
    return detail::crc32c_impl_name();
}

} // namespace TracEon
