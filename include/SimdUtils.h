#pragma once
/**
 * @file SimdUtils.h
 * @brief Cross-platform SIMD character-search with runtime dispatch (AVX2/NEON/scalar)
 *
 * This header provides three concrete implementations and a small runtime
 * dispatcher:
 *   - simd_find_char_avx2  (compiled with AVX2 target attribute when available)
 *   - simd_find_char_neon  (NEON implementation for AArch64)
 *   - simd_find_char_scalar (std::memchr fallback)
 *
 * The dispatcher simd_find_char() uses __builtin_cpu_supports("avx2") on
 * GCC/Clang x86-64 hosts to select the AVX2 implementation at runtime when
 * available.  CMake sets TRACEON_HAS_AVX2/TRACEON_HAS_NEON when the platform
 * and compiler allow compiling the specialized variants.
 */

#include <cstring>
#include <cstdint>
#include <cstddef>

// Intrinsics (conditionally included)
#if defined(__x86_64__) || defined(_M_X64)
#  include <immintrin.h>
#endif

#if defined(__aarch64__) || defined(__ARM_NEON) || defined(TRACEON_HAS_NEON)
#  include <arm_neon.h>
#endif

#if defined(_MSC_VER)
#  include <intrin.h>
#endif

// Portable count-trailing-zeros (undefined if x == 0)
namespace TracEon::detail {

#if defined(_MSC_VER)
#  pragma intrinsic(_BitScanForward)
    inline int ctz32(unsigned int x) noexcept {
        unsigned long i; _BitScanForward(&i, x); return static_cast<int>(i);
    }
#  if defined(_WIN64)
#    pragma intrinsic(_BitScanForward64)
    inline int ctz64(unsigned long long x) noexcept {
        unsigned long i; _BitScanForward64(&i, x); return static_cast<int>(i);
    }
#  else
    inline int ctz64(unsigned long long x) noexcept {
        unsigned lo = static_cast<unsigned>(x);
        if (lo) return ctz32(lo);
        return 32 + ctz32(static_cast<unsigned>(x >> 32));
    }
#  endif
#else
    inline int ctz32(unsigned int x) noexcept { return __builtin_ctz(x); }
    inline int ctz64(unsigned long long x) noexcept { return __builtin_ctzll(x); }
#endif

} // namespace TracEon::detail

namespace TracEon {

// Scalar fallback using libc's memchr (typically auto-vectorised)
[[nodiscard]] inline const char* simd_find_char_scalar(const char* ptr, const char* end, char c) noexcept {
    const void* found = std::memchr(ptr, static_cast<unsigned char>(c), static_cast<size_t>(end - ptr));
    return found ? static_cast<const char*>(found) : end;
}

// NEON implementation (AArch64)
#if defined(__aarch64__) || defined(__ARM_NEON) || defined(TRACEON_HAS_NEON)
[[nodiscard]] inline const char* simd_find_char_neon(const char* ptr, const char* end, char c) noexcept {
    const uint8x16_t needle = vdupq_n_u8(static_cast<uint8_t>(c));
    while (ptr + 16 <= end) {
        uint8x16_t chunk = vld1q_u8(reinterpret_cast<const uint8_t*>(ptr));
        uint8x16_t eq = vceqq_u8(chunk, needle);
        uint64x2_t eq64 = vreinterpretq_u64_u8(eq);
        uint64_t lo = vgetq_lane_u64(eq64, 0);
        uint64_t hi = vgetq_lane_u64(eq64, 1);
        if (lo | hi) {
            if (lo) return ptr + (detail::ctz64(lo) >> 3);
            else    return ptr + 8 + (detail::ctz64(hi) >> 3);
        }
        ptr += 16;
    }
    while (ptr < end) { if (*ptr == c) return ptr; ++ptr; }
    return end;
}
#endif

// AVX2 implementation compiled with target("avx2") on GCC/Clang when available
#if defined(TRACEON_HAS_AVX2) && (defined(__GNUC__) || defined(__clang__))
[[nodiscard]] __attribute__((target("avx2"))) inline const char* simd_find_char_avx2(const char* ptr, const char* end, char c) noexcept {
    const __m256i needle = _mm256_set1_epi8(c);
    while (ptr + 32 <= end) {
        __m256i chunk = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(ptr));
        __m256i eq = _mm256_cmpeq_epi8(chunk, needle);
        unsigned mask = static_cast<unsigned>(_mm256_movemask_epi8(eq));
        if (mask) return ptr + detail::ctz32(mask);
        ptr += 32;
    }
    while (ptr < end) { if (*ptr == c) return ptr; ++ptr; }
    return end;
}
#else
// If AVX2 variant was not compiled, provide a cheap alias to the scalar path
[[nodiscard]] inline const char* simd_find_char_avx2(const char* ptr, const char* end, char c) noexcept {
    (void)ptr; (void)end; (void)c; return simd_find_char_scalar(ptr, end, c);
}
#endif

// Optional MSVC runtime AVX2 detection (used only if TRACEON_HAS_AVX2 is set for MSVC)
#if defined(_MSC_VER) && defined(TRACEON_HAS_AVX2)
static inline bool cpu_has_avx2_msvc() noexcept {
    // Check CPUID for AVX/OSXSAVE and extended features for AVX2
    int cpuInfo[4] = {0,0,0,0};
    __cpuid(cpuInfo, 0);
    int nIds = cpuInfo[0];
    if (nIds < 1) return false;
    __cpuid(cpuInfo, 1);
    bool osxsave = (cpuInfo[2] & (1 << 27)) != 0;
    bool avx     = (cpuInfo[2] & (1 << 28)) != 0;
    if (!osxsave || !avx) return false;
    unsigned long long xcr = _xgetbv(0);
    if ((xcr & 0x6ULL) != 0x6ULL) return false;
    if (nIds >= 7) {
        __cpuidex(cpuInfo, 7, 0);
        return (cpuInfo[1] & (1 << 5)) != 0; // AVX2 bit in EBX
    }
    return false;
}
#endif

/**
 * @brief Runtime-dispatching simd_find_char.
 *
 * On AArch64 this calls the NEON path.  On x86-64 (GCC/Clang) this checks
 * __builtin_cpu_supports("avx2") and calls the AVX2 variant when available;
 * otherwise the scalar fallback is used.  The AVX2 function is compiled only
 * when TRACEON_HAS_AVX2 was defined by CMake.
 */
[[nodiscard]] inline const char* simd_find_char(const char* ptr, const char* end, char c) noexcept {
#if defined(__aarch64__) || defined(TRACEON_HAS_NEON)
    return simd_find_char_neon(ptr, end, c);
#elif defined(__x86_64__) || defined(_M_X64)
    // Prefer compiler builtin on GCC/Clang which also checks OS support
#  if defined(TRACEON_HAS_AVX2) && (defined(__GNUC__) || defined(__clang__))
    if (__builtin_cpu_supports("avx2")) return simd_find_char_avx2(ptr, end, c);
    return simd_find_char_scalar(ptr, end, c);
#  elif defined(_MSC_VER) && defined(TRACEON_HAS_AVX2)
    if (cpu_has_avx2_msvc()) return simd_find_char_avx2(ptr, end, c);
    return simd_find_char_scalar(ptr, end, c);
#  else
    return simd_find_char_scalar(ptr, end, c);
#  endif
#else
    return simd_find_char_scalar(ptr, end, c);
#endif
}

} // namespace TracEon
