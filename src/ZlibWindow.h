#ifndef TRACEON_ZLIB_WINDOW_H
#define TRACEON_ZLIB_WINDOW_H

#include <atomic>
#include <cstddef>
#include <limits>

namespace TracEon::detail {

// zlib's avail_in/avail_out use uInt even when the surrounding file and
// buffers are addressed with size_t. Feed at most one uInt-sized window.
constexpr unsigned int zlibWindowSize(size_t remaining) noexcept {
    constexpr size_t max_window = std::numeric_limits<unsigned int>::max();
    return static_cast<unsigned int>(remaining > max_window ? max_window : remaining);
}

// Atomic admission: competing workers cannot spend the same remaining bytes.
inline bool reserveGzipBytes(std::atomic<size_t>& reserved, size_t bytes, size_t cap) {
    size_t current = reserved.load(std::memory_order_relaxed);
    do {
        if (bytes > cap || current > cap - bytes) return false;
    } while (!reserved.compare_exchange_weak(current, current + bytes,
                                             std::memory_order_relaxed));
    return true;
}

// Checked geometric growth used by parallel GZIP decoding. Returns false if
// doubling would overflow, fail to make progress, or exceed the caller's cap.
constexpr bool nextGzipCapacity(size_t current, size_t cap, size_t& next) noexcept {
    if (current == 0) {
        next = cap == 0 ? 0 : 1;
        return next != 0;
    }
    if (current > std::numeric_limits<size_t>::max() / 2) return false;
    next = current * 2;
    return next > current && next <= cap;
}

} // namespace TracEon::detail

#endif
