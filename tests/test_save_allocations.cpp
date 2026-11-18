// Separate process: observe C++ buffer allocations without changing the library allocator.
#include "Cache.h"
#include <atomic>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <new>

static std::atomic<size_t> largest{0};
void* operator new(size_t n) {
    size_t old = largest.load();
    while (n > old && !largest.compare_exchange_weak(old, n)) {}
    if (void* p = std::malloc(n ? n : 1)) return p;
    throw std::bad_alloc();
}
void operator delete(void* p) noexcept { std::free(p); }
void operator delete(void* p, size_t) noexcept { std::free(p); }

int main() {
    const std::string path = "audit_save_allocations.traceon";
    struct Cleanup {
        const std::string& path;
        ~Cleanup() { std::error_code ec; std::filesystem::remove(path, ec); }
    } cleanup{path};
    try {
        TracEon::Cache cache;
        const std::string sequence(32 * 1024 * 1024, 'A');
        cache.set("chr1", sequence);
        for (auto storage : {TracEon::CacheStorage::Lz4, TracEon::CacheStorage::Mapped}) {
            largest = 0;
            cache.save(path, storage);
            const size_t maximum = largest.load();
            if (maximum > 2 * 1024 * 1024) {
                std::cerr << "unbounded save allocation: " << maximum << '\n';
                return 1;
            }
            TracEon::Cache restored;
            restored.restore(path);
            if (restored.getView("chr1") != sequence) return 1;
        }
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
