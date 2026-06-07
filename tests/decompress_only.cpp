#include <cstdint>
#include <fstream>
#include <iostream>
#include <thread>
#include <atomic>
#include <chrono>
#include <string>
#include <vector>
#include <memory>
#include <cstring>
#include <zlib.h>
#include <filesystem>
#include <unistd.h>
#include <sys/resource.h>

namespace fs = std::filesystem;

static long get_vmrss_kb() {
    std::ifstream f("/proc/self/status");
    if (!f) return -1;
    std::string line;
    while (std::getline(f, line))
        if (line.find("VmRSS:") == 0) {
            long v = 0;
            std::sscanf(line.c_str(), "%*s %ld", &v);
            return v;
        }
    return -1;
}

static long get_vmhwm_kb() {
    std::ifstream f("/proc/self/status");
    if (!f) return -1;
    std::string line;
    while (std::getline(f, line))
        if (line.find("VmHWM:") == 0) {
            long v = 0;
            std::sscanf(line.c_str(), "%*s %ld", &v);
            return v;
        }
    return -1;
}

static long get_vmpeak_kb() {
    std::ifstream f("/proc/self/status");
    if (!f) return -1;
    std::string line;
    while (std::getline(f, line))
        if (line.find("VmPeak:") == 0) {
            long v = 0;
            std::sscanf(line.c_str(), "%*s %ld", &v);
            return v;
        }
    return -1;
}

static long get_vmdata_kb() {
    std::ifstream f("/proc/self/status");
    if (!f) return -1;
    std::string line;
    while (std::getline(f, line))
        if (line.find("VmData:") == 0) {
            long v = 0;
            std::sscanf(line.c_str(), "%*s %ld", &v);
            return v;
        }
    return -1;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <gzip-file>\n";
        return 1;
    }

    const std::string filepath = argv[1];
    auto csize = fs::file_size(filepath);

    std::cout << "=== Decompression-Only Memory Profiler ===\n";
    std::cout << "File: " << filepath << "\n";
    std::cout << "Compressed: " << csize << " B (" << (csize / 1048576.0) << " MB)\n";

    const size_t CHUNK_SIZE = 1024 * 1024;
    std::vector<char> arena;

    long bl_vmrss  = get_vmrss_kb();
    long bl_vmhwm  = get_vmhwm_kb();
    long bl_vmpeak = get_vmpeak_kb();
    long bl_vmdata = get_vmdata_kb();

    std::cout << "\n--- Baseline ---\n";
    std::cout << "VmRSS="  << bl_vmrss  << " kB  VmHWM="  << bl_vmhwm
              << " kB  VmPeak=" << bl_vmpeak << " kB  VmData=" << bl_vmdata << " kB\n";

    // Pre-size: file_size * 3
    {
        const size_t estimated = static_cast<size_t>(csize * 3);
        arena.reserve(estimated);
        long ar = get_vmrss_kb();
        long ap = get_vmpeak_kb();
        std::cout << "\nAfter 3x pre-size reserve (" << (estimated / 1048576.0)
                  << " MB):\n  VmRSS=" << ar << " kB (delta +" << (ar - bl_vmrss)
                  << " kB)  VmPeak=" << ap << " kB\n";
    }

    // Background RSS poller (5ms interval)
    std::atomic<bool> done{false};
    std::atomic<long> peak_rss{0};
    std::thread poller([&]() {
        while (!done.load(std::memory_order_acquire)) {
            long cur = get_vmrss_kb();
            long p = peak_rss.load(std::memory_order_relaxed);
            while (cur > p) {
                if (peak_rss.compare_exchange_weak(p, cur,
                        std::memory_order_relaxed, std::memory_order_relaxed))
                    break;
            }
            usleep(5000);
        }
        long cur = get_vmrss_kb();
        long p = peak_rss.load(std::memory_order_relaxed);
        while (cur > p) {
            if (peak_rss.compare_exchange_weak(p, cur,
                    std::memory_order_relaxed, std::memory_order_relaxed))
                break;
        }
    });

    gzFile gz = gzopen(filepath.c_str(), "rb");
    if (!gz) {
        std::cerr << "Failed to open GZIP file\n";
        return 1;
    }

    auto chunk = std::make_unique<char[]>(CHUNK_SIZE);
    size_t write_pos = 0;

    auto t0 = std::chrono::steady_clock::now();

    std::cout << "\n--- Decompression loop ---\n";
    while (true) {
        size_t required = write_pos + CHUNK_SIZE;
        if (required > arena.capacity()) {
            size_t old_cap = arena.capacity();
            size_t new_cap = std::max(old_cap * 2, required);
            arena.reserve(new_cap);
            long cur_rss  = get_vmrss_kb();
            long cur_peak = get_vmpeak_kb();
            long cur_data = get_vmdata_kb();
            std::cout << "  reserve(" << (new_cap / 1048576.0)
                      << " MB) VmRSS=" << cur_rss << " kB  VmPeak=" << cur_peak
                      << " kB  VmData=" << cur_data << " kB"
                      << "  (delta: VmRSS +" << (cur_rss - bl_vmrss)
                      << " kB, VmPeak +" << (cur_peak - bl_vmpeak) << " kB)\n";
        }
        arena.resize(required);
        int bytes_read = gzread(gz, chunk.get(), CHUNK_SIZE);
        if (bytes_read < 0) {
            std::cerr << "GZIP read error\n";
            return 1;
        }
        if (bytes_read == 0) break;
        std::memcpy(arena.data() + write_pos, chunk.get(), static_cast<size_t>(bytes_read));
        write_pos += static_cast<size_t>(bytes_read);
    }
    gzclose(gz);

    // Trim
    arena.resize(write_pos);
    arena.shrink_to_fit();

    auto t1 = std::chrono::steady_clock::now();
    done.store(true, std::memory_order_release);
    poller.join();

    double elapsed = std::chrono::duration<double>(t1 - t0).count();
    long peak_polled  = peak_rss.load(std::memory_order_relaxed);
    long po_vmrss     = get_vmrss_kb();
    long po_vmhwm     = get_vmhwm_kb();
    long po_vmpeak    = get_vmpeak_kb();
    long po_vmdata    = get_vmdata_kb();

    std::cout << "\n=== DECOMPRESSION-ONLY RESULTS ===\n";
    std::cout << "Time: " << elapsed << " s\n";
    std::cout << "Final data size:     " << arena.size() << " B ("
              << (arena.size() / 1048576.0) << " MB)\n";
    std::cout << "Final capacity:      " << arena.capacity() << " B ("
              << (arena.capacity() / 1048576.0) << " MB)\n";
    std::cout << "Excess capacity:     " << (arena.capacity() - arena.size())
              << " B (" << ((arena.capacity() - arena.size()) / 1048576.0) << " MB)\n";
    std::cout << "Chunk buffer:        " << CHUNK_SIZE << " B (1 MB)\n";

    std::cout << "\nMemory readings:\n";
    std::cout << "  Baseline:  VmRSS="  << bl_vmrss  << " kB  VmHWM="  << bl_vmhwm
              << " kB  VmPeak=" << bl_vmpeak << " kB  VmData=" << bl_vmdata << " kB\n";
    std::cout << "  Post-decompress: VmRSS="  << po_vmrss  << " kB  VmHWM="  << po_vmhwm
              << " kB  VmPeak=" << po_vmpeak << " kB  VmData=" << po_vmdata << " kB\n";

    std::cout << "\nMemory deltas (peak - baseline):\n";
    std::cout << "  VmRSS  delta: "  << (po_vmrss  - bl_vmrss)  << " kB ("
              << ((po_vmrss - bl_vmrss) / 1024.0) << " MB)\n";
    std::cout << "  VmHWM  delta: "  << (po_vmhwm  - bl_vmhwm)  << " kB ("
              << ((po_vmhwm - bl_vmhwm) / 1024.0) << " MB)\n";
    std::cout << "  VmPeak delta: "  << (po_vmpeak - bl_vmpeak) << " kB ("
              << ((po_vmpeak - bl_vmpeak) / 1024.0) << " MB)\n";
    std::cout << "  VmData delta: "  << (po_vmdata - bl_vmdata) << " kB ("
              << ((po_vmdata - bl_vmdata) / 1024.0) << " MB)\n";
    std::cout << "  Polled RSS  peak delta: " << (peak_polled - bl_vmrss) << " kB ("
              << ((peak_polled - bl_vmrss) / 1024.0) << " MB)\n";

    std::cout << "\n--- VERDICT ---\n";
    std::cout << "Target: <= 5 MB peak RSS increase during decompression\n";
    long rss_delta = po_vmrss - bl_vmrss;
    long polled_delta = peak_polled - bl_vmrss;
    long worst = (polled_delta > rss_delta) ? polled_delta : rss_delta;

    std::cout << "  Final RSS delta:  " << (rss_delta / 1024.0) << " MB\n";
    std::cout << "  Peak (polled) RSS delta: " << (polled_delta / 1024.0) << " MB\n";
    std::cout << "  Worst observed:  " << (worst / 1024.0) << " MB\n";
    if (worst <= 5 * 1024) {
        std::cout << "  PASS\n";
    } else {
        std::cout << "  FAIL (exceeds 5 MB by " << ((worst - 5*1024) / 1024.0) << " MB)\n";
    }

    std::cout << "\n--- Note ---\n";
    std::cout << "The 5 MB target appears infeasible for " << (arena.size() / 1048576.0)
              << " MB uncompressed data.\n";
    std::cout << "Pre-size at 3x (" << (csize*3/1048576.0) << " MB) guarantees at least "
              << (csize*3/1048576.0) << " MB VM usage.\n";
    std::cout << "Geometric growth peak capacity reaches ~240 MB (one doubling).\n";
    std::cout << "Actual peak RSS reflects the uncompressed data written to RAM.\n";

    return (worst <= 5 * 1024) ? 0 : 1;
}
