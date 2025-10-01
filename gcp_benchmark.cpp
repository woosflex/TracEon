#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <vector>
#include <filesystem>
#include <algorithm>
#include <sstream>
#include <ctime>
#include <thread>
#include <zlib.h>
#include "Cache.h"

#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#include <psapi.h>
#else
#include <unistd.h>
#include <sys/resource.h>
#endif

namespace fs = std::filesystem;

struct SystemInfo {
    size_t total_memory_mb, available_memory_mb;
    std::string cpu_info, hostname;
    unsigned int num_cores;

    SystemInfo() {
        num_cores = std::thread::hardware_concurrency();
#ifdef _WIN32
        MEMORYSTATUSEX memInfo;
        memInfo.dwLength = sizeof(MEMORYSTATUSEX);
        GlobalMemoryStatusEx(&memInfo);
        total_memory_mb = memInfo.ullTotalPhys / (1024 * 1024);
        available_memory_mb = memInfo.ullAvailPhys / (1024 * 1024);
        cpu_info = "Windows CPU";
        char hostname_buf[256];
        DWORD size = sizeof(hostname_buf);
        GetComputerNameA(hostname_buf, &size);
        hostname = hostname_buf;
#else
        std::ifstream meminfo("/proc/meminfo");
        std::string line;
        total_memory_mb = available_memory_mb = 0;
        while (std::getline(meminfo, line)) {
            if (line.rfind("MemTotal:", 0) == 0) total_memory_mb = std::stoul(line.substr(10)) / 1024;
            if (line.rfind("MemAvailable:", 0) == 0) available_memory_mb = std::stoul(line.substr(14)) / 1024;
        }
        std::ifstream cpuinfo("/proc/cpuinfo");
        while (std::getline(cpuinfo, line)) {
            if (line.rfind("model name", 0) == 0) {
                cpu_info = line.substr(line.find(':') + 2);
                break;
            }
        }
        char hostname_buf[256];
        gethostname(hostname_buf, sizeof(hostname_buf));
        hostname = hostname_buf;
#endif
    }
};

struct MemoryUsage {
    size_t rss_kb;
    MemoryUsage() {
#ifdef _WIN32
        PROCESS_MEMORY_COUNTERS pmc;
        GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc));
        rss_kb = pmc.WorkingSetSize / 1024;
#else
        std::ifstream status("/proc/self/status");
        std::string line;
        rss_kb = 0;
        while (std::getline(status, line)) {
            if (line.rfind("VmRSS:", 0) == 0) {
                rss_kb = std::stoul(line.substr(7));
                break;
            }
        }
#endif
    }
    double rss_mb() const { return rss_kb / 1024.0; }
};

class Timer {
    std::chrono::high_resolution_clock::time_point start;
public:
    Timer() : start(std::chrono::high_resolution_clock::now()) {}
    double elapsed_ms() const {
        auto end = std::chrono::high_resolution_clock::now();
        return std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0;
    }
};

void printSection(const std::string& title) {
    std::cout << "\n" << std::string(80, '=') << "\n" << title << "\n" << std::string(80, '=') << "\n";
}

void printSubSection(const std::string& title) {
    std::cout << "\n" << std::string(80, '-') << "\n" << title << "\n" << std::string(80, '-') << "\n";
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_file>\n";
        return 1;
    }

    std::string filepath = argv[1];
    std::cout << std::fixed << std::setprecision(4);

    try {
        std::cout << "\n╔" << std::string(78, '=') << "╗\n"
                  << "║" << std::string(20, ' ') << "TracEon NGS Benchmark - Production Grade" << std::string(17, ' ') << "║\n"
                  << "╚" << std::string(78, '=') << "╝\n";

        SystemInfo sysinfo;
        auto t = std::time(nullptr);
        std::cout << "\n[SYSTEM INFORMATION]\n"
                  << "  Timestamp: " << std::put_time(std::localtime(&t), "%Y-%m-%d %H:%M:%S") << "\n"
                  << "  Hostname: " << sysinfo.hostname << "\n"
                  << "  CPU: " << sysinfo.cpu_info << "\n"
                  << "  CPU Cores: " << sysinfo.num_cores << "\n"
                  << "  Total RAM: " << sysinfo.total_memory_mb << " MB\n"
                  << "  Available RAM: " << sysinfo.available_memory_mb << " MB\n"
                  << "  Input File: " << filepath << "\n";

        bool is_gzipped = (filepath.size() > 3 && filepath.substr(filepath.size() - 3) == ".gz");
        bool is_fastq = (filepath.find(".fastq") != std::string::npos || filepath.find(".fq") != std::string::npos);
        double file_size_mb = fs::file_size(filepath) / (1024.0 * 1024.0);

        std::cout << "  File Type: " << (is_fastq ? "FASTQ" : "FASTA") << (is_gzipped ? " (gzipped)" : " (plain text)") << "\n"
                  << "  File Size: " << file_size_mb << " MB (" << (file_size_mb / 1024.0) << " GB)\n";

        printSection("PHASE 1: LOADING DATA INTO TRACEON CACHE");
        MemoryUsage mem_before;
        std::cout << "  Initial Memory Usage: " << mem_before.rss_mb() << " MB (RSS)\n"
                  << "  Loading file (multithreaded if > 1MB, single-threaded otherwise)...\n\n";

        Timer load_timer;
        TracEon::Cache cache;
        cache.loadFile(filepath);
        double load_time = load_timer.elapsed_ms();
        MemoryUsage mem_after;

        size_t num_sequences = cache.size();
        double mem_increase = mem_after.rss_mb() - mem_before.rss_mb();
        double load_throughput = file_size_mb / (load_time / 1000.0);

        std::cout << "[LOAD RESULTS]\n"
                  << "  ✓ Sequences Loaded: " << num_sequences << "\n"
                  << "  ✓ Load Time: " << load_time << " ms (" << (load_time / 1000.0) << " seconds)\n"
                  << "  ✓ Load Throughput: " << load_throughput << " MB/s\n"
                  << "  ✓ Memory After Load: " << mem_after.rss_mb() << " MB (RSS)\n"
                  << "  ✓ Memory Increase: " << mem_increase << " MB\n"
                  << "  ✓ Memory Efficiency: " << (file_size_mb / mem_increase) << "× (original/memory ratio)\n"
                  << "  ✓ Bytes per Sequence: " << (mem_increase * 1024 * 1024 / num_sequences) << " bytes\n";

        printSection("PHASE 2: BINARY CACHE PERSISTENCE");
        std::string cache_file = "traceon_cache.bin";

        std::cout << "  Saving cache to disk (format: SMRT binary)...\n";
        Timer save_timer;
        cache.save(cache_file);
        double save_time = save_timer.elapsed_ms();

        if (!fs::exists(cache_file)) {
            std::cerr << "\n  ✗ ERROR: Cache file was not created!\n"
                      << "  This usually means the save() method needs to call SmartStrategy::saveBinary()\n";
            return 1;
        }

        double cache_size_mb = fs::file_size(cache_file) / (1024.0 * 1024.0);
        double compression_ratio = ((file_size_mb - cache_size_mb) / file_size_mb) * 100.0;
        double save_throughput = cache_size_mb / (save_time / 1000.0);

        std::cout << "\n[SAVE RESULTS]\n"
                  << "  ✓ Cache File: " << cache_file << "\n"
                  << "  ✓ Cache Size: " << cache_size_mb << " MB (" << (cache_size_mb / 1024.0) << " GB)\n"
                  << "  ✓ Compression: " << compression_ratio << "% reduction\n"
                  << "  ✓ Save Time: " << save_time << " ms\n"
                  << "  ✓ Write Throughput: " << save_throughput << " MB/s\n"
                  << "  ✓ Space Savings: " << (file_size_mb - cache_size_mb) << " MB\n";

        std::cout << "\n  Restoring cache from disk...\n";
        TracEon::Cache cache2;
        Timer restore_timer;
        cache2.restore(cache_file);
        double restore_time = restore_timer.elapsed_ms();
        size_t restored_sequences = cache2.size();
        double restore_throughput = cache_size_mb / (restore_time / 1000.0);

        std::cout << "\n[RESTORE RESULTS]\n"
                  << "  ✓ Sequences Restored: " << restored_sequences << "\n"
                  << "  ✓ Restore Time: " << restore_time << " ms\n"
                  << "  ✓ Read Throughput: " << restore_throughput << " MB/s\n"
                  << "  ✓ Integrity Check: " << (restored_sequences == num_sequences ? "PASS" : "FAIL") << "\n";

        if (restored_sequences != num_sequences) {
            std::cerr << "\n  ✗ WARNING: Restored " << restored_sequences << " sequences but expected " << num_sequences << "!\n";
        }

        printSection("PHASE 3: RANDOM ACCESS PERFORMANCE");
        const int num_lookups = 10000;

        std::cout << "  Running " << num_lookups << " random memory lookups...\n";
        Timer mem_timer;
        for (int i = 0; i < num_lookups; ++i) {
            std::string key = "chr" + std::to_string(i % num_sequences);
            volatile auto result = cache.get(key);
        }
        double mem_total_time = mem_timer.elapsed_ms();
        double mem_avg_time = mem_total_time / num_lookups;
        double mem_throughput = num_lookups / (mem_total_time / 1000.0);

        std::cout << "\n[MEMORY ACCESS PERFORMANCE]\n"
                  << "  ✓ Total Lookups: " << num_lookups << "\n"
                  << "  ✓ Total Time: " << mem_total_time << " ms\n"
                  << "  ✓ Average Time: " << mem_avg_time << " ms/lookup\n"
                  << "  ✓ Throughput: " << std::scientific << mem_throughput << std::fixed << " lookups/second\n";

        printSection("PHASE 4: COMPARATIVE ANALYSIS");
        double load_vs_restore_speedup = load_time / restore_time;
        double cache_vs_original_ratio = cache_size_mb / file_size_mb;

        std::cout << "\n[KEY PERFORMANCE INDICATORS]\n"
                  << "  • Original File Size: " << file_size_mb << " MB\n"
                  << "  • Cached File Size: " << cache_size_mb << " MB\n"
                  << "  • Size Ratio: " << cache_vs_original_ratio << "× (lower is better)\n"
                  << "  • Load Time: " << (load_time / 1000.0) << " seconds\n"
                  << "  • Restore Time: " << (restore_time / 1000.0) << " seconds\n"
                  << "  • Restore Speedup: " << load_vs_restore_speedup << "× faster than initial load\n"
                  << "  • Memory Footprint: " << mem_increase << " MB for " << num_sequences << " sequences\n"
                  << "  • Lookup Performance: " << (mem_avg_time * 1000) << " microseconds average\n"
                  << "\n[RECOMMENDATION]\n"
                  << "  Break-even: After " << (int)(load_time / 100) << " lookups, caching becomes beneficial\n"
                  << "  Use Case: " << (compression_ratio > 50 ? "Excellent" : compression_ratio > 30 ? "Good" : "Moderate")
                  << " space savings with " << (load_vs_restore_speedup > 5 ? "exceptional" : "good") << " reload speed\n";

        fs::remove(cache_file);

        std::cout << "\n╔" << std::string(78, '=') << "╗\n"
                  << "║" << std::string(25, ' ') << "BENCHMARK COMPLETE" << std::string(35, ' ') << "║\n"
                  << "╚" << std::string(78, '=') << "╝\n\n";

    } catch (const std::exception& e) {
        std::cerr << "\n✗ ERROR: " << e.what() << "\n";
        return 1;
    }

    return 0;
}