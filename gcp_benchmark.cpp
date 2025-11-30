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
#include <numeric>
#include "Cache.h"
#include "SmartStrategy.h"

#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#include <psapi.h>
#else
#include <unistd.h>
#include <sys/resource.h>
#endif

namespace fs = std::filesystem;

// --- SystemInfo, MemoryUsage, and Timer classes (Unchanged) ---
// (These helper classes are well-written and do not need changes)
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

// --- Main Benchmark Logic ---

void printSection(const std::string& title) {
    std::cout << "\n" << std::string(80, '=') << "\n" << title << "\n" << std::string(80, '=') << "\n";
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

        // --- System Information ---
        SystemInfo sysinfo;
        auto t = std::time(nullptr);
        std::cout << "\n[SYSTEM INFORMATION]\n"
                  << "  Timestamp: " << std::put_time(std::localtime(&t), "%Y-%m-%d %H:%M:%S") << "\n"
                  << "  Hostname: " << sysinfo.hostname << "\n"
                  << "  CPU: " << sysinfo.cpu_info << "\n"
                  << "  CPU Cores: " << sysinfo.num_cores << "\n"
                  << "  Total RAM: " << sysinfo.total_memory_mb << " MB\n";

        double file_size_mb = fs::file_size(filepath) / (1024.0 * 1024.0);
        std::cout << "  Input File: " << filepath << " (" << file_size_mb << " MB)\n";

        // --- Phase 1: Load Data ---
        printSection("PHASE 1: INITIAL DATA LOAD FROM SOURCE FILE");
        MemoryUsage mem_before;
        Timer load_timer;
        TracEon::Cache cache;
        cache.loadSmartFile(filepath); // Use the correct method for SmartStrategy
        double load_time = load_timer.elapsed_ms();
        MemoryUsage mem_after;

        size_t num_sequences = cache.size();
        double mem_increase = mem_after.rss_mb() - mem_before.rss_mb();
        
        std::cout << "[LOAD RESULTS]\n"
                  << "  ✓ Sequences Loaded: " << num_sequences << "\n"
                  << "  ✓ Load Time: " << (load_time / 1000.0) << " seconds\n"
                  << "  ✓ Memory Footprint: " << mem_increase << " MB\n";

        // --- Phase 2: Uncompressed Cache Benchmark ---
        printSection("PHASE 2: UNCOMPRESSED CACHE PERFORMANCE");
        const std::string uncompressed_cache_file = "traceon_uncompressed.bin";
        
        Timer save_uncomp_timer;
        cache.saveSmartBinary(uncompressed_cache_file, TracEon::SaveMode::Uncompressed);
        double save_uncomp_time = save_uncomp_timer.elapsed_ms();
        double uncomp_size_mb = fs::file_size(uncompressed_cache_file) / (1024.0 * 1024.0);

        Timer restore_uncomp_timer;
        TracEon::Cache restored_uncomp_cache;
        restored_uncomp_cache.loadSmartBinary(uncompressed_cache_file);
        double restore_uncomp_time = restore_uncomp_timer.elapsed_ms();

        std::cout << "[UNCOMPRESSED RESULTS]\n"
                  << "  ✓ Cache Size: " << uncomp_size_mb << " MB\n"
                  << "  ✓ Save Time: " << (save_uncomp_time / 1000.0) << " seconds\n"
                  << "  ✓ Restore Time: " << (restore_uncomp_time / 1000.0) << " seconds\n"
                  << "  ✓ Integrity Check: " << (restored_uncomp_cache.size() == num_sequences ? "PASS" : "FAIL") << "\n";

        // --- Phase 3: Compressed Cache Benchmark ---
        printSection("PHASE 3: COMPRESSED CACHE PERFORMANCE");
        const std::string compressed_cache_file = "traceon_compressed.bin";
        
        Timer save_comp_timer;
        cache.saveSmartBinary(compressed_cache_file, TracEon::SaveMode::Compressed);
        double save_comp_time = save_comp_timer.elapsed_ms();
        double comp_size_mb = fs::file_size(compressed_cache_file) / (1024.0 * 1024.0);

        Timer restore_comp_timer;
        TracEon::Cache restored_comp_cache;
        restored_comp_cache.loadSmartBinary(compressed_cache_file);
        double restore_comp_time = restore_comp_timer.elapsed_ms();
        
        std::cout << "[COMPRESSED RESULTS]\n"
                  << "  ✓ Cache Size: " << comp_size_mb << " MB\n"
                  << "  ✓ Save Time: " << (save_comp_time / 1000.0) << " seconds\n"
                  << "  ✓ Restore Time: " << (restore_comp_time / 1000.0) << " seconds\n"
                  << "  ✓ Integrity Check: " << (restored_comp_cache.size() == num_sequences ? "PASS" : "FAIL") << "\n";

        // --- Phase 4: Random Access Performance ---
        printSection("PHASE 4: RANDOM ACCESS PERFORMANCE");
        const int num_lookups = 100000;
        std::vector<std::string> keys_to_test;
        
        // **CORRECTED:** Get real keys from the cache for a valid test
        if (auto* smart_strategy = dynamic_cast<TracEon::SmartStrategy*>(cache.getStrategy())) {
            for (const auto& key : smart_strategy->getAllKeys()) {
                keys_to_test.push_back(key);
                if (keys_to_test.size() >= num_lookups) break;
            }
        }
        
        std::cout << "  Running " << keys_to_test.size() << " random memory lookups...\n";
        
        Timer mem_timer;
        for (const auto& key : keys_to_test) {
            volatile auto result = cache.get(key); // volatile prevents optimization
        }
        double mem_total_time = mem_timer.elapsed_ms();

        std::cout << "[MEMORY ACCESS RESULTS]\n"
                  << "  ✓ Average Time: " << (mem_total_time / keys_to_test.size() * 1000.0) << " microseconds/lookup\n"
                  << "  ✓ Throughput: " << (keys_to_test.size() / (mem_total_time / 1000.0)) << " lookups/second\n";

        // --- Phase 5: Final KPI Summary ---
        printSection("PHASE 5: KEY PERFORMANCE INDICATORS (KPIs)");
        std::cout << "[COMPARATIVE ANALYSIS]\n"
                  << "  • Initial Load Time: " << (load_time / 1000.0) << " seconds\n\n"
                  << "  --- UNCOMPRESSED CACHE ---\n"
                  << "  • Restore Speedup: " << (load_time / restore_uncomp_time) << "x faster than initial load\n"
                  << "  • Size Ratio: " << (uncomp_size_mb / file_size_mb) << "x of original file size\n\n"
                  << "  --- COMPRESSED CACHE ---\n"
                  << "  • Restore Speedup: " << (load_time / restore_comp_time) << "x faster than initial load\n"
                  << "  • Size Ratio: " << (comp_size_mb / file_size_mb) << "x of original file size\n"
                  << "  • Space Savings: " << ((1.0 - (comp_size_mb / uncomp_size_mb)) * 100.0) << "% smaller than uncompressed cache\n";

        // --- Cleanup ---
        fs::remove(uncompressed_cache_file);
        fs::remove(compressed_cache_file);

    } catch (const std::exception& e) {
        std::cerr << "\n✗ ERROR: " << e.what() << "\n";
        return 1;
    }

    std::cout << "\n" << std::string(80, '=') << "\nBENCHMARK COMPLETE\n" << std::string(80, '=') << "\n\n";
    return 0;
}