#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <vector>
#include <filesystem>
#include <algorithm>
#include <sstream>
#include <ctime>
#include <sys/resource.h>
#include "Cache.h"

namespace fs = std::filesystem;

// System info utilities
struct SystemInfo {
    size_t total_memory_mb;
    size_t available_memory_mb;
    std::string cpu_info;
    std::string hostname;
    
    SystemInfo() {
        // Get memory info
        std::ifstream meminfo("/proc/meminfo");
        std::string line;
        while (std::getline(meminfo, line)) {
            if (line.substr(0, 9) == "MemTotal:") {
                total_memory_mb = std::stoul(line.substr(10)) / 1024;
            }
            if (line.substr(0, 13) == "MemAvailable:") {
                available_memory_mb = std::stoul(line.substr(14)) / 1024;
            }
        }
        
        // Get CPU info
        std::ifstream cpuinfo("/proc/cpuinfo");
        while (std::getline(cpuinfo, line)) {
            if (line.substr(0, 10) == "model name") {
                cpu_info = line.substr(line.find(':') + 2);
                break;
            }
        }
        
        // Get hostname
        char hostname_buf[256];
        gethostname(hostname_buf, sizeof(hostname_buf));
        hostname = hostname_buf;
    }
};

struct MemoryUsage {
    size_t rss_kb;  // Resident Set Size
    size_t vm_kb;   // Virtual Memory
    
    MemoryUsage() {
        std::ifstream status("/proc/self/status");
        std::string line;
        while (std::getline(status, line)) {
            if (line.substr(0, 6) == "VmRSS:") {
                rss_kb = std::stoul(line.substr(7));
            }
            if (line.substr(0, 6) == "VmSize:") {
                vm_kb = std::stoul(line.substr(7));
            }
        }
    }
    
    double rss_mb() const { return rss_kb / 1024.0; }
    double vm_mb() const { return vm_kb / 1024.0; }
};

class DetailedTimer {
    std::chrono::high_resolution_clock::time_point start_time;
    MemoryUsage start_mem;
    
public:
    struct Result {
        double elapsed_ms;
        double memory_delta_mb;
    };
    
    DetailedTimer() {
        start_mem = MemoryUsage();
        start_time = std::chrono::high_resolution_clock::now();
    }
    
    Result stop() {
        auto end_time = std::chrono::high_resolution_clock::now();
        MemoryUsage end_mem;
        
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
        
        Result r;
        r.elapsed_ms = duration.count() / 1000.0;
        r.memory_delta_mb = end_mem.rss_mb() - start_mem.rss_mb();
        return r;
    }
};

struct FileStats {
    size_t num_sequences;
    size_t total_bases;
    size_t min_length;
    size_t max_length;
    double avg_length;
    double file_size_mb;
    
    FileStats() : num_sequences(0), total_bases(0), min_length(SIZE_MAX), 
                  max_length(0), avg_length(0), file_size_mb(0) {}
};

// Analyze FASTA file without loading into memory
FileStats analyzeFASTA(const std::string& filepath) {
    FileStats stats;
    stats.file_size_mb = fs::file_size(filepath) / (1024.0 * 1024.0);
    
    std::ifstream file(filepath);
    std::string line;
    size_t current_seq_len = 0;
    
    while (std::getline(file, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        
        if (line[0] == '>') {
            if (current_seq_len > 0) {
                stats.num_sequences++;
                stats.total_bases += current_seq_len;
                stats.min_length = std::min(stats.min_length, current_seq_len);
                stats.max_length = std::max(stats.max_length, current_seq_len);
                current_seq_len = 0;
            }
        } else {
            current_seq_len += line.length();
        }
    }
    
    // Handle last sequence
    if (current_seq_len > 0) {
        stats.num_sequences++;
        stats.total_bases += current_seq_len;
        stats.min_length = std::min(stats.min_length, current_seq_len);
        stats.max_length = std::max(stats.max_length, current_seq_len);
    }
    
    if (stats.num_sequences > 0) {
        stats.avg_length = static_cast<double>(stats.total_bases) / stats.num_sequences;
    }
    
    return stats;
}

// Analyze FASTQ file
FileStats analyzeFASTQ(const std::string& filepath) {
    FileStats stats;
    stats.file_size_mb = fs::file_size(filepath) / (1024.0 * 1024.0);
    
    std::ifstream file(filepath);
    std::string header, sequence, plus, quality;
    
    while (std::getline(file, header) &&
           std::getline(file, sequence) &&
           std::getline(file, plus) &&
           std::getline(file, quality)) {
        
        if (!sequence.empty() && sequence.back() == '\r') sequence.pop_back();
        
        size_t len = sequence.length();
        stats.num_sequences++;
        stats.total_bases += len;
        stats.min_length = std::min(stats.min_length, len);
        stats.max_length = std::max(stats.max_length, len);
    }
    
    if (stats.num_sequences > 0) {
        stats.avg_length = static_cast<double>(stats.total_bases) / stats.num_sequences;
    }
    
    return stats;
}

// Disk-based sequential search
double benchmarkDiskSearch(const std::string& filepath, const std::string& target_id, bool is_fastq) {
    DetailedTimer timer;
    
    std::ifstream file(filepath);
    std::string line, current_id, result;
    bool found = false;
    
    if (is_fastq) {
        std::string header, sequence, plus, quality;
        while (std::getline(file, header) && !found) {
            std::getline(file, sequence);
            std::getline(file, plus);
            std::getline(file, quality);
            
            if (!header.empty() && header[0] == '@') {
                size_t space = header.find(' ');
                current_id = header.substr(1, space - 1);
                if (current_id == target_id) {
                    result = sequence;
                    found = true;
                }
            }
        }
    } else {
        std::string sequence;
        while (std::getline(file, line) && !found) {
            if (!line.empty() && line[0] == '>') {
                if (!current_id.empty() && current_id == target_id) {
                    result = sequence;
                    found = true;
                    break;
                }
                size_t space = line.find(' ');
                current_id = line.substr(1, space - 1);
                sequence.clear();
            } else {
                sequence += line;
            }
        }
    }
    
    return timer.stop().elapsed_ms;
}

void printHeader(const std::string& title) {
    std::cout << "\n" << std::string(80, '=') << "\n";
    std::cout << title << "\n";
    std::cout << std::string(80, '=') << "\n";
}

void printSubHeader(const std::string& title) {
    std::cout << "\n" << std::string(80, '-') << "\n";
    std::cout << title << "\n";
    std::cout << std::string(80, '-') << "\n";
}

std::string getCurrentTimestamp() {
    auto now = std::time(nullptr);
    std::stringstream ss;
    ss << std::put_time(std::localtime(&now), "%Y-%m-%d %H:%M:%S");
    return ss.str();
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_file> [output_report.txt]\n";
        std::cerr << "\nSupported formats: FASTA (.fasta, .fa) and FASTQ (.fastq, .fq)\n";
        std::cerr << "\nExample:\n";
        std::cerr << "  " << argv[0] << " /path/to/genome.fasta report.txt\n";
        return 1;
    }
    
    std::string input_file = argv[1];
    std::string output_report = (argc >= 3) ? argv[2] : "traceon_benchmark_report.txt";
    
    if (!fs::exists(input_file)) {
        std::cerr << "Error: File not found: " << input_file << "\n";
        return 1;
    }
    
    // Detect file type
    std::string ext = fs::path(input_file).extension().string();
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
    bool is_fastq = (ext == ".fastq" || ext == ".fq");
    bool is_fasta = (ext == ".fasta" || ext == ".fa");
    
    if (!is_fasta && !is_fastq) {
        std::cerr << "Error: Unsupported file format. Use .fasta, .fa, .fastq, or .fq\n";
        return 1;
    }
    
    // Redirect output to both console and file
    std::ofstream report(output_report);
    auto cout_buf = std::cout.rdbuf();
    
    // Print to both console and file
    #define DUAL_PRINT(x) do { std::cout << x; report << x; } while(0)
    
    DUAL_PRINT("\n╔════════════════════════════════════════════════════════════════════════════╗\n");
    DUAL_PRINT("║                  TracEon NGS Benchmark - Real-World Data                  ║\n");
    DUAL_PRINT("╚════════════════════════════════════════════════════════════════════════════╝\n");
    
    // System information
    SystemInfo sysinfo;
    DUAL_PRINT("\nBenchmark Information:\n");
    DUAL_PRINT("  Timestamp: " << getCurrentTimestamp() << "\n");
    DUAL_PRINT("  Hostname: " << sysinfo.hostname << "\n");
    DUAL_PRINT("  CPU: " << sysinfo.cpu_info << "\n");
    DUAL_PRINT("  Total Memory: " << sysinfo.total_memory_mb << " MB\n");
    DUAL_PRINT("  Available Memory: " << sysinfo.available_memory_mb << " MB\n");
    DUAL_PRINT("  Input File: " << input_file << "\n");
    DUAL_PRINT("  File Type: " << (is_fastq ? "FASTQ" : "FASTA") << "\n");
    
    // Analyze input file
    printHeader("Phase 1: File Analysis");
    DUAL_PRINT("Analyzing input file structure...\n");
    
    FileStats stats = is_fastq ? analyzeFASTQ(input_file) : analyzeFASTA(input_file);
    
    DUAL_PRINT("\nFile Statistics:\n");
    DUAL_PRINT("  File Size: " << std::fixed << std::setprecision(2) << stats.file_size_mb << " MB\n");
    DUAL_PRINT("  Number of Sequences: " << stats.num_sequences << "\n");
    DUAL_PRINT("  Total Bases: " << stats.total_bases << "\n");
    DUAL_PRINT("  Sequence Length (min/avg/max): " 
                << stats.min_length << " / " 
                << std::fixed << std::setprecision(1) << stats.avg_length << " / " 
                << stats.max_length << " bp\n");
    
    // Load into TracEon cache
    printHeader("Phase 2: Loading into TracEon Cache");
    
    MemoryUsage mem_before;
    DUAL_PRINT("Memory before load: " << std::fixed << std::setprecision(2) 
                << mem_before.rss_mb() << " MB (RSS)\n");
    
    TracEon::Cache cache;
    DetailedTimer load_timer;
    DUAL_PRINT("\nLoading file into cache...\n");
    cache.loadFile(input_file);
    auto load_result = load_timer.stop();
    
    MemoryUsage mem_after;
    DUAL_PRINT("\nLoad Complete:\n");
    DUAL_PRINT("  Time: " << std::fixed << std::setprecision(2) << load_result.elapsed_ms << " ms\n");
    DUAL_PRINT("  Sequences loaded: " << cache.size() << "\n");
    DUAL_PRINT("  Memory after load: " << mem_after.rss_mb() << " MB (RSS)\n");
    DUAL_PRINT("  Memory increase: " << load_result.memory_delta_mb << " MB\n");
    DUAL_PRINT("  Memory efficiency: " << std::fixed << std::setprecision(2) 
                << (stats.file_size_mb / load_result.memory_delta_mb) << "× (file/memory ratio)\n");
    
    // Save to binary cache
    printHeader("Phase 3: Binary Cache Persistence");
    
    std::string cache_file = "traceon_cache.bin";
    DUAL_PRINT("Saving to binary cache file...\n");
    
    DetailedTimer save_timer;
    cache.save(cache_file);
    auto save_result = save_timer.stop();
    
    double cache_size_mb = fs::file_size(cache_file) / (1024.0 * 1024.0);
    double compression_ratio = (1.0 - cache_size_mb / stats.file_size_mb) * 100.0;
    
    DUAL_PRINT("\nSave Complete:\n");
    DUAL_PRINT("  Time: " << std::fixed << std::setprecision(2) << save_result.elapsed_ms << " ms\n");
    DUAL_PRINT("  Cache file: " << cache_file << "\n");
    DUAL_PRINT("  Cache size: " << cache_size_mb << " MB\n");
    DUAL_PRINT("  Compression ratio: " << std::fixed << std::setprecision(1) 
                << compression_ratio << "%\n");
    DUAL_PRINT("  Write speed: " << std::fixed << std::setprecision(2) 
                << (cache_size_mb / (save_result.elapsed_ms / 1000.0)) << " MB/s\n");
    
    // Restore from binary cache
    DUAL_PRINT("\nTesting cache restore...\n");
    TracEon::Cache cache2;
    DetailedTimer restore_timer;
    cache2.restore(cache_file);
    auto restore_result = restore_timer.stop();
    
    DUAL_PRINT("Restore Complete:\n");
    DUAL_PRINT("  Time: " << std::fixed << std::setprecision(2) << restore_result.elapsed_ms << " ms\n");
    DUAL_PRINT("  Sequences restored: " << cache2.size() << "\n");
    DUAL_PRINT("  Read speed: " << std::fixed << std::setprecision(2) 
                << (cache_size_mb / (restore_result.elapsed_ms / 1000.0)) << " MB/s\n");
    
    // Random access performance
    printHeader("Phase 4: Random Access Performance");
    
    const size_t NUM_LOOKUPS = std::min(size_t(1000), stats.num_sequences);
    std::vector<std::string> test_ids;
    
    // Get sample sequence IDs
    DUAL_PRINT("Collecting sample sequence IDs...\n");
    std::ifstream id_file(input_file);
    std::string line;
    size_t count = 0;
    size_t step = std::max(size_t(1), stats.num_sequences / NUM_LOOKUPS);
    
    while (std::getline(id_file, line) && test_ids.size() < NUM_LOOKUPS) {
        if (!line.empty() && (line[0] == '>' || line[0] == '@')) {
            if (count % step == 0) {
                size_t space = line.find(' ');
                test_ids.push_back(line.substr(1, space - 1));
            }
            count++;
        }
    }
    
    DUAL_PRINT("Testing " << test_ids.size() << " random lookups...\n\n");
    
    // Memory-based lookups
    std::vector<double> memory_times;
    for (const auto& id : test_ids) {
        DetailedTimer t;
        volatile std::string result = cache.get(id);
        memory_times.push_back(t.stop().elapsed_ms);
    }
    
    // Disk-based lookups (sample a few to avoid taking too long)
    std::vector<double> disk_times;
    size_t disk_samples = std::min(size_t(10), test_ids.size());
    DUAL_PRINT("Comparing with disk-based access (" << disk_samples << " samples)...\n");
    
    for (size_t i = 0; i < disk_samples; i++) {
        size_t idx = i * (test_ids.size() / disk_samples);
        double time = benchmarkDiskSearch(input_file, test_ids[idx], is_fastq);
        disk_times.push_back(time);
    }
    
    // Calculate statistics
    auto calc_stats = [](const std::vector<double>& times) {
        double sum = 0, min_t = times[0], max_t = times[0];
        for (double t : times) {
            sum += t;
            min_t = std::min(min_t, t);
            max_t = std::max(max_t, t);
        }
        return std::make_tuple(min_t, sum / times.size(), max_t);
    };
    
    auto [mem_min, mem_avg, mem_max] = calc_stats(memory_times);
    auto [disk_min, disk_avg, disk_max] = calc_stats(disk_times);
    
    DUAL_PRINT("\nMemory-Based Access (TracEon):\n");
    DUAL_PRINT("  Lookups: " << memory_times.size() << "\n");
    DUAL_PRINT("  Min time: " << std::fixed << std::setprecision(4) << mem_min << " ms\n");
    DUAL_PRINT("  Avg time: " << mem_avg << " ms\n");
    DUAL_PRINT("  Max time: " << mem_max << " ms\n");
    DUAL_PRINT("  Throughput: " << std::fixed << std::setprecision(0) 
                << (1000.0 / mem_avg) << " lookups/second\n");
    
    DUAL_PRINT("\nDisk-Based Access:\n");
    DUAL_PRINT("  Lookups: " << disk_times.size() << "\n");
    DUAL_PRINT("  Min time: " << std::fixed << std::setprecision(2) << disk_min << " ms\n");
    DUAL_PRINT("  Avg time: " << disk_avg << " ms\n");
    DUAL_PRINT("  Max time: " << disk_max << " ms\n");
    
    DUAL_PRINT("\nSpeedup: " << std::fixed << std::setprecision(1) 
                << (disk_avg / mem_avg) << "× faster\n");
    
    // Summary
    printHeader("Benchmark Summary");
    
    DUAL_PRINT("\nKey Performance Metrics:\n");
    DUAL_PRINT("  • File size: " << std::fixed << std::setprecision(2) << stats.file_size_mb << " MB\n");
    DUAL_PRINT("  • Sequences: " << stats.num_sequences << "\n");
    DUAL_PRINT("  • Load time: " << load_result.elapsed_ms << " ms\n");
    DUAL_PRINT("  • Memory overhead: " << load_result.memory_delta_mb << " MB\n");
    DUAL_PRINT("  • Compression: " << std::fixed << std::setprecision(1) << compression_ratio << "%\n");
    DUAL_PRINT("  • Avg lookup (memory): " << std::fixed << std::setprecision(4) << mem_avg << " ms\n");
    DUAL_PRINT("  • Avg lookup (disk): " << std::fixed << std::setprecision(2) << disk_avg << " ms\n");
    DUAL_PRINT("  • Speedup: " << std::fixed << std::setprecision(1) << (disk_avg / mem_avg) << "×\n");
    
    // Cost-benefit analysis
    double break_even = load_result.elapsed_ms / (disk_avg - mem_avg);
    DUAL_PRINT("\nCost-Benefit Analysis:\n");
    DUAL_PRINT("  • One-time load cost: " << std::fixed << std::setprecision(2) 
                << load_result.elapsed_ms << " ms\n");
    DUAL_PRINT("  • Time saved per lookup: " << (disk_avg - mem_avg) << " ms\n");
    DUAL_PRINT("  • Break-even point: " << std::fixed << std::setprecision(0) 
                << break_even << " lookups\n");
    DUAL_PRINT("  • After 1000 lookups: " << std::fixed << std::setprecision(2)
                << ((disk_avg * 1000 - load_result.elapsed_ms - mem_avg * 1000) / 1000.0) 
                << " seconds saved\n");
    
    DUAL_PRINT("\n" << std::string(80, '=') << "\n");
    DUAL_PRINT("Benchmark complete! Report saved to: " << output_report << "\n");
    DUAL_PRINT(std::string(80, '=') << "\n\n");
    
    // Cleanup
    fs::remove(cache_file);
    
    report.close();
    return 0;
}
