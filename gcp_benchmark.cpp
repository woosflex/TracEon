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
#include <zlib.h>
#include "Cache.h"

namespace fs = std::filesystem;

// Wrapper for reading both regular and gzipped files
class FileReader {
    gzFile gz_file;
    std::ifstream regular_file;
    bool is_gzipped;
    char buffer[8192];

public:
    FileReader(const std::string& filepath) {
        is_gzipped = (filepath.size() > 3 && filepath.substr(filepath.size() - 3) == ".gz");

        if (is_gzipped) {
            gz_file = gzopen(filepath.c_str(), "rb");
            if (!gz_file) {
                throw std::runtime_error("Cannot open gzipped file: " + filepath);
            }
        } else {
            regular_file.open(filepath);
            if (!regular_file.is_open()) {
                throw std::runtime_error("Cannot open file: " + filepath);
            }
        }
    }

    ~FileReader() {
        if (is_gzipped && gz_file) {
            gzclose(gz_file);
        }
    }

    bool getline(std::string& line) {
        line.clear();

        if (is_gzipped) {
            char* result = gzgets(gz_file, buffer, sizeof(buffer));
            if (!result) return false;
            line = buffer;
            // Remove newline
            if (!line.empty() && line.back() == '\n') line.pop_back();
            if (!line.empty() && line.back() == '\r') line.pop_back();
            return true;
        } else {
            return static_cast<bool>(std::getline(regular_file, line));
        }
    }

    bool is_open() const {
        return is_gzipped ? (gz_file != nullptr) : regular_file.is_open();
    }
};

// System info utilities
struct SystemInfo {
    size_t total_memory_mb;
    size_t available_memory_mb;
    std::string cpu_info;
    std::string hostname;

    SystemInfo() {
#ifdef _WIN32
        // Windows implementation
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
        // Unix/Linux/macOS implementation
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
#endif
    }
};

struct MemoryUsage {
    size_t rss_kb;  // Resident Set Size
    size_t vm_kb;   // Virtual Memory

    MemoryUsage() {
#ifdef _WIN32
        PROCESS_MEMORY_COUNTERS pmc;
        GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc));
        rss_kb = pmc.WorkingSetSize / 1024;
        vm_kb = pmc.PagefileUsage / 1024;
#else
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
#endif
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
    bool is_compressed;

    FileStats() : num_sequences(0), total_bases(0), min_length(SIZE_MAX),
                  max_length(0), avg_length(0), file_size_mb(0), is_compressed(false) {}
};

// Analyze FASTA file without loading into memory
FileStats analyzeFASTA(const std::string& filepath) {
    FileStats stats;
    stats.file_size_mb = fs::file_size(filepath) / (1024.0 * 1024.0);
    stats.is_compressed = (filepath.size() > 3 && filepath.substr(filepath.size() - 3) == ".gz");

    FileReader file(filepath);
    std::string line;
    size_t current_seq_len = 0;

    while (file.getline(line)) {
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
    stats.is_compressed = (filepath.size() > 3 && filepath.substr(filepath.size() - 3) == ".gz");

    FileReader file(filepath);
    std::string header, sequence, plus, quality;

    while (file.getline(header) &&
           file.getline(sequence) &&
           file.getline(plus) &&
           file.getline(quality)) {

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

// Decompress .gz file to temporary file for TracEon to load
std::string decompressFile(const std::string& gz_filepath) {
    std::string temp_file = "temp_decompressed_" +
                           fs::path(gz_filepath).stem().string();

    gzFile gz_in = gzopen(gz_filepath.c_str(), "rb");
    if (!gz_in) {
        throw std::runtime_error("Cannot open gzipped file: " + gz_filepath);
    }

    std::ofstream out(temp_file, std::ios::binary);
    if (!out.is_open()) {
        gzclose(gz_in);
        throw std::runtime_error("Cannot create temporary file");
    }

    char buffer[8192];
    int bytes_read;
    while ((bytes_read = gzread(gz_in, buffer, sizeof(buffer))) > 0) {
        out.write(buffer, bytes_read);
    }

    gzclose(gz_in);
    out.close();

    return temp_file;
}

// Disk-based sequential search
double benchmarkDiskSearch(const std::string& filepath, const std::string& target_id, bool is_fastq) {
    DetailedTimer timer;

    FileReader file(filepath);
    std::string line, current_id, result;
    bool found = false;