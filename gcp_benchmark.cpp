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

// Conditional includes for Windows/Linux
#ifdef _WIN32
#include <windows.h>
#include <psapi.h>
#else
#include <unistd.h>
#endif


namespace fs = std::filesystem;

// Wrapper for reading both regular and gzipped files
class FileReader {
    gzFile gz_file;
    std::ifstream regular_file;
    bool is_gzipped;
    char buffer[8192];

public:
    FileReader(const std::string& filepath) : gz_file(nullptr) {
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
            // Remove newline characters
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
        total_memory_mb = 0;
        available_memory_mb = 0;
        while (std::getline(meminfo, line)) {
            if (line.rfind("MemTotal:", 0) == 0) {
                total_memory_mb = std::stoul(line.substr(10)) / 1024;
            }
            if (line.rfind("MemAvailable:", 0) == 0) {
                available_memory_mb = std::stoul(line.substr(14)) / 1024;
            }
        }

        // Get CPU info
        std::ifstream cpuinfo("/proc/cpuinfo");
        while (std::getline(cpuinfo, line)) {
            if (line.rfind("model name", 0) == 0) {
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
    size_t rss_kb; // Resident Set Size

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

class DetailedTimer {
    std::chrono::high_resolution_clock::time_point start_time;
public:
    struct Result {
        double elapsed_ms;
    };

    DetailedTimer() {
        start_time = std::chrono::high_resolution_clock::now();
    }

    Result stop() {
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
        Result r;
        r.elapsed_ms = duration.count() / 1000.0;
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
    std::string line, result;
    bool found = false;

    if (is_fastq) {
        std::string header, sequence, plus, quality;
        while (file.getline(header) && !found) {
            file.getline(sequence);
            file.getline(plus);
            file.getline(quality);

            if (!header.empty() && header[0] == '@') {
                size_t space = header.find(' ');
                std::string current_id = header.substr(1, space != std::string::npos ? space - 1 : std::string::npos);
                if (current_id == target_id) {
                    result = sequence;
                    found = true;
                }
            }
        }
    } else { // FASTA
        std::string current_id;
        std::stringstream sequence_stream;
        while(file.getline(line) && !found) {
            if(!line.empty() && line[0] == '>') {
                if(!current_id.empty()) {
                    if(current_id == target_id) {
                        result = sequence_stream.str();
                        found = true;
                        break;
                    }
                }
                size_t space = line.find(' ');
                current_id = line.substr(1, space != std::string::npos ? space - 1 : std::string::npos);
                sequence_stream.str(""); // clear
            } else {
                sequence_stream << line;
            }
        }
        // check last sequence
        if(!found && !current_id.empty() && current_id == target_id) {
            result = sequence_stream.str();
        }
    }

    return timer.stop().elapsed_ms;
}

// Function to get sample IDs from a file
std::vector<std::string> getSampleIds(const std::string& filepath, bool is_fastq, int count) {
    std::vector<std::string> ids;
    FileReader file(filepath);
    std::string line;

    while (ids.size() < count && file.getline(line)) {
         if ((is_fastq && line[0] == '@') || (!is_fastq && line[0] == '>')) {
            size_t space = line.find(' ');
            ids.push_back(line.substr(1, space != std::string::npos ? space - 1 : std::string::npos));
            if(is_fastq) {
                // Skip the next 3 lines of the FASTQ record
                for(int i=0; i<3; ++i) file.getline(line);
            }
        }
    }
    return ids;
}


void printHeader(const std::string& title) {
    std::string border(80, '=');
    std::cout << border << std::endl;
    std::cout << title << std::endl;
    std::cout << border << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <path_to_fasta_or_fastq_file>" << std::endl;
        return 1;
    }

    std::string filepath = argv[1];

    try {
        std::cout << std::fixed << std::setprecision(4);

        // --- Header ---
        std::cout << "╔════════════════════════════════════════════════════════════════════════════╗" << std::endl;
        std::cout << "║                  TracEon NGS Benchmark - Real-World Data                  ║" << std::endl;
        std::cout << "╚════════════════════════════════════════════════════════════════════════════╝" << std::endl << std::endl;

        // --- System Info ---
        SystemInfo sys_info;
        auto t = std::time(nullptr);
        auto tm = *std::localtime(&t);
        std::cout << "Benchmark Information:" << std::endl;
        std::cout << "  Timestamp: " << std::put_time(&tm, "%Y-%m-%d %H:%M:%S") << std::endl;
        std::cout << "  Hostname: " << sys_info.hostname << std::endl;
        std::cout << "  CPU: " << sys_info.cpu_info << std::endl;
        std::cout << "  Total Memory: " << sys_info.total_memory_mb << " MB" << std::endl;
        std::cout << "  Available Memory: " << sys_info.available_memory_mb << " MB" << std::endl;
        std::cout << "  Input File: " << filepath << std::endl;

        // --- Phase 1: Analysis ---
        printHeader("Phase 1: File Analysis");
        std::cout << "Analyzing input file structure..." << std::endl << std::endl;

        bool is_fastq = (filepath.find(".fastq") != std::string::npos || filepath.find(".fq") != std::string::npos);
        std::cout << "  File Type: " << (is_fastq ? "FASTQ" : "FASTA") << std::endl;

        FileStats stats = is_fastq ? analyzeFASTQ(filepath) : analyzeFASTA(filepath);

        std::cout << "File Statistics:" << std::endl;
        std::cout << "  File Size: " << stats.file_size_mb << " MB" << std::endl;
        std::cout << "  Number of Sequences: " << stats.num_sequences << std::endl;
        std::cout << "  Total Bases: " << stats.total_bases << std::endl;
        std::cout << "  Sequence Length (min/avg/max): " << stats.min_length << " / " << std::fixed << std::setprecision(2) << stats.avg_length << " / " << stats.max_length << " bp" << std::endl;

        // --- Phase 2: Loading ---
        printHeader("Phase 2: Loading into TracEon Cache");

        std::string file_to_load = filepath;
        std::string temp_file_to_delete = "";
        if (stats.is_compressed) {
            std::cout << "Decompressing file for loading..." << std::endl;
            temp_file_to_delete = decompressFile(filepath);
            file_to_load = temp_file_to_delete;
        }

        MemoryUsage mem_before = MemoryUsage();
        std::cout << "Memory before load: " << mem_before.rss_mb() << " MB (RSS)" << std::endl << std::endl;

        std::cout << "Loading file into cache..." << std::endl;
        DetailedTimer load_timer;
        TracEon::Cache cache;
        cache.loadFile(file_to_load);
        auto load_result = load_timer.stop();
        MemoryUsage mem_after = MemoryUsage();

        std::cout << "\nLoad Complete:" << std::endl;
        std::cout << "  Time: " << load_result.elapsed_ms << " ms" << std::endl;
        std::cout << "  Sequences loaded: " << cache.size() << std::endl;
        std::cout << "  Memory after load: " << mem_after.rss_mb() << " MB (RSS)" << std::endl;
        double mem_increase = mem_after.rss_mb() - mem_before.rss_mb();
        std::cout << "  Memory increase: " << mem_increase << " MB" << std::endl;
        std::cout << "  Memory efficiency: " << (stats.file_size_mb / mem_increase) << "× (file/memory ratio)" << std::endl;

        // --- Phase 3: Persistence ---
        printHeader("Phase 3: Binary Cache Persistence");

        std::string cache_path = "traceon_cache.bin";
        std::cout << "Saving to binary cache file..." << std::endl;
        DetailedTimer save_timer;
        cache.save(cache_path);
        auto save_result = save_timer.stop();

        double cache_size_mb = fs::file_size(cache_path) / (1024.0 * 1024.0);

        std::cout << "\nSave Complete:" << std::endl;
        std::cout << "  Time: " << save_result.elapsed_ms << " ms" << std::endl;
        std::cout << "  Cache file: " << cache_path << std::endl;
        std::cout << "  Cache size: " << cache_size_mb << " MB" << std::endl;
        std::cout << "  Compression ratio: " << (1.0 - (cache_size_mb / stats.file_size_mb)) * 100.0 << "%" << std::endl;
        std::cout << "  Write speed: " << (cache_size_mb / (save_result.elapsed_ms / 1000.0)) << " MB/s" << std::endl;

        std::cout << "\nTesting cache restore..." << std::endl;
        TracEon::Cache restore_cache;
        DetailedTimer restore_timer;
        restore_cache.restore(cache_path);
        auto restore_result = restore_timer.stop();

        std::cout << "Restore Complete:" << std::endl;
        std::cout << "  Time: " << restore_result.elapsed_ms << " ms" << std::endl;
        std::cout << "  Sequences restored: " << restore_cache.size() << std::endl;
        std::cout << "  Read speed: " << (cache_size_mb / (restore_result.elapsed_ms / 1000.0)) << " MB/s" << std::endl;

        // --- Phase 4: Random Access ---
        printHeader("Phase 4: Random Access Performance");

        std::cout << "Collecting sample sequence IDs..." << std::endl;
        std::vector<std::string> sample_ids = getSampleIds(file_to_load, is_fastq, 10);

        if (sample_ids.empty()) {
            throw std::runtime_error("Could not find any sequence IDs in the input file to use for lookup test.");
        }

        // --- Memory Access ---
        std::cout << "Testing 1000 random lookups..." << std::endl << std::endl;
        std::cout << "Memory-Based Access (TracEon):" << std::endl;
        const int num_lookups = 1000;
        double min_mem_time = 1e9, max_mem_time = 0;
        DetailedTimer mem_lookup_timer;
        for (int i = 0; i < num_lookups; ++i) {
            DetailedTimer single_lookup_timer;
            cache.get(sample_ids[i % sample_ids.size()]);
            auto single_res = single_lookup_timer.stop();
            min_mem_time = std::min(min_mem_time, single_res.elapsed_ms);
            max_mem_time = std::max(max_mem_time, single_res.elapsed_ms);
        }
        auto mem_lookup_result = mem_lookup_timer.stop();
        double avg_mem_time = mem_lookup_result.elapsed_ms / num_lookups;

        std::cout << "  Lookups: " << num_lookups << std::endl;
        std::cout << "  Min time: " << min_mem_time << " ms" << std::endl;
        std::cout << "  Avg time: " << avg_mem_time << " ms" << std::endl;
        std::cout << "  Max time: " << max_mem_time << " ms" << std::endl;
        std::cout << "  Throughput: " << num_lookups / (mem_lookup_result.elapsed_ms / 1000.0) << " lookups/second" << std::endl;

        // --- Disk Access ---
        std::cout << "\nDisk-Based Access:" << std::endl;
        const int num_disk_lookups = 10;
        double total_disk_time = 0;
        for (int i = 0; i < num_disk_lookups; ++i) {
             total_disk_time += benchmarkDiskSearch(file_to_load, sample_ids[i % sample_ids.size()], is_fastq);
        }
        double avg_disk_time = num_disk_lookups > 0 ? total_disk_time / num_disk_lookups : 0;
        std::cout << "  Lookups: " << num_disk_lookups << std::endl;
        std::cout << "  Avg time: " << avg_disk_time << " ms" << std::endl;

        std::cout << "\nSpeedup: " << avg_disk_time / avg_mem_time << "× faster" << std::endl;

        // --- Summary ---
        printHeader("Benchmark Summary");
        std::cout << "Key Performance Metrics:" << std::endl;
        std::cout << "  • File size: " << stats.file_size_mb << " MB" << std::endl;
        std::cout << "  • Sequences: " << cache.size() << std::endl;
        std::cout << "  • Load time: " << load_result.elapsed_ms << " ms" << std::endl;
        std::cout << "  • Memory overhead: " << mem_increase << " MB" << std::endl;
        std::cout << "  • Compression: " << (1.0 - (cache_size_mb / stats.file_size_mb)) * 100.0 << "%" << std::endl;
        std::cout << "  • Avg lookup (memory): " << avg_mem_time << " ms" << std::endl;
        std::cout << "  • Avg lookup (disk): " << avg_disk_time << " ms" << std::endl;
        std::cout << "  • Speedup: " << avg_disk_time / avg_mem_time << "×" << std::endl;

        // --- Cleanup ---
        std::cout << "\nCleaning up temporary files..." << std::endl;
        if (!temp_file_to_delete.empty()) fs::remove(temp_file_to_delete);
        fs::remove(cache_path);

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}

