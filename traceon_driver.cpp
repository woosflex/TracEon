#include <iostream>
#include <chrono>
#include <string>
#include <vector>
#include <cstdlib>
#include <random>
#include <iomanip> // For std::fixed
#include "Cache.h"
#include "SmartStrategy.h" 

volatile size_t global_sink = 0;

class Timer {
    std::chrono::high_resolution_clock::time_point start;
public:
    Timer() : start(std::chrono::high_resolution_clock::now()) {}
    double elapsed_sec() {
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff = end - start;
        return diff.count();
    }
};

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <mode> <file> [args]" << std::endl;
        return 1;
    }

    std::string mode = argv[1];
    std::string filepath = argv[2];
    TracEon::Cache cache;

    try {
        if (mode == "load") {
            Timer t;
            cache.loadFile(filepath);
            std::cout << "Load_Time_s: " << t.elapsed_sec() << std::endl;
            std::cout << "Sequences: " << cache.size() << std::endl;
        } 
        else if (mode == "save") {
            cache.loadFile(filepath);
            if (argc < 4) return 1;
            Timer t;
            cache.save(argv[3]);
            std::cout << "Save_Time_s: " << t.elapsed_sec() << std::endl;
        }
        else if (mode == "restore") {
            Timer t;
            cache.restore(filepath);
            std::cout << "Restore_Time_s: " << t.elapsed_sec() << std::endl;
            std::cout << "Sequences: " << cache.size() << std::endl;
        }
        else if (mode == "lookup") {
            if (argc < 6) return 1;
            int num_lookups = std::stoi(argv[3]);
            std::string prefix = argv[4];
            int max_index = std::stoi(argv[5]);

            cache.restore(filepath);
            if (cache.size() == 0) {
                std::cerr << "Error: Cache is empty after restore!" << std::endl;
                return 1;
            }

            auto strategy = (TracEon::SmartStrategy*)cache.getStrategy();
            std::vector<std::string> all_keys = strategy->getAllKeys();
            
            if (all_keys.empty()) {
                std::cerr << "Error: No keys found." << std::endl;
                return 1;
            }

            std::vector<std::string> keys;
            keys.reserve(num_lookups);
            std::mt19937 rng(42);
            std::uniform_int_distribution<size_t> dist(0, all_keys.size() - 1);
            
            std::cout << "Generating " << num_lookups << " lookup keys..." << std::endl;
            for(int i=0; i<num_lookups; ++i) {
                keys.push_back(all_keys[dist(rng)]);
            }

            std::cout << "Validating keys..." << std::endl;
            size_t misses = 0;
            for (const auto& k : keys) {
                if (!cache.hasSequence(k)) misses++; 
            }
            std::cout << "Key misses: " << misses << "/" << keys.size() << "\n";
            if (misses > 0) {
                std::cerr << "CRITICAL WARNING: Benchmark will fail due to key mismatches!" << std::endl;
            }

            std::cout << "Benchmarking..." << std::endl;
            volatile size_t bytes = 0;
            Timer t;
            for (const auto& key : keys) {
                std::string_view val = cache.getView(key);
                size_t s = val.size();
                global_sink += s; 
                bytes += s; // <--- FIX 1: Correctly count bytes locally for reporting
            }
            double sec = t.elapsed_sec();

            // DEBUG: Print raw values
            std::cout << "DEBUG: Lookups=" << num_lookups << ", Time=" << std::fixed << std::setprecision(6) << sec << ", Sink=" << global_sink << std::endl;

            if (bytes == 0) {
                 std::cerr << "WARNING: Total bytes retrieved is 0. Benchmark may be invalid." << std::endl;
            }
            
            std::cout << "Lookup_Time_s: " << sec << std::endl;
            std::cout << "Data_Retrieved_MB: " << (bytes / 1024.0 / 1024.0) << std::endl;
            
            if (sec > 1e-9) {
                double ops = num_lookups / sec;
                // FIX 2: REMOVE " OPS/s" suffix so Python can parse the float
                std::cout << "Throughput: " << std::fixed << std::setprecision(2) << ops << std::endl;
            } else {
                std::cout << "Throughput: 0" << std::endl; // Return 0 instead of Inf string
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}