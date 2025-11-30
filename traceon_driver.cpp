#include <iostream>
#include <chrono>
#include <string>
#include <vector>
#include <cstdlib>
#include "Cache.h" 

// Simple timer class
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
    // Usage: ./traceon_driver <mode> <file> [args...]
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <mode> <file> [args]" << std::endl;
        return 1;
    }

    std::string mode = argv[1];
    std::string filepath = argv[2];
    TracEon::Cache cache;

    if (mode == "load") {
        Timer t;
        cache.loadFile(filepath);
        std::cout << "Load_Time_s: " << t.elapsed_sec() << std::endl;
        std::cout << "Sequences: " << cache.size() << std::endl;
    } 
    else if (mode == "save") {
        // Load first (untimed), then time the save
        cache.loadFile(filepath);
        if (argc < 4) return 1;
        std::string outfile = argv[3];
        
        Timer t;
        cache.save(outfile);
        std::cout << "Save_Time_s: " << t.elapsed_sec() << std::endl;
    }
    else if (mode == "restore") {
        Timer t;
        cache.restore(filepath); // filepath is the .bin file
        std::cout << "Restore_Time_s: " << t.elapsed_sec() << std::endl;
        std::cout << "Sequences: " << cache.size() << std::endl;
    }
    else if (mode == "lookup") {
        // Command: lookup <cache_bin> <num_lookups> <key_prefix> <max_index>
        if (argc < 6) return 1;
        int num_lookups = std::stoi(argv[3]);
        std::string prefix = argv[4];
        int max_index = std::stoi(argv[5]);

        // Restore first (setup cost, not counted in throughput)
        cache.restore(filepath);
        
        // Benchmark Random Access
        std::srand(42); 
        volatile size_t bytes = 0; // Prevent compiler optimization
        
        Timer t;
        for (int i = 0; i < num_lookups; ++i) {
            int id = std::rand() % max_index;
            std::string key = prefix + std::to_string(id);
            std::string val = cache.get(key);
            bytes += val.size();
        }
        double sec = t.elapsed_sec();
        
        std::cout << "Lookup_Time_s: " << sec << std::endl;
        std::cout << "Throughput: " << (num_lookups / sec) << std::endl;
    }
    else {
        std::cerr << "Unknown mode: " << mode << std::endl;
        return 1;
    }

    return 0;
}