//
// Created by Adnan Raza on 9/27/2025.
//
#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include "Cache.h"

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <fasta_file> <sequence_id>\n";
        return 1;
    }
    std::string filepath = argv[1];
    std::string target_id = argv[2];

    TracEon::Cache cache;

    // --- Phase 1: One-time setup cost ---
    auto load_start = std::chrono::high_resolution_clock::now();
    cache.loadFasta(filepath);
    auto load_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> load_duration = load_end - load_start;
    std::cout << "TracEon load time: " << load_duration.count() << " ms\n";

    // --- Phase 2: Retrieval from memory ---
    auto get_start = std::chrono::high_resolution_clock::now();
    std::string sequence = cache.get(target_id);
    auto get_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> get_duration = get_end - get_start;

    std::cout << "TracEon in-memory retrieval for '" << target_id << "' took: "
              << get_duration.count() << " ms\n";

    return 0;
}