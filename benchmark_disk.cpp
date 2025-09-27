//
// Created by Adnan Raza on 9/27/2025.
//
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <chrono>

// A simple parser to find one sequence
std::string find_sequence_on_disk(const std::string& filepath, const std::string& target_id) {
    std::ifstream file(filepath);
    std::string line, current_id;
    bool reading_target = false;
    std::string sequence;

    while (std::getline(file, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            reading_target = false; // Reset on new header
            size_t first_space = line.find(' ');
            current_id = line.substr(1, first_space - 1);
            if (current_id == target_id) {
                reading_target = true;
            }
        } else if (reading_target) {
            sequence += line;
        }
    }
    return sequence;
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <fasta_file> <sequence_id>\n";
        return 1;
    }
    std::string filepath = argv[1];
    std::string target_id = argv[2];

    auto start = std::chrono::high_resolution_clock::now();

    std::string sequence = find_sequence_on_disk(filepath, target_id);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = end - start;

    std::cout << "Disk-based retrieval for '" << target_id << "' took: "
              << duration.count() << " ms\n";

    return 0;
}