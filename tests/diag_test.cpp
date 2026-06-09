#include <cstdio>
#include <cstdint>
#include <fstream>
#include <string>
#include <cstring>
#include "SmartStrategy.h"

int main() {
    // Test: FASTA last record without trailing newline
    std::string test_file = "/tmp/diag_no_trail.fasta";
    {
        std::ofstream out(test_file, std::ios::binary);
        out << ">seq1\nACGTACGTACGTACGT\n>seq2\nGGGGCCCC";
    }

    TracEon::SmartStrategy strategy;
    strategy.loadFile(test_file);

    // Get seq1 and seq2
    auto v1 = strategy.getView("seq1");
    auto v2 = strategy.getView("seq2");

    printf("seq1 SIZE=%zu HEX:", v1.size());
    for (size_t i = 0; i < v1.size(); i++) printf(" %02x", (unsigned char)v1[i]);
    printf("\n");

    printf("seq2 SIZE=%zu HEX:", v2.size());
    for (size_t i = 0; i < v2.size(); i++) printf(" %02x", (unsigned char)v2[i]);
    printf("\n");

    // Check the actual bytes around seq2's data to see what's before and after
    printf("\nBytes around seq2 in arena:\n");
    // Walk backwards/forwards from seq2 start to find boundaries
    const char* seq2_start = v2.data();
    printf("seq2_start points to arena offset %td\n", seq2_start - v1.data());
    
    // Print some context bytes
    printf("Bytes before seq2 (offset -5 to -1):");
    for (int i = -5; i < 0; i++) printf(" %02x", (unsigned char)seq2_start[i]);
    printf("\n");
    
    printf("seq2 itself (offset 0 to %zu):", v2.size());
    for (size_t i = 0; i < v2.size(); i++) printf(" %02x", (unsigned char)seq2_start[i]);
    printf("\n");
    
    printf("Bytes after seq2 (offset +0 to +5):");
    for (int i = 0; i < 6; i++) printf(" %02x", (unsigned char)seq2_start[v2.size() + i]);
    printf("\n");

    std::remove(test_file.c_str());
    
    bool pass = (v2.size() == 8);
    for (size_t i = 0; i < v2.size(); i++) if (v2[i] == '\0') pass = false;
    printf("\n%s\n", pass ? "PASS" : "FAIL");
    return pass ? 0 : 1;
}
