#include <cstdio>
#include <vector>
int main() {
    // Simulate what normalizeFastaArena does:
    // 1. reserve(original_size + 1) for capacity
    // 2. Write past size via raw pointer
    // 3. resize() to include the written bytes
    
    std::vector<char> v;
    v.resize(5, 'A');       // size=5, "AAAAA"
    
    const size_t original_size = v.size();
    v.reserve(original_size + 1);  // capacity >= 6
    char* write = v.data() + original_size;
    *write++ = 'B';  // write past end: v.data()[5] = 'B'
    
    printf("Before resize: size=%zu, capacity=%zu\n", v.size(), v.capacity());
    printf("data()[5] via raw ptr: %02x\n", (unsigned char)*(v.data() + 5));
    
    v.resize(6);
    
    printf("After resize(6): size=%zu\n", v.size());
    printf("v[5] = %02x '%c'\n", (unsigned char)v[5], v[5]);
    
    // Check if resize zeroed the byte we wrote
    if (v[5] == 'B') printf("PASS: byte preserved\n");
    else if (v[5] == '\0') printf("FAIL: resize zeroed the byte!\n");
    else printf("UNEXPECTED: v[5]=%02x\n", (unsigned char)v[5]);
    
    return 0;
}
