// tools/traceon_make_cache.cpp — build a v4 binary cache from a FASTA/FASTQ.
//
// Small utility used by the Docker remote-access testbed (docker/remote-test)
// and by the manual verification flow:
//
//   traceon_make_cache <input.fasta|fastq[.gz]> <output.traceon>
//
// Usage in the testbed: the server container converts the shared test FASTA
// into a v4 cache once; the server then serves that cache over TCP while the
// client container benchmarks remote lookups.

#include "Cache.h"

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>

int main(int argc, char** argv) {
    if (argc != 3) {
        std::fprintf(stderr,
                     "Usage: %s <input.fasta|fastq[.gz]> <output.traceon>\n", argv[0]);
        return 2;
    }
    const std::string input = argv[1];
    const std::string output = argv[2];

    try {
        TracEon::Cache cache;
        cache.loadFile(input);
        if (cache.size() == 0) {
            std::fprintf(stderr, "traceon_make_cache: input parsed to 0 records\n");
            return 1;
        }
        cache.save(output);
        std::printf("wrote %s: %zu records, format=%d, arena_bytes=%zu\n",
                    output.c_str(), cache.size(),
                    static_cast<int>(cache.getDetectedFormat()),
                    cache.getArenaBytes());
        return 0;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "traceon_make_cache: error: %s\n", e.what());
        return 1;
    }
}