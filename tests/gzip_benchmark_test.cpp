/**
 * @file gzip_benchmark_test.cpp
 * @brief Synthetic benchmark instrumenting decompress, parse, and total load
 *        times using std::chrono for each component.
 *
 * Generates FASTA records, writes .gz via zlib-ng, then measures:
 *   1. Decompression (gzread into arena)
 *   2. Parsing (parseArena into map)
 *   3. Total (both combined)
 *
 * Reports JSON-per-line results.
 *
 * Compile: g++ -std=c++20 -O2 -I../include -I../third_party/lz4/lib \
 *         -L../build -o gzip_benchmark_test gzip_benchmark_test.cpp \
 *         -ltraceon_core -lz -llz4
 */

#include <chrono>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <cstring>
#include <memory>
#include <random>
#include <filesystem>
#include <zlib.h>
#include <cstdio>

// ============================================================
// Wall clock timer using C++20 high_resolution_clock (chrono)
// ============================================================
class BenchTimer {
    using Clock = std::chrono::high_resolution_clock;
    Clock::time_point start_;
public:
    BenchTimer() : start_(Clock::now()) {}
    double elapsed_sec() const {
        return std::chrono::duration<double>(Clock::now() - start_).count();
    }
    void reset() { start_ = Clock::now(); }
};

// ============================================================
// Phase 1: Pure decompression via gzread
// ============================================================
struct DecompTiming {
    double seconds;    ///< wall time for decompression alone
    size_t raw_bytes;  ///< uncompressed byte count
};

static DecompTiming time_decompression(const std::string& gz_path,
                                       size_t chunk = 1u << 20) {
    gzFile f = gzopen(gz_path.c_str(), "rb");
    if (!f) { std::cerr << "FAIL: gzopen " << gz_path << "\n"; return {0,0}; }
    auto buf = std::make_unique<char[]>(chunk);
    size_t total = 0;
    BenchTimer t;
    int n;
    while ((n = gzread(f, buf.get(), chunk)) > 0) total += static_cast<size_t>(n);
    double sec = t.elapsed_sec();
    if (n < 0) { int e; std::cerr << "gzread err: " << gzerror(f, &e) << "\n"; }
    gzclose(f);
    return {sec, total};
}

// ============================================================
// Phase 2: Full load via traceon_core (decompress + parse)
// ============================================================
#include "Cache.h"

struct FullTiming {
    double seconds;
    size_t records;
    size_t raw_bytes;
};

static FullTiming time_full_load(const std::string& path) {
    // Pre-measure decompressed size
    size_t raw_est = 0;
    {
        gzFile gz = gzopen(path.c_str(), "rb");
        if (gz) { char b[65536]; int n; while ((n = gzread(gz, b, sizeof(b))) > 0) raw_est += n; gzclose(gz); }
    }

    TracEon::Cache cache;
    BenchTimer t;
    cache.loadFile(path);
    double sec = t.elapsed_sec();

    return {sec, cache.size(), raw_est};
}

// ============================================================
// Synthetic FASTA generator
// ============================================================
static void generate_fasta_gz(const std::string& outpath,
                               size_t target_bytes,
                               size_t seq_len = 1000) {
    std::mt19937 rng(42);
    static const char bases[] = "ACGT";
    std::uniform_int_distribution<int> dist(0, 3);

    // Header overhead: ">seqN\n" + "\n" after each sequence
    const size_t bytes_per_seq = seq_len + 8;  // approx
    size_t num_seqs = std::max<size_t>(1, target_bytes / bytes_per_seq);

    // Buffer for sequence generation
    std::string seq;
    seq.reserve(seq_len);
    for (size_t i = 0; i < seq_len; ++i)
        seq.push_back(bases[dist(rng)]);

    // Determine resulting uncompressed size
    std::string uncompressed;
    for (size_t i = 0; i < num_seqs; ++i) {
        uncompressed += ">seq" + std::to_string(i) + "\n";
        uncompressed += seq;
        uncompressed += "\n";
    }

    // Write gzip
    gzFile gf = gzopen(outpath.c_str(), "wb9");
    if (!gf) { std::cerr << "FAIL: gzopen write " << outpath << "\n"; return; }
    gzwrite(gf, uncompressed.data(), uncompressed.size());
    gzclose(gf);
}

// ============================================================
int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <fasta.gz> [--generate <target_bytes>]\n";
        return 1;
    }

    std::string path = argv[1];
    bool gen_mode = (argc >= 4 && std::string(argv[2]) == "--generate");
    if (gen_mode) {
        size_t target = static_cast<size_t>(std::stoll(argv[3]));
        generate_fasta_gz(path, target);
        std::cout << "Generated: " << path << " (~" << (target >> 20) << " MB)\n";
        return 0;
    }

    if (!std::filesystem::exists(path)) {
        std::cerr << "File not found: " << path << "\n";
        return 1;
    }

    size_t comp_size = std::filesystem::file_size(path);

    // --- DECOMPRESSION ONLY ---
    DecompTiming dec = time_decompression(path);

    // --- FULL LOAD (decompress + parse) ---
    FullTiming full = time_full_load(path);

    // --- DERIVE PARSE TIME ---
    double parse_s = full.seconds - dec.seconds;
    if (parse_s < 0.0) parse_s = 0.0;

    // --- JSON OUTPUT ---
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "{"
              << "\"file\":\""             << std::filesystem::path(path).filename().string() << "\","
              << "\"compressed_bytes\":"   << comp_size << ","
              << "\"raw_bytes\":"          << full.raw_bytes << ","
              << "\"records\":"            << full.records << ","
              << "\"decompress_s\":"       << dec.seconds << ","
              << "\"parse_s\":"            << parse_s << ","
              << "\"total_s\":"            << full.seconds
              << "}\n";

    return 0;
}
