// 80k/80k multithreaded-FASTQ record-preservation repro driver.
//
// Verifies the strict-4-line chunk-boundary classifier (Bug 1 fix, P0):
// a file whose quality lines all START with '+' (Phred+33 Q10) must load
// ALL 80k records — the old backward-looking classifier dropped every
// chunk after chunk 0 (80k -> 4,468). Also verifies the control file
// (normal quality) loads 80k/80k.
//
// Build: g++ -O2 -std=c++20 -I../include -I../build/_deps/ankerl_unordered_dense-src/include \
//        repro_driver.cpp -o repro_driver ../build/libtraceon_core.a -lz -llz4 -pthread
// Run (from outputs/): ./repro_driver
#include "Cache.h"
#include "SmartStrategy.h"
#include <cstdio>
#include <string>

int main() {
    constexpr size_t EXPECT = 80000;
    int failures = 0;

    for (const char* path : {"repro_ctrl.fastq", "repro_plus.fastq"}) {
        TracEon::Cache cache;
        cache.loadFile(path);
        const size_t n = cache.size();
        const bool ok = (n == EXPECT);
        printf("%-18s loaded %zu (expected %zu) -> %s\n", path, n, EXPECT,
               ok ? "PASS" : "FAIL");
        if (!ok) ++failures;
    }

    // Spot-check the last record of the '+'-quality file through the API.
    TracEon::Cache plus;
    plus.loadFile("repro_plus.fastq");
    if (!plus.hasSequence("read79999")) { printf("read79999 missing -> FAIL\n"); ++failures; }
    else {
        auto rec = plus.getFastqRecord("read79999");
        const bool have = rec.has_value();
        const size_t sl = have ? rec->sequence.size() : 0;
        const size_t ql = have ? rec->quality.size() : 0;
        const char q0 = (have && !rec->quality.empty()) ? rec->quality[0] : '?';
        printf("read79999 seq_len=%zu qual_len=%zu qual[0]=%c -> %s\n", sl, ql, q0,
               (have && sl == 100 && ql == 100 && q0 == '+') ? "PASS" : "FAIL");
        if (!have || sl != 100 || ql != 100 || q0 != '+') ++failures;
    }

    printf(failures == 0 ? "REPRO RESULT: 80000/80000 PASS\n" : "REPRO RESULT: %d FAILURES\n", failures);
    return failures == 0 ? 0 : 1;
}
