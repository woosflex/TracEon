// fuzz_kmer_api.cpp — libFuzzer target for the kmerindex C API
// (kmerindex_c_api.h / src/kmerindex_c_api.cpp).
//
// Contract under test (vuln-0003): no C++ exception may ever cross the C
// ABI boundary and no call may exhibit UB — with arbitrary keys/values, a
// NULL handle, NULL output pointers, and the freeze/pointer-stability
// state machine (get() blocks later inserts). Every operation here returns
// a status code; the harness drives a full create → reserve → insert/get →
// destroy lifecycle per input with bounds derived from the bytes. Any
// crash / leak / UB / thrown exception escaping the boundary is a finding.

#include "kmerindex_c_api.h"
#include "fuzz_common.h"

#include <cstdint>
#include <cstddef>

extern "C" int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
    size_t pos = 0;
    auto next_byte = [&]() -> uint8_t {
        return pos < size ? data[pos++] : 0x5A;
    };
    auto next64 = [&]() -> uint64_t {
        uint64_t v = 0;
        for (int i = 0; i < 8; ++i) v = (v << 8) | next_byte();
        return v;
    };

    kmerindex_t* h = kmerindex_create();
    if (!h) return 0; // allocation failure is a legal (null) result

    // Pre-reserve a bounded derived capacity; failure is legal (returns 0).
    (void)kmerindex_reserve(h, static_cast<size_t>(next64() % 65536));

    size_t ops = (size == 0) ? 4 : (1 + size / 16);
    if (ops > 2048) ops = 2048;
    for (size_t i = 0; i < ops; ++i) {
        const uint64_t key = next64();
        const uint64_t val = next64();
        switch (i & 3) {
        case 0: // get with the key just derived (get() freezes the index)
        {
            uint64_t matched = 0;
            const uint64_t* p = kmerindex_get(h, key, &matched);
            if (p) { volatile uint64_t sink = *p; (void)sink; }
            break;
        }
        case 1: // NULL output-pointer guard (legal error, never a crash)
            (void)kmerindex_get(h, key, nullptr);
            break;
        default: // plain insert (may fail once frozen — legal)
            (void)kmerindex_insert(h, key, val);
            break;
        }
        (void)kmerindex_last_error(h); // thread-local diagnostic, never NULL
    }

    (void)kmerindex_size(h);

    // NULL-handle coverage: every function must return its failure value.
    (void)kmerindex_reserve(nullptr, 1);
    (void)kmerindex_size(nullptr);
    (void)kmerindex_insert(nullptr, 1, 2);
    (void)kmerindex_get(nullptr, 1, nullptr);
    (void)kmerindex_last_error(nullptr);

    kmerindex_destroy(h);
    kmerindex_destroy(nullptr); // NULL-safe
    return 0;
}
