// fuzz_fasta.cpp — libFuzzer target for FASTA text parsing
// (SmartStrategy::loadFile with format detection, in-place multi-line
// normalization in normalizeFastaArena(), and both single-threaded and
// multithreaded (>10MB) parsers).
//
// Contract under test: arbitrary bytes as FASTA text must parse or be
// rejected — never a crash, an OOB read/write (normalizeFastaArena compacts
// in place), an infinite loop, or a leak. Interesting corners: multi-line
// sequences, header-only records, sequences abutting the next header (no
// newline), CRLF, missing trailing newline, no leading '>'.

#include "SmartStrategy.h"
#include "Cache.h"
#include "fuzz_common.h"

#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace {

void exercise_reads(TracEon::SmartStrategy& s) {
    volatile size_t sink = 0;
    const auto keys = s.getAllKeys();
    for (const auto& k : keys) {
        const std::string_view v = s.getView(k);
        for (char c : v) sink += static_cast<unsigned char>(c);
        sink += s.getSequence(k).size();
    }
    (void)sink;
}

} // namespace

extern "C" int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
    const std::string path = traceon_fuzz::write_input_to_tmp(data, size);
    // Alternate GENOME / NGS index modes (both have their own parse path).
    const bool ngs_mode = (size > 0) && ((data[0] & 1) != 0);
    try {
        TracEon::SmartStrategy s(ngs_mode ? TracEon::IndexMode::NGS
                                          : TracEon::IndexMode::GENOME);
        s.loadFile(path);
        exercise_reads(s);
    } catch (const std::exception&) {
        // Expected for malformed/truncated input.
    }
    return 0;
}
