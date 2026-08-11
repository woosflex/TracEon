// fuzz_fastq.cpp — libFuzzer target for FASTQ text parsing
// (SmartStrategy::loadFile with format detection + both the single-threaded
// and multithreaded (>10MB) strict-4-line parsers).
//
// Contract under test: arbitrary bytes as FASTQ text must parse or be
// rejected — never a crash, an OOB read, an infinite loop, or a leak. The
// strict 4-line framing rules (empty seq/qual lines, '+'-leading quality,
// '@' inside quality, CRLF, missing trailing newline) are the interesting
// corners; the seed corpus carries crafted variants of each.
//
// Note: inputs starting with gzip magic (0x1f 0x8b) are routed by loadFile's
// detection to the gzip path — real-world behavior, still under test here.

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
        sink += s.getQuality(k).size();
    }
    (void)sink;
}

} // namespace

extern "C" int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
    const std::string path = traceon_fuzz::write_input_to_tmp(data, size);
    // Alternate GENOME / NGS index modes: both map types have their own
    // parse + read paths.
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
