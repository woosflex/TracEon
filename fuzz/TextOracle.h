#pragma once

#include "SmartStrategy.h"
#include <cstdlib>
#include <map>
#include <sstream>
#include <string>

// Independent line-oriented oracle for accepted, uncompressed text. The core
// parser uses pointer spans/SIMD and parallel chunks; compare every ID/payload.
inline void verifyText(const uint8_t* data, size_t size, TracEon::SmartStrategy& cache) {
    if (size >= 2 && data[0] == 0x1f && data[1] == 0x8b) return;
    std::map<std::string, std::pair<std::string, std::string>> expected;
    std::istringstream input(std::string(reinterpret_cast<const char*>(data), size));
    auto line = [&](std::string& out) {
        if (!std::getline(input, out)) return false;
        if (!out.empty() && out.back() == '\r') out.pop_back();
        return true;
    };
    auto id = [](const std::string& header) {
        const auto end = header.find_first_of(" \t\r", 1);
        return header.substr(1, end == std::string::npos ? end : end - 1);
    };
    std::string header;
    while (line(header) && header.empty()) {}
    if (!header.empty() && header[0] == '@') {
        do {
            std::string sequence, plus, quality;
            if (!line(sequence) || !line(plus) || !line(quality)) std::abort();
            expected.emplace(id(header), std::make_pair(sequence, quality));
        } while (line(header));
    } else if (!header.empty() && header[0] == '>') {
        std::string key = id(header), sequence, next;
        auto record = [&] {
            while (!sequence.empty() && sequence.back() == ' ')
                sequence.pop_back();
            expected.emplace(key, std::make_pair(sequence, std::string{}));
        };
        while (line(next)) {
            if (!next.empty() && next[0] == '>') {
                record();
                key = id(next);
                sequence.clear();
            } else {
                sequence += next;
            }
        }
        record();
    }
    if (cache.getFileCacheSize() != expected.size()) std::abort();
    for (const auto& [key, payload] : expected) {
        if (!cache.hasSequence(key) || cache.getView(key) != payload.first ||
            cache.getQuality(key) != payload.second) std::abort();
    }
}
