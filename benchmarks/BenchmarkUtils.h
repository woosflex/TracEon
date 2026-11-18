#pragma once

#include <charconv>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <string_view>

inline uint64_t benchmarkUnsigned(std::string_view value, std::string_view option) {
    uint64_t result = 0;
    const auto [end, error] = std::from_chars(value.data(), value.data() + value.size(), result);
    if (error != std::errc{} || end != value.data() + value.size())
        throw std::invalid_argument(std::string(option) + " requires an unsigned integer");
    return result;
}
