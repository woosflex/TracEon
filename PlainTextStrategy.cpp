#include "PlainTextStrategy.h"

// This is our simplest strategy, it just handles the basic conversion
// from a string to a byte vector for storage.
std::vector<unsigned char> PlainTextStrategy::encode(const std::string& data) const {
    // Just convert the string into a vector of its raw bytes.
    return {data.begin(), data.end()};
}

// And this just converts the byte vector right back into a string.
std::string PlainTextStrategy::decode(const std::vector<unsigned char>& data) const {
    return {data.begin(), data.end()};
}