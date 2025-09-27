#include "PlainTextStrategy.h"

PlainTextStrategy::~PlainTextStrategy() = default;

std::vector<unsigned char> PlainTextStrategy::encode(const std::string& data) const {
    return {data.begin(), data.end()};
}

std::string PlainTextStrategy::decode(const std::vector<unsigned char>& data) const {
    return {data.begin(), data.end()};
}