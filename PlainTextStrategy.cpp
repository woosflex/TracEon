//
// Created by Adnan Raza on 9/26/2025.
//

#include "PlainTextStrategy.h"

std::vector<unsigned char> PlainTextStrategy::encode(const std::string& data) const {
    // Convert the string to a vector of its raw byte data.
    return {data.begin(), data.end()};
}

std::string PlainTextStrategy::decode(const std::vector<unsigned char>& data) const {
    // Convert the vector of bytes back into a string.
    return {data.begin(), data.end()};
}