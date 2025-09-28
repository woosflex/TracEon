#include "QualityScoreRleStrategy.h"

QualityScoreRleStrategy::~QualityScoreRleStrategy() = default;

// Compresses a string using Run-Length Encoding.
std::vector<unsigned char> QualityScoreRleStrategy::encode(const std::string& data) const {
    if (data.empty()) {
        return {};
    }

    std::vector<unsigned char> encoded_data;
    for (size_t i = 0; i < data.length(); ++i) {
        unsigned char current_char = data[i];
        unsigned char count = 1;
        // Count how many times the character repeats.
        while (i + 1 < data.length() && data[i + 1] == current_char && count < 255) {
            count++;
            i++;
        }
        // Store the count, then the character.
        encoded_data.push_back(count);
        encoded_data.push_back(current_char);
    }
    return encoded_data;
}

// Decompresses a string using Run-Length Encoding.
std::string QualityScoreRleStrategy::decode(const std::vector<unsigned char>& data) const {
    std::string decoded_string;
    // Data is stored in pairs: [count, char, count, char, ...]
    for (size_t i = 0; i < data.size(); i += 2) {
        unsigned char count = data[i];
        unsigned char character = data[i + 1];
        decoded_string.append(count, character);
    }
    return decoded_string;
}