#include "TwoBitDnaStrategy.h"
#include <cmath>
#include <cstdint>

TwoBitDnaStrategy::~TwoBitDnaStrategy() = default;

namespace { // Use an anonymous namespace for internal helper functions
    // Maps a DNA base to its 2-bit value.
    inline unsigned char baseToBits(char base) {
        switch (base) {
            case 'A': return 0b00;
            case 'C': return 0b01;
            case 'G': return 0b10;
            case 'T': return 0b11;
            default:  return 0b00; // Let 'N' and others default to 'A'
        }
    }

    // Maps a 2-bit value back to a DNA base.
    inline char bitsToBase(unsigned char bits) {
        switch (bits) {
            case 0b00: return 'A';
            case 0b01: return 'C';
            case 0b10: return 'G';
            case 0b11: return 'T';
            default:   return 'N';
        }
    }
}


std::vector<unsigned char> TwoBitDnaStrategy::encode(const std::string& data) const {
    if (data.empty()) return {};

    uint32_t original_length = data.length();
    // We need space for 4 bytes of length, plus the packed data.
    size_t packed_data_size = std::ceil(static_cast<double>(original_length) / 4.0);
    std::vector<unsigned char> encoded_data(4 + packed_data_size, 0);

    // Manually write the 32-bit length into the first 4 bytes of the vector.
    encoded_data[0] = (original_length >> 24) & 0xFF;
    encoded_data[1] = (original_length >> 16) & 0xFF;
    encoded_data[2] = (original_length >> 8) & 0xFF;
    encoded_data[3] = original_length & 0xFF;

    // Pack 4 bases into each subsequent byte.
    for (size_t i = 0; i < original_length; ++i) {
        unsigned char two_bit_value = baseToBits(data[i]);
        size_t byte_index = 4 + (i / 4);
        // First base goes into the most significant bits (7-6), second into (5-4), etc.
        int bit_shift = (3 - (i % 4)) * 2;
        encoded_data[byte_index] |= (two_bit_value << bit_shift);
    }
    return encoded_data;
}

std::string TwoBitDnaStrategy::decode(const std::vector<unsigned char>& data) const {
    if (data.size() < 4) return ""; // Not enough data for the length header

    // Reconstruct the original 32-bit length from the first 4 bytes.
    uint32_t original_length = (static_cast<uint32_t>(data[0]) << 24) |
                               (static_cast<uint32_t>(data[1]) << 16) |
                               (static_cast<uint32_t>(data[2]) << 8)  |
                               static_cast<uint32_t>(data[3]);

    std::string decoded_dna;
    decoded_dna.reserve(original_length);

    // Unpack exactly 'original_length' bases.
    for (size_t i = 0; i < original_length; ++i) {
        size_t byte_index = 4 + (i / 4);
        int bit_shift = (3 - (i % 4)) * 2;
        // Shift the byte to align the 2 bits we want, then mask them off.
        unsigned char two_bit_value = (data[byte_index] >> bit_shift) & 0b11;
        decoded_dna += bitsToBase(two_bit_value);
    }
    return decoded_dna;
}