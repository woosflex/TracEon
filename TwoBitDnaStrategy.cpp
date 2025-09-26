//
// Created by Adnan Raza on 9/26/2025.
//

#include "TwoBitDnaStrategy.h"
#include <stdexcept> // For std::runtime_error
#include <cmath>     // For std::ceil

// Helper function to convert a base to its 2-bit representation.
inline unsigned char baseToBits(char base) {
    switch (base) {
        case 'A': return 0b00; // 0
        case 'C': return 0b01; // 1
        case 'G': return 0b10; // 2
        case 'T': return 0b11; // 3
        default:  return 0b00; // Default for 'N' etc.
    }
}

// Helper function to convert 2 bits back to a base character.
inline char bitsToBase(unsigned char bits) {
    switch (bits) {
        case 0b00: return 'A';
        case 0b01: return 'C';
        case 0b10: return 'G';
        case 0b11: return 'T';
        default:   return 'N';
    }
}


std::vector<unsigned char> TwoBitDnaStrategy::encode(const std::string& data) const {
    if (data.empty()) {
        return {};
    }

    // We'll store the original length in the first 4 bytes.
    uint32_t original_length = data.length();

    // Calculate how many bytes we need for the packed data. Each byte holds 4 bases.
    size_t packed_data_size = std::ceil(static_cast<double>(original_length) / 4.0);

    // Total size is 4 bytes for the length + the packed data size.
    std::vector<unsigned char> encoded_data(4 + packed_data_size, 0);

    // Write the length into the first 4 bytes.
    encoded_data[0] = (original_length >> 24) & 0xFF;
    encoded_data[1] = (original_length >> 16) & 0xFF;
    encoded_data[2] = (original_length >> 8) & 0xFF;
    encoded_data[3] = original_length & 0xFF;

    // Now, pack the DNA sequence.
    for (size_t i = 0; i < original_length; ++i) {
        unsigned char two_bit_value = baseToBits(data[i]);
        size_t byte_index = 4 + (i / 4);
        // We pack from left to right (e.g., first base goes into bits 7-6)
        int bit_shift = (3 - (i % 4)) * 2;
        encoded_data[byte_index] |= (two_bit_value << bit_shift);
    }

    return encoded_data;
}


std::string TwoBitDnaStrategy::decode(const std::vector<unsigned char>& data) const {
    if (data.size() < 4) {
        // Not enough data to even have the length, something is wrong.
        return "";
    }

    // First, read the original length from the first 4 bytes.
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
        // Shift the bits to the right, then mask to get only the 2 we care about.
        unsigned char two_bit_value = (data[byte_index] >> bit_shift) & 0b11;
        decoded_dna += bitsToBase(two_bit_value);
    }

    return decoded_dna;
}