//
// Created by Adnan Raza on 9/26/2025.
//

#ifndef TRACEON_RECORDTYPES_H
#define TRACEON_RECORDTYPES_H

#include <string>
#include <vector>

struct FastqRecord {
    std::string identifier;
    std::vector<unsigned char> compressed_sequence;
    std::vector<unsigned char> compressed_quality;
};

#endif //TRACEON_RECORDTYPES_H